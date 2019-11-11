import numpy as np
import tarfile, os, urllib, gzip, re, time, requests
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup

defaultAgent = {'User-Agent': 'Mozilla/5.0 (compatible; Googlebot/2.1; +http://www.google.com/bot.html)'}

numberDict = {
    'zero': 0,
    'one': 1,
    'two': 2,
    'three': 3,
    'four': 4,
    'five': 5,
    'six': 6,
    'seven': 7,
    'eight': 8,
    'nine': 9,
    'ten': 10,
    'eleven': 11,
    'twelve': 12,
    'thirteen': 13,
    'fourteen': 14,
    'fifteen': 15,
    'sixteen': 16,
    'seventeen': 17,
    'eighteen': 18,
    'nineteen': 19,
    'twenty': 20,
    'thirty': 30,
    'forty': 40,
    'fifty': 50,
    'sixty': 60,
    'seventy': 70,
    'eighty': 80,
    'ninety': 90,
    'hundred': 100,
    'thousand': 1000,
    'million': 1000000,
    'billion': 1000000000,
    'point': '.'
}


def tar_gz_extracter(url):
    """
    This function will take a FTP URL, download then extract the tar-gzipped file (.tgz).
    The first try: except: loop will catch the rare `ContentTooShortError`. I believe this appears in very big runs but its so infrequent, I haven't been able to fully diagnosis it.
    The second try: except: loop will trigger when there isn't an .xml file in the .tgz file. Currently it prints the files as I was debugging but this probably isn't necessary. Would be good to mark it somehow..?

        Args:
            `url` - Str: Series FTP URL leading to a .tgz file

        Returns:
            `xml_name` - Str: Name of .xml file that was extracted.
    """

    try:
        f_url = urllib.request.urlretrieve(url, filename=None)[0]
    except urllib.error.ContentTooShortError as e:
        print(f'URL too short: {url}')
        return None
    except urllib.error.URLError as e:
        print(f'waiting 60 seconds and retrying')
        time.sleep(61)
        try:
            f_url = urllib.request.urlretrieve(url, filename=None)[0]
        except:
            print('tgz extraction failed after waiting.')
            return None
    base_name = os.path.basename(url)

    file_name, file_extension = os.path.splitext(base_name)
    tar = tarfile.open(f_url)

    try:
        xml_name = [i for i in tar.getnames() if re.match('.*.xml', str(i))][0]
    except:
        xml_name = ''
        print(base_name)
        print(tar.getnames())
        return None
    tar.extract(xml_name, path = f'output/xml/')

    tar.close()

    return xml_name



def geoAgeExtract(sample_dict, checkCell = True,
    parseAgeIDs = ['characteristics', 'description', 'treatment_protocol', 'growth_protocol'], convertTo = 'week',
    nullReturn = 'n/a', checkConverts = ['day', 'week', 'month', 'year'],
    tryAgeStudy = True, parseStudyIDs = ['series_summary', 'series_design'],
    tryAgePMID = True, pmidSection = 'methods'):
    """ Extract age from GEO sample (GSM) text values in the `sample_dict`. Either/or the 'characteristics'
        or protocols/descriptions entries. See numericTimeConvert() docstring for
        more detail on the approach. Return null on cell detect. If the same age
        is detected in both 'Characteristics' and descriptions/protocols, then
        assume that they are duplicated, and only add non-matching ages in the
        descriptions/protocols. NumberDict is an english word to integer
        conversion dictionary pre-defined in globals().

        Args:
            sample_dict: Dictionary - dictionary of a sample's data
            parseAgeIDs: List - Entries within url text
            convertTo: Str - Converted time units ['day', 'week', 'month', 'year']
            nullReturn: Str - Return if no numbers found
            checkConverts: List - Time units to convert to convertTo unit
            tryAgeStudy: Bool - If no age is found, attempt to find age in the
                GSE study page
            parseStudyIDs: List - Entries within url text if GSE study is checked
            tryAgePMID: Bool - Resort to full text extraction if no age can be
                extracted some GSE nor GSM
            pmidSection: Str - Full-text paper section in which to expect an age
                value. Default to 'methods'
        Return:
            age: Float - Estimated age
            source: Str - Source of age scrape('Sample', 'Study', 'Text')
    """
    if (sample_dict['cells'] is True) and (checkCell == True):
        return nullReturn, nullReturn

    ageCounter, charAge = [], nullReturn

    for i in parseAgeIDs:
        if (i == 'characteristics') and (sample_dict['sample_age'] != ''):
            charAge = numericTimeConvert(text = sample_dict['sample_age'], convertTo = convertTo, nullReturn = nullReturn)
        elif (i == 'characteristics') and (sample_dict['sample_age'] == ''):
            charAge = nullReturn
        elif i != 'characteristics':
            cleanText = re.sub(r'<.*?>|\\n', ' ', sample_dict[i])
            ageCounter.append(numericTimeConvert(text = cleanText, convertTo = convertTo, nullReturn = nullReturn))

    ageCounter = [x for x in ageCounter if x != nullReturn]

    if charAge != nullReturn:
        if charAge in ageCounter:
            return sum(ageCounter), 'Sample'
        elif charAge not in ageCounter and len(ageCounter) > 0:
            return sum([charAge, sum(ageCounter)]), 'Sample'
        else:
            return charAge, 'Sample'
    elif charAge == nullReturn and len(ageCounter) > 0:
        return sum(ageCounter), 'Sample'
    else: ### charAge == nullReturn
        if tryAgeStudy is True:
            gseAttempt, flagged = gseAgeExtract(sample_dict = sample_dict, convertTo = convertTo,
                checkConverts = checkConverts, parseIDs = parseStudyIDs,
                nullReturn = nullReturn)
            if gseAttempt == nullReturn and tryAgePMID is True:
                pmidAttempt = pmidAgeExtract(sample_dict = sample_dict, sectionID = pmidSection,
                    convertTo = convertTo, checkConverts = checkConverts,
                    nullReturn =  nullReturn)
                return pmidAttempt, 'Text'

            elif gseAttempt == nullReturn and tryAgePMID is False:
                return nullReturn, 'Study'
            else:
                return gseAttempt, 'Study'
        else:
            return nullReturn, 'Sample'

def numericTimeConvert(text, convertTo = 'week', checkConverts = ['day', 'week',
    'month', 'year'], nullReturn = 'n/a'):
    """ Convert text into a numerical time. Some general rules: Check will be
        made for a range of numbers, e.g. 12-13, and averaged. If no match, then
        check for single numbers, e.g. 12. However, if "for" is in front of this
        number, then it needs to be added to other numeric matches (e.g. "12
        weeks old for 2 weeks" == 14). This syntax is ideal, however one could
        possibly check that the 12 != 2, with the expectation that the
        likelihood is low that it would read "7 weeks old for 7 weeks".
        These checks can be added in the future. Currently only able to process
        days, weeks, months, and years (shorter time units are typically reserved
        for exposures, such as a drug, rather than the animal's age).
        Args:
            text: Str - text to be parsed
            convertTo: Str - Converted time units ['day', 'week', 'month', 'year']
            nullReturn: Str - Return if no numbers found
            checkConverts: List - Time units to convert to convertTo unit

        Return:
            nums - Float - Converted number in the text (to weeks)

    TODO: add sentence split, excl terms (cell-stuff)
    """
    if text == nullReturn:
        return nullReturn


    checkConverts = [re.sub('s$', '', x.lower()) for x in checkConverts]
    convertTo = re.sub('s$', '', convertTo.lower())

    if (any([x for x in checkConverts if x not in ['day', 'week', 'month', 'year']]) or
        convertTo not in ['day', 'week', 'month', 'year']):
        raise ValueError("Times must be in ['day', 'week', 'month', 'year']")

    text = re.sub('(\D)\-(\D)', '\\1 \\2', text) #keep '7-8', convert 'seven-week'
    text = re.sub('(\d+)\-(\D)', '\\1 \\2', text) #keep '7-8', convert '7-week'

    splitWords = text.lower().strip().split()
    strToNumConvert = []
    for word in splitWords:
        if word in numberDict:
            strToNumConvert.append(str(numberDict[word]))
        else:
            strToNumConvert.append(word)
    strToNumConvert = ' '.join(strToNumConvert)

    nums = 0
    for convertFrom in checkConverts:
        nums = addAgeStrings(nums, strToNumConvert = strToNumConvert, convertFrom = convertFrom, convertTo = convertTo)

    if nums == 0:
        return nullReturn
    else:
        return nums


def addAgeStrings(nums, strToNumConvert, convertFrom = 'day', convertTo = 'week'):
    """ Find instances of time, and convert to a quantitative value (perhaps of
        another time unit). Add time matches to counter, converted. Checks made
        for durations, e.g. "for 2 weeks", and if found, these are added to nums.
        For other matches, checks are made for multiple numbers (12-13 weeks),
        and if found, these are averaged. Assume no ranges are reported for the
        duration (e.g. no instances of "adminstered for 2-3 weeks"). Arguments
        for convertFrom and convertTo are checked in numericTimeConvert().
        Args:
            nums: Float - counter (in weeks)
            strToNumConvert: Str - string, presumably converted words to integers
            convertFrom: Str - Any of ['day', 'week', 'month', 'year]
            convertTo: Str - Any of ['day', 'week', 'month', 'year]
        Returns:
            Updated nums
    """
    convertCoef = timeConversions(convertFrom = convertFrom, convertTo = convertTo)

    if convertFrom == 'day':
        re1, re2 = '[ \-]day', r'[ \-]d[ \.s\,]'
    if convertFrom == 'week':
        re1, re2 = '[ \-]week', r'[ \-]wk[ \.s\,]'
    if convertFrom == 'month':
        re1, re2 = '[ \-]month', r'[ \-]mo[ \.s\,]'
    if convertFrom == 'year':
        re1, re2 = '[ \-]year', r'[ \-]yr[ \.s\,]'

    timeMatches = re.compile(f'(?<!for )(\d+\-?\d+?){re1}|(?<!for )(\d+\-?\d+?){re2}|(?<!for )(\d+){re1}|(?<!for )(\d+){re2}|(?<!for )(\d+\.?\d+?){re1}|(?<!for )(\d+\.?\d+?){re2}')
    timeMatchesDur = re.compile(f'for (\d+\-?\d+?){re1}|for (\d+\-?\d+?){re2}|for (\d+){re1}|for (\d+){re2}|for (\d+\.?\d+?){re1}|for (\d+\.?\d+?){re2}')

    times = set(re.findall(timeMatches, strToNumConvert))
    if re.findall(timeMatchesDur, strToNumConvert):
        timeDurs = set(re.findall(timeMatchesDur, strToNumConvert)[0])
    else:
        timeDurs = []

    for match in times:
        for x in match:
            if len(x.split('-')) == 2:
                try:
                    nums += np.mean([float(n) for n in x.split('-') if float(n) not in timeDurs])*convertCoef
                except ValueError:
                    continue
                except Exception as err:
                    print(str(err))
            else:
                try:
                    if x not in timeDurs:
                        nums += float(x)*convertCoef
                except ValueError:
                    continue
                except Exception as err:
                    print(str(err))

    for match in timeDurs:
        for x in match:
            if len(x.split('-')) == 2:
                try:
                    nums += np.mean([float(n) for n in x.split('-')])*convertCoef
                except ValueError:
                    continue
                except Exception as err:
                    print('found me')
                    print(str(err))
            else:
                try:
                    nums += int(x)*convertCoef
                except ValueError:
                    continue
                except Exception as err:
                    print(str(err))

    return nums


def timeConversions(convertFrom, convertTo):
    """ Return coefficients of converting from e.g. days to weeks. Inherits from
        addAgeStrings. Some matches are approximate (30 days per month) """

    if convertFrom == convertTo:
        return 1
    if convertFrom == 'day' and convertTo == 'week':
        return (1/7)
    if convertFrom == 'day' and convertTo == 'month':
        return (1/30)
    if convertFrom == 'day' and convertTo == 'year':
        return (1/365)
    if convertFrom == 'week' and convertTo == 'month':
        return (1/4)
    if convertFrom == 'week' and convertTo == 'year':
        return (1/52)
    if convertFrom == 'week' and convertTo == 'day':
        return 7
    if convertFrom == 'month' and convertTo == 'day':
        return 30
    if convertFrom == 'year' and convertTo == 'day':
        return 365
    if convertFrom == 'month' and convertTo == 'week':
        return 4
    if convertFrom == 'year' and convertTo == 'week':
        return 52



def gseAgeExtract(sample_dict, convertTo = 'week',
    checkConverts = ['day', 'week', 'month', 'year'],
    parseIDs = ['series_summary', 'series_design'], nullReturn = 'n/a',
    flagRange = True):
    """ Akin to extracting age from individual sample descriptions, extract age
        from GSE study information.
        Return:
            age: Float - Estimated age
            flagged: Bool - Wide age range flag, if called for by flagRange arg
    TODO: Different ages may exist for different experimental groups! Check
        GSE7191 as an example. Currently this function (incorrectly) combines
        and adds the ages of each group

    """
    flagged = False
    ageCounter, charAge, flags = [], nullReturn, []
    for i in parseIDs:
        #try:
        matches = sample_dict[i]
        if len(matches) != '':
            iAge = numericTimeConvert(text = sample_dict[i],
                convertTo = convertTo, nullReturn = nullReturn)
            ageCounter.append(iAge)
            #flags.append(iRange)
        #except KeyError:
        #    pass

    ageCounter = [x for x in ageCounter if x != nullReturn]
    if flagRange is True:
        flagged = any(flags)
    else:
        flagged = False
    if charAge != nullReturn:
        if charAge in ageCounter:
            return sum(ageCounter), flagged
        elif charAge not in ageCounter and len(ageCounter) > 0:
            return sum([charAge, sum(ageCounter)]), flagged
        else:
            return charAge, flagged
    elif charAge == nullReturn and len(ageCounter) > 0:
        return sum(ageCounter), flagged
    else:
        return nullReturn, flagged


def pmidAgeExtract(sample_dict, sectionID = 'methods', convertTo = 'week',
    checkConverts = ['day', 'week', 'month', 'year'], nullReturn = 'n/a',
    flagRange = True):
    """ Attempt age extraction from PMID ID link """

    if len(sample_dict['pmid']) == 0:
        print('No PMIDs found from GEO GSE metadata. Returning')
        return nullReturn

    pubContent = requests.get(sample_dict['pmid']).content

    sectionText = pmidSectionExtraction(pmidContent = pubContent, sectionID = sectionID)

    age = numericTimeConvert(text = sectionText, convertTo = convertTo,
                            nullReturn = nullReturn)

    #if flagRange is False:
    #    flagged = False

    return age


def pmidSectionExtraction(pmidContent, sectionID, possibleSections = ['abstract',
    'introduction', 'figures', 'materials and methods', 'methods',
    'experimental procedures', 'results', 'discussion', 'method summary',
    'supplementary material', 'acknowledgements', 'references', 'conclusions',
    'supporting information', 'funding'], abstractTags = ['abstract'],
    introTags = ['intro', 'introduction'], methodTags = ['methods', 'procedures'],
    resultTags = ['results'], conclusionTags = ['discussion', 'conclusions'],
    soupAttempts = ['h2', 'h1', 'h6'], pmcPreference = True, nullReturn = 'n/a'):
    """ From the text of a pmid link, extract the subsection of the paper text.
        A set of possible section IDs are scanned for a given BeautifulSoup
        parsing mechanism (e.g. the sections may fall under 'h2' headers). If
        no sections are found or no sectionID-specific labels are found, then
        alternative parsing will be attempted (null returned otherwise). A
        preference may be given to PMC texts, if multiple links are available,
        as they have a more predictable format.
        Args:
            pmidText: Str - Text of a PubMed PMID abstact page, with links to
                possible full-text papers
            sectionID: Str - Desired sub-section of interest (e.g. 'methods')
            possibleSections: List - Possible section header names
            abstractTags: List - Set of possible abstract partial name matches
            introTags: List - Set of possible intro partial name matches
            methodTags: List - Set of possible methods partial name matches
            resultTags: List - Set of possible results partial name matches
            conclusionTags: List - Set of possible conclusion partial name matches
            soupAttempts: List - Header IDs for soup.find_all(). On empirical
                evidence, there may be 'div' tag and CSS class 'h6', instead of
                solely h1/h2 tags
            pmcPreference: Bool - If multiple paper texts exist, choose a PMC
                link (more predictable format)

        Returns:
            sectionText: Str - sectionID text
    """
    sectionID = sectionID.lower()
    allTags = abstractTags + introTags + methodTags + resultTags + conclusionTags
    if sectionID not in allTags:
        raise ValueError('{0} is not a valid section ID. Try one of '
                            '{1}'.format(sectionID, allTags))

    root = ET.fromstring(pmidContent)
    for elem in root.getiterator():
        if not hasattr(elem.tag, 'find'): continue
        i = elem.tag.find('}')
        if i >= 0:
            elem.tag = elem.tag[i+1:]

    if sectionID in abstractTags:
        tags = abstractTags
    elif sectionID in introTags:
        tags = introTags
    elif sectionID in methodTags:
        tags = methodTags
    elif sectionID in resultTags:
        tags = resultTags
    elif sectionID in conclusionTags:
        tags = resultTags

    links = [i.get('href') for i in root.findall('.//div[@class="supplemental col three_col last"]//a')]
    if len(links) == 0:
        print('No links to full text papers found on PubMed PMID site. Returning')
        return nullReturn

    if len(links) > 1:
        if pmcPreference is True:
            paperURL = [x for x in links if ('/pmc/articles/pmid/' in x.lower())]
            paperURL = [x for x in links if ('http' in x.lower())]
            if len(paperURL) > 1:
                raise Exception('Multiple PMC links found for a single paper?')
            elif len(paperURL) == 1:
                links = paperURL
            else:
                pass # links = links[0]

    winLink, winText, winDiv = maxSectionMatch(links = links,
                possibleSections = possibleSections, soupAttempts = soupAttempts,
                nullReturn = nullReturn)

    sectionText = findSectionText(fullText = winText, sectionTags = tags,
                soupDiv = winDiv, nullReturn = nullReturn)

    return sectionText



def findTextInSectionTexts(links, possibleSections, nullReturn):
    MAX_RETRIES = 20
    adapter = requests.adapters.HTTPAdapter(max_retries=MAX_RETRIES)

    for link in links:
        linkResults[link] = dict()

        session = requests.Session()
        session.mount('https://', adapter)
        session.mount('http://', adapter)

        linkGet = session.get(link).text.lower()
        linkResults[link]['text'] = linkGet

        linkResults[link]['divs'] = dict()



def findSectionText(fullText, sectionTags, soupDiv, nullReturn = 'n/a'):
    """ Within the text of a paper (from maxSectionMatch()), and using the Soup
        div argument (from maxSectionMatch()), look for the section headers
        defined by sectionTags. If the sections are found but cannot be located
        by an index-based re.findall(), an attempt is made to remove extra
        newlines, otherwise the section cannot be found and null returned """

    if fullText == nullReturn:
        return nullReturn

    updated = False
    methodSplit = []
    countSplit, splitIndex = 0, 0
    soup = BeautifulSoup(fullText, features = 'html.parser')
    if soupDiv in ['h1', 'h2']:
        for link in soup.find_all(soupDiv):
            methodSplit.append(str(link))
            countSplit += 1
            for tag in sectionTags:
                if tag in link.get_text().lower() and updated is False:
                    updated = True
                    splitIndex = countSplit

    elif soupDiv in ['h6']:
        for link in soup.find_all('div', class_ = soupDiv):
            methodSplit.append(str(link))
            countSplit += 1
            for tag in sectionTags:
                if tag in link.get_text().lower() and updated is False:
                    updated = True
                    splitIndex = countSplit

    if splitIndex == 0:
        print('No section headings tagged with {0} found, returning '
                'null'.format(sectionTags))

        return nullReturn

    findMe = '{0}(.*){1}'.format(methodSplit[splitIndex - 1],
                                    methodSplit[splitIndex])

    findText = re.findall(findMe, fullText)
    if len(findText) == 0:
        fullText = re.sub(r'\n', '', fullText)
        findText = re.findall(findMe, fullText)
        if len(findText) == 0:
            print('No text found with {0} section headers, returning '
                    'null'.format(sectionTags))

            return nullReturn

    if len(findText) > 1:
        print('Warning, multiple method sections matched. Returning first')

    return findText[0]

def maxSectionMatch(links, possibleSections = ['abstract', 'introduction',
    'figures', 'materials and methods', 'methods', 'experimental procedures',
    'results', 'discussion', 'method summary', 'supplementary material',
    'acknowledgements', 'references', 'conclusions', 'supporting information',
    'funding'], soupAttempts = ['h2', 'h1', 'h6'], nullReturn = 'n/a'):
    """ Over a list of links, scan through a set of html headers and choose that
        which has the maximal number of matches to possible section IDs. Then
        return the url and the get() text of the link which wins (max-max) and
        the Soup div. """

    MAX_RETRIES = 20
    adapter = requests.adapters.HTTPAdapter(max_retries=MAX_RETRIES)

    missings = [x for x in soupAttempts if x not in ['h1', 'h2', 'h6']]
    if len(missings) > 0:
        print('Warning. Do not know how to handle Soup parsers {0}. These will '
            'be ignored'.format(missings))

    soupAttempts = [x for x in soupAttempts if x not in missings]

    linkResults = dict()
    for link in links:
        linkResults[link] = dict()
        session = requests.Session()
        session.mount('https://', adapter)
        session.mount('http://', adapter)

        linkGet = session.get(link).text.lower()
        linkResults[link]['text'] = linkGet

        linkResults[link]['divs'] = dict()

        for attempt in soupAttempts:
            linkCount = 0
            soup = BeautifulSoup(linkGet, features = 'html.parser')
            ### Count the number of h1/h2 headings
            if attempt in ['h1', 'h2']:
                for sec in soup.find_all(attempt):
                    for x in possibleSections:
                        if x in sec:
                            linkCount += 1
            ### Count the number of div containers with h6 class
            elif attempt in ['h6']:
                for sec in soup.find_all('div', class_ = attempt):
                    for x in possibleSections:
                        if x in sec:
                            linkCount += 1
            ### linkResults = {link : {'divs' : {attempt : linkCount}}}
            linkResults[link]['divs'][attempt] = linkCount

    maxes = dict()
    for link in linkResults:
        maxes[link] = max(linkResults[link]['divs'].values())

    maxVal = max(maxes.values())
    if maxVal == 0:
        print(f'No section matches found for any link. Returning null for {links}')
        return nullReturn, nullReturn, nullReturn

    maxLink = [x for x in linkResults if maxes[x] == maxVal]
    if len(maxLink) > 1:
        print('Tie for maximum number of matches, will return first link')

    maxLink = maxLink[0]
    maxDivVal = max(linkResults[maxLink]['divs'].values())
    maxDiv = [x for x in linkResults[maxLink]['divs'] if linkResults[maxLink]['divs'][x] == maxDivVal]
    if len(maxDiv) > 1:
        print('Tie for maximum div heading, will return first')


    return maxLink[0], linkResults[maxLink]['text'], maxDiv[0]


def geoSampleCellCheck(sample_dict,
                    parseLocations = ['treatment_protocol', 'growth_protocol'],
                    checkLocations = ['sample_cell_line', 'sample_cell_type'],
                    cellDetectKWs = ['DMEM', 'FBS', 'bovine serum', 'passage']):

    """
    Function used to create the 'cells' field in the output dictinary.
    Function will check locations in `parseLocations` using a keyword match. If it finds any keywords listed in `cellDetectKWs` argument then the functin returns True.
    The function will also return True if a field in `checkLocations` is not blank. These fields are characteristic fields that are left blank in non-cell experiments.

        Args:
            `sample_dict` - Dict: Contains {sampleID : data} key-value pairs for each element in the selected sample.
            `parseLocations` - List: Locations to use keywords to parse whether there are cells or not
            `checkLocations` - List: Location to check for any content at all. No cells means no content in these places.
            `cellDetectKWs` - List:

        Returns:
            True/False depending on the input data.
    """
    for parseLocation in parseLocations:
        for cellDetectKW in cellDetectKWs:
            if cellDetectKW in sample_dict[parseLocation]:
                return True
    for checkLocation in checkLocations:
        try:
            if sample_dict[checkLocation] != '':
                return True
        except:
            continue

    return False
