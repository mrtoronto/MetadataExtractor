import re, requests, os, random, sys 
import numpy as np
from bs4 import BeautifulSoup
defaultAgent = {'User-Agent': 'SomeAgent 11.0'}

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


def extractGEOSampleInfo(sampleID, organism = 'Mus musculus', extracts = ['ID', 
    'Study', 'Organism', 'Sample type', 'Extracted molecule', 'Age', 'Gender', 
    'Expression', 'Cells'], keepCells = True, cellDetectChar = 'cell lines?\:', 
    cellDetectProt = ['DMEM', 'FBS', 'bovine serum', 'passage'],  
    parseCellIDs = ['Treatment', 'Growth'], keepMultiChannel = False, 
    parseAgeIDs = ['Characteristics', 'Description', 'Treatment protocol', 
    'Growth protocol'], convertAgeTo = 'week', checkAgeConverts = ['day', 'week',
    'month', 'year'], nullReturn = 'n/a', tryAgeStudy = True,
    parseStudyIDs = ['Summary', 'Overall design']):
    """ Extract metadata from a GEO GSM ID. Options to keep detected
        in vitro samples, and to exclude multichannel expression assays (e.g. 
        microarrays). Not all structured entries need to be parsed, however
        characteristics and protocols contain most salient information. Each
        major item, e.g. "Age", may require a unique set of functions. Cells may 
        have some utility, however multi-channel microarray reports are of very 
        limited use, thus its suggested to remove these samples. Main args will be
        checked here, however age and cell-specific parameters will be evaluated
        within those functions. This function should house the major set of data
        extraction methods, starting with age and gender, but eventual expansion
        to study design, gene knockouts, doses, drugs, etc.. Note these regular
        expressions are critically dependent on the format of the GEO's sample
        page texts.

        Args:
            sampleID - Str: GSM ID
            organism - Str: Organism name. Organism is especially critical for
                age/gender extraction
            keepCells - Bool: Keep in vitro samples
            cellDetectChar - Str: see geoSampleCellCheck()
            cellDetectProt - List: see geoSampleCellCheck()
            keepMultiChannel - Bool: Keep microarray multichannel samples
            extracts - List: Items to parse/infer
            parseAgeIDs - List: see sampleAgeExtract()
            parseCellIDs - List: see geoSampleCellCheck()
            checkAgeConverts - List: see sampleAgeExtract()
            nullReturn - Str: Null return value
            tryAgeStudy: Bool - If no age is found, attempt to find age in the 
                GSE study page
            parseStudyIDs: List geoAgeExtract()

        Returns:
            meta - Dict: Items in extracts for sampleID
    """
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={0}'.format(sampleID)
    urlGetText = requests.get(url).text

    if 'Could not find a public or private accession' in urlGetText:
        raise ValueError('Is {0} a valid GEO sample?'.format(sampleID))
    
    IDs = [x for x in re.findall('\<td nowrap\>(.*)\<\/td\>', urlGetText)]
    misMatches = [x for x in extracts if x not in IDs + ['ID', 'Study', 'Age', 
        'Gender', 'Expression', 'Cells']]
    if len(misMatches) > 0:
        print('Warning, {0} IDs not found. Will be returned as null'.format(misMatches))

    cleanText = re.sub(r'<.*?>|\\n', ' ', urlGetText) 

    if keepMultiChannel == False:
        if 'Channel 1' in cleanText and 'Channel 2' in cleanText:
            raise AttributeError('Multi-channel not allowed!')
    
    if keepCells == False:
        if geoSampleCellCheck(urlText = cleanText, cellDetectChar = cellDetectChar,
            cellDetectProt = cellDetectProt, protocolEntries = parseCellIDs):
            raise AttributeError('Cells not allowed!')

    GEOOrganism = re.findall(r'Organism \n  (.*)  \n', cleanText)[0]
    if GEOOrganism != organism:
        raise ValueError('Organism mismatch for sample {0}'.format(sampleID))

    sampType = re.findall(r'Sample type \n (.*) \n', cleanText)[0]
    moleExtract = re.findall(r'Extracted molecule \n (.*) \n', cleanText)[0]

    meta = dict()
    for extract in extracts:
        if extract == 'ID':
            meta[extract] = sampleID
        elif extract == 'Study':
            meta[extract] = re.findall(r'Series \(\d+\) \n     (.*)  \n', 
                cleanText)[0]
        elif extract == 'Sample type':
            meta[extract] = sampType
        elif extract == 'Extracted molecule':
            meta[extract] = moleExtract
        elif extract == 'Organism':
            meta[extract] = GEOOrganism
        elif extract == 'Age':
            meta[extract] = geoAgeExtract(urlText = urlGetText, 
                    parseAgeIDs = parseAgeIDs, convertTo = convertAgeTo, 
                    nullReturn = nullReturn, checkConverts = checkAgeConverts,
                    tryAgeStudy = tryAgeStudy, parseStudyIDs = parseStudyIDs)
        elif extract == 'Gender':
            meta[extract] = geoGenderExtract(urlGetText)
        elif extract == 'Expression':
            if ('RNA' in sampType or 'RNA' in moleExtract):
                meta[extract] = True
            else:
                meta[extract] = False
        elif extract == 'Cells':
            meta[extract] = geoSampleCellCheck(urlText = cleanText, 
                cellDetectChar = cellDetectChar, cellDetectProt = cellDetectProt, 
                protocolEntries = parseCellIDs)
        else:
            meta[extract] = nullReturn
    
    return meta


def geoSampleCellCheck(urlText, cellDetectChar = 'cell lines?\:', 
    cellDetectProt = ['DMEM', 'FBS', 'bovine serum', 'passage'], 
    protocolEntries = ['Treatment', 'Growth'], nullReturn = False):
    """ Boolean check for presence of cells in a GEO sample (GSM) in either the 
        'Characteristics' or protocols entries. nullReturn is in the event that
        cells are not detected in characteristics, nor are the protocol entries
        found in the sample text (missing info). Note that these terms are 
        designed for in vitro studies, and do not necessarily exclude cell-sorting.
        Sorted cells can come from an animal of a certain age, however cultured
        cells does not likely correspond to age in a meaningful way.

        Args:
            urlText - Str: requests.get(sampleID).text
            cellDetectChar - Str: Regex for cell detection under "Characteristics"
            cellDetectProt - List: Common keywords to indicate in vitro assays 
            protocolEntries - List: Protocol entries to check for cellDetectProt
                items
        Return:
            bool
    """
    cleanText = re.sub(r'<.*?>|\\n', ' ', urlText)
    if re.findall(r'Characteristics \n(.*)\n', cleanText):
        if re.findall(cellDetectChar, re.findall(r'Characteristics \n(.*)\n', 
            cleanText)[0]):
            return True 

    for prot in protocolEntries:
        if re.findall('{0} protocol \n(.*)\n'.format(prot), cleanText):
            if any([x for x in cellDetectProt if x in re.findall('{0} protocol \n(.*)\n'.format(prot), 
                cleanText)[0]]):
                return True

    return False


def geoAgeExtract(urlText, checkCell = True, parseAgeIDs = ['Characteristics',
     'Description', 'Treatment protocol', 'Growth protocol'], convertTo = 'week',
    nullReturn = 'n/a', checkConverts = ['day', 'week', 'month', 'year'],
    tryAgeStudy = True, parseStudyIDs = ['Summary', 'Overall design']):
    """ Extract age from GEO sample (GSM) text. Either/or the 'Characteristics'
        or protocols/descriptions entries. See numericTimeConvert() docstring for
        more detail on the approach. Return null on cell detect. If the same age
        is detected in both 'Characteristics' and descriptions/protocols, then 
        assume that they are duplicated, and only add non-matching ages in the
        descriptions/protocols. NumberDict is an english word to integer 
        conversion dictionary pre-defined in globals().
            
        Args:   
            urlText: Str - requests.get(sampleID).text
            parseAgeIDs: List - Entries within url text
            convertTo: Str - Converted time units ['day', 'week', 'month', 'year']
            nullReturn: Str - Return if no numbers found
            checkConverts: List - Time units to convert to convertTo unit
            tryAgeStudy: Bool - If no age is found, attempt to find age in the 
                GSE study page
            parseStudyIDs: List - Entries within url text if GSE study is checked

        Return:
            age: Float - Estimated age
    """
    if geoSampleCellCheck(urlText) is True and checkCell == True:
        return nullReturn

    ageCounter, charAge = [], nullReturn
    
    cleanText = re.sub(r'<.*?>|\\n', ' ', urlText)
    for i in parseAgeIDs:
        if re.findall('({0}) \n (.*)'.format(i), cleanText):
            matches = list(re.findall('({0}) \n (.*)'.format(i), cleanText))[0]
            if len(matches) > 1:
                if i == 'Characteristics' and 'age:' in matches[1]:
                    ageMatch = re.findall('\<br\>age\:(.*)\<br\>', urlText)[0]
                    charAge = numericTimeConvert(text = ageMatch, 
                        convertTo = convertTo, nullReturn = nullReturn)

                elif i != 'Characteristics':
                    ageCounter.append(numericTimeConvert(text = matches[1], 
                        convertTo = convertTo, nullReturn = nullReturn))

    ageCounter = [x for x in ageCounter if x != nullReturn]

    if charAge != nullReturn:
        if charAge in ageCounter:
            return sum(ageCounter)
        elif charAge not in ageCounter and len(ageCounter) > 0:
            return sum([charAge, sum(ageCounter)])
        else:
            return charAge
    elif charAge == nullReturn and len(ageCounter) > 0:
        return sum(ageCounter)
    else:
        if tryAgeStudy is True:
            gseAttempt = gseAgeExtract(urlText = urlText, convertTo = convertTo, 
                checkConverts = checkConverts, parseIDs = parseStudyIDs,
                nullReturn = nullReturn)
            if gseAttempt == nullReturn:
                return nullReturn #New function
            else:
                return gseAttempt
        else:
            return nullReturn


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
    text = re.sub('(\d+)\-(\D)', '\\1 \\2', text) #keep '7-8', convert 'seven-week'
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
        nums = addAgeStrings(nums, strToNumConvert = strToNumConvert, 
            convertFrom = convertFrom, convertTo = convertTo)

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

    timeMatches = re.compile('(?<!for )(\d+\-?\d+?){0}|(?<!for )(\d+\-?\d+?){1}|(?<!for )(\d+){0}|(?<!for )(\d+){1}|(?<!for )(\d+\.?\d+?){0}|(?<!for )(\d+\.?\d+?){1}'.format(re1, re2))
    timeMatchesDur = re.compile('for (\d+\-?\d+?){0}|for (\d+\-?\d+?){1}|for (\d+){0}|for (\d+){1}|for (\d+\.?\d+?){0}|for (\d+\.?\d+?){1}'.format(re1, re2))

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


def gseAgeExtract(urlText, convertTo = 'week', 
    checkConverts = ['day', 'week', 'month', 'year'], 
    parseIDs = ['Summary', 'Overall design'], nullReturn = 'n/a'):
    """ Akin to extracting age from individual sample descriptions, extract age
        from GSE study information.

        Return:
            age: Float - Estimated age

    TODO: Different ages may exist for different experimental groups! Check 
        GSE7191 as an example. Currently this function (incorrectly) combines
        and adds the ages of each group
    
    """
    ageCounter, charAge = [], nullReturn

    newText = sampToGSEText(urlText)
    cleanText = re.sub(r'<.*?>|\\n', ' ', newText)
    for i in parseIDs:
        if re.findall('({0}) \n (.*)'.format(i), cleanText):
            matches = list(re.findall('({0}) \n (.*)'.format(i), cleanText))[0]
            if len(matches) > 1:
                ageCounter.append(numericTimeConvert(text = matches[1], 
                        convertTo = convertTo, nullReturn = nullReturn))

    ageCounter = [x for x in ageCounter if x != nullReturn]

    if charAge != nullReturn:
        if charAge in ageCounter:
            return sum(ageCounter)    
        elif charAge not in ageCounter and len(ageCounter) > 0:
            return sum([charAge, sum(ageCounter)])
        else:
            return charAge
    elif charAge == nullReturn and len(ageCounter) > 0:
        return sum(ageCounter)
    else:
        return nullReturn


def sampToGSEText(urlText):
    """ Extract the GEO study ID from a GSM/sample request get """

    sampleID = set(re.findall(r'acc\=(GSM\d+)\"', urlText))
    if len(sampleID) > 1:
        raise Exception('Multiple sample IDs found under GEO GSM page')
    elif len(sampleID) == 0:
        raise Exception('Sample ID extraction from GEO GSM page failed')
    else:
        if list(sampleID)[0] == ['']:
            raise Exception('Sample ID extraction from GEO GSM page failed')
        sampleID = list(sampleID)[0]

    cleanText = re.sub(r'<.*?>|\\n', ' ', urlText) 
    GSEExtract = re.findall(r'Series \(\d+\) \n     (.*)  \n', cleanText)

    if len(GSEExtract) == 0:
        raise Exception('No GSE study IDs found from GEO GSM page')

    elif len(GSEExtract) > 1:
        print('Warning, multiple GSEs found for a samples url text. An attempt '
            'will be made to select the study with the sample in the sample list')
        
        matched = False 
        for GSE in GSEExtract:
            url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={0}'.format(GSE)
            newText = requests.get(url).text
            if sampleID in newText:
                matched = True 
                break
        
        if matched is False:
            raise Exception('Original sample ID not found in any of the potential '
                'GSE studies')
        
    elif len(GSEExtract) == 1: 
        url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={0}'.format(GSEExtract[0])
        newText = requests.get(url).text

    return newText 


def geoGenderExtract(urlText, protocolEntries = ['Treatment', 'Growth'], 
    nullReturn = 'n/a'):
    """ Male/female check for a GEO sample (GSM) in either the 'Characteristics' 
        or protocols entries. May need adjustment for non-mammals (e.g. 
        hermaphrodites in C elegans)

        Args:
            urlText - Str: requests.get(sampleID).text
            cellDetectChar - Str: Regex for cell detection under "Characteristics"
            cellDetectProt - List: Common keywords to indicate in vitro assays 

        Return:
            'Male'/'Female'/nullReturn
    """
    cleanText = re.sub(r'<.*?>|\\n', ' ', urlText)
    charCheck = re.findall(r'Characteristics \n(.*)\n', cleanText)
    if len(charCheck) > 0:
        if re.findall('female[ s]', charCheck[0]):
            return 'Female' 
        if re.findall('(?<!fe)male[ s]', charCheck[0]):
            return 'Male' 

    for prot in protocolEntries:
        protCheck = re.findall('{0} protocol \n(.*)\n'.format(prot), cleanText)
        if len(protCheck) > 0:
            if re.findall('female[ s]', protCheck[0]):
                return 'Female' 
            if re.findall('(?<!fe)male[ s]', protCheck[0]):
                return 'Male' 

    return nullReturn



def pmidAgeExtract(urlText, sectionID = 'methods', convertTo = 'week', 
    checkConverts = ['day', 'week', 'month', 'year'], nullReturn = 'n/a'):
    """ Attempt age extraction from PMID ID link """

    newText = sampToGSEText(urlText)
    pmidExtract = re.findall(r'\/pubmed\/(\d+)', newText)

    if len(pmidExtract) == 0 or 'Citation missing' in newText:
        print('No PMIDs found from GEO GSE study page. Returning')
        return nullReturn

    elif len(pmidExtract) > 1:
        print('Warning, multiple PMIDs found from a GEO GSE study page. An attempt '
            'to extract age on the first publication will be made')
    
    pmid = pmidExtract[0]
    url = 'https://www.ncbi.nlm.nih.gov/pubmed/{0}'.format(pmid)
    pubText = requests.get(url).text

    sectionText = pmidSectionExtraction(pmidText = pubText, sectionID = sectionID)

    age = numericTimeConvert(text = sectionText, convertTo = convertTo, 
                            nullReturn = nullReturn)
                        
    return age


def pmidSectionExtraction(pmidText, sectionID, possibleSections = ['abstract', 
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
    
    linkExtract = re.findall(r'(Full text links.*)\n', pmidText)
    if len(linkExtract) == 0:
        print('No links to full text papers found on PubMed PMID site. Returning')

        return nullReturn 

    linkExtract = linkExtract[0]
    
    soup = BeautifulSoup(linkExtract, features = 'html.parser')
    links = []
    for link in soup.find_all('a', href = True):
        links.append(link['href'])

    if len(links) > 1:
        if pmcPreference is True:
            paperURL = [x for x in links if '/pmc/' in x.lower()]
            if len(paperURL) > 1:
                raise Exception('Multiple PMC links found for a single paper?')

            elif len(paperURL) == 1:
                links = paperURL

    winLink, winText, winDiv = maxSectionMatch(links = links, 
                possibleSections = possibleSections, soupAttempts = soupAttempts,
                nullReturn = nullReturn)

    sectionText = findSectionText(fullText = winText, sectionTags = tags, 
                soupDiv = winDiv, nullReturn = nullReturn)

    return sectionText


def maxSectionMatch(links, possibleSections = ['abstract', 'introduction', 
    'figures', 'materials and methods', 'methods', 'experimental procedures', 
    'results', 'discussion', 'method summary', 'supplementary material', 
    'acknowledgements', 'references', 'conclusions', 'supporting information', 
    'funding'], soupAttempts = ['h2', 'h1', 'h6'], nullReturn = 'n/a'):
    """ Over a list of links, scan through a set of bs4 parsers and choose that
        which has the maximal number of matches to possible section IDs. Then 
        return the url and the get() text of the link which wins (max-max) and
        the Soup div. """

    missings = [x for x in soupAttempts if x not in ['h1', 'h2', 'h6']]
    if len(missings) > 0:
        print('Warning. Do not know how to handle Soup parsers {0}. These will '
            'be ignored'.format(missings))
    
    soupAttempts = [x for x in soupAttempts if x not in missings]

    linkResults = dict() 
    for link in links:
        linkResults[link] = dict()
        linkGet = requests.get(link, headers = defaultAgent).text.lower()
        linkResults[link]['text'] = linkGet 
        
        linkResults[link]['divs'] = dict()
        for attempt in soupAttempts:
            linkCount = 0
            soup = BeautifulSoup(linkGet, features = 'html.parser')
            if attempt in ['h1', 'h2']:
                for sec in soup.find_all(attempt):
                    for x in possibleSections:
                        if x in sec:
                            linkCount += 1
                        
            elif attempt in ['h6']:
                for sec in soup.find_all('div', class_ = attempt):
                    for x in possibleSections:
                        if x in sec:
                            linkCount += 1

            linkResults[link]['divs'][attempt] = linkCount 

    maxes = dict()
    for link in linkResults:
        maxes[link] = max(linkResults[link]['divs'].values())

    maxVal = max(maxes.values())
    if maxVal == 0:
        print('No section matches found for any link. Returning null')
        return nullReturn, nullReturn, nullReturn
    
    maxLink = [x for x in linkResults if maxes[x] == maxVal]
    if len(maxLink) > 1:
        print('Tie for maximum number of matches, will return first link')

    maxLink = maxLink[0]
    maxDivVal = max(linkResults[maxLink]['divs'].values())
    maxDiv = [x for x in linkResults[maxLink]['divs'] if linkResults[maxLink]['divs'][x] == maxDivVal]
    if len(maxDiv) > 1:
        print('Tie for maximum div heading, will return first')


    return maxLink, linkResults[maxLink]['text'], maxDiv[0]


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