from lxml import etree, objectify
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import tarfile, os, urllib, gzip, re

def tar_gz_extracter(url):
    """
    This function will take a FTP URL and download then extract the tar-gzipped file (.tgz).
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
    base_name = os.path.basename(url)

    file_name, file_extension = os.path.splitext(base_name)
    tar = tarfile.open(f_url)

    try:
        xml_name = [i for i in tar.getnames() if re.match('.*.xml', str(i))][0]
    except:
        xml_name = ''
        print(tar.getnames())
        return None
    tar.extract(xml_name, path = f'output/xml/')

    tar.close()

    return xml_name


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
        if sample_dict[checkLocation] != '':
            return True

    return False

def xml_parser(url, parse_platforms = False, DEBUG = 0, multichannel = False, keep_files = [None]):
    """
    This function will take in the URL for a Series FTP .tgz file, download the file, extract the .xml file from the compressed directory and then parse the .xml.

    Try: except lxml.etree.XMLSyntaxError as e: section is to prevent illegal XML characters from breaking the script. Ideally, we'd replace all illegal characters before parsing but I couldn't find good code to do that.

    First the ./Series tags are looped through and the corresponding pubmed ID is taken from each series.

    Then the ./Samples tags are looped through and all relevant information is captured for each sample.

    If multichannel == False and the parser finds more than 1 channel-count tag within a sample, it'll append an essentially empty data dict knowing it'll be filtered later.
        Args:
            `url` - Str: Series FTP URL leading to a .tgz file
            `parse_platforms` - Bool: If set to its default value, False, a platform URL will be rejected instead of parsed poorly.
            `DEBUG` - Int: Sets depth of debug messages. Mostly deprecated.
            `multichannel` - Bool: If set to False, multichannel files will be filtered out of the final result.
            `keep_files` - List: File types to keep during a full run. If 'xml' is in `keep_files` then the xml file won't be deleted.

        Returns:
            `series_pmid_dict` - Dict: Mappings from Series Accession numbers to PMIDs
            `samples_dict` - Dict: {sampleID : data} key-value pairs for all samples in input list.
                - `data` is a dictionary of the following information:
                    {'sample_id' : sample_id,
                        'sample_source' : sample_source,
                        'sample_title' : sample_title,
                        'molecule' : molecule,
                        'organism' : organism,
                        'treatment_protocol' : treatment_protocol,
                        'extract_protocol' : extract_protocol,
                        'growth_protocol' : growth_protocol,
                        'description' : sample_description,
                        'sample_cell_type' : sample_cell_type,
                        'sample_type' : sample_type,
                        'sample_sex' : sample_sex,
                        'sample_tissue' : sample_tissue,
                        'sample_age' : sample_age,
                        'sample_indication' : sample_indication,
                        'sample_genotype' : sample_genotype,
                        'sample_cell_line' : sample_cell_line,
                        'expression' : expression,
                        'cells' : cells
                        }
    """


    n=1000
    if DEBUG >= 2:
        print(f'Parsing: {url}')
    if parse_platforms == False and re.search('.*GPL.*', url):
        print('URL is to a platform and platform flag is off.')
        column_list = ['sample_id', 'source', 'title', 'molecule', 'treatment_protocol', 'extract_protocol', 'growth_protocol', 'characteristic_list']
        return dict()
    ### Files are gzipped so open the url then extract only the XML file from the tar file
    xml_filename = tar_gz_extracter(url)
    ### If you can't extract the XML, return a blank dictionary
    if xml_filename is None:
        return {}, {}

    parser = etree.XMLParser(recover=True)
    try:
        tree_ab = ET.parse(f'output/xml/{xml_filename}', parser=parser)
    except etree.XMLSyntaxError as e:
        print(e)
        return {}, {}
    root = tree_ab.getroot()
    if 'xml' not in keep_files:
        os.unlink('output/xml/' + xml_filename)

    #### I stole this code from the internet but it will replace the prefix that gets added
    #### There's probably a better way to do this
    for elem in root.getiterator():
        if not hasattr(elem.tag, 'find'): continue
        i = elem.tag.find('}')
        if i >= 0:
            elem.tag = elem.tag[i+1:]
    ####


    for series in root.findall('./Series'):
        pmids = [i.text for i in series.findall('./Pubmed-ID')]
        series_acc_numbers = [i.text for i in series.findall('./Accession')]
        series_pmid_dict = {key:value for (key,value) in list(zip(series_acc_numbers, pmids))}

    sample_row_list = []
    samples_dict = {}
    series_pubmed_id_dict = {}
    for samples in root.findall('./Sample'):
        sample_id = samples.attrib['iid']
        channel_num = np.max([int(i.text) for i in samples.findall('./Channel-Count')])
        if multichannel == False:
            if int(channel_num) > 1:
                #print('Multichannel set to false, returning blank meta-data for multichannel sample')
                sample_dict = {'sample_id' : sample_id,
                                'sample_source' : '',
                                'sample_title' : '',
                                'molecule' : '',
                                'organism' : '',
                                'treatment_protocol' : '',
                                'extract_protocol' : '',
                                'growth_protocol' : '',
                                'description' : 'multi',
                                'sample_cell_type' : '',
                                'sample_type' : '',
                                'sample_sex' : '',
                                'sample_tissue' : '',
                                'sample_age' : '',
                                'sample_indication' : '',
                                'sample_genotype' : '',
                                'sample_cell_line' : '',
                                'expression' : '',
                                'cells' : ''
                                }
                samples_dict[sample_id] = sample_dict

                sample_row_list.append(sample_dict)
                continue
        sample_source = ''
        organism = ''
        sample_title = ''
        molecule = ''
        sample_indication = ''
        treatment_protocol = ''
        sample_description = ''
        extract_protocol = ''
        growth_protocol = ''
        sample_type = ''
        sample_cell_line = ''
        sample_cell_type = ''
        sample_sex = ''
        sample_tissue = ''
        sample_age = ''
        sample_genotype = ''

        ### Loop through elements of the sample
        for sample_element in samples:
            ### Potentially more than one channel?
            ### Adding this for later but not used yet
            channel_count = 1
            if sample_element.tag == 'Title':
                sample_title = sample_element.text

            elif sample_element.tag == 'Accession':
                accession_number = sample_element.text
                accession_db = sample_element.attrib['database']

            elif sample_element.tag == 'Channel-Count':
                channel_count = sample_element.text

            elif sample_element.tag == 'Description':
                sample_description = sample_element.text.strip()
            elif sample_element.tag == 'Channel':
                channel_num = sample_element.attrib['position']
                for channel_value in sample_element:
                    if channel_value.tag == 'Source':
                        sample_source = channel_value.text

                    if channel_value.tag == 'Treatment-Protocol':
                        treatment_protocol = channel_value.text.strip()

                    if channel_value.tag == 'Molecule':
                        molecule = channel_value.text.strip()

                    if channel_value.tag == 'Organism':
                        organism = channel_value.text.strip()

                    if channel_value.tag == 'Extract-Protocol':
                        extract_protocol = channel_value.text.strip()

                    if channel_value.tag == 'Growth-Protocol':
                        growth_protocol = channel_value.text.strip()

                for char_value in sample_element.findall('./Characteristics'):
                    try:
                        if lower(char_value.attrib['tag']) == 'cell type':
                            sample_cell_type = char_value.text.strip()
                        if lower(char_value.attrib['tag']) == 'cell line':
                            sample_cell_line = char_value.text.strip()
                        if lower(char_value.attrib['tag']) == 'sample type':
                            sample_type = char_value.text.strip()
                        if lower(char_value.attrib['tag']) == 'sex':
                            sample_sex = char_value.text.strip()
                        if lower(char_value.attrib['tag']) == 'tissue':
                            sample_tissue = char_value.text.strip()
                        if lower(char_value.attrib['tag']) == 'age':
                            sample_age = char_value.text.strip()
                        if lower(char_value.attrib['tag']) == 'indication':
                            sample_indication = char_value.text.strip()
                        if lower(char_value.attrib['tag']) == 'genotype':
                            sample_genotype = char_value.text.strip()
                    except:
                        pass

        if 'RNA' in sample_source or 'RNA' in molecule:
            expression = True
        else:
            expression = False



        sample_dict = {'sample_id' : sample_id,
                        'sample_source' : sample_source,
                        'sample_title' : sample_title,
                        'molecule' : molecule,
                        'organism' : organism,
                        'treatment_protocol' : treatment_protocol,
                        'extract_protocol' : extract_protocol,
                        'growth_protocol' : growth_protocol,
                        'description' : sample_description,
                        'sample_cell_type' : sample_cell_type,
                        'sample_type' : sample_type,
                        'sample_sex' : sample_sex,
                        'sample_tissue' : sample_tissue,
                        'sample_age' : sample_age,
                        'sample_indication' : sample_indication,
                        'sample_genotype' : sample_genotype,
                        'sample_cell_line' : sample_cell_line,
                        'expression' : expression,
                        'cells' : ''
                        }
        sample_dict['cells'] = geoSampleCellCheck(sample_dict)
        samples_dict[sample_id] = sample_dict
        sample_row_list.append(sample_dict)

    return series_pmid_dict, samples_dict
