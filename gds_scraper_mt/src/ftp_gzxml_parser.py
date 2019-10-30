import tarfile
import os
import urllib.request
import urllib
import gzip
from lxml import etree, objectify
import xml.etree.ElementTree as ET
import re
import pandas as pd

### Digging through some of the layers of the XML
### Look for <Sample *> section

def tar_gz_extracter(url):

    try:
        f_url = urllib.request.urlretrieve(url, filename=None)[0]
    except urllib.error.ContentTooShortError as e:
        ### This happens super rarely but sometimes the URL is too short..?
        ### Would be good to formally log these URLs and see what's up
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

    #f_url.close()
    tar.close()

    return xml_name


def geoSampleCellCheck(sample_dict,
                    parseLocations = ['treatment_protocol', 'growth_protocol'],
                    checkLocations = ['sample_cell_line', 'sample_cell_type'],
                    cellDetectKWs = ['DMEM', 'FBS', 'bovine serum', 'passage']):
    for parseLocation in parseLocations:
        for cellDetectKW in cellDetectKWs:
            if cellDetectKW in sample_dict[parseLocation]:
                return True
    for checkLocation in checkLocations:
        if sample_dict[checkLocation] != '':
            return True

    return False

def xml_parser(url, parse_platforms = False, DEBUG = 0, multichannel = False, keep_files = [None]):
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
        return [], {}

    parser = etree.XMLParser(recover=True)
    tree_ab = ET.parse(f'output/xml/{xml_filename}', parser=parser)
    root = tree_ab.getroot()
    if 'xml' not in keep_files:
        os.unlink('output/xml/' + xml_filename)

    #### I stole this code from the internet but it will replace the prefix that gets added
    #### There's probably a better way but doing this let me get to the data"""
    for elem in root.getiterator():
        if not hasattr(elem.tag, 'find'): continue  # (1)
        i = elem.tag.find('}')
        if i >= 0:
            elem.tag = elem.tag[i+1:]
    ####


    for series in root.findall('./Series'):
        pmids = [i.text for i in series.findall('./Pubmed-ID')]
        series_acc_numbers = [i.text for i in series.findall('./Accession')]
        series_pmid_dict = {key:value for (key,value) in list(zip(series_acc_numbers, pmids))}

    sample_row_list = []
    series_pubmed_id_dict = {}
    for samples in root.findall('./Sample'):
        if multichannel == False:
            channel_count = samples.find('./Channel-Count').text
            if int(channel_count) > 1:
                #print('Multichannel set to false, returning blank meta-data for multichannel sample')
                sample_row_list.append(
                    {'sample_id' : samples.attrib['iid'],
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
                        'sample_cell_line' : '',
                        'sample_genotype' : '',
                        'expression' : '',
                        'cells' : ''})
                continue
        sample_id = ''
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

        sample_id = samples.attrib['iid']
        if DEBUG >= 3:
            print(sample_id)
        ### Loop through elements of the sample
        for sample_element in samples:
            ### Potentially more than one channel?
            ### Adding this for later but not used yet
            channel_count = 1
            if sample_element.tag == 'Title':
                sample_title = sample_element.text
                if DEBUG >= 3:
                    print('Title:', sample_title)
            elif sample_element.tag == 'Accession':
                accession_number = sample_element.text
                accession_db = sample_element.attrib['database']
                if DEBUG >= 3:
                    print(f'Accession number (db): {accession_number} ({accession_db})')
            elif sample_element.tag == 'Channel-Count':
                channel_count = sample_element.text
                if DEBUG >= 3:
                    print(f'Channel count: {channel_count}')
            ### Channel section has a bunch of cool stuff
            elif sample_element.tag == 'Description':
                sample_description = sample_element.text.strip()
            elif sample_element.tag == 'Channel':
                channel_num = sample_element.attrib['position']
                for channel_value in sample_element:
                    if channel_value.tag == 'Source':
                        sample_source = channel_value.text
                        if DEBUG >= 3:
                            print(f'Sample source: {sample_source}')
                    ### This is cool, we could grab some metadata from keywords maybe?
                    if channel_value.tag == 'Treatment-Protocol':
                        treatment_protocol = channel_value.text.strip()
                        if DEBUG >= 3:
                            print(f'Treatment protocol: {treatment_protocol[:50]}')
                    ### I saw you had this on github I think
                    if channel_value.tag == 'Molecule':
                        molecule = channel_value.text.strip()
                        if DEBUG >= 3:
                            print(f'Molecule: {molecule}')
                    if channel_value.tag == 'Organism':
                        organism = channel_value.text.strip()
                        if DEBUG >= 3:
                            print(f'organism: {organism}')
                    ### Not sure if its important
                    if channel_value.tag == 'Extract-Protocol':
                        extract_protocol = channel_value.text.strip()
                        if DEBUG >= 3:
                            print(f'Extract-Protocol: {extract_protocol[:50]}')

                    if channel_value.tag == 'Growth-Protocol':
                        growth_protocol = channel_value.text.strip()
                        if DEBUG >= 3:
                            print(f'Growth-Protocol: {growth_protocol[:50]}')

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

            else:
                if DEBUG >= 3:
                    print('Tag:', sample_element.tag)
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

        sample_row_list.append(sample_dict)
    if int(channel_count) > 1:
        print('Multichannel set to false, returning blank meta-data for multichannel sample')
    return sample_row_list, series_pmid_dict

"""if __name__ == 'main':
    url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133304/miniml/GSE133304_family.xml.tgz'

    column_list = ['sample_id', 'sample_source', 'sample_title', 'molecule', 'treatment_protocol', 'extract_protocol', 'growth_protocol', 'characteristic_list']

    fetch_df = pd.DataFrame(xml_parser(url), columns = column_list)

    print(fetch_df.head())"""
