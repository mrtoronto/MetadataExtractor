import tarfile
import os
import urllib.request
import gzip
from lxml import etree, objectify
import xml.etree.ElementTree as ET
import re
import pandas as pd

### Digging through some of the layers of the XML
### Look for <Sample *> section

def tar_gz_extracter(url):

    f_url = urllib.request.urlretrieve(url, filename=None)[0]
    base_name = os.path.basename(url)

    file_name, file_extension = os.path.splitext(base_name)
    tar = tarfile.open(f_url)

    try:
        xml_name = [i for i in tar.getnames() if re.match('.*.xml', str(i))][0]
    except:
        xml_name = ''
        print(tar.getnames())
        return None
    tar.extract(xml_name)

    #f_url.close()
    tar.close()

    return xml_name

def xml_parser(url, parse_platforms = False, DEBUG = 0):
    blocksize = 1 << 16
    n=1000
    if DEBUG >= 2:
        print(f'Parsing: {url}')
    if parse_platforms == False and re.search('.*GPL.*', url):
        print('URL is to a platform and platform flag is off.')
        column_list = ['sample_id', 'source', 'title', 'molecule', 'treatment_protocol', 'extract_protocol', 'growth_protocol', 'characteristic_list']
        return pd.DataFrame(columns=column_list)
    ### Files are gzipped so open the url then extract only the XML file from the tar file
    xml_filename = tar_gz_extracter(url)
    ### If you can't extract the XML, return a blank dictionary
    if xml_filename is None:
        return dict()

    tree_ab = ET.parse(xml_filename)
    root = tree_ab.getroot()
    os.unlink(xml_filename)

    #### I stole this code from the internet but it will replace the prefix that gets added
    #### There's probably a better way but doing this let me get to the data"""
    try:
        for elem in root.getiterator():
            if not hasattr(elem.tag, 'find'): continue  # (1)
            i = elem.tag.find('}')
            if i >= 0:
                elem.tag = elem.tag[i+1:]
    ### Used to fail here often but recently fixed
    except Exception as e:
        if DEBUG >= 1:
            print(f'First {n} chars of XML')
            print(file_content[:n])
            print(e)
            print(f'Failed on url: {url}')
        return []
    ####



    sample_row_list = []
    for samples in root.findall('./Sample'):
        sample_id = ''
        sample_source = ''
        sample_title = ''
        molecule = ''
        treatment_protocol = ''
        extract_protocol = ''
        growth_protocol = ''
        sample_cell_type = ''
        sample_sex = ''
        sample_tissue = ''
        sample_age = ''

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
                        if char_value.attrib['tag'] == 'cell type':
                            sample_cell_type = char_value.text.strip()
                        if char_value.attrib['tag'] == 'sex':
                            sample_sex = char_value.text.strip()
                        if char_value.attrib['tag'] == 'tissue':
                            sample_tissue = char_value.text.strip()
                        if char_value.attrib['tag'] == 'age':
                            sample_age = char_value.text.strip()
                    except:
                        pass

            else:
                if DEBUG >= 3:
                    print('Tag:', sample_element.tag)

        sample_dict = {'sample_id' : sample_id,
                        'sample_source' : sample_source,
                        'sample_title' : sample_title,
                        'molecule' : molecule,
                        'treatment_protocol' : treatment_protocol,
                        'extract_protocol' : extract_protocol,
                        'growth_protocol' : growth_protocol,
                        'sample_cell_type' : sample_cell_type,
                        'sample_sex' : sample_sex,
                        'sample_tissue' : sample_tissue,
                        'sample_age' : sample_age
                        }
        sample_row_list.append(sample_dict)

    return sample_row_list

"""if __name__ == 'main':
    url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133304/miniml/GSE133304_family.xml.tgz'

    column_list = ['sample_id', 'sample_source', 'sample_title', 'molecule', 'treatment_protocol', 'extract_protocol', 'growth_protocol', 'characteristic_list']

    fetch_df = pd.DataFrame(xml_parser(url), columns = column_list)

    print(fetch_df.head())"""
