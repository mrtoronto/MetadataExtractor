from tqdm import tqdm
from lxml import etree, objectify
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import tarfile, os, urllib, gzip, re
from src.geo_extraction_funcs import *

def xml_parser(url = "", filename = "", sample_list = [], parse_platforms = False, DEBUG = 0, multichannel = False, keep_files = [None]):
    """
    This function will take in the URL for a Series FTP .tgz file, download the file, extract the .xml file from the compressed directory and then parse the .xml.

    Try: except lxml.etree.XMLSyntaxError as e: section is to prevent illegal XML characters from breaking the script. Ideally, we'd replace all illegal characters before parsing but I couldn't find good code to do that.

    First the ./Series tags are looped through and the corresponding pubmed ID is taken from each series.

    Then the ./Samples tags are looped through and all relevant information is captured for each sample.

    If multichannel == False and the parser finds more than 1 channel-count tag within a sample, it'll append an essentially empty data dict knowing it'll be filtered later.
        Args:
            `url` - Str: Series FTP URL leading to a .tgz file. Blank if the file is already downloaded.
            `filename` - Str: File location of XML file to be parsed. Used when program finds a relevant cached file.
            `parse_platforms` - Bool: If set to its default value, False, a platform URL will be rejected instead of parsed poorly.
            `DEBUG` - Int: Sets depth of debug messages. Mostly deprecated.
            `multichannel` - Bool: If set to False, multichannel files will be filtered out of the final result.
            `keep_files` - List: File types to keep during a full run. If 'xml' is in `keep_files` then the xml file won't be deleted.

        Returns:
            `series_pmid_dict` - Dict: Mappings from Series Accession numbers to PMIDs
            `samples_dict` - Dict: {sampleID : data} key-value pairs for all samples in input list.
                - `data` is a dictionary of the following information:
                    ***Maybe fill in these slots as fields become useful for specific things.***
                    {'sample_id' : GSM sample ID,
                        'sample_source' : ,
                        'sample_title' : ,
                        'molecule' : ,
                        'organism' : organism,
                        'treatment_protocol' : ,
                        'extract_protocol' : ,
                        'growth_protocol' : ,
                        'description' : ,
                        'sample_cell_type' : Can be used as a check to filter cell experiments,
                        'sample_type' : ,
                        'sample_sex' : ,
                        'sample_tissue' : ,
                        'sample_age' : Sometimes has age but needs to be parsed and converted.
                        'sample_indication' : ,
                        'sample_genotype' : ,
                        'sample_cell_line' : Can be used as a check to filter cell experiments,
                        'expression' : ,
                        'cells' : Flag for any positive cell experiment check
                        }
    """


    n=1000
    if DEBUG >= 2:
        print(f'Parsing: {url}')
    if parse_platforms == False and re.search('.*GPL.*', url):
        print(f'{url} is to a platform and platform flag is off.')
        return dict()

    ### Files are gzipped so open the url then extract only the XML file from the tar file
    if filename == "":
        xml_filename = tar_gz_extracter(url)
    elif url == "":
        xml_filename = filename

    ### If you can't extract the XML, return a blank dictionary
    if xml_filename is None:
        return {}, {}

    parser = etree.XMLParser(recover=True, encoding="utf-8")
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

    pmids, series_acc_numbers, series_summary = [], [], []
    series_overalldesign, series_zip = [], []
    for series in root.findall('./Series'):
        pmids += [i.text for i in series.findall('./Pubmed-ID')]
        series_acc_numbers += [i.text for i in series.findall('./Accession')]
        series_summary += [i.text.strip() for i in series.findall('./Summary')]
        series_overalldesign += [i.text.strip() for i in series.findall('./Overall-Design')]
        series_zip += list(zip(series_acc_numbers, pmids, series_summary, series_overalldesign))

    series_pmid_dict = {key: {'pmid':pmid, 'series_summary': series_summary, 'series_design':series_design} for (key, pmid, series_summary, series_design) in series_zip}


    sample_row_list = []
    samples_dict = {}
    series_pubmed_id_dict = {}
    for samples in root.findall('./Sample'):
        sample_id = samples.attrib['iid']
        if sample_id not in sample_list:
            continue
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
                                'cells' : '',
                                'age_func' : ''
                                }
                samples_dict[sample_id] = sample_dict

                #sample_row_list.append(sample_dict)
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
                        if char_value.attrib['tag'].lower() == 'cell type':
                            sample_cell_type = char_value.text.strip()
                        if char_value.attrib['tag'].lower() == 'cell line':
                            sample_cell_line = char_value.text.strip()
                        if char_value.attrib['tag'].lower() == 'sample type':
                            sample_type = char_value.text.strip()
                        if char_value.attrib['tag'].lower() == 'sex':
                            sample_sex = char_value.text.strip()
                        if char_value.attrib['tag'].lower() == 'tissue':
                            sample_tissue = char_value.text.strip()
                        if char_value.attrib['tag'].lower() == 'age':
                            sample_age = char_value.text.strip()
                        if char_value.attrib['tag'].lower() == 'indication':
                            sample_indication = char_value.text.strip()
                        if char_value.attrib['tag'].lower() == 'genotype':
                            sample_genotype = char_value.text.strip()
                    except Exception as e:
                        pass

        if ('RNA' in sample_source) or ('RNA' in molecule):
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
                        'age_func' : ''
                        }
        samples_dict[sample_id] = sample_dict

    return series_pmid_dict, samples_dict
