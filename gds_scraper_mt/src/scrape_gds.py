import xml.etree.ElementTree as ET
import json, os, time, sys, re, requests
from tqdm import tqdm
import pandas as pd
from src.search_samples import get_sample_data
from src.geo_parser import geo_txt_parse
from src.ftp_gzxml_parser import xml_parser
from src.geo_extraction_funcs import *

def scrape_gds(query_terms,
                api_key,
                DEBUG = 0,
                multichannel = False,
                cells_flag = False,
                metadata_filter = False,
                out_path = '/output/test',
                keep_files = [None],
                out_types = ['json', 'csv']):

    """
    Function to take in an iterator of sample IDs and output either a .json or .csv file of the sample, series and platform data and the metadata associated with the sample's series.
    In `src.search_samples.get_sample_data()`, the NCBI API is used to download text files of samples using the GSM ID. This will return the filename of the downloaded data.
    The list of filenames and associated sampleIDs is fed to `src.geo_parser.geo_txt_parse()` which will extract the sample's sample, series and platform data. This function returns a dictionary with GSM ID keys and data dictionary values.
    The series FTP links are fed individual to `src.ftp_gzxml_parser.xml_parser()` to gather metadata on samples in each series. This function returns a dictionary of samples from that series.
    Before the data is returned, a final parsing loop is run that will update the samples with the sample's metadata, run regex checks over certain text blocks, and run the age-extraction function.
    The data is returned in the format specified by `out_types` at the location specified in `out_path`.

        Args:
            `query_terms` - Series/list: Terms to feed into the API
            `api_key` - Str: If you have an API key, call limit goes up to 10/s so you could lower time.sleep(1) to time.sleep(.5) which may be worth it for pulling 10,000+ records at once
            `DEBUG` - Int: Adds some DEBUG print options. Being phased out. Soon to deprecate.
            `multichannel` - Boolean: Flag telling the program to skip samples with multichannels found in meta-data.
            `cells_flag` - Boolean: Flag telling the program to skip samples with 'cells' data found in meta-data.
            `metadata_filter` - Boolean: Flag telling the program to skip samples without any series metadata
            `out_path` - Str: Output path, does not include extension as the correct extension will be added. Defaults to /output/test
            `keep_files` - List: Flags to keep certain files used during the process. Options include 'txt' and 'xml'.
            `out_types` - List: Output file type. Options include 'json' and 'csv'

        Returns:
            `text_file_dict` - Dict: Contains {sampleID : data} key-value pairs for all the requested samples.

            This function will also create output files.
                - File type is decided by files in `out_type`. By default, create both CSV and JSON.
                - Files will be named based on `out_path` with the corresponding extension added.
                - Files have sampleID as key/index and then associated data as values/columns.
    """

    for file_type in keep_files:
        if (file_type != 'txt') and (file_type != 'xml') and (file_type is not None):
            raise ValueError(f'{file_type} is not a valid option for `keep_files`. Please enter either `xml` and/or `csv`.')
    for file_type in out_types:
        if (file_type != 'json') and (file_type != 'csv') and (file_type is not None):
            raise ValueError(f'{file_type} is not a valid option for `out_type`. Please enter either `csv` and/or `json`.')

    ### `parse_list` elements will be lists containing [filename of sample API data, sample ID]
    ### [GSM4667_19-10-31-1954_fetch.txt, GSM4667]
    parse_list=[]
    print('Creating sample .txt files')
    query_terms = list(set(query_terms))
    for query_term in tqdm(query_terms, total=len(query_terms)):
        search_ids_file = get_sample_data(query=query_term, api_key=api_key)
        parse_list.append([search_ids_file, query_term])
        time.sleep(.5)

    """
    `txt_file_dict` elements are lists containing the following:
        Include description and cells so later you can filter on them without an Error
        txt_file_dict[sample_id] = {'sample_id' : query_name,
                                'series' : series,
                                'series_accession' : series_accession,
                                'series_ftp' : series_ftp,
                                'platform' : platform,
                                'platform_accession' : platform_accession,
                                'platform_ftp' : platform_ftp,
                                'sample' : sample,
                                'contents' : contents,
                                'description' : '',
                                'cells' : ''}

    """
    print('Parsing the .txt files.')
    text_file_dict = geo_txt_parse(parse_list, keep_files = keep_files)

    ### Lists for later
    samples_series_accession_numbers = [i['series_accession'] for i in text_file_dict.values()]
    samples_series_ftp_links = [i['series_ftp'] for i in text_file_dict.values()]

    """
    Prepare a slot in the metadata dictionary for all the requested samples.
        - Add blank sample_age to each one to use it for the geoAgeExtract()
    Also created a series -> PMID data dict using all the requested samples series' number
        - Also add blank PMID, Series Summary and Design
    """
    samples_metadata_dict = {key:{'sample_age' : ''} for key in query_terms}
    series_pmid_dict = {key: {'pmid':'', 'series_summary': '', 'series_design':''} for key in samples_series_accession_numbers}

    print('Downloading and collecting meta-data')
    for url in tqdm(list(set(samples_series_ftp_links)), total=len(list(set(samples_series_ftp_links)))):
        if url == None:
            ### Empty dict if broken URL
            series_meta_dict = dict()
            xml_pmid_list = []
        else:
            """
            `xml_parser` will take in a Series FTP Download URL and return,
                -   `xml_pmid_list` is a dictionary of series -> PMID data mappings
                -   `series_meta_dict` is a dictionary of the series meta data
            """
            xml_pmid_list, series_meta_dict = xml_parser(url, sample_list=query_terms, multichannel = multichannel, keep_files = keep_files)
        ### Update relevant dictionaries
        series_pmid_dict.update(xml_pmid_list)
        samples_metadata_dict.update(series_meta_dict)

    ### Run through all the requested samples, add metadata and PMID to `text_file_dict`
    print('Final parsing loop')
    text_file_dict_copy = text_file_dict.copy()
    for sample_id in tqdm(text_file_dict_copy.keys(), total=len(text_file_dict_copy.keys())):

        text_file_dict[sample_id].update(samples_metadata_dict[sample_id])
        ### If multichannel == False, `xml_parser` will mark the multichannel samples in the metadata,
        ### this actually removes them.
        if (multichannel == False) and (text_file_dict[sample_id]['description'] == 'multi'):
            text_file_dict.pop(sample_id, None)
            continue

        ### Filters samples without metadata
        if (metadata_filter == True) and (text_file_dict[sample_id]['description'] == ''):
            text_file_dict.pop(sample_id, None)
            continue

        ### Filters cell samples
        if (cells_flag == False) and (text_file_dict[sample_id]['cells'] == True):
            text_file_dict.pop(sample_id, None)
            continue

        sample_series_acc_num = text_file_dict[sample_id]['series_accession']
        series_metadict = series_pmid_dict[sample_series_acc_num]

        ### Hoping to create a high level split for different types of samples but didn't get very far into it
        """if re.findall('(wild\s?type)', text_file_dict[sample_id]['sample'].lower()):
            text_file_dict[sample_id]['wild_type'] = True
        else:
            text_file_dict[sample_id]['wild_type'] = False

        if re.findall('(treatment[a-z\s]*with)|(treated[a-z\s]with)', text_file_dict[sample_id]['sample'].lower()):
            text_file_dict[sample_id]['molecule_bool'] = True
        else:
            text_file_dict[sample_id]['molecule_bool'] = False

        if re.findall('([+-]/[+-])', text_file_dict[sample_id]['sample'].lower()):
            text_file_dict[sample_id]['ko_bool'] = True
        else:
            text_file_dict[sample_id]['ko_bool'] = False"""

        ### If the sample's series has a PMID,
        if series_metadict["pmid"] != '':
            text_file_dict[sample_id]['pmid'] = f'https://www.ncbi.nlm.nih.gov/pubmed/{series_metadict["pmid"]}'

            text_file_dict[sample_id]['series_summary'] = series_metadict['series_summary']
            text_file_dict[sample_id]['series_design'] = series_metadict['series_design']

            ### May be worth doing running this section on the set of PMID links rather than every sample
            root = ET.fromstring(requests.get(text_file_dict[sample_id]['pmid']).content)
            for elem in root.getiterator():
                if not hasattr(elem.tag, 'find'): continue
                i = elem.tag.find('}')
                if i >= 0:
                    elem.tag = elem.tag[i+1:]
            try:
                ### Index [0] of findall results is the title <h3>Abstract</h3>
                text_file_dict[sample_id]['abstract_article'] = list(root.findall('.//div[@class="rprt_all"]//div[@class="abstr"]/')[1].itertext())[0]
            except:
                text_file_dict[sample_id]['article_abstract'] = ''
            try:
                text_file_dict[sample_id]['article_title'] = root.findall('.//div[@class="rprt_all"]//h1')[0].text
            except:
                text_file_dict[sample_id]['article_title'] = ''
        ### If no PMID,
        else:
            text_file_dict[sample_id]['pmid'] = ''
            text_file_dict[sample_id]['series_summary'] = ''
            text_file_dict[sample_id]['series_design'] = ''
            text_file_dict[sample_id]['abstract_article'] = ''
            text_file_dict[sample_id]['article_title'] = ''

        text_file_dict[sample_id]['age_func'] = geoAgeExtract(text_file_dict[sample_id])

    ### Export data
    if 'json' in out_types:
        with open(f'{out_path}.json', 'w') as fout:
            json.dump(text_file_dict, fout, indent = 4)
    if 'csv' in out_types:
        pd.DataFrame.from_dict(text_file_dict, orient='index').to_csv(f'{out_path}.csv')

    return text_file_dict
