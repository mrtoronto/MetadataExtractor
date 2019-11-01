import json
from tqdm import tqdm
import time
import pandas as pd
from src.search_samples import get_sample_data
from src.geo_parser import geo_txt_parse
import sys
import os
from src.ftp_gzxml_parser import xml_parser

"""
Parameters:
    `query_terms`     - Series/list of terms to feed into the API
    `api_key`         - If you have an API key, call limit goes up to 10/s so you could lower time.sleep(1) to time.sleep(.5) which may be worth it for pulling 10,000+ records at once
    `DEBUG`           - Adds some DEBUG print options. Being phased out. Soon to deprecate.
    `multichannel`    - Flag telling the program to skip samples with multichannels found in meta-data.
    `out_path`        - Output path, does not include extension as the correct extension will be added. Defaults to /output/test
    `keep_files`      - Option to keep certain files used during the process. Options include 'txt' and 'xml'.
    `out_types`       - File type of output files. Options include 'json' and 'csv'

Output:
    - File type is decided by files in `out_type`. By default, create both CSV and JSON.
    - Files will be named based on `out_path` with the corresponding extension added.
    - Files have sampleID as key/index and then associated data as values/columns.


"""

def scrape_gds(query_terms,
                api_key,
                DEBUG = 0,
                multichannel = False,
                out_path = '/output/test',
                keep_files = [None],
                out_types = ['json', 'csv']):

    for file_type in keep_files:
        if (file_type != 'txt') and (file_type != 'xml'):
            raise ValueError(f'{file_type} is not a valid option for `keep_files`. Please enter either `xml` and/or `csv`.')
    for file_type in out_types:
        if (file_type != 'json') and (file_type != 'csv'):
            raise ValueError(f'{file_type} is not a valid option for `out_type`. Please enter either `csv` and/or `json`.')

    ### `parse_list` elements will be lists containing [filename of sample API data, sample ID]
    ### [GSM4667_19-10-31-1954_fetch.txt, GSM4667]
    parse_list=[]
    print('Creating sample .txt files')
    for query_term in tqdm(query_terms, total=len(query_terms)):
        search_ids_list = get_sample_data(query=query_term, api_key=api_key)
        parse_list.append(search_ids_list)
        time.sleep(1)

    """
    `txt_file_dict` elements are lists containing the following:
        txt_file_dict[sample_id] = {'sample_id' : sample_id,
                                'series' : series,
                                'series_accession' : series_accession,
                                'series_ftp' : series_ftp,
                                'platform' : platform,
                                'platform_accession' : platform_accession,
                                'platform_ftp' : platform_ftp,
                                'sample' : sample,
                                'contents' : contents}

    """
    print('Parsing the .txt files.')
    text_file_dict = geo_txt_parse(parse_list, keep_files = keep_files)

    samples_series_accession_numbers = [i['series_accession'] for i in text_file_dict.values()]
    samples_series_ftp_links = [i['series_ftp'] for i in text_file_dict.values()]

    """
    Prepare a slot in the metadata dictionary for all the requested samples.
    Also created a series -> PMID dictionary using all the requested samples series' number
    """
    samples_metadata_dict = {key:dict() for key in query_terms}
    series_pmid_dict = {key:'' for key in samples_series_accession_numbers}

    print('Downloading and collecting meta-data')
    for url in tqdm(list(set(samples_series_ftp_links)), total=len(list(set(samples_series_ftp_links)))):
        if url == None:
            ### Empty dict if broken URL
            series_meta_dict = dict()
            xml_pmid_list = []
        else:
            """
            `xml_parser` will take in a Series FTP Download URL and return,
                -   `xml_pmid_list` is a dictionary of series -> PMID mappings from the file
                -   `series_meta_dict`
            """
            xml_pmid_list, series_meta_dict = xml_parser(url, multichannel = multichannel, keep_files = keep_files)
        ### Update relevant dictionaries
        series_pmid_dict.update(xml_pmid_list)
        samples_metadata_dict.update(series_meta_dict)


    ### Run through all the requested samples, add metadata and PMID to `text_file_dict`
    for sample_id in text_file_dict.keys():
        text_file_dict[sample_id].update(samples_metadata_dict[sample_id])
        text_file_dict[sample_id]['pmid'] = series_pmid_dict[text_file_dict[sample_id]['series_accession']]

    ### If multichannel == False, `xml_parser` will mark the multichannel samples in the metadata,
    ### this dictionary comprehension actually removes them.
    if multichannel == False:
        text_file_dict = {k:v for (k,v) in text_file_dict.items() if v['description'] != 'multi'}

    ### Export data
    if 'json' in out_types:
        with open(f'{out_path}.json', 'w') as fout:
            json.dump(text_file_dict, fout, indent = 4)
    if 'csv' in out_types:
        pd.DataFrame.from_dict(text_file_dict, orient='index').to_csv(f'{out_path}.csv')

    return text_file_dict
