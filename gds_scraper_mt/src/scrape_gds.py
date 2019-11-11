import xml.etree.ElementTree as ET
import json, os, time, sys, re, requests
from tqdm import tqdm
import pandas as pd
from src.search_samples import get_sample_data
from src.geo_parser import geo_txt_parse
from src.ftp_gzxml_parser import xml_parser
from src.geo_extraction_funcs import *
from src.processing_data import final_processing_loop

def scrape_gds(query_terms,
                api_key,
                DEBUG = 0,
                multichannel = False,
                cells_flag = False,
                metadata_filter = False,
                out_path = '/output/test',
                keep_files = [None],
                run_type = 'append',
                out_types = ['json', 'csv']):

    """
    Function to take in an iterator of sample IDs and output either a .json or .csv file of the sample, series and platform data and the metadata associated with the sample's series.
    In `src.search_samples.get_sample_data()`, the NCBI API is used to download text files of samples using the GSM ID. This will return the filename of the downloaded data.
    The list of filenames and associated sampleIDs is fed to `src.geo_parser.geo_txt_parse()` which will extract the sample's sample, series and platform data. This function returns a dictionary with GSM ID keys and data dictionary values.
    The series FTP links are fed individual to `src.ftp_gzxml_parser.xml_parser()` to gather metadata on samples in each series. This function returns a dictionary of samples from that series.
    Before the data is returned, a final parsing loop is run that will update the samples with the sample's metadata, run regex checks over certain text blocks, and run the age-extraction function.
    The data is returned in the format specified by `out_types` at the location specified in `out_path`. If `run_type` is set to `append`, progam will check for a pre-existing file to append to and create one if there isn't an existing option. Setting `run_type` to `new` will overwrite any existing file at `out_path`.

        Args:
            `query_terms` - Series/list: Terms to feed into the API
            `api_key` - Str: If you have an API key, call limit goes up to 10/s so you could lower time.sleep(1) to time.sleep(.5) which may be worth it for pulling 10,000+ records at once
            `DEBUG` - Int: Adds some DEBUG print options. Being phased out. Soon to deprecate.
            `multichannel` - Boolean: Flag telling the program to skip samples with multichannels found in meta-data.
            `cells_flag` - Boolean: Flag telling the program to skip samples with 'cells' data found in meta-data.
            `metadata_filter` - Boolean: Flag telling the program to skip samples without any series metadata
            `out_path` - Str: Output path, does not include extension as the correct extension will be added. Defaults to /output/test
            `keep_files` - List: Flags to keep certain files used during the process. Options include 'txt' and 'xml'.
            `run_type` - Str: Set to `new` to create a new file at `out_path`. Set to `append` to append to an existing file at `out_path`.
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
    if (run_type != 'new') and (run_type != 'append'):
        raise ValueError(f'{run_type} is not a valid option for `run_type`. Please enter either `new` or `append`.')
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

    ### Do final processing on samples
    ### Start by adding metadata then adding derived fields
    print('Final parsing loop')
    text_file_dict = final_processing_loop(text_file_dict, samples_metadata_dict, series_pmid_dict, multichannel, metadata_filter, cells_flag)

    ### Export data
    if run_type == 'new':
        if 'json' in out_types:
            with open(f'{out_path}.json', 'w') as fout:
                json.dump(text_file_dict, fout, indent = 4)
        if 'csv' in out_types:
            pd.DataFrame.from_dict(text_file_dict, orient='index').to_csv(f'{out_path}.csv')
    elif run_type == 'append':
        if 'json' in out_types:
            try:
                with open(f'{out_path}.json', 'r') as f:
                    pre_existing_data = json.load(f)
            except:
                print('No pre-existing data.')
                pre_existing_data = dict()
            text_file_dict = {k : v for (k,v) in text_file_dict.items() if k not in pre_existing_data.keys()}
            pre_existing_data.update(text_file_dict)
            with open(f'{out_path}.json', 'w') as fout:
                json.dump(pre_existing_data, fout, indent = 4)
        if 'csv' in out_types:
            pre_existing_data = pd.read_csv(f'{out_path}.csv', header=0, index_col=0).to_dict(orient='index')
            text_file_dict = {k : v for (k,v) in text_file_dict.items() if k not in pre_existing_data.keys()}
            pre_existing_data.update(text_file_dict)

            pd.DataFrame.from_dict(pre_existing_data, orient='index').to_csv(f'{out_path}.csv')
    return text_file_dict
