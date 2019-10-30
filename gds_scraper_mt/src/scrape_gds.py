from tqdm import tqdm
import time
import pandas as pd
from src.search_samples import get_sample_data
from src.geo_parser import geo_txt_parse
import sys
import os
from src.ftp_gzxml_parser import xml_parser

"""
Arguments:
    `query_terms`     - Series/list of terms to feed into the API
    `res_per_query`   - Max number of results one query will pull. Defaults to 50 and shouldn't need to be touched.
    `api_key`         - If you have an API key, call limit goes up to 10/s so you could lower time.sleep(1) to time.sleep(.5) which may be worth it for pulling 10,000+ records at once
    `out_path`        - Output path, includes filename, defaults to /output/test.csv


Output:
    files created:
        `out_path`            - This file will have the final table with Samples as the index and associated metadata merged

"""


def scrape_gds_to_csv(query_terms,
                        api_key,
                        res_per_query = 50,
                        DEBUG = 0,
                        multichannel = False,
                        out_path = '/output/test.csv',
                        keep_files = [None]):

    if DEBUG >= 1:
        print(f'length of query_terms: {len(query_terms)}.')

    ### List elements will be lists containing [filename of sample API data, sample ID]
    parse_list=[]

    print('Creating sample .txt files')
    for query_term in tqdm(query_terms, total=len(query_terms)):
        search_ids_list = get_sample_data(
                    have_ids=False,
                    query=query_term,
                    retmax=res_per_query,
                    sort='relevance',
                    filename=None,
                    api_key=api_key,
                    DEBUG=0
                )
        parse_list.append(search_ids_list)
        time.sleep(1)

    """
    `txt_file_list` elements are lists containing the following:
        [query_name, series, series_accession, series_ftp, platform, platform_accession, platform_ftp, sample, contents]
    """
    print('Parsing the .txt files.')
    txt_file_list = geo_txt_parse(parse_list, DEBUG=False, keep_files = keep_files)

    samples_column_list = ['query_name', 'series_text', 'series_accession', 'series_ftp', 'platform_text', 'platform_accession', 'platform_ftp', 'sample_text', 'contents']
    samples = pd.DataFrame(txt_file_list, columns=samples_column_list).set_index('query_name')

    if DEBUG >= 2:
        print(f'samples shape: {samples.shape}')


    """
    `series_ftp_list` elements are dictionaries containing the following:
        {'sample_id' : sample_id,
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
    """
    series_ftp_list=[]
    samples['pmid'] = ''
    series_pmid_dict = {key:value for key,value in zip(samples['series_accession'],samples['pmid'])}

    if DEBUG >= 1:
        ftp_link_len = len(samples['series_ftp'])
        uni_ftp_link_len = len(samples['series_ftp'].unique())
        print(f'Number of FTP Series links: {ftp_link_len}     //  Number of unique FTP Series links: {uni_ftp_link_len}')

    print('Downloading and collecting meta-data')
    for url in tqdm(samples['series_ftp'].unique(), total=len(samples['series_ftp'].unique())):
        if url == None:
            ### Empty dict if broken URL
            url_meta = dict()
            xml_pmid_list = []
        else:
            """
            Returns list of dictionaries
                 Each dict == 1 Sample in this Series
                 Could give it a list of already parsed samples so its not re-parsing anything but I skipped this for now
            """
            url_meta, xml_pmid_list = xml_parser(url, parse_platforms = False, DEBUG=0, multichannel = multichannel, keep_files = keep_files)
        ### Series_ftp_list == list of list of dictionaries
        series_ftp_list += url_meta
        series_pmid_dict.update(xml_pmid_list)
    ### Flattening to list of dictionaries
    #series_ftp_list = [series_dict for sublist in series_ftp_list for series_dict in sublist]
    #series_pmid_list = [pmid_dict for sublist in series_pmid_list for pmid_dict in sublist]
    #metadata_columns = ['sample_id', 'sample_source', 'sample_title', 'molecule', 'treatment_protocol', 'extract_protocol', 'growth_protocol', 'sample_cell_type', 'sample_sex', 'sample_tissue', 'sample_age']
    samples_metadata = pd.DataFrame(series_ftp_list, columns=series_ftp_list[0].keys()).set_index('sample_id')

    if DEBUG >= 2:
        print(f'samples_metadata shape: {samples_metadata.shape}')

    merg_samples = samples.merge(samples_metadata, left_index = True, right_index = True, how='left')
    if multichannel == False:
        merg_samples = merg_samples[merg_samples['description'] != 'multi']
    merg_samples['pmid'] = [series_pmid_dict[i] for i in merg_samples['series_accession']]
    ### This requires no lists in the DataFrame
    ### Not normally an issue to turn it off but I like it better with it on
    merg_samples = merg_samples.drop_duplicates()

    if DEBUG >= 1:
        print(f'samples shape post-merge: {merg_samples.shape}')
    ### Set up with date to save multiple runs, I don't because I'm doing lots of runs
    if out_path is not None:
        merg_samples.to_csv(out_path, index=True)
    if DEBUG >= 2:
        meta_csv_filename = 'output/meta_debug.csv'
        presample_csv_filename = 'output/pre_debug.csv'
        samples_metadata.to_csv(meta_csv_filename, index=True)
        samples.to_csv(presample_csv_filename, index=True)
    return merg_samples
