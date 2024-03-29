import time, os, requests, datetime, csv, re
import pandas as pd
import xml.etree.ElementTree as ET

def get_sample_data(query, retmax=50, sort='relevance', api_key="", DEBUG=0):
    """
    This function will query the NCBI e-utils API using the `query` parameter, write the results to a text file and return the file name.
    Works in two stages:
    - E-search will send the query to the API and return a QueryKey and WebEnv.
    - The QueryKey and WebEnv are fed into the GDS E-fetch which will return a text file with Sample,

    There are try: except: loops in each of the request.get sections are these will sometimes fail on big runs. These loops are currently untested as of 11/01/19.

        Args:
            `query` - query phrase to send to the API
            `retmax` - Maximum results to pull. `Samples` in GEO seem to have 3 max (Series, Platform, sample)
            `sort` - Sort order for results. Not relevant for GEO searches
            `api_key` - Supply API key for the API request

        Returns:
            `file_name_fetch` - Filename of text file downloaded from e-fetch
    """

    now = datetime.datetime.now()
    if DEBUG >= 1:
        print(query)

    query_id = query
    file_name_fetch = 'output/txt/' + query + '_' + now.strftime("%y-%m-%d-%H%M") + '_fetch.txt'
    esearch_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

    ### More DB options here : https://www.ncbi.nlm.nih.gov/books/NBK3837/
    db = '?db=' + 'gds'
    query_kw = '&term=' + query
    hist_api = "&usehistory=y"
    ret_max = '&retmax=' + str(retmax)
    sort = '&sort=' + sort

    ### Get the webpage with the IDs for the articles you'll want to fetch
    if api_key != "":
        url_search = esearch_base + db + query_kw + hist_api + ret_max + sort + f'&api_key={api_key}'
    else:
        url_search = esearch_base + db + query_kw + hist_api + ret_max + sort

    if DEBUG >= 2:
        print(f'Search URL : {url_search}')

    try:
        docsearch_resp = requests.get(url_search)
    except requests.exceptions.ConnectionError as e:
        print(e)
        print('Sleeping for 1 min then retrying')
        time.sleep(61)
        docsearch_resp = requests.get(url_search)


    ### Search the results
    root_search = ET.fromstring(docsearch_resp.content)
    QK = "&query_key=" + root_search.findall('./QueryKey')[0].text
    WE = "&WebEnv=" + root_search.findall('./WebEnv')[0].text

    ### Get Abstracts with efetch
    efetch_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    rettype_mode = "&rettype=abstract&retmode=xml"

    if api_key != "":
        url_ab = efetch_base + '?db=gds' + QK + WE + rettype_mode + ret_max  + f'&api_key={api_key}'
    else:
        url_ab = efetch_base + '?db=gds' + QK + WE + rettype_mode + ret_max
    try:
        docsab_resp = requests.get(url_ab)
    except requests.exceptions.ConnectionError as e:
        print(e)
        print('Sleeping for 1 min then retrying')
        time.sleep(61)
        docsab_resp = requests.get(url_ab)

    with open(file_name_fetch, 'wb') as f:
        f.write(docsab_resp.content)

    return file_name_fetch
