import os
import pandas as pd
import xml.etree.ElementTree as ET
import requests
import datetime
import csv
import re

"""

Arguments:
    query       - query phrase to send to the API
    filename    - If user has a list of UIDs in a file, filename would go here. (DEPRECATED)
    retmax      - Maximum results to pull. `Samples` in GEO seem to have 3 max (Series, Platform, sample)
    sort        - Sort order for results. Not relevant for GEO searches
    have_ids    - Boolean if user has a list of file IDs
    api_key     - Supply API key for the API request


"""

def get_sample_data(query, filename, retmax, sort, have_ids=False, api_key="", DEBUG=0):

    now = datetime.datetime.now()
    if DEBUG >= 1:
        print(query)
    if os.path.isdir("output") == False:
        os.mkdir('output')
    if os.path.isdir("output/txt") == False:
        os.mkdir('output/txt')

    ### If `have_ids` == TRUE, user has a list of UIDs saved at `filename` already
    if (have_ids == True):

        epost_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=gds&id='

        ### Read in the file IDs
        with open(filename, 'r') as f:
            uid_list = list(csv.reader(f))

        ### Create the URL to query
        url_post = epost_base
        if DEBUG >= 2:
            print(url_post)

        ### Get webpage with QK and WebEnv values for fetching articles
        docsearch_resp = requests.get(url_post)

        ### Set to blank because they don't apply with ID list
        ret_max, query_str = '', ''

        ### Build filename for the _fetch .XML file
        file_name_fetch = 'output/xml/' + filename[:-4] + '_idlist_' + now.strftime("%y%m%d_%H%M") + '_fetch.txt'

    ### If have_ids != TRUE, then user wants to query using a search term/phrase
    else:
        query_id = query
        # file_name_search = query + now.strftime("%y%m%d_%H%M") + '_search.xml'
        file_name_fetch = 'output/txt/' + query + '_' + now.strftime("%y-%m-%d-%H%M") + '_fetch.txt'
        esearch_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

        ### More DB options here : https://www.ncbi.nlm.nih.gov/books/NBK3837/
        db = '?db=' + 'gds'
        query = '&term=' + query
        hist_api = "&usehistory=y"
        ret_max = '&retmax=' + str(retmax)
        sort = '&sort=' + sort

        ### Get the webpage with the IDs for the articles you'll want to fetch
        url_search = esearch_base + db + query + hist_api + ret_max + sort
        if DEBUG >= 2:
            print(f'Search URL : {url_search}')

        docsearch_resp = requests.get(url_search)



    ### Search the results
    root_search = ET.fromstring(docsearch_resp.content)
    #root_search = tree_search.getroot()
    QK = "&query_key=" + root_search.findall('./QueryKey')[0].text
    WE = "&WebEnv=" + root_search.findall('./WebEnv')[0].text


    ### Get Abstracts with efetch
    efetch_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    rettype_mode = "&rettype=abstract&retmode=xml"

    url_ab = efetch_base + '?db=gds' + QK + WE + rettype_mode + ret_max
    if DEBUG >= 2:
        print(f'Data URL : {url_ab}')

    docsab_resp = requests.get(url_ab)
    with open(file_name_fetch, 'wb') as f:
        f.write(docsab_resp.content)

    ret_list = [file_name_fetch, query_id]
    return ret_list
