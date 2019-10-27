import itertools
import os
import re
import pandas as pd
import xml.etree.ElementTree as ET
import requests
import datetime

"""

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSM3855772&usehistory=y&retmax=20&sort=relevance

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&query_key=1&WebEnv=NCID_1_43256172_130.14.18.97_9001_1571760482_60187628_0MetA0_S_MegaStore&rettype=abstract&retmode=xml&retmax=20

Arguments:

    - list_of_ids_querys = [ids_fetch_file, query_term]
        ids_fetch_file =

"""

def geo_txt_parse(list_of_ids_querys, DEBUG = False):

    now = datetime.datetime.now()
    row_list = []

    for filename, query_name in list_of_ids_querys:
        print('1', query_name)
        series_text_list = []
        sample_text_list = []
        platform_text_list = []
        platform_text = ''
        sample_text = ''
        if DEBUG == True:
            print(filename, query_name)

        with open(filename, 'r') as f:
            contents = f.read()

        contents = contents.strip()

        if DEBUG == True:
            print(contents)

        results_list = contents.split("\n\n")

        for result in results_list:
            rows_in_result = result.split('\n')

            if re.search('Series\s*Accession:', result):
                series_text_list.append(result)
                series_text = result

            if re.search('Platform\s*Accession:', result):
                platform_text = result
                platform_text_list.append(result)

            if re.search('Sample\s*Accession:', result):
                sample_text_list.append(result)

            title = rows_in_result[0]



        zipped_lists = list(itertools.product(series_text_list, platform_text_list, sample_text_list))
        ### Parse out specific info from the series, platforms and sample text blocks
        for series, platform, sample in zipped_lists:

            series_accession = ''
            series_ftp = ''
            platform_accession = ''
            platform_ftp = ''

            series_ftp = re.search('FTP download: .*GSE.*', series)
            series_accession = re.search('Series\s*Accession: GSE.*', series)
            try:
                series_accession = series_accession.group()
                series_accession = re.sub('.* GSE', 'GSE', series_accession)
                series_accession = series_accession.split(': ')[0][:-2].strip()
            except:
                pass
            try:
                # ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL17nnn/GPL17021/miniml/GPL17021_family.xml.tgz
                series_ftp = series_ftp.group()
                series_ftp = re.sub('.* ftp://', 'ftp://', series_ftp)
                series_ftp = series_ftp + f'miniml/{series_accession}_family.xml.tgz'
            except:
                pass

            platform_ftp = re.search('FTP download: .*GPL.*', platform)
            try:
                platform_ftp = platform_ftp.group()
                platform_ftp = re.sub('.* ftp://', 'ftp://', platform_ftp)
            except:
                pass

            platform_accession = re.search('Platform\s*Accession: GPL.*', platform)
            try:
                platform_accession = platform_accession.group()
                platform_accession = re.sub('.* GPL', 'GPL', platform_accession)
                platform_accession = platform_accession.split(': ')[0][:-2].strip()
            except:
                pass

            row = [query_name, series, series_accession, series_ftp, platform, platform_accession, platform_ftp, sample, contents]
            row_list.append(row)

        os.unlink(filename)

    return row_list
