import itertools, os, re, requests, datetime
import pandas as pd
import xml.etree.ElementTree as ET

def geo_txt_parse(loc_query_list, keep_files = [None]):
    """
    This function will take in a list of 2-element lists. Each internal list will contain a filename and a sampleID.

    Looping through this top level list, each file will be read and data will be parsed from the individual sections of the file.

    Key data on the sample, series and platform in these files includes a summary-esque block of text for each, a series and platform FTP download link and series and platform accession numbers (GSE and GPL respectively).

    To manually check the input file, submit this link, copy the WebEnv string, and paste it into the lower URL.
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSM3855772&usehistory=y&retmax=20&sort=relevance
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&query_key=1&WebEnv=NCID_1_125842368_130.14.18.97_9001_1572652184_589518123_0MetA0_S_MegaStore&rettype=abstract&retmode=xml&retmax=20

        Args:
            `loc_query_list` - List: List of Two element lists. Elements created by `src.search_samples.get_sample_data`
                - [location of text file, query_term]
                - [GSM4667_19-10-31-1954_fetch.txt, GSM4667]
            `keep_files` - List: File types to keep during a full run. If 'txt' is in `keep_files` then the txt file won't be deleted.

        Returns:
            `rows_dict` - Dict: {sampleID : data} key-value pairs for all samples in input list.
                - `data` is a dictionary of the following information:
                    {'sample_id' : query_name,
                        'series' : series,
                        'series_accession' : series_accession,
                        'series_ftp' : series_ftp,
                        'platform' : platform,
                        'platform_accession' : platform_accession,
                        'platform_ftp' : platform_ftp,
                        'sample' : sample,
                        'contents' : contents}
    """
    
    row_list = []
    rows_dict = {}

    for filename, query_name in loc_query_list:
        series_text_list = []
        sample_text_list = []
        platform_text_list = []
        platform_text = ''
        sample_text = ''

        with open(filename, 'r') as f:
            contents = f.read()

        contents = contents.strip()

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
            rows_dict[query_name] = {'sample_id' : query_name,
                                    'series' : series,
                                    'series_accession' : series_accession,
                                    'series_ftp' : series_ftp,
                                    'platform' : platform,
                                    'platform_accession' : platform_accession,
                                    'platform_ftp' : platform_ftp,
                                    'sample' : sample,
                                    'contents' : contents}
            row_list.append(row)

        if 'txt' not in keep_files:
            os.unlink(filename)

    return rows_dict
