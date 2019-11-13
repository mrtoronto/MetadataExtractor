from src.geo_extraction_funcs import *
from tqdm import tqdm
import requests
import xml.etree.ElementTree as ET

def final_processing_loop(text_file_dict, samples_metadata_dict, series_pmid_dict, multichannel, metadata_filter, cells_flag):

    text_file_dict_copy = text_file_dict.copy()
    for sample_id in tqdm(text_file_dict_copy.keys(), total=len(text_file_dict_copy.keys())):

        text_file_dict[sample_id].update(samples_metadata_dict[sample_id])

        text_file_dict[sample_id]['cells'] = geoSampleCellCheck(text_file_dict[sample_id])
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

        ### Split for different categories of samples
        wildtype_bool_pattern = '(((wild)|(Wild)|([Ww]))[ -_]?((type)|(Type)|([Tt])))'
        if (re.findall(wildtype_bool_pattern, text_file_dict[sample_id]['sample'].lower())) or (re.findall(wildtype_bool_pattern, text_file_dict[sample_id]['description'].lower())):
            text_file_dict[sample_id]['wild_type'] = True
        else:
            text_file_dict[sample_id]['wild_type'] = False

        molecule_bool_pattern = '(treatment[a-z\s]*with)|(treated[a-z\s]with)'
        if re.findall(molecule_bool_pattern, text_file_dict[sample_id]['sample'].lower()):
            text_file_dict[sample_id]['molecule_bool'] = True
        else:
            text_file_dict[sample_id]['molecule_bool'] = False

        ko_bool_pattern = '(([Kk][-_ ]?[Oo])|(((knock)|(Knock))[- _]?((Out)|(out))))'
        if re.search(ko_bool_pattern, text_file_dict[sample_id]['sample'].lower()):
            text_file_dict[sample_id]['ko_bool'] = True
        else:
            text_file_dict[sample_id]['ko_bool'] = False

        if text_file_dict[sample_id]['ko_bool'] == True:
            ### pattern built from `sample` field
            ko_pattern3 = "(((.*)[ -_]?){1,3})(([kK])|(knock)|(Knock))[- _]?((out)|(Out)|([Oo]))"
            ### Group 1 grabs everything before the knockout/ko/Knock-Out group
            for pattern, key, group_val in [(ko_pattern3, 'ko_gene3', 1)]:
                if re.search(pattern, text_file_dict[sample_id]['sample']):
                    text_file_dict[sample_id][key] = re.search(pattern, text_file_dict[sample_id]['sample']).group(group_val)

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

    return text_file_dict
