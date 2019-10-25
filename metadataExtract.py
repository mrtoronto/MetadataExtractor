""" Start small: age, gender, cell/tissue types in Mouse studies off GEO 
    All samples downloaded off https://www.ncbi.nlm.nih.gov/gds/?term=mus+musculus%5Borganism%5D

    Author: Seth Rhoades

Expected values:
GSM400641 - 4.28 - check
GSM32928 - 7.5 - check
GSM1402452 - n/a (however ages exist in the journal article but very hard to know which
                ... age corresponds to the sample. However cell type is in description so 
                    ... probably should be skipped anyways)
GSM1418737 - 13.5 - check
GSM1394796 - skipped - check (for dual channel)
GSM172972 - n/a (however ages exist in the GSE series description, should be 8 weeks. Currently 
                    two groups are sampled at different times (4-7 weeks and 8 weeks). Thus 
                    13.5 is not correct, however WT is not explicityly mentioned either. Without
                    access to article, should this be left blank or all groups just set at 13.5?)
GSM1338039 - 9.42 - check
GSM218336 - n/a - check (5hr is not valid, but treatment durations exist in the GSE description 
                ... and sample names, and ages exist in the journal article.. should be 8 weeks)
GSM1282831 - 16 - check
GSM937915 - n/a - check (for cell == True)
GSM557138 - 6 - check (however, induction began at 6 weeks and mice were only analyzed upon death, 
                ... which ranged over months, this may be edge case where the true age cannot be extracted)
"""

import re, requests, os, random, sys, json, random
import pandas as pd
sys.path.append('./src')
import setup_metadataExtract as util
refDirectory = 'refFiles'


organism = 'Mus musculus'
queryIDs = False
queryNum = 10000
overwrite = False

if queryIDs == True:
    if organism == 'Mus musculus':
        ids = pd.read_csv('{0}/GEO_{1}Samples.csv'.format(refDirectory, re.sub(' ', 
            '', organism)))['0'].values.tolist()
        ids = list(set(ids))
        sampleIDs = random.sample(ids, queryNum)
else:
    if organism == 'Mus musculus':
        sampleIDs = ['GSM400641', 'GSM32928', 'GSM1402452', 'GSM1418737', 'GSM1394796',
            'GSM172972', 'GSM1338039', 'GSM218336', 'GSM1282831', 'GSM937915', 
            'GSM557138']

if 'GEO_{0}Metadata.json'.format(re.sub(' ', '', organism)) in os.listdir(refDirectory):
    with open('{0}/GEO_{1}Metadata.json'.format(refDirectory, re.sub(' ', '',   
        organism)), 'r') as fin:
        metaDict = json.load(fin)

def main(sampleIDs, organism, metaDict):
    
    for sample in sampleIDs:
        if sample not in metaDict or overwrite is True:
            try:
                metaData = util.extractGEOSampleInfo(sampleID = sample, 
                    organism = organism)
                metaDict[sample] = metaData
            except AttributeError as err:
                print('Will continue, but take a look at this: ', str(err))
            except Exception as err:
                print('Actual error: ', str(err))

    if overwrite is True:
        with open('{0}/GEO_{1}Metadata.json'.format(refDirectory, re.sub(' ', '', 
            organism)), 'w') as fout:
            json.dump(metaDict, fout, indent = 4)
    

if __name__ == "__main__":

    main(sampleIDs, organism, metaDict)