"Basic queries of GEO studies --> sample lists"

import csv, re, requests
import pandas as pd

geoURL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
refDirectory = 'refFiles'

def main(url, directory):
        
    studyList = []
    with open('{0}/GEO_MusmusculusStudies.txt'.format(directory), newline = '') as fin:
        csvRead = csv.reader(fin, delimiter=' ', quotechar='|')
        for row in csvRead:
            rowJoin = ''.join(row)
            if 'Accession:GSE' in rowJoin:
                studyList += re.findall(r'Accession\:(GSE\d+)', rowJoin)

    sampleList = []
    counter = 0
    for geo in studyList:
        urlGetText = requests.get(url + geo).text
        sampleList += list(set(re.findall(r'acc\=(GSM\d+)\"', urlGetText)))
        counter += 1
        if counter % 1000 == 0:
            print('1k studies done')

    pd.DataFrame(sampleList).to_csv('{0}/GEO_MusmusculusSamples.csv'.format(directory))

if __name__ == "__main__":
    main(geoURL, refDirectory)