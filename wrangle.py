""" Main wrangle function for gene lists, tissues, organisms, diets, etc.. """

import json, os, re, sys 
sys.path.append('./src')
import setup_wrangler as wrangler
import pandas as pd 

refDirectory = 'refFiles'
organism = 'Mus musculus'

geneFileTag = 'gene_result_'
geneOrgFile = geneFileTag + re.sub(' ', '_', organism).lower()
valsToAlphaNumLower = True
overwrite = True

def main(geneOrgFile, refDirectory, organism, overwrite, valsToAlphaNumLower):

    if geneOrgFile + '.txt' not in os.listdir(refDirectory):
        raise FileNotFoundError('Gene list for {0} not found'.format(organism))

    if geneOrgFile + '.json' not in os.listdir(refDirectory) or overwrite is True:

        df = pd.read_csv('{0}/{1}.txt'.format(refDirectory, geneOrgFile),
            sep = '\t')[['Org_name', 'GeneID', 'Symbol', 'Aliases']].drop_duplicates()
    
        geneDict = wrangler.geneWrangler(df = df, organism = organism, 
            valsToAlphaNumLower = valsToAlphaNumLower)

        with open('{0}/{1}.json'.format(refDirectory, geneOrgFile), 'w') as fout:
            json.dump(geneDict, fout, indent = 4)

if __name__ == "__main__":
    
    main(geneOrgFile, refDirectory, organism, overwrite, valsToAlphaNumLower)