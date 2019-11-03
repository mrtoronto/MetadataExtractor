""" Utility functions for aggregating gene lists, diets, tissues, etc

    Gene lists Downloaded off https://www.ncbi.nlm.nih.gov/gene/?term=(organism) 
    "Send to:" --> File --> "Create File"

    Basic aggregation of gene IDs. More formality and programmatic queries would
    be nice, but this kikely does not need to be run very often anyways, and thus
    a fixed df structure and column is assumed from the 
"""
import json, os, re
import pandas as pd 

def geneWrangler(df, organism, valsToAlphaNumLower = True):
    """ Wrangle gene .txt dataframes from www.ncbi.nim.nih.gov. Return dataframe
        converted to a dict with gene IDs as keys, and symbols and list of 
        aliases as values """

    if valsToAlphaNumLower is True:
        for col in df.columns:
            df[col] = df[col].astype(str).replace(r'[^0-9a-zA-Z\,\ +]+', '', 
                regex = True).str.lower()
        
        organism = organism.lower()

    geneDict = dict()
    for i, j, k, l in zip(df['Org_name'], df['GeneID'], df['Symbol'], df['Aliases']):
        if i == organism:     
            if j not in geneDict:
                geneDict[j] = dict() 
                geneDict[j]['Aliases'] = []
            
            geneDict[j]['Symbol'] = k 
            alias = [re.sub(r'^ | $', '', x) for x in re.split(r'\,', l)]
            geneDict[j]['Aliases'] += [x for x in alias if x != '']

    return geneDict

