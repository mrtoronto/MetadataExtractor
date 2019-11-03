""" Utility functions for aggregating gene lists, diets, tissues, etc

    Gene lists Downloaded off https://www.ncbi.nlm.nih.gov/gene/?term=(organism) 
    "Send to:" --> File --> "Create File"

    Basic aggregation of gene IDs. More formality and programmatic queries would
    be nice, but this kikely does not need to be run very often anyways, and thus
    a fixed df structure and column is assumed from the 
"""
import json, os, re, time, requests
import pandas as pd 
from bs4 import BeautifulSoup
defaultAgent = {'User-Agent': 'SomeAgent 11.0'}

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


def mouseDietWrangler(baseUrl = 'https://www.researchdiets.com',
    dietLinks = ['/opensource-diets/in-stock-diets', 
    '/opensource-diets/dio-series-diets']):
    """ Extract macronutrient informatin from diets off researchdiets.com. These
        diets are most commonly used in research. Macro numbers are in percentages
        of total calories. """

    dietDict = dict()

    for link in dietLinks:
        urlGet = requests.get(baseUrl+link, headers = defaultAgent)
        soup = BeautifulSoup(urlGet.text, features = 'html.parser')
        for subLink in soup.find_all('a', href=True):
            if 'formula' in subLink['href']:
                dietID = re.sub('\/formulas\/', '', subLink['href'])
                dietDict[dietID] = dict()
                dietDict[dietID]['url'] = baseUrl + subLink['href']

        for dietID in dietDict:
            time.sleep(0.25)
            urlGetText = requests.get(dietDict[dietID]['url'], 
                headers = defaultAgent).text
            
            prot = re.findall('Protein\:.*\n(.*)\n', urlGetText)[0]
            fat = re.findall('Fat\:.*\n(.*)\n', urlGetText)[0]
            carb = re.findall('Carbohydrate\:.*\n(.*)\n', urlGetText)[0]
            density = re.findall('Energy Density\:.*\n(.*)\n', urlGetText)
            if len(density) == 0:
                density = 'n/a'
            else:
                density = density[0] 
            
            dietDict[dietID]['Protein'] = re.findall('\>(\d+)\<', prot)[0]
            dietDict[dietID]['Fat'] = re.findall('\>(\d+)\<', fat)[0]
            dietDict[dietID]['Carbohydrate'] = re.findall('\>(\d+)\<', carb)[0]

            if density != 'n/a':
                dietDict[dietID]['Energy Density'] = re.findall('\>(\d+\.?\d+?)\<', density)[0]
            else:
                dietDict[dietID]['Energy Density'] = 'n/a'
    
    return dietDict 
