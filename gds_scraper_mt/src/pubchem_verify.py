import requests, urllib
import xml.etree.ElementTree as ET


def pubChem_verify(chemical_name):
    """
    Function will take in a chemical name, ping it to pubchem and then return a boolean based on the response page.
    This will return True for a partial match to a valid chemical which may cause issues.
        Args:
            chemical_name - Str: Name of chemical to be sent to PubChem

        Returns:
            Boolean: Value based on whether the chemical is valid according to PubChem
    """
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/xml'
    root = ET.fromstring(requests.get(url).text)

    if 'Fault' in root.tag:
        return False
    else:
        print(url)
        return True
