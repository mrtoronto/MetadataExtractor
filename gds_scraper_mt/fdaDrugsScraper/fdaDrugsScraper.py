def fdaDrugScraper(export=False):
    """
    Scrapes https://www.accessdata.fda.gov for drugs and associated data.
    """
    drugs_dict = {}
    for letter in letters:
        url = f'https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=browseByLetter.page&productLetter={letter}'
        root = lxml.html.document_fromstring(requests.get(url).text)
        drug_list = root.findall('.//tbody/tr/td')
        for drug in drug_list:
            drug_cols = [re.sub('[\r\n\t]', '', i) for i in drug.itertext() if re.sub('[\r\n\t ]', '', i) != '']
            drug_name = drug_cols[0]
            active_ingred = drug_cols[1].split('(')[0]#[:-1]
            anda_number = drug_cols[2].split('|')[1].split('#')[1]
            dosage_form = drug_cols[2].split('|')[2]
            company = drug_cols[2].split('|')[3]

            drug_dict = {'drug_name' : drug_name,
                        'ANDA' : anda_number,
                        'active_ingred' : active_ingred,
                        'dosage_form' : dosage_form,
                        'company' : company
                        }
            drugs_dict[drug_name] = drug_dict
    if export == True:
        with open('fdaDrug_data.json', 'w') as fout:
            json.dump(drugs_dict, fout, indent = 4)
    return drugs_dict


fdaDrugScraper(export=True)
