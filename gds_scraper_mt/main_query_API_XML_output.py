from src.scrape_gds import scrape_gds_to_csv
import pandas as pd
import random

num_samples = 20
### Set to false for `num_samples` random samples, else pull `num_samples` samples in a row
samples_in_order = False


query_terms = pd.read_csv('data/GEO_MusmusculusSamples0.csv', index_col=[0], names=['sample_ID'], header=0)

if samples_in_order == True:
    rand_num = random.randint(0, len(query_terms) - num_samples)
    query_terms = list(query_terms['sample_ID'][rand_num:rand_num+num_samples].unique())
else:
    query_terms = list(query_terms['sample_ID'].sample(num_samples))

#query_terms = query_terms + ['GSM1394796']
out_path = 'output/test_small.csv'

csv_name = scrape_gds_to_csv(query_terms = query_terms,
                                api_key = "",
                                multichannel = False, ### In progress
                                DEBUG=1,
                                out_path = out_path,
                                keep_files = [None],
                                )

print(out_path)
