from src.scrape_gds import scrape_gds_to_csv
import pandas as pd

num_samples = 10

query_terms = pd.read_csv('data/GEO_MusmusculusSamples0.csv', index_col=[0], names=['sample_ID'], header=0)
query_terms = query_terms.sample(num_samples)['sample_ID'].unique()

csv_name = scrape_gds_to_csv(query_terms = query_terms, api_key = "", DEBUG=1)

print(csv_name)
