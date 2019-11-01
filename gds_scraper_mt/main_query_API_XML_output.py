from src.scrape_gds import scrape_gds
import pandas as pd
import random

def main(sampleIDs, organism, out_path, multichannel_flag, keep_files_list, out_types_list):

    scrape_dict = scrape_gds(query_terms = sampleIDs,
                    api_key = "",
                    multichannel = multichannel_flag, ### In progress
                    DEBUG=1,
                    out_path = out_path,
                    keep_files = keep_files_list,
                    out_types = out_types_list
                    )
    print(f'Scrape successful. File saved to {out_path}')
    return scrape_dict


#######################################################

### SampleIDs parameters
num_samples = 200
### Set to False for random sample, else pull `num_samples` samples in a row starting at a random sample.
samples_in_order = False

sampleIDs = pd.read_csv('data/GEO_MusmusculusSamples.csv', index_col=[0], names=['sample_ID'], header=0)

if samples_in_order == True:
    rand_num = random.randint(0, len(sampleIDs) - num_samples)
    sampleIDs = list(sampleIDs['sample_ID'][rand_num:rand_num+num_samples].unique())
else:
    sampleIDs = list(sampleIDs['sample_ID'].sample(num_samples))

organism = 'Mus musculus'
out_path = 'output/test_200'
multichannel_flag = True
keep_files_list = ['xml', 'txt']
out_types_list = ['json', 'csv']


if __name__ == "__main__":

    main(sampleIDs, organism, out_path, multichannel_flag, keep_files_list, out_types_list)
