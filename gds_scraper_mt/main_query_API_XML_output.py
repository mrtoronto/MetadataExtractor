import random
import pandas as pd
from src.scrape_gds import scrape_gds

def main(sampleIDs, organism, out_path, multichannel_flag, keep_files_list, out_types_list):
    """
    This function will take in a list of sampleIDs and create the specified output files based on the sample's data and associated meta-data.

        Args:
            `query_terms` - List: List of file names generated from the `data/GEO_MusmusculusSamples.csv` data file. Samples can be selected randomly or in order starting from a random spot.
            `api_key` - Str: NCBI e-utils API key
            `multichannel` - Boolean: Flag to filter out multichannel samples
            `DEBUG` - Int: Number to determine the depth of debug print output (Deprecated).
            `out_path` - Str: Filename of the output file. Excludes extensions as that'll be added depending on the type of the output.
            `keep_files` - List: Intermediate files to keep. Options currently are 'xml' or 'txt'. Both may be selected at once.
            `out_types` - List: File type of the output. Options currently include 'json' and 'csv'. Both may be selected at once.
    """
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
### Set to False for random sample, else pull `num_samples` samples consecutively starting at a random sample.
samples_in_order = False

sampleIDs = pd.read_csv('data/GEO_MusmusculusSamples.csv', index_col=[0], names=['sample_ID'], header=0)

if samples_in_order == True:
    rand_num = random.randint(0, len(sampleIDs) - num_samples)
    sampleIDs = list(sampleIDs['sample_ID'][rand_num:rand_num+num_samples].unique())
else:
    sampleIDs = list(sampleIDs['sample_ID'].sample(num_samples))

organism = 'Mus musculus'
out_path = 'output/200_sample_test_run'
multichannel_flag = True
keep_files_list = ['xml', 'txt']
out_types_list = ['json', 'csv']


if __name__ == "__main__":

    main(sampleIDs, organism, out_path, multichannel_flag, keep_files_list, out_types_list)
