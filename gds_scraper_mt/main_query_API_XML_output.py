import numpy as np
import random
import pandas as pd
from src.scrape_gds import scrape_gds

def main(sampleIDs, organism, out_path, multichannel_flag, cells_flag, metadata_filter, keep_files_list, out_types_list):
    """
    This function will take in a list of sampleIDs and create the specified output files based on the sample's data and associated meta-data.

        Args:
            `query_terms` - List: List of file names generated from the `data/GEO_MusmusculusSamples.csv` data file. Samples can be selected randomly or in order starting from a random spot.
            `api_key` - Str: NCBI e-utils API key
            `multichannel` - Boolean: Flag to filter out multichannel samples
            `cells_flag` - Boolean: Flag to filter out samples marked as Cells (See `geoSampleCellCheck()`)
            `metadata_filter` - Boolean: Flag to filter out samples that fail to join to any downloaded metadata.
            `DEBUG` - Int: Number to determine the depth of debug print output (Deprecated).
            `out_path` - Str: Filename of the output file. Excludes extensions as that'll be added depending on the type of the output.
            `keep_files` - List: Intermediate files to keep. Options currently are 'xml' or 'txt'. Both may be selected at once.
            `out_types` - List: File type of the output. Options currently include 'json' and 'csv'. Both may be selected at once.
    """
    scrape_dict = scrape_gds(query_terms = sampleIDs,
                    api_key = "",
                    multichannel = multichannel_flag, ### In progress
                    cells_flag = cells_flag,
                    metadata_filter = metadata_filter,
                    DEBUG=1,
                    out_path = out_path,
                    keep_files = keep_files_list,
                    out_types = out_types_list
                    )
    print(f'Scrape successful. File saved to {out_path}')
    return scrape_dict


#######################################################

### SampleIDs parameters
num_samples = 15
### Set to False for random sample, else pull `num_samples` samples consecutively starting at a random sample.
samples_in_order = False

sampleIDs = pd.read_csv('data/GEO_MusmusculusSamples.csv', index_col=[0], names=['sample_ID'], header=0)['sample_ID'].unique()

if samples_in_order == True:
    rand_num = random.randint(0, len(sampleIDs) - num_samples)
    sampleIDs = list(sampleIDs[rand_num:rand_num+num_samples])
else:
    sampleIDs = list(np.random.choice(sampleIDs, num_samples))

samples_with_age_char = ['GSM1369233', 'GSM1255384', 'GSM1369233']
sampleIDs = sampleIDs + samples_with_age_char

organism = 'Mus musculus'
out_path = 'output/15_run1'
multichannel_flag = False
cells_flag = False
metadata_filter = False
keep_files_list = [None]
out_types_list = ['json', 'csv']


if __name__ == "__main__":

    main(sampleIDs, organism, out_path, multichannel_flag, cells_flag, metadata_filter, keep_files_list, out_types_list)
