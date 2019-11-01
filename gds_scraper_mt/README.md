# GEO_scraper

Copying and updating what I wrote for the Pubmed scraper.

## Program Goal

Allow user to scrap data from GDS samples in an organized and efficient way.

## How to Run

- Run `main_query_API_XML_output.py` from the terminal to start the process. Variables to control results can be edited within the script.


## Output File

- The output will have a filename like `output/test_output.csv` or `output/test_output.json` and can be edited from the `main_query_API_XML_output.py` script.


## Program Flow

Scripts run in the following order:

1. `main_query_API_XML_output.py`
2. `scrape_gds.py`
3. Running once per query, `src/search_sample.py`
4. Running once on a list of filenames, `src/geo_parser.py`
5. Running once for each unique series FTP link found in the above parser's output, `src/ftp_gzxml_parser.py`
6. Data organizing.

### `Individual Scripts`

I'd like to do little writeups about each of the scripts like I did for the PubMed parser but I haven't gotten a chance yet.
