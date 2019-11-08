#/bin/bash
python3.6 -m cProfile -o metadataSpeed.prof metadataExtract.py
snakeviz metadataSpeed.prof
