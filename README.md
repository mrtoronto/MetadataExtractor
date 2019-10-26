# Metadata Extractor

## Description
The NIH's gene expression omnibus (GEO) is a large repository for gene expression data. While the metadata for these experiments is largely structured, the entries are human generated, and thus prone to erroneous or inconsistent syntax. These inconsistent naming mechanisms may apply to gene names, treatments, diet summary, and animal/human age. An initial step in this project is to extract the ages of the animals/humans in these studies, using regular expressions, contextualization, shallow parsing, etc... This small piece of metadata extraction will be part of a larger aim to programmatically extract experimental metadata within and across such databases

### Objectives
- Query study-specific text regarding protocols and experimental designs
- Convert numeric and textual descriptors of age, treatment durations, and other time-based metrics to uniform and numerical metrics of days, weeks, or months
- Generate distribution of animal ages, organism-specific