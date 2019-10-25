rule all:
	input:
        "refFiles/GEO_MusmusculusStudies.txt",
        "Figures/MusMusculusAgeHistogram.png"

rule sampleListGen:
    input:
        "refFiles/GEO_MusmusculusStudies.txt"
    output:
        "refFiles/GEO_MusmusculusSamples.csv"
    shell:
        "python sampleListGen.py {input} > {output}"

rule metaExtract:
    input:
        "src/setup_metadataExtract.py",
        "refFiles/GEO_MusmusculusSamples.csv"
    output:
        "refFiles/GEO_MusmusculusMetadata.json"
    shell:
        "python metadataExtract.py {input} > {output}"

rule ageHisto:
    input:
        "refFiles/GEO_MusmusculusMetadata.json"
	output:
        "Figures/MusMusculusAgeHistogram.png"
    shell:
		"Rscript studyAgeViz.R {input} > {output}"