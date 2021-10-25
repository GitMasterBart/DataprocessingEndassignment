"""
Main script.
Includes all necessary Snakefiles needed
to execute the pipeline.
"""
configfile: "config.yaml"

include: "fasta.snakefile"
include: "sortedbam.snakefile"
include: 'vcf.snakefile'
include: "createshist.snakefile"

rule all:
    input:
        'results/histograms.pdf'
onsuccess:
    print("The pipeline had a successful walk trough, no error")
onerror:
    print("An error occurred, please check the log files.")


