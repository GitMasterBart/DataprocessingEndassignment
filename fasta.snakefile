"""
All rules create fastq files are in this file.
#twoBitToFa
#index_fastqFile
"""

rule twoBitToFa:
    input: config["2bitfile"]
    output: 'fasta_files/mm10.fa'
    log: 'log_files/fasta_files/mm10.log'
    message: '{input} is transformed to {output}!'
    shell: 'bit_files/./twoBitToFa {input} {output} 2> {log}'

rule index_fastqFile:
    input: 'fasta_files/mm10.fa'
    output: 'fasta_files/mm10.fa.fai'
    log: 'log_files/fasta_files/mm10.fai.log'
    message: 'fastq is indexed'
    shell: 'samtools faidx {input} 2> {log}'