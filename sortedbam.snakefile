"""
Sort al .bam files
"""

rule sorted_bamFiles:
    input: 'bam_files/{sample}.bam'
    output: 'sorted_read/{sample}.bam'
    log: 'log_files/bam_files/{sample}.log'
    threads: 8
    message: 'All bam files are sorted! with {threads}, to the file {output}'
    shell:
        'samtools sort {input} > {output} --threads {threads} 2> {log}'