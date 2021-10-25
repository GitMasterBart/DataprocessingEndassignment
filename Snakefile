SAMPLES = ['KO1','KO2','KO3','WT']


rule all:
    input:
        'results/graph.pdf'

       # expand('vcf_inputs/{sample}.q.d.vcf.gz.tbi',sample=SAMPLES)

rule twoBitToFa:
    input: 'bit_files/mm10.2bit'
    output: 'fasta_files/mm10.fa'
    message: 'bit files is transformed to fa file!!'
    shell: 'bit_files/./twoBitToFa {input} {output}'

rule IndexFastaFile:
    input: 'fasta_files/mm10.fa'
    output: 'fasta_files/mm10.fa.fai'
    shell: 'samtools faidx {input}'

rule sortBam:
    input: 'bam_files/{sample}.bam'
    output: 'sorted_read/{sample}.bam'
    shell:
        'samtools sort {input} > {output} --threads 6'

rule Callvarians:
    input:
        refFasta = 'fasta_files/mm10.fa',
        reffastaindex ='fasta_files/mm10.fa.fai',
        bamfiles= 'sorted_read/{sample}.bam',
    output: 'vcf_inputs/{sample}.raw.bcf'
    shell:
        'bcftools mpileup -f {input.refFasta} {input.bamfiles} | bcftools call -mv -Ou -o {output}'

rule bcftovcf:
    input:'vcf_inputs/{sample}.raw.bcf'
    output:'vcf_inputs/{sample}.vcf'
    shell:
         'bcftools view {input} > {output}'

rule vcf_filterQuality:
    input: 'vcf_inputs/{sample}.vcf'
    output: 'vcf_inputs/{sample}.q.vcf'
    shell: 'python3 scripts/vcffilter.py Quality -i {input} -q 30 -o {output}'
#
rule vcf_filterdepth:
    input: 'vcf_inputs/{sample}.q.vcf'
    output: 'vcf_inputs/{sample}.q.d.vcf'
    shell: 'python3 scripts/vcffilter.py Depth -i {input} -d 10 -o {output}'
#
rule bzip:
    input: 'vcf_inputs/{sample}.q.d.vcf'
    output: 'vcf_inputs/{sample}.q.d.vcf.gz'
    shell : 'bgzip -i {input} '

rule indexx:
    input: 'vcf_inputs/{sample}.q.d.vcf.gz'
    output: 'vcf_inputs/{sample}.q.d.vcf.gz.tbi'
    shell: 'tabix {input}'

rule MergeKOfiles:
    input:
        file = 'vcf_inputs/KO1.q.d.vcf.gz',
        file2 = 'vcf_inputs/KO2.q.d.vcf.gz',
        file3 = 'vcf_inputs/KO3.q.d.vcf.gz',
        id_file1 = 'vcf_inputs/KO1.q.d.vcf.gz.tbi',
        id_file2 = 'vcf_inputs/KO2.q.d.vcf.gz.tbi',
        id_file3 = 'vcf_inputs/KO3.q.d.vcf.gz.tbi'
    output : 'vcf_inputs/MergedKO.qd.vcf'
    shell: 'bcftools merge --force-samples {input.file} {input.file2} {input.file3} > {output}'

rule make_histogram:
    """ rule that creates histogram from gene expression counts"""
    input:
         WTfile = 'vcf_inputs/WT.vcf',
         KOfile ='vcf_inputs/MergedKO.qd.vcf'
    output:
         'results/graph.pdf'
    shell:
        "Rscript scripts/plot.R {input.WTfile} {input.KOfile} {output}"


rule cleanall:
    shell: 'rm -r dag1.png'

# snakemake --forceall --dag | dot -Tpng > dag1.png
#           snakemake --dag vcf_files/{KO1,KO2,KO3,WT}.qdd.vcf.gz.tbi | dot -Tsvg > dag.svg

