SAMPLES = ['KO1','KO2','KO3','WT']


rule all:
    input:
        expand('vcf_inputs/{sample}.qdd.vcf.gz.tbi', sample = SAMPLES)

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
        bamfiles= 'sorted_read/{sample}.bam',
    output: 'vcf_inputs/{sample}.raw.bcf'
    shell:
        'bcftools mpileup -f {input.refFasta} {input.bamfiles} | bcftools call -mv -Ou -o {output}'

rule bcftovcf:
    input:'vcf_inputs/{sample}.raw.bcf'
    output:'vcf_inputs/{sample}.vcf'
    shell:
         'bcftools view {input}> {output}'

rule vcf_filterQuality:
    input: 'vcf_inputs/{sample}.vcf'
    output: 'vcf_inputs/{sample}.q.vcf'
    shell: 'python3 scripts/vcffilter.py Quality -i {input} -q 30 -o {output}'

rule vcf_filterdepth:
    input: 'vcf_inputs/{sample}.q.vcf'
    output: 'vcf_inputs/{sample}.qd.vcf'
    shell: 'python3 scripts/vcffilter.py Depth -i {input} -d 10 -o {output}'

rule vcf_filter_decompose:
    input: 'vcf_inputs/{sample}.qd.vcf'
    output: 'vcf_inputs/{sample}.qdd.vcf'
    shell: 'python3 scripts/vcffilter.py Decompose -i {input} -o {output}'

rule bzip:
    input: 'vcf_inputs/{sample}.qdd.vcf'
    output: 'vcf_inputs/{sample}.qdd.vcf.gz'
    shell : 'bgzip -i {input} '

rule indexx:
    input: 'vcf_inputs/{sample}.qdd.vcf.gz'
    output: 'vcf_inputs/{sample}.qdd.vcf.gz.tbi'
    shell: 'tabix {input}'

rule MergeKOfiles:
    output : 'MergedKO.vcf'
    shell: 'bcftools merge --force-samples vcf_inputs/KO1.qdd.vcf.gz vcf_inputs/KO2.qdd.vcf.gz vcf_inputs/KO3.qdd.vcf.gz > {output}'


rule cleanall:
    shell: 'rm -r dag1.png'

# snakemake --forceall --dag | dot -Tpng > dag1.png
#           snakemake --dag vcf_files/{KO1,KO2,KO3,WT}.qdd.vcf.gz.tbi | dot -Tsvg > dag.svg

