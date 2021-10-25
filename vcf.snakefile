"""
Cals and filters all variants and merge it to one file called mergedKO.qd.vcf
#call_Variants
#bcfToVcf
#vcf_filterQuality
#vcf_filterDepth
#zip_vcf
#index_VcfFile
#mergeKOFiles
"""

rule call_Variants:
    input:
        refFasta = 'fasta_files/mm10.fa',
        reffastaindex ='fasta_files/mm10.fa.fai',
        bamfiles= 'sorted_read/{sample}.bam',
    output: 'vcf_inputs/{sample}.raw.bcf'
    log: 'log_files/vcf_files/{sample}.raw.bcf.log'
    threads: 8
    message: 'Variants are called with {threads} and saved in {output}'
    shell:
        'bcftools mpileup --threads {threads} -f {input.refFasta} {input.bamfiles} | bcftools call -mv -Ou -o {output} 2> log'

rule bcfToVcf:
    input:'vcf_inputs/{sample}.raw.bcf'
    output:'vcf_inputs/{sample}.vcf'
    log: 'log_files/vcf_files/{sample}vcf.log'
    message: 'bcf is transferred tot vcf '
    shell:
         'bcftools view {input} > {output} 2> {log}'

rule vcf_filterQuality:
    input: 'vcf_inputs/{sample}.vcf'
    output: 'vcf_inputs/{sample}.q.vcf'
    log: "log_files/vcf_files/{sample}.q.d.vcf.log"
    message: 'Vcf file is filtered on Quality'
    shell: 'python3 scripts/vcffilter.py Quality -i {input} -q 30 -o {output} 2> {log}'
#
rule vcf_filterDepth:
    input: 'vcf_inputs/{sample}.q.vcf'
    output: 'vcf_inputs/{sample}.q.d.vcf'
    log: 'log_files/vcf_files/{sample}.q.d.vcf.log'
    message: 'Vcf file is filtered on Depth'
    shell: 'python3 scripts/vcffilter.py Depth -i {input} -d 10 -o {output} 2> {log}'
#
rule zip_vcf:
    input: 'vcf_inputs/{sample}.q.d.vcf'
    output: 'vcf_inputs/{sample}.q.d.vcf.gz'
    log: 'log_files/vcf_files/{sample}.q.d.vcf.gz.log'
    message: 'vcf file is zipped and ready almost for merge'
    shell : 'bgzip -i {input} 2> {log}'

rule index_VcfFile:
    input: 'vcf_inputs/{sample}.q.d.vcf.gz'
    output: 'vcf_inputs/{sample}.q.d.vcf.gz.tbi'
    log: 'log_files/vcf_files/{sample}.q.d.vcf.gz.tbi.log'
    message: 'vcf file is indexed and ready for merge'
    shell: 'tabix {input} 2> {log}'

rule mergeKOFiles:
    input:
        file = 'vcf_inputs/KO1.q.d.vcf.gz',
        file2 = 'vcf_inputs/KO2.q.d.vcf.gz',
        file3 = 'vcf_inputs/KO3.q.d.vcf.gz',
        id_file1 = 'vcf_inputs/KO1.q.d.vcf.gz.tbi',
        id_file2 = 'vcf_inputs/KO2.q.d.vcf.gz.tbi',
        id_file3 = 'vcf_inputs/KO3.q.d.vcf.gz.tbi'
    threads: 8
    output : 'vcf_inputs/MergedKO.qd.vcf'
    log: 'log_files/vcf_files/mergedKO.log'
    message: 'Vcf KO files are merged to {output}'
    shell: 'bcftools merge --force-samples {input.file} {input.file2} {input.file3} > {output}'



