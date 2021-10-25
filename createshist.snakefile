"""
This creat the graph with histograms
"""
rule make_histograms:
    """ rule that creates histogram from gene expression counts"""
    input:
         WTfile = 'vcf_inputs/WT.q.d.vcf',
         KOfile ='vcf_inputs/MergedKO.qd.vcf'
    output:
         'results/histograms.pdf'
    log: 'log_files/results/histograms.log'
    message: 'Graph is created and added to results file'
    shell:
        "Rscript scripts/plot.R {input.WTfile} {input.KOfile} {output} 2> {log}"