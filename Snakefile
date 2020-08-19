import os
import glob
import re

workdir: os.getcwd()


# Query the final the trees and the benchmark
rule final:
    input: 
        "output/benchmark/fig1.png",
        expand("output/graph/ologram_result_tree_{cell_line}.pdf", cell_line = ['mcf7', 'mcf7_filtered','mcf7_manual','artificial'])


# ---------------------------------------------------------------------------- #
#                                 Data retrieval                               #
# ---------------------------------------------------------------------------- #

rule prepare_incl:
    """
    Restrict the shuffling to "pseudo-CRMs", obtained by merging all TF peaks.
    Then move the query (FOXA1 for MCF7) to prevent them from being used by MODL.
    """
    output:
        query = "input/foxa1_mcf7.bed",
        incl = "input/crm_mcf7.bed"

    shell: """
        # Unzip all
        gunzip input/mcf7/*.gz

        # Create incl file by merging all BEDs
        cat input/mcf7/*.bed | bedtools sort | bedtools merge > {output.incl}

        # Once done, move future queries
        mv input/mcf7/foxa1.bed {output.query}
    """


rule prepare_artificial:
    """
    Prepare artificial demonstration data.
    It consits of random regions for the query, compared against (a) the same regions as the query, (b) half of the query, (c) the other half of the query and (d) random regions.
    Shuffling is restricted to a slop of 2x query size around the query itself.
    """
    output:
        query = "output/artificial_data/query.bed", incl = "output/artificial_data/incl.bed",
        a =  "output/artificial_data/same.bed", b =  "output/artificial_data/half.bed", c =  "output/artificial_data/other_half.bed", d =  "output/artificial_data/neg_control.bed",

    params:
        size=1000, length=20000
    
    shell: """
    mkdir -p output/artificial_data

    # Generate A
    bedtools random -n $((2*{params.size})) -l $(({params.length})) -seed 123456 -g input/hg38.genome > {output.query}
    cp {output.query} {output.a}

    # Generate B and C
    head -n $(({params.size})) {output.a} > {output.b}
    tail -n $(({params.size})) {output.b} > {output.c}

    # Generate incl
    # so files cover roughly 1/3 of shuffling region
    bedtools slop -b $(({params.size})) -i {output.query} -g input/hg38.genome > {output.incl}

    # Generate D
    bedtools shuffle -i {output.query} -g input/hg38.genome -incl {output.incl} -noOverlapping -seed 654321 > {output.d}
    """



# ---------------------------------------------------------------------------- #
#                                 Data processing                              #
# ---------------------------------------------------------------------------- #

# Get the list of files in the various input directory
def get_peaks_mcf7(wildcards):
    file_list = glob.glob('input/mcf7/*bed')
    return " ".join(file_list)

def get_peaks_artificial(wildcards):
    file_list = glob.glob('output/artificial_data/*bed')
    return " ".join(file_list)

def how_many_peaks_mcf7(wildcards):
    file_list = glob.glob('input/mcf7/*bed')
    return len(file_list)

rule compute_combi_enrichment_mcf7:
    """
    For a the MCF7 cell line, compute the enrichment in n-wise TF combinations using OLOGRAM-MODL.
    The query is FOXA1.
    """
    input: 
        query = 'input/foxa1_mcf7.bed',
        incl = "input/crm_mcf7.bed"
    params:
        trs = get_peaks_mcf7,
        minibatch_number = 10, minibatch_size = 5, threads = 8,

    output: 'output/ologram_result_mcf7/00_ologram_stats.tsv', 

    shell: """
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.trs} \
        -o output/ologram_result_mcf7 --force-chrom-peak --force-chrom-more-bed  \
        -V 3 -k {params.threads} -mn {params.minibatch_number} -ms {params.minibatch_size} \
        --more-bed-multiple-overlap --bed-incl {input.incl} --no-date
    """


rule compute_mcf7_modl_selection:
    """
    Run MCF7 on a single shuffle.
    This is NOT designed to get results, only to show which combinations will be selected by MODL.
    Indeed, MODL works only on the matrix of true overlaps and does not care for the shuffles.
    """
    input: 
        query = 'input/foxa1_mcf7.bed',
        incl = "input/crm_mcf7.bed"

    params:
        trs = get_peaks_mcf7,
        minibatch_number = 1, minibatch_size = 1, threads = 4,
        max_combis = 20
    
    output:'output/ologram_result_mcf7_filtered/00_ologram_stats.tsv'

    shell: """
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.trs}\
        -o output/ologram_result_mcf7_filtered --force-chrom-peak --force-chrom-more-bed --no-date \
        -k {params.threads} -mn {params.minibatch_number} -ms {params.minibatch_size} -V 3 \
        --more-bed-multiple-overlap --bed-incl {input.incl} \
        --multiple-overlap-max-number-of-combinations {params.max_combis} \
        #--multiple-overlap-target-combi-size XX                 
    """


rule mcf7_manual_filtering:
    """
    Manual filtering of the displayed combinations for MCF7
    """
    input:
        res = 'output/ologram_result_mcf7/00_ologram_stats.tsv',
        filtered = 'input/desired_combis.txt'
    output:
        'output/ologram_result_mcf7_manual/00_ologram_stats.tsv'
    shell: """
    head -n 1 {input.res} > {output}                    # Keep the header
    grep -w -F -f {input.filtered} {input.res} >> {output}
    """

# ---------------------------------------------------------------------------- #
#                              Benchmarking                                    #
# ---------------------------------------------------------------------------- #

rule produce_modl_comparison:
    """
    Benchmarking MODL on an artificial data test case (AB,ABCD,EF, and poisson noise)
    and comparing it to Apriori.
    """
    output: 
        fig1 = "output/benchmark/fig1.png",
        fig2 = "output/benchmark/fig2.png",
        fig3 = "output/benchmark/fig3.png"

    log:
        err= "output/benchmark/comparison_benchmark_ERROR_LOG.txt",
        out= "output/benchmark/comparison_benchmark.txt"

    shell:"""
    mkdir -p output/benchmark



    # python scripts/modl_comparison.py 2> {log.err} 1> {log.out}


    touch {output.fig1}
    touch {output.fig2}
    touch {output.fig3}
    """

rule run_artificial:
    """
    Run OLOGRAM-MODL on artificial data
    """
    input:
        query = "output/artificial_data/query.bed", incl = "output/artificial_data/incl.bed",
        a =  "output/artificial_data/same.bed", b =  "output/artificial_data/half.bed", c =  "output/artificial_data/other_half.bed", d =  "output/artificial_data/neg_control.bed"

    output:
        'output/ologram_result_artificial/00_ologram_stats.tsv'

    params:
        minibatch_number = 20, minibatch_size = 10, threads = 4,
        peaks = get_peaks_artificial

    shell:"""
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.peaks} --bed-incl {input.incl} \
        -o output/ologram_result_artificial --force-chrom-peak --force-chrom-more-bed \
        -V 3 -k {params.threads} -mn {params.minibatch_number} -ms {params.minibatch_size} \
        --more-bed-multiple-overlap --no-date -K output/ologram_result_artificial_TEMP
    """

# ---------------------------------------------------------------------------- #
#                                 Visual                                       #
# ---------------------------------------------------------------------------- #

# TODO Hardcoded for now. Of course, this is no longer apropos when using different queries and cell lines.
def get_queryname(wildcards):
    run = wildcards.cell_line
    if "mcf7" in run : return "foxa1"
    if run == "artificial": return "Query"
    return " ".join(file_list)


rule treeify:
    """
    For visual purposes.
    Turns an OLOGRAM-MODL result file into a graph of the found combinations.
    """
    input: "output/ologram_result_{cell_line}/00_ologram_stats.tsv"
    output: "output/graph/ologram_result_tree_{cell_line}.pdf"
    params:
        queryname = get_queryname
    shell:"""
    gtftk ologram_modl_treeify -i {input} -o {output} -l {params.queryname}
    """
