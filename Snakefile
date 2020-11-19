import os
import glob
import re

workdir: os.getcwd()


# Query the final trees and the benchmark
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
    It consits of random regions for the query, compared against (a) a third of 
    the query, (b) another third, and (c) a negative control.
    """
    output:
        query = "output/artificial_data/query.bed",
        a = "output/artificial_data/data/third.bed", a2 = "output/artificial_data/data/third_bis.bed", b = "output/artificial_data/data/other_third.bed", c = "output/artificial_data/data/neg_control.bed"

    params:
        size=1000, length=200000
    
    shell: """
    mkdir -p output/artificial_data

    # Generate query
    bedtools random -n $((3*{params.size})) -l $(({params.length})) -seed 123456 -g input/hg38.genome > {output.query}

    # Generate A and B
    head -n $(({params.size})) {output.query} > {output.a}
    bedtools shuffle -i {output.a} -incl {output.query} -excl {output.a} -g input/hg38.genome -noOverlapping -seed 654321 > {output.b}
    # Do not put all regions of query into A and B

    # Ensure there is at least one common line in A and B and C so they are displayed
    echo "chr1	1	2	0	1	-" >> {output.a}
    echo "chr1	1	2	0	1	-" >> {output.b}
    echo "chr1	1	2	0	1	-" >> {output.query}

    cp {output.a} {output.a2}

    # Generate D
    bedtools random -n $((2*{params.size})) -l $(({params.length})) -seed 654321 -g input/hg38.genome > {output.c}
    """



# ---------------------------------------------------------------------------- #
#                                 Data processing                              #
# ---------------------------------------------------------------------------- #

# Get the list of files in the various input directory
def get_peaks_mcf7(wildcards):
    file_list = sorted(glob.glob('input/mcf7/*bed'))
    return " ".join(file_list)

def get_peaks_artificial(wildcards):
    file_list = sorted(glob.glob('output/artificial_data/data/*bed'))
    return " ".join(file_list)

def how_many_peaks_mcf7(wildcards):
    file_list = sorted(glob.glob('input/mcf7/*bed'))
    return len(file_list)

rule compute_combi_enrichment_mcf7:
    """
    For a the MCF7 cell line, compute the enrichment in n-wise TF combinations 
    using OLOGRAM-MODL. The query is FOXA1.
    """
    input: 
        query = 'input/foxa1_mcf7.bed',
        incl = "input/crm_mcf7.bed"
    params:
        trs = get_peaks_mcf7,
        minibatch_number = 10, minibatch_size = 5#, threads = 8,
    threads: 8

    output: 'output/ologram_result_mcf7/00_ologram_stats.tsv', 

    shell: """
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.trs} \
        -o output/ologram_result_mcf7 --force-chrom-peak --force-chrom-more-bed  \
        -V 3 -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} \
        --more-bed-multiple-overlap --bed-incl {input.incl} --no-date
    """


rule compute_mcf7_modl_selection:
    """
    Run MCF7 on a single shuffle.
    This is NOT designed to get results, only to show which combinations will be
    selected by MODL. Indeed, MODL works only on the matrix of true overlaps and
    does not care for the shuffles.
    """
    input: 
        query = 'input/foxa1_mcf7.bed',
        incl = "input/crm_mcf7.bed"

    params:
        trs = get_peaks_mcf7,
        minibatch_number = 1, minibatch_size = 1,#, threads = 4,
        max_combis = 20
    threads: 4

    output:'output/ologram_result_mcf7_filtered/00_ologram_stats.tsv'

    shell: """
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.trs}\
        -o output/ologram_result_mcf7_filtered --force-chrom-peak --force-chrom-more-bed --no-date \
        -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} -V 3 \
        --more-bed-multiple-overlap --bed-incl {input.incl} \
        --multiple-overlap-max-number-of-combinations {params.max_combis}               
    """


rule mcf7_manual_filtering:
    """
    Manual filtering of the displayed combinations for MCF7.
    Remember to keep the header.
    """
    input:
        res = 'output/ologram_result_mcf7/00_ologram_stats.tsv',
        filtered = 'input/desired_combis.txt'
    output:
        'output/ologram_result_mcf7_manual/00_ologram_stats.tsv'
    shell: """
    head -n 1 {input.res} > {output}
    grep -w -f {input.filtered} {input.res} >> {output}
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
        fig3 = "output/benchmark/fig3.png",
        fig4 = "output/benchmark/fig4.png"

    log:
        err= "output/benchmark/comparison_benchmark_ERROR_LOG.txt",
        out= "output/benchmark/comparison_benchmark.txt"

    shell:"""
    mkdir -p output/benchmark
    python scripts/modl_comparison.py 2> {log.err} 1> {log.out}
    """


rule run_artificial:
    """
    Run OLOGRAM-MODL on artificial data
    """
    input:
        query = "output/artificial_data/query.bed"

    output:
        'output/ologram_result_artificial/00_ologram_stats.tsv'

    params:
        minibatch_number = 20, minibatch_size = 10, threads = 4,
        peaks = get_peaks_artificial

    shell:"""
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.peaks} \
        -o output/ologram_result_artificial --force-chrom-peak --force-chrom-more-bed \
        -V 3 -k {params.threads} -mn {params.minibatch_number} -ms {params.minibatch_size} \
        --more-bed-multiple-overlap --no-date -K output/ologram_result_artificial_TEMP
    """

# ---------------------------------------------------------------------------- #
#                                 Visual                                       #
# ---------------------------------------------------------------------------- #

# TODO Hardcoded for now. Of course, this is no longer apropos when using 
# different queries and cell lines.
def get_queryname(wildcards):
    run = wildcards.cell_line
    if "mcf7" in run : return "foxa1"
    if run == "artificial": return "Query"


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
