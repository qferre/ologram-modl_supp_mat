import os
import glob
import re

workdir: os.getcwd()


# Query the final heatmap, the trees and the benchmark
rule final:
    input: 
        "output/supp_fig/comparison_benchmark.txt",
        expand("output/ologram/ologram_result_tree_{cell_line}.pdf", cell_line = ['mcf7', 'mcf7_filtered','artificial'])



# ---------------------------------------------------------------------------- #
#                                 Data retrieval                               #
# ---------------------------------------------------------------------------- #


rule prepare_incl:
    """
    Get H3K27ac ChIP-seq peaks to act as background for FOXA1 in MCF7, where we restrict the shuffling to those regions only.
    Then move the query (FOXA1 for MCF7) to prevent them from being used by MODL.
    """
    output:
        hist = "output/H3K27ac_mcf7.bed",
        query = "input/foxa1_mcf7.bed"

    shell: """
        # Get histone
        wget https://www.encodeproject.org/files/ENCFF784QFH/@@download/ENCFF784QFH.bed.gz -O H3K27ac_mcf7_ENCFF784QFH.bed.gz
        gunzip -f H3K27ac_mcf7_ENCFF784QFH.bed.gz
        cat H3K27ac_mcf7_ENCFF784QFH.bed | perl -ne 'print if(/^chr[0-9XY]+\t/)' > {output.hist}
        rm -f H3K27ac_mcf7_ENCFF784QFH.bed

        # Unzip all
        gunzip input/mcf7/*.gz

        # Once done, move future queries
        mv input/mcf7/foxa1.bed {output.query}
    """


rule prepare_artificial:
    """
    Prepare artificial demonstration data.
    """

    output:
        query = "output/artificial/query.bed", incl = "output/artificial/incl.bed",
        a =  "output/artificial/A.bed", b =  "output/artificial/B.bed", c =  "output/artificial/C.bed", d =  "output/artificial/D.bed",
    
    shell: """
    # TODO To add some noise to that, remove at random a few regions from all those files ? Not needed for proof of concept
    mkdir output; mkdir output/artificial

    # Generate A
    bedtools random -n $((4*SIZE)) -l $((LENGTH)) -seed 123456 -g input/hg38.genome > {output.query}
    cp {output.query} {output.a}

    # Generate B and C
    tail -n $((2*SIZE)) {output.a} > {output.b}
    tail -n $((1*SIZE)) {output.b} > {output.c}

    # Generate C D
    bedtools shuffle -i {output.query} -g input/hg38.genome -excl {output.query} -noOverlapping -seed 265159 > {output.d}

    # Generate incl
    bedtools slop -b 500 -i {output.query} -g input/hg38.genome > {output.incl}
    """



# ---------------------------------------------------------------------------- #
#                                 Data processing                              #
# ---------------------------------------------------------------------------- #

# Get the list of files in the `input/mcf7` directory


def get_peaks_mcf7(wildcards):
    file_list = glob.glob('input/mcf7/*bed')
    return " ".join(file_list)

def how_many_peaks_mcf7(wildcards):
    file_list = glob.glob('input/mcf7/*bed')
    return len(file_list)

# TODO : do not pass the queries (FOXA1 for mcf7 and FOS for K562 !)

rule compute_combi_enrichment_mcf7:
    """
    For a the MCF7 cell line, compute the enrichment in n-wise TF combinations using OLOGRAM-MODL
    The query is FOXA1.
    """
    input: 
        query = 'input/foxa1_mcf7.bed',
        incl = "output/H3K27ac_mcf7.bed"
    params:
        trs = get_peaks_mcf7,
        minibatch_number = 6, minibatch_size = 5, threads = 4,

    output: 'output/ologram_result_mcf7/00_ologram_stats.tsv', 

    shell: """
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.trs} |\
        -o ologram_results_mcf7 --force-chrom-peak --force-chrom-more-bed  |\
        -V 3 -k {params.threads} -mn {params.minibatch_number} -ms {params.minibatch_size} |\
        --more-bed-multiple-overlap --bed-incl {input.incl} --no-date
    """

rule compute_combi_enrichment_mcf7_filtered:
    """
    Run MCF7 on a single shuffle

    This is NOT designed to get results. Only to show which combinations will be selected by MODL.

    Indeed, MODL works only on the matrix of true overlaps and does not care for the shuffles.
    """

    input: 
        query = 'input/foxa1_mcf7.bed',
        incl = "output/H3K27ac_mcf7.bed"

    params:
        trs = get_peaks_mcf7,
        minibatch_number = 1, minibatch_size = 1, threads = 4,
        max_combis = 20
        #size = 4,                       # NOT SURE... if the dict learning learns very big words they will be discarded, not truncated ! Maybe don't and see what happens
                                        # Maybe simply do exact like I originally intended
        


    output:'output/ologram_result_mcf7_filtered/00_ologram_stats.tsv'

    shell: """
    # Do the same but WITH MODL and CAPPED AT 4. Likely not display it or only in supp mat
    gtftk ologram -z -c hg38 -p {input.query}  --more-bed {params.trs} |\
        -o ologram_results_mcf7_filtered --force-chrom-peak --force-chrom-more-bed --no-date |\
        -V 3 -k {params.threads} -mn {params.minibatch_number} -ms {params.minibatch_size} |\    
        --more-bed-multiple-overlap --bed-incl {input.incl} |\
        --multiple-overlap-max-number-of-combinations {params.max_combis} #|\
        #--multiple-overlap-target-combi-size XXX                 
    """

# ---------------------------------------------------------------------------- #
#                              Benchmarking                                    #
# ---------------------------------------------------------------------------- #

rule produce_modl_comparison:
    """
    Benchmarking MODL on an artificial data test case (AB,ABCD,EF, and poisson noise)
    and comparing it to Apriori.
    """
    output: "output/supp_fig/comparison_benchmark.txt"
    script:
        "scripts/modl_comparison.py"



rule run_artificial:
    """
    Run on artificial data
    """
    input:
        query = "output/artificial/query.bed", incl = "output/artificial/incl.bed",
        a = "output/artificial/A.bed", b = "output/artificial/B.bed", c = "output/artificial/C.bed", d =  "output/artificial/D.bed",

    output:
        'output/ologram_result_artificial/00_ologram_stats.tsv'

    params:
        minibatch_number = 8, minibatch_size = 3, threads = 4,

    shell:"""
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {input.a} {input.b} {input.c} {input.d} -o output/ologram_result_artificial |\
     --force-chrom-peak --force-chrom-more-bed -V 3 -k ${params.threads} -mn ${params.minibatch_number} -ms ${params.minibatch_size} |\
     --more-bed-multiple-overlap --bed-incl {input.incl} --no-date
    """

# ---------------------------------------------------------------------------- #
#                                 Visual                                       #
# ---------------------------------------------------------------------------- #

# TODO Hardcoded for now. Of course, this is no longer germane when using different queries and cell lines.
def get_queryname(wildcards):
    run = wildcards.cell_line
    if run == "mcf7" or "mcf7_filtered" : return "FOXA1"
    if run == "artificial": return "Query"
    return " ".join(file_list)


rule treeify:
    """
    For visual purposes.
    Turns an OLOGRAM-MODL result file into a graph of the found combinations.
    """
    input: "output/ologram_result_{cell_line}/00_ologram_stats.tsv"
    output: "output/ologram/ologram_result_tree_{cell_line}.pdf"
    params:
        queryname = get_queryname
    shell:"""
    gtftk ologram_modl_treeify -i {input} -o {output} -l {params.queryname}
    """
