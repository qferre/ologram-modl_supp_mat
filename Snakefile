import os
import glob
import re

workdir: os.getcwd()


# ---------------------------------------------------------------------------- #
#                                 Parameters                                   #
# ---------------------------------------------------------------------------- #

# When making many runs and merging them with ologram_merge_runs, how many 
# should we make ?
N_RUNS_TO_MERGE = 16 

# How many threads to use for...
THREADS_SIMPLE = 8                  # Jobs with small files (and small RAM cost)
THREADS_SIMPLE_HIGH_COMPUTE = 8     # Simple jobs but with demanding algorithms (ie. MODL)
THREADS_DEMANDING = 4               # Jobs with a potentially high RAM cost per thread

# ---------------------------------------------------------------------------- #

## Query the final trees and the benchmarks
rule final:
    input: 
        ## OLOGRAM results
        # Artificial data and calibration
        expand("output/tree_results/ologram_result_tree_{testing_set}.pdf",
                testing_set = ['artificial','artificial_calibrate']), 
        "output/multovl_result_artificial_calibrate/calibrated.txt", # MULTOVL test
        # MCF7 - FOXA1 as query
        expand("output/tree_results/ologram_result_tree_{testing_set}.pdf",
            testing_set = ['mcf7', 'mcf7_filtered','mcf7_manual']),
        # MCF7 - full DHS as query
        "output/tree_results/ologram_result_tree_mcf7_full_dhs.pdf",
        # sc-ATAC-Seq and combination entropy
        "output/ologram_result_scatacseq_pbmc/done", 
        # Murine promoters
        "output/murine_result/murine_fig.png", "output/murine_result_restricted/murine_fig.png",
        # Comparison with GINOM
        expand("output/tree_results/ologram_result_tree_{testing_set}.pdf",
                testing_set = ['ginom','ginom_filtered']),
        ## MODL benchmarks
        "output/benchmark/comparison/done",
        expand("output/benchmark/scaling/fig{n}.png", n = [1,2]),
        # Elementary benchmarks
        expand("output/benchmark/scaling/fig{n}.png", n = [3,4,5])
    # shell: """
    # # Produce a summary graph
    # snakemake --forceall --dag | dot -Tsvg > output/dag.svg
    # snakemake --forceall --rulegraph | dot -Tsvg > output/rulegraph.svg
    # """

# ---------------------------------------------------------------------------- #
#                        Data retrieval and preparation                        #
# ---------------------------------------------------------------------------- #

rule prepare_incl:
    """
    Restrict the shuffling to "pseudo-CRMs", obtained by merging all TF peaks.
    Then move the query (FOXA1 for MCF7) to prevent them from being used by MODL.

    Also do that for the additional MCF7 TFs
    """
    input:
        "input/mcf7/foxa1.bed" # Ensure files have been uncompressed
    output:
        query = "input/foxa1_mcf7.bed",
        incl = "input/crm_mcf7.bed",
        incl_extended = "input/crm_mcf7_extended.bed"

    shell: """
        # Create incl file by merging all BEDs
        cat input/mcf7/*.bed | bedtools sort | bedtools merge > {output.incl}

        # And add another one
        cat input/mcf7_additional/*.bed > input/extended_all_mcf7.bed
        awk -F"\t" 'NF{{NF-=3}};3' < input/extended_all_mcf7.bed > input/extended_all_mcf7_3col.bed
        cat {output.incl} input/extended_all_mcf7_3col.bed | sed -e 's/ /\t/g' | bedtools sort | bedtools merge > {output.incl_extended}

        # Once done, move future query
        mv input/mcf7/foxa1.bed {output.query}
    """


rule uncompress:
    """
    Uncompress certain relevant files and signal to Snakemake that we now have them.
    """
    output:
        "input/desired_combis.txt",
        "input/DNAse-seq_MCF7_ENCSR000EPH_rep1_1_se_bwa_biorep_filtered_peaks_aka_ENCFF961ZCT.bed",
        "input/murine_chipatlas.txt",
        "input/as_ginom/mappability_human.bed",
        "input/as_ginom/query_som_trans_lung.bed",
        "input/mcf7/foxa1.bed",
        "input/hg38.genome",
        "input/artificial_simple.genome",

    shell: """
        gunzip input/*.gz
        gunzip input/*/*.gz    
    """


rule prepare_artificial:
    """
    Prepare artificial demonstration data.
    It consits of random regions for the query, compared against (a) a third of 
    the query, (b) another third, and (c) a negative control.
    """
    input:
        "input/hg38.genome"
    output:
        query = "output/artificial_data/query.bed",
        a = "output/artificial_data/data/third.bed", 
        a2 = "output/artificial_data/data/third_bis.bed", 
        b = "output/artificial_data/data/other_third.bed", 
        c = "output/artificial_data/data/neg_control.bed"

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



rule prepare_artificial_calibrate:
    """
    Artificial data with precisely calibrated regions.

    Show that MULTOVL precision will bottom to its number of shuffles, but our precision will still be good.
    These are an "easy" case where the true p-vals may be computed by MULTOVL so we can compare. 
    Also run MULTOVL on original artifical data, the neg control will serve as further confirmation.
    """
    input:
        "input/artificial_simple.genome",
    output:
        query = "output/artificial_data_calibrate/query.bed",
        filler1 = "output/artificial_data_calibrate/filler1.bed",
        filler2 = "output/artificial_data_calibrate/filler2.bed",
        a = "output/artificial_data_calibrate/data/21pc.bed", 
        b = "output/artificial_data_calibrate/data/28pc.bed"

    params:
        size=10, length=250   # Note that size is multiplied by 100 later

    shell: """
   
    mkdir -p output/artificial_data_calibrate/data

    # Prepare query
    bedtools random -n $((100*{params.size})) -l $(({params.length})) -seed 123456 -g input/artificial_simple.genome > {output.query}

    # And filler regions that do NOT overlap with query
    # Use a different set each time
    bedtools shuffle -i {output.query} -excl {output.query} -g input/artificial_simple.genome -noOverlapping -seed 654321 > {output.filler1}
    bedtools shuffle -i {output.query} -excl {output.query} -g input/artificial_simple.genome -noOverlapping -seed 789789 > {output.filler2}


    # Calibrating : get a third overlapping, a fifth, a twentieth, etc.
    # and fill the rest with random regions that do NOT overlap with the query

    head -n $((21*{params.size})) {output.query} > {output.a}
    head -n $((79*{params.size})) {output.filler1} >> {output.a}

    tail -n $((26*{params.size})) {output.query} > {output.b}
    tail -n $((74*{params.size})) {output.filler2} >> {output.b}
    """



rule prepare_sc_atac_seq_bedfiles:
    """
    From the Signac vignette, download and prepare the data.
    For all the cells, turn their sc-ATAC-seq sites into BED files, one per cell.
    """
    output:
        translate = "output/sc_atac_seq_pbmc_data/cell_to_class.tsv",
        allreg = "output/sc_atac_seq_pbmc_data/all_pbmc_regions.bed",
        allreg_slopped = "output/sc_atac_seq_pbmc_data/all_pbmc_regions_slopped_250_both.bed"

    shell:"""
    mkdir -p output/sc_atac_seq_pbmc_data/
    Rscript scripts/pbmc_download_data.R

    # Preparing the incl files
    cat output/sc_atac_seq_pbmc_data/bed/*.bed | bedtools sort | bedtools merge > {output.allreg}
    bedtools slop -i {output.allreg} -g input/hg19.genome -b 250 > {output.allreg_slopped}
    """




rule download_murine_chip:
    """
    Download murine ChIP-Seq and prepare bed-incl file.

    Also slop TSS by 1kb to the left to estimate promoters, if required.
    """
    input:
        sources = "input/murine_chipatlas.txt"
    output: 
        incl = "output/murine_data/murine_incl.bed"
    params:
        incl_prom_left_slop = 1000
    
    shell:"""
    mkdir -p output/murine_data/selected_tf/

    # Download data
    for i in `cat {input.sources} | grep mm10 | awk '{{print $0,"\t","bla"}}' | sed 's/\"//g' | cut -f1,6,7 | perl -npe 's/\t/-/g'| perl -npe 's/ +/_/g'`
    do
        n=`echo $i | sed 's/\-.*//'`
        wget http://dbarchive.biosciencedbc.jp/kyushu-u/mm10/eachData/bed10/$n.10.bed -O output/murine_data/selected_tf/$i.bed --cipher 'DEFAULT:!DH'
    done


    # Prepare incl
    wget reftss.clst.riken.jp/datafiles/current/mouse/refTSS_v3.1_mouse_coordinate.mm10.bed.gz -O output/murine_data/refTSS_v3.1_mouse_coordinate.mm10.bed.gz --cipher 'DEFAULT:!DH'
    gunzip output/murine_data/refTSS_v3.1_mouse_coordinate.mm10.bed.gz

    cat output/murine_data/refTSS_v3.1_mouse_coordinate.mm10.bed | bedtools slop -l {params.incl_prom_left_slop} -r 1 -s -g input/mm10.genome | bedtools sort | bedtools merge > {output.incl}
    """


## Get the list of files in the various input directories

def get_peaks_artificial(wildcards):
    file_list = sorted(glob.glob('output/artificial_data/data/*bed'))
    return " ".join(file_list)

def get_peaks_artificial_calibrate(wildcards):
    file_list = sorted(glob.glob('output/artificial_data_calibrate/data/*bed'))
    return " ".join(file_list)

def get_peaks_artificial_calibrate_as_list(wildcards): 
    res = sorted(glob.glob('output/artificial_data_calibrate/data/*bed'))
    if res != []: return res
    else: return [None] * 1000

def get_peaks_ginom(wildcards):
    file_list = sorted(glob.glob('input/as_ginom/reference_*bed'))
    return " ".join(file_list)

def get_peaks_mcf7(wildcards):
    file_list = sorted(glob.glob('input/mcf7/*bed'))
    return " ".join(file_list)

def how_many_peaks_mcf7(wildcards):
    file_list = sorted(glob.glob('input/mcf7/*bed'))
    return len(file_list)

def get_peaks_mcf7_with_additional(wildcards):
    file_list_main = []
    file_list_additional = sorted(glob.glob('input/mcf7_additional/*bed'))
    return " ".join(file_list_main + file_list_additional)

def get_peaks_sc_atac_seq(wildcards): return sorted(glob.glob('output/sc_atac_seq_pbmc_data/bed_selected/*bed'))


# ---------------------------------------------------------------------------- #
#                                 Data processing                              #
# ---------------------------------------------------------------------------- #

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
        minibatch_number = 10, minibatch_size = 5
    threads: THREADS_SIMPLE

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
        minibatch_number = 1, minibatch_size = 1,
        max_combis = 20

    threads: THREADS_SIMPLE_HIGH_COMPUTE

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


rule mcf7_full_dhs:
    """
    For a the MCF7 cell line, compute the enrichment in n-wise TF combinations 
    using OLOGRAM-MODL. The query is not FOXA1, but the DHS regions in MCF7. 
    We restrict to those regions as well.
    This time however, I will use ReMap 2018 full data for 40 TFs, not something
    curated in advance.
    """
    input: 
        query = "input/DNAse-seq_MCF7_ENCSR000EPH_rep1_1_se_bwa_biorep_filtered_peaks_aka_ENCFF961ZCT.bed",
        incl = "input/crm_mcf7.bed"
    params:
        trs = get_peaks_mcf7_with_additional,
        minibatch_number = 20, minibatch_size = 1
    threads: THREADS_DEMANDING 

    output: 'output/ologram_result_mcf7_full_dhs/00_ologram_stats.tsv', 

    shell: """
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.trs} \
        -o output/ologram_result_mcf7_full_dhs --force-chrom-peak --force-chrom-more-bed  \
        -V 3 -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} \
        --more-bed-multiple-overlap --bed-incl {input.incl} --no-date
    """



# ---------------------------------------------------------------------------- #
#                              Benchmarking  MODL                              #
# ---------------------------------------------------------------------------- #

rule produce_modl_comparison:
    """
    Benchmarking MODL on an artificial data test case (AB,ABCD,EF, and poisson noise)
    and comparing it to Apriori, FP-Growth and others.
    """
    output: "output/benchmark/comparison/done",
    threads: THREADS_SIMPLE_HIGH_COMPUTE
    log:
        err= "output/benchmark/comparison/comparison_benchmark_ERROR_LOG.txt",
        out= "output/benchmark/comparison/comparison_benchmark.txt"

    shell:"""
    mkdir -p output/benchmark
    python scripts/modl_comparison.py 2> {log.err} 1> {log.out}

    # Signal we are done
    touch {output}
    """



"""
Rules to evaluate MODL time scaling, and scaling of the elementary operations, in parallel
"""
rule produce_modl_scaling:
    output: 
        fig1 = "output/benchmark/scaling/fig1.png",
        fig2 = "output/benchmark/scaling/fig2.png",
    threads: THREADS_SIMPLE_HIGH_COMPUTE
    log:
        err= "output/benchmark/scaling/scaling_benchmark_1_2_ERROR_LOG.txt",
        out= "output/benchmark/scaling/scaling_benchmark_1_2.txt"

    shell:"""  
    mkdir -p output/benchmark/scaling

    # Print information about the CPU
    lscpu > output/cpu_info.txt

    python scripts/modl_scaling_1_2.py 2> {log.err} 1> {log.out}
    """

rule modl_elementary_scaling_1:
    output: "output/benchmark/scaling/fig3.png"
    threads: 1
    log:
        err= "output/benchmark/scaling/scaling_benchmark_3_ERROR_LOG.txt",
        out= "output/benchmark/scaling/scaling_benchmark_3.txt"

    shell:"""
    mkdir -p output/benchmark/scaling
    python scripts/modl_scaling_3.py 2> {log.err} 1> {log.out}
    """

rule modl_elementary_scaling_2:
    output: "output/benchmark/scaling/fig4.png"
    threads: 1
    log:
        err= "output/benchmark/scaling/scaling_benchmark_4_ERROR_LOG.txt",
        out= "output/benchmark/scaling/scaling_benchmark_4.txt"

    shell:"""
    mkdir -p output/benchmark/scaling
    python scripts/modl_scaling_4.py 2> {log.err} 1> {log.out}
    """

rule modl_elementary_scaling_3:
    output: "output/benchmark/scaling/fig5.png"
    threads: 1
    log:
        err= "output/benchmark/scaling/scaling_benchmark_5_ERROR_LOG.txt",
        out= "output/benchmark/scaling/scaling_benchmark_5.txt"

    shell:"""
    mkdir -p output/benchmark/scaling
    python scripts/modl_scaling_5.py 2> {log.err} 1> {log.out}
    """





rule run_on_ginom_data:
    """
    For comparison with GINOM, run on their example data files instead and 
    compare the results.

    NOTE: here is how I converted their data.

    # For the xslx, in bash
    for filename in "$DATA_DIR"/*.xlsx; do
        libreoffice --calc --headless --convert-to csv "$filename" --outdir "{params.outdir}"
    done
    for filename in "{params.outdir}"/*.csv; do
        tr '\t' ',' < "$filename" > "$filename".bed
    done

    # For the mat files, in Octave
    load filename.mat
    # Then simply use the variable editor and copy paste into a libreoffice
    """

    input: 
        query = "input/as_ginom/query_som_trans_lung.bed",
        incl = "input/as_ginom/mappability_human.bed"

    params:
        trs = get_peaks_ginom,
        minibatch_number = 10, minibatch_size = 10,
        max_combis = 7 # As in the GINOM example

    threads: THREADS_SIMPLE

    output:
        with_modl = 'output/ologram_result_ginom_filtered/00_ologram_stats.tsv',
        without_modl = 'output/ologram_result_ginom/00_ologram_stats.tsv'

    shell: """
    mkdir -p output/ologram_result_ginom
    
    # With MODL
    gtftk ologram -z -c hg19 -p {input.query} --more-bed {params.trs} \
        -o output/ologram_result_ginom_filtered --force-chrom-peak --force-chrom-more-bed --no-date \
        -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} -V 3 \
        --more-bed-multiple-overlap --bed-incl {input.incl} \
        --multiple-overlap-max-number-of-combinations {params.max_combis}   

    # Without MODL
    gtftk ologram -z -c hg19 -p {input.query} --more-bed {params.trs}\
        -o output/ologram_result_ginom --force-chrom-peak --force-chrom-more-bed --no-date \
        -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} -V 3 \
        --more-bed-multiple-overlap --bed-incl {input.incl}      
    """


rule run_ologram_murine:
    """
    Run OLOGRAM on the murine data. All at once, those are not big files.
    """
    input: 
        incl = "output/murine_data/murine_incl.bed"
    output: expand("output/murine_{kind}/{run}/00_ologram_stats.tsv", run=["SRX1815531-Ctcf","SRX1583885-Irf1","ERX1633247-Nanog"], kind=["result","result_restricted"])
    params: 
        minibatch_number = 20, minibatch_size = 10,
    threads: THREADS_SIMPLE

    shell:"""
    mkdir -p output/murine_result
    mkdir -p output/murine_result_restricted

    # Without restricting to estimated promoters only

    gtftk ologram -z -f -w -q -c mm10 -p output/murine_data/selected_tf/SRX1815531-Ctcf-Blood.bed \
        --more-bed `ls output/murine_data/selected_tf/* | grep -vi ctcf` --more-bed-multiple-overlap --no-date \
        -k {threads} -V 3 -j summed_bp_overlaps_pvalue -a summed_bp_overlaps_pvalue -g 0.05 -o output/murine_result/SRX1815531-Ctcf \
        -mn {params.minibatch_number} -ms {params.minibatch_size}

    gtftk ologram -z -f -w -q -c mm10 -p output/murine_data/selected_tf/SRX1583885-Irf1-Blood.bed \
        --more-bed `ls output/murine_data/selected_tf/* | grep -vi Irf1` --more-bed-multiple-overlap --no-date \
        -k {threads} -V 3  -j summed_bp_overlaps_pvalue -a summed_bp_overlaps_pvalue -g 0.05 -o output/murine_result/SRX1583885-Irf1 \
        -mn {params.minibatch_number} -ms {params.minibatch_size}

    gtftk ologram -z -f -w -q -c mm10 -p output/murine_data/selected_tf/ERX1633247-Nanog-Pluripotent_stem_cell.bed \
        --more-bed `ls output/murine_data/selected_tf/* | grep -vi nanog` --more-bed-multiple-overlap --no-date \
        -k {threads} -V 3  -j summed_bp_overlaps_pvalue -a summed_bp_overlaps_pvalue -g 0.05 -o output/murine_result/ERX1633247-Nanog \
        -mn {params.minibatch_number} -ms {params.minibatch_size}



    # With the restriction

    gtftk ologram -z -f -w -q -c mm10 -p output/murine_data/selected_tf/SRX1815531-Ctcf-Blood.bed \
        --more-bed `ls output/murine_data/selected_tf/* | grep -vi ctcf` --more-bed-multiple-overlap --no-date \
        -k {threads} -V 3 -j summed_bp_overlaps_pvalue -a summed_bp_overlaps_pvalue -g 0.05 -o output/murine_result_restricted/SRX1815531-Ctcf \
        -mn {params.minibatch_number} -ms {params.minibatch_size} --bed-incl {input.incl}

    gtftk ologram -z -f -w -q -c mm10 -p output/murine_data/selected_tf/SRX1583885-Irf1-Blood.bed \
        --more-bed `ls output/murine_data/selected_tf/* | grep -vi Irf1` --more-bed-multiple-overlap --no-date \
        -k {threads} -V 3  -j summed_bp_overlaps_pvalue -a summed_bp_overlaps_pvalue -g 0.05 -o output/murine_result_restricted/SRX1583885-Irf1 \
        -mn {params.minibatch_number} -ms {params.minibatch_size} --bed-incl {input.incl}

    gtftk ologram -z -f -w -q -c mm10 -p output/murine_data/selected_tf/ERX1633247-Nanog-Pluripotent_stem_cell.bed \
        --more-bed `ls output/murine_data/selected_tf/* | grep -vi nanog` --more-bed-multiple-overlap --no-date \
        -k {threads} -V 3  -j summed_bp_overlaps_pvalue -a summed_bp_overlaps_pvalue -g 0.05 -o output/murine_result_restricted/ERX1633247-Nanog \
        -mn {params.minibatch_number} -ms {params.minibatch_size} --bed-incl {input.incl}
    """



# ---------------------------------------------------------------------------- #
#                       Benchmarking  OLOGRAM (p-values)                       #
# ---------------------------------------------------------------------------- #


rule run_artificial:
    """
    Run OLOGRAM-MODL on artificial data.

    Calibrated data is an easier class of data, where p-values are higher, for 
    comparision with tools using empirical p-valuess.
    """
    input:
        query = "output/artificial_data/query.bed",
        query_calibrate = "output/artificial_data_calibrate/query.bed"

    output:
        regular = 'output/ologram_result_artificial/00_ologram_stats.tsv',
        calibrated = 'output/ologram_result_artificial_calibrate/00_ologram_stats.tsv'

    params:
        minibatch_number = 20, minibatch_size = 10,
        peaks = get_peaks_artificial,
        peaks_calibrate = get_peaks_artificial_calibrate
    threads: THREADS_SIMPLE
    shell:"""

    # Run on artificial data
    gtftk ologram -z -c hg38 -p {input.query} --more-bed {params.peaks} \
        -o output/ologram_result_artificial --force-chrom-peak --force-chrom-more-bed \
        -V 3 -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} \
        --more-bed-multiple-overlap --no-date -K output/TMP_ologram_result_artificial

    # Run on calibrated artificial data
    gtftk ologram -z -c input/artificial_simple.genome -p {input.query_calibrate} --more-bed {params.peaks_calibrate} \
        -o output/ologram_result_artificial_calibrate --force-chrom-peak --force-chrom-more-bed \
        -V 3 -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} \
        --more-bed-multiple-overlap --no-date -K output/TMP_ologram_result_artificial_calibrate
    """



rule run_multovl_on_artificial:
    """
    Run MULTOVL to confirm the p-values on the artificial calibrated data. 
    Also run on the regular artificial data to show its precision bottoms out.

    Comparably to OLOGRAM, the free regions here (equivalent of bed-incl) will
    be made of all the genome.
    """
    input:
        query = "output/artificial_data/query.bed",
        query_calibrate = "output/artificial_data_calibrate/query.bed"

    output:
        artificial = 'output/multovl_result_artificial_calibrate/artificial.txt',
        calibrated = 'output/multovl_result_artificial_calibrate/calibrated.txt'

    params:
        shuffle_numbers = 200,
        peaks = get_peaks_artificial,
        peaks_calibrate = get_peaks_artificial_calibrate,
        peaks_calibrate_as_list = get_peaks_artificial_calibrate_as_list,
        multovl_exec = "ext/multovl-1.3/bin/multovlprob"

    shell: """
    mkdir -p output/multovl_result_artificial_calibrate

    # On full artificial data
    {params.multovl_exec} {input.query} {params.peaks} --nointrack \
        --reshufflings {params.shuffle_numbers} --free input/hg38.genome.bed \
        > {output.artificial}
    

    # On simpler, calibrated artificial data
    {params.multovl_exec} {input.query_calibrate} {params.peaks_calibrate} --nointrack \
        --reshufflings {params.shuffle_numbers} --free input/artificial_simple.genome.bed \
        > output/multovl_result_artificial_calibrate/partial1.txt

    {params.multovl_exec} {input.query_calibrate} {params.peaks_calibrate_as_list[0]} --nointrack \
        --reshufflings {params.shuffle_numbers} --free input/artificial_simple.genome.bed \
        > output/multovl_result_artificial_calibrate/partial2.txt

    {params.multovl_exec} {input.query_calibrate} {params.peaks_calibrate_as_list[1]} --nointrack \
        --reshufflings {params.shuffle_numbers} --free input/artificial_simple.genome.bed \
        > output/multovl_result_artificial_calibrate/partial3.txt


    # Concatenate results for calibrated, with new line between each
    awk 'FNR==1{{print ""}}1' output/multovl_result_artificial_calibrate/partial* > {output.calibrated}
    """





rule run_ologram_sc_atac_seq:
    """
    Enrichment depending on the entropy of the combinations.
    Elements are grouped by supercluster of cells (CD4+8, CD14, B)

    The idea is that each of these runs will make a few shuffles due to RAM cost, then they will all be merged.
    """
    input: "output/sc_atac_seq_pbmc_data/all_pbmc_regions_slopped_250_both.bed"
    output: "output/ologram_result_scatacseq_pbmc/run_{i}/00_ologram_stats.tsv"
    params:
        minibatch_number = 2, minibatch_size = 2, 
        peaks = get_peaks_sc_atac_seq
    threads: THREADS_DEMANDING

    shell:"""
    mkdir -p output/ologram_result_scatacseq_pbmc/run_{wildcards.i}

    gtftk ologram -z -c hg19 -p {input} --more-bed {params.peaks} --no-date \
        -o output/ologram_result_scatacseq_pbmc/run_{wildcards.i} \
        --force-chrom-peak --force-chrom-more-bed \
        -k {threads} -mn {params.minibatch_number} -ms {params.minibatch_size} -V 3 \
        --more-bed-multiple-overlap --bed-incl {input}
    """



rule ologram_sc_atac_seq_merge_runs:
    """
    Merge the various runs using ologram_merge_runs to recalculate the statistics.
    """
    input: expand("output/ologram_result_scatacseq_pbmc/run_{runid}/00_ologram_stats.tsv", runid = range(N_RUNS_TO_MERGE))
    output: "output/ologram_result_scatacseq_pbmc/merged_batches_result.tsv"
    shell: """
    gtftk ologram_merge_runs --inputfiles `ls output/ologram_result_scatacseq_pbmc/*/*.tsv` -o {output} -V 3
    """


rule ologram_sc_atac_seq_analysis:
    """
    Python script to draw the figure and perform the analysis for sc-ATAC-seq data.
    """
    input: "output/ologram_result_scatacseq_pbmc/merged_batches_result.tsv"
    output: "output/ologram_result_scatacseq_pbmc/done"
    shell: """
    mkdir -p output/ologram_result_scatacseq_pbmc/entropy_graph/

    # Run the Python script to produce the figures
    python scripts/combi_entropy_analysis.py

    # Signal we are done
    touch {output}
    """



rule process_ologram_results_murine:
    """
    Grab the headers, and run the R script for murine data.
    """
    input: expand("output/murine_{kind}/{run}/00_ologram_stats.tsv", run=["SRX1815531-Ctcf","SRX1583885-Irf1","ERX1633247-Nanog"], kind=["result","result_restricted"])
    output: "output/murine_{kind}/murine_fig.png"
    shell: """
    
    ## Merging OLOGRAM results
    # Print the data and a supplementary column for file names
    cat output/murine_{wildcards.kind}/*X*/*tsv | head -n 1 | awk 'BEGIN{{FS=OFS="\t"}}{{print $0,"query"}}'> output/murine_{wildcards.kind}/all_ologram_results.txt
    # Print data without header
    for i in `ls output/murine_{wildcards.kind}/*X*/*tsv`; do  awk -v f=$i 'BEGIN{{FS=OFS="\t"}}{{print $0,f}}' $i ; done | grep -v nb_intersections_expectation_shuffled >> output/murine_{wildcards.kind}/all_ologram_results.txt

    ## Now call the R script, with the point of execution in argument
    Rscript scripts/murine_analysis.R "./output/murine_{wildcards.kind}/"
    """



# ---------------------------------------------------------------------------- #
#                                 Visual                                       #
# ---------------------------------------------------------------------------- #

# TODO Hardcoded for now. Of course, this is no longer apropos when using 
# different queries and cell lines.
def get_queryname(wildcards):
    run = wildcards.cell_line

    if "mcf7" in run:
        if "full" in run: return "DHS"
        else: return "foxa1"
    if "artificial" in run: return "Query"


rule treeify:
    """
    For visual purposes.
    Turns an OLOGRAM-MODL result file into a graph of the found combinations.
    """
    input: "output/ologram_result_{cell_line}/00_ologram_stats.tsv"
    output: "output/tree_results/ologram_result_tree_{cell_line}.pdf"
    params:
        queryname = get_queryname
    shell: """
    gtftk ologram_modl_treeify -i {input} -o {output} -l {params.queryname}
    """