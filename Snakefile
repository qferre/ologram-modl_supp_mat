import os
import glob
from snakemake.utils import R
import re

workdir: os.getcwd()


REGIONS = ["exon-11-363", "exon-1-2",
           "exon-2-3", "exon-3-6",
           "exon-6-11", 'introns',
           'promoters', 'exons',
           'random']

H3K4_TYPE = ['H3K4me3_ENCFF616DLO_midpoint',
             'H3K4me3_ENCFF616DLO']
# Genome size
# See https://tinyurl.com/y6j367hv

HG38_SIZE = 2913022398

rule all:
    conda: "env/conda.yaml"
    input:  "output/exons_classes/ologram_exon-11-363_pygtftk.bed"
            #expand("output/supp_table2/binomial_test/{h3k4type}/binomial_test_{region}.txt", region=REGIONS, h3k4type=H3K4_TYPE)

#--------------------------------------------------------------
# Retrieve a set of peaks (GRCh38, H3K4me3, Homo sapiens K562,
# ENCFF616DLO). Keep only conventional chromosomes
#--------------------------------------------------------------

rule get_h3k4me3:
    output: "input/peaks/H3K4me3_ENCFF616DLO.bed"
    conda: "env/conda.yaml"
    shell: '''
    wget https://www.encodeproject.org/files/ENCFF616DLO/@@download/ENCFF616DLO.bed.gz -O H3K4me3_ENCFF616DLO.bed.gz
    gunzip -f H3K4me3_ENCFF616DLO.bed.gz
    cat H3K4me3_ENCFF616DLO.bed | perl -ne 'print if(/^chr[0-9XY]+\t/)' > {output}
    rm -f H3K4me3_ENCFF616DLO.bed
    '''

#--------------------------------------------------------------
# Retrieve a GTF from ensembl (Homo sapiens, release 92)
# add 'chr' prefix (-C) and select conventional chromosomes
#--------------------------------------------------------------

rule get_gtf_from_ensembl:
    conda: "env/conda.yaml"
    output: "input/Homo_sapiens_GRCh38_92_chr.gtf"
    shell: '''
    gtftk  retrieve -V 1 -Ccd -r 92 -s homo_sapiens | \
    gtftk select_by_regexp -V 1 -k chrom -r '^chr[0-9XY]+$'  -o {output}
    '''

#--------------------------------------------------------------
# Retrieve chromosome sizes subsequently used by Bedtools
#--------------------------------------------------------------

rule get_chrom_size:
    conda: "env/conda.yaml"
    output: "input/hg38.chromInfo"
    shell: '''
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
       -e "select chrom, size from hg38.chromInfo" | \
       perl -ne 'print if(/^chr[0-9XY]+\t/)' > {output}
    '''

#--------------------------------------------------------------
# Supplementary file 1:
# We calculate the significance of combinaison of epigenetic marks
# falling in GTF-defined numbered exons. 
#--------------------------------------------------------------

rule get_exon_classes:
    conda: "env/conda.yaml"
    input: gtf="input/Homo_sapiens_GRCh38_92_chr.gtf", peak="input/peaks/H3K4me3_ENCFF616DLO.bed"
    output: pdf="output/supp_fig2/supplfig2.pdf", \
            bed1="output/exons_classes/ologram_exon-11-363_pygtftk.bed", \
            bed2="output/exons_classes/ologram_exon-1-2_pygtftk.bed", \
            bed3="output/exons_classes/ologram_exon-2-3_pygtftk.bed", \
            bed4="output/exons_classes/ologram_exon-3-6_pygtftk.bed", \
            bed5="output/exons_classes/ologram_exon-6-11_pygtftk.bed"
    shell: """
       mkdir -p output/exons_classes
       gtftk add_exon_nb -k exon_nbr -i {input.gtf} | \
       gtftk discretize_key -p -d exon_nbr_cat -n 5 -k exon_nbr | \
       gtftk ologram -p {input.peak} -c hg38 -D -n -y -m exon_nbr_cat \
       -pf {output.pdf}  -k 8 -V 3 \
       -j summed_bp_overlaps_true -k 8 -D -o output/exons_classes \
       -K output/exons_classes/tmp
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__11_0_363_0_*.bed {output.bed1}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__1_0_2_0_*.bed {output.bed2}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__2_0_3_0_*.bed {output.bed3}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__3_0_6_0_*.bed {output.bed4}
    mv output/supp_fig2/tmp/ologram_exon_nbr_cat__6_0_11_0_*.bed {output.bed5}
    """

def get_label(wildcards):
    return re.sub('\W+', '_', wildcards.region)

rule supp_figure_1:
    conda: "env/conda.yaml"
    input:  bed="input/peaks/H3K4me3_ENCFF616DLO.bed", region="output/exons_classes/ologram_{region}_pygtftk.bed"
    output: pdf="output/supp_table2/ologram_{region}/ologram_{region}_regions.pdf"
    params: label=get_label
    shell: '''
    mkdir -p output/supp_table2/ologram/tmp
    gtftk ologram -z -y -V 2 -c hg38 -p {input.bed} -k 8 -o output/supp_table2/ologram_{wildcards.region} -D \
    -b {input.region} -K output/supp_table2/ologram_{wildcards.region}/tmp -pf {output.pdf} -l {params.label}
    '''

