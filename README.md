# OLOGRAM-MODL Supplementary material

This is the code accompanying the paper describing OLOGRAM-MODL. It was used to create supplementary files and the paper's figures.

To run it, create first a conda environment containing *pygtftk* with `conda env create -f ologram_modl_env.yaml`. This will download the latest pygtftk release from bioconda and install Snakemake, as well as all requires dependencies.

Then, to run the workflow itself, use these commands:

```{bash}
# Full run, where -j is the number of jobs in parallel
snakemake -j4

# Dry run (nothing is done)
snakemake -n -p --reason

# Example: full run on a cluster, with qsub, on the "tagc" queue
snakemake -c 'qsub -q tagc -V -l nodes=1:ppn={threads}' -j32

```

To clean the working directory:

```{bash}
rm -rf output
rm -rf input
git reset --hard
```

This code is available under GNU general public license v3.