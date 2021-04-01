# OLOGRAM-MODL Supplementary material

This is the code accompanying the paper describing OLOGRAM-MODL. It was used to create supplementary files and the paper's figures.

To run it, create first a conda environment containing pygtftk with `conda env create -f ologram_modl_env.yaml`. This will download the latest pygtftk release from bioconda and install Snakemake, as well as all requires dependencies.

To run the workflow itself, use these commands:

```{bash}
# Full run, where -j is the number of jobs
snakemake -j4
```

Here are some auxiliary commands:

```{bash}
# Dry run: nothing will be done, but instead print the shell commands and their reasons.
snakemake -n -p --reason


# Example: running it on a cluster with qsub, with resources per node
snakemake -c 'qsub -q tagc -V -l nodes=1:ppn={threads}' -j25

```

To clean the working directory:

```{bash}
rm -rf output
rm -rf input
git reset --hard
```

This code is available under GNU general public license v3.