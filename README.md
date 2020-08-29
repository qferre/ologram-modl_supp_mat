# OLOGRAM-MODL Supplementary material

This is the code accompanying the paper describing OLOGRAM-MODL. It was used to create supplementary files and the paper's figures.

To run it, create first a conda environment containing pygtftk with `conda env create -f ologram_modl_env.yaml`. This will download the latest pygtftk release from bioconda and install snakemake.

To run the workflow itself, use these commands:

```{bash}
# Full run
snakemake

# Dry run
snakemake -n

# Don't run, produce a graph
snakemake --dag | dot -Tsvg > dag.svg
```

To clean the working directory:

```{bash}
rm -r output
find . -name '*.bed' -delete
git reset --hard
```

Data was selected from ReMap 2018. This code is available under GNU general public license v3.
