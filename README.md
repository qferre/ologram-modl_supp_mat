# OLOGRAM Supplementary material

This is the code accompanying the paper describing OLOGRAM-MODL. This code was used to create supplementary files and the paper's figures.

Here are the three possible run commands

```{bash}
# Full run
snakemake

# Dry run
    snakemake -n

# Don't run, produce a graph
snakemake --dag | dot -Tsvg > dag.svg
```

To clean the working directory :

```{bash}
rm -r output
git checkout # As some files were moved
```

TODO Remove the conda stuff and just say you need to install pygtftk first via conda and give command line

Say in readme data was taken from ReMap 2018 + REF

Available under GNU general public license v3