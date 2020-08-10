# OLOGRAM Supplementary material

This is the code accompanying the paper describing OLOGRAM-MODL. This code was used to create supplementary files and the paper's figures.

To run :

```{bash}
dry run :
    snakemake -n

graph :
    snakemake --dag | dot -Tsvg > dag.svg
```



To clean run :

```{bash}
rm -r output
git checkout # Cause some files were moved
```






TODO Remove the conda stuff and just say you need to install pygtftk first via conda and give command line


To check : I think I am now the owner of this repository, Denis transfered it to me. Careful with the links in the paper.







Say in readme data was taken from ReMap 2018 + REF and ENCODE + REF









Available under GNU general public license v3