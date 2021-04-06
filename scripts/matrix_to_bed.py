# Read snakemake parameters ?
# input = snakemake.input[0]

# Quick way to convert an intersection matrix (one line per intersection/position, one column per set) into a BED file for each set
import csv
import numpy as np

def intersection_matrix_to_bedfiles(X, output_root):
    """
    One bed file per column
    """
    bedfiles = [[]]*X.shape[1] # One BED file per region set

    # Consider increments of 100 base pairs for an abstraction, and all on the (fake) chr1
    pos = 0
    step = 100
    chrom = "chr1"

    # For each line in the matrix
    for line in X:
        # If the flag for this set was not zero, add a region to the BED file
        for j in len(line):
            if line[j] != 0:
                bedfiles[j] += [(chrom, pos, pos+step)]
        
        pos =+ step	# Increment
        
    # Finally, print the BED files
    filenames = [output_root+str(k)+".bed" for k in range(len(bedfiles))]
    for k in range(len(filenames)):
        with open(filenames[k], 'w', newline='') as f_output:
            tsv_output = csv.writer(f_output, delimiter='\t')
            tsv_output.writerows(bedfiles[k])

    return filenames