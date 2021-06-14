"""
A perspective application of MODL using a supervised loss, based on the 
performance of a Naive Bayes classifier.
This will select the combinations (of reference sets) that best predict the
query set using a Naive Bayes classifier.

This is currently in beta and requires some tuning of the parameters, so is not
yet available as option in the main command line, but is an interesting perspective.

This file is a fully functional example. To use it, simply replace the filepaths
at the beginning with the paths to your own files, adjust the parameters if 
desired, and run the script. 
"""

## Paths to BED files of regions

QUERY_PATH  = "input/as_ginom/query_som_trans_lung.bed"        # The query set
INCL_PATH   = "input/as_ginom/mappability_human.bed"           # Inclusion: the analysis will be limited to the sub-genome designated in this files. Regions outside of it will be discarded.
GENOME_PATH = "input/hg19.genome"                              # Chromosome sizes

# Reference sets
MORE_PATHS  = [
    "input/as_ginom/reference_1_hg19_Refseq_promoters_2000bp.bed",
    "input/as_ginom/reference_2_hg19_Refseq_gene_bodies.bed",
    "input/as_ginom/reference_3_Broad_A549_H3k04me3_Dex.bed",
    "input/as_ginom/reference_4_Broad_A549_H3k36me3_Dex.bed",
    "input/as_ginom/reference_5_Broad_A549_H3k79me2_Dex.bed"
]
# The meaning of the combinations  returnedwill always be : query, followed by the reference sets in the order they are given
# For example, '[101010]' means "Query + Reference 2 (gene bodies) + Reference 4 (H3K36me3)"


## Parameters

# Subsampling parameters
KEEP_N_TIMES_QUERY = 100          # For each '1' line with the query, how many '0' lines without it do we keep?
QUERY_WEIGHT = 100                # In the Naive Bayes classifier, how much more do we weigh the presence of the query?

# MODL parameters
DESIRED_ITEMSETS = 7              # How many interesting combinations should MODL return?

RANDOM_SEED = 42

# ---------------------------------------------------------------------------- #

def presence_vector(row, combis):
    """
    Given a list of combinations, returns a vector saying whether they are
    present or absent in this row. The order in the vector is the same as given
    in the combis list.

    This works exactly like an inexact (transitive) count, which is OLOGRAM-MODL's default operating mode.

    Example:

    >>> row = [1,0,1,1,0,1]
    >>> combis = [1,0,1,0,0,0],[0,1,1,0,0,0]
    >>> res = presence_vector(row, combis)
    >>> assert list(res) == [1,0]

    """

    result = []
    for combi in combis:
        mask = row - combi

        # If any element from the combi is missing, it's absent, so
        # add a 0. Otherwise add a 1.
        if min(mask) < 0 : result +=[0]
        else: result += [1]

    return result


import pandas as pd
import numpy as np
import sklearn
import pybedtools
from sklearn.naive_bayes import GaussianNB

from pygtftk.stats.intersect.overlap_stats_compute import compute_true_intersection
from pygtftk.stats.intersect.overlap.overlap_regions import does_combi_match_query
from pygtftk.stats.intersect.modl.dict_learning import Modl
from pygtftk import utils
from pygtftk.stats.intersect.read_bed import read_bed_as_list as read_bed
from pygtftk.utils import chrom_info_as_dict
from pygtftk import arg_formatter

np.random.seed(RANDOM_SEED)  # Random seed for reproducibility



## Prepare files
bedA = pybedtools.BedTool(QUERY_PATH).sort().merge()
bedsB = [pybedtools.BedTool(bedfilepath).sort().merge() for bedfilepath in MORE_PATHS]


# Do the exclusion manually Generate a fake bed for the entire genome, using the chromsizes
bed_incl = pybedtools.BedTool(INCL_PATH)
chrom_len = chrom_info_as_dict(open(GENOME_PATH, 'r'))

full_genome_bed = [str(chrom) + '\t' + '0' + '\t' + str(chrom_len[chrom]) + '\n' for chrom in chrom_len if chrom != 'all_chrom']
full_genome_bed = pybedtools.BedTool(full_genome_bed)
bed_excl = full_genome_bed.subtract(bed_incl)

bedA = read_bed.exclude_concatenate(bedA, bed_excl)
bedsB = [read_bed.exclude_concatenate(bedB, bed_excl) for bedB in bedsB]
full_genome_bed_after_excl = read_bed.exclude_concatenate(full_genome_bed, bed_excl)


# Note that by definition, in this intersections' matrix only regions where at 
# least two sets are open are given. For example {4} alone is not found. 
# To fix it, use a fake full genome bed as query, so there is always one file 
# open, then truncate it to get the flags_matrix. 
#true_intersection = compute_true_intersection(bedA, bedsB)
true_intersection = compute_true_intersection(full_genome_bed_after_excl, [bedA] + bedsB)
flags_matrix = np.array([i[3] for i in true_intersection])

"""
# NOTE : The length of each intersection is also accessible, in case you want
# the final matrix to have one row per X base pairs instead of 1 row per intersection/event.
# For reference, the true_intersection object looks like this: [('chr1', 150, 180, np.array([1,1,0])), ('chr1', 180, 200, np.array([1,1,1])), ('chr1', 350, 380, np.array([1,1,0]))]
for i in true_intersection: # For each observed intersection...
    inter_chr, inter_start, inter_stop = i[0], i[1], i[2]   # Get the chromosome, start and end of the observed intersection
    intersection_length = inter_stop - inter_start          # Deduce the length of the intersection
    intersection_flags = i[3]                               # Which combination was observed here?
    # Now you can process it.
"""

flags_matrix = flags_matrix[:,1:] # Remove first column that represents the fake, full genome BED.


# Print unique rows with transitive (inexact) counting
uniqueValues, occurCount = np.unique(flags_matrix, return_counts=True, axis=0)
for uniqueRow in uniqueValues:
    count = 0

    for row in flags_matrix:
        if does_combi_match_query(
            tuple(row), 
            tuple(uniqueRow), 
            exact = False):
            count += 1
    print(uniqueRow, ' - freq = ', count)
print("------------------")


# Perform subsampling of the majority class ('0', meaning query is absent)
flags_matrix_with_query = flags_matrix[flags_matrix[:, 0] != 0]
draw_this_many = int(np.around(KEEP_N_TIMES_QUERY * len(flags_matrix_with_query)))
random_rows = flags_matrix[np.random.randint(flags_matrix.shape[0], size=draw_this_many),:]
flags_matrix = np.concatenate((flags_matrix_with_query,random_rows))



# This is the custom error function that we will use
def custom_error_function_bnb(X_true, X_rebuilt, encoded, dictionary):
    """
    This is a non-pure function that only cares about the dictionary and discards X_true in 
    favor of the previous flags_matrix.
    """

    # Always assume that the first column is the query, so we split it.
    # Truncate the first column, with the query, since that is what we want to predict
    Y = flags_matrix[:,0]
    X = flags_matrix[:,1:]

    new_dictionary = []
    for atom in dictionary:
        new_dictionary += [atom[1:]]
    dictionary = np.array(new_dictionary)


    # Convert to a binary matrix giving, for each atom in the dictionary, 
    # whether it is present at this row
    X_bydict_bin = [presence_vector(row, dictionary) for row in X]
    X_bydict_bin = np.array(X_bydict_bin)

    # Use Gaussian Naive Bayes on the data and retrieve accuracy
    bnb = GaussianNB()

    # Potential sample weights
    sample_weights = []
    for y in Y:
        if y == 1 : sample_weights += [QUERY_WEIGHT]
        else: sample_weights += [1]
    

    bnb.fit(X_bydict_bin, Y, sample_weight= sample_weights)
    Y_pred = bnb.predict(X_bydict_bin)

    # Compute F1-score
    score = sklearn.metrics.f1_score(Y, Y_pred)
    
    # print(pd.crosstab(Y,Y_pred))
    
    return -score # We want an error, higher means worse



## Run MODL with custom error function based on error of BNB
utils.VERBOSITY = 3 # Ensure DEBUG messages are shown


# Keep only combis WITH query (first element is not 0)
flags_matrix_with_query = flags_matrix[flags_matrix[:, 0] != 0]

combi_miner = Modl(flags_matrix_with_query,
    multiple_overlap_max_number_of_combinations = DESIRED_ITEMSETS,        # How many words to find ?
    error_function = custom_error_function_bnb)                            # Custom error function in step 2
interesting_combis = combi_miner.find_interesting_combinations()