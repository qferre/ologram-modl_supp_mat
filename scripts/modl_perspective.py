"""
A perspective application of MODL using a supervised loss, based on the 
performance of a Naive Bayes classifier.


This currently is not stable enough and requires heavy tuning of the parameters,
so it can't be used yet in regular code, but is an interesting perspective.
"""

# Paths to BED files of regions
QUERY_PATH  = "input/as_ginom/query_som_trans_lung.bed"
EXCL_PATH   = "input/as_ginom/mappability_human.bed"
GENOME_PATH = "input/hg19.genome"
MORE_PATHS  = [
    "input/as_ginom/reference_1_hg19_Refseq_promoters_2000bp.bed",
    "input/as_ginom/reference_2_hg19_Refseq_gene_bodies.bed",
    "input/as_ginom/reference_3_Broad_A549_H3k04me3_Dex.bed",
    "input/as_ginom/reference_4_Broad_A549_H3k36me3_Dex.bed",
    "input/as_ginom/reference_5_Broad_A549_H3k79me2_Dex.bed"
]


# Subsampling parameters
KEEP_N_TIMES_QUERY = 30          # For each '1' line with the query, how many '0' lines without it do we keep ?
QUERY_WEIGHT = 30                # In the Naive Bayes classifier, how much more do we weigh the presence of the query ?


# MODL parameters
DESIRED_ITEMSETS = 7



# ---------------------------------------------------------------------------- #

def presence_vector(row, combis):
    """
    Given a list of combinations, returns a vector saying whether they are
    present or absent in this row. The order in the vector is the same as given
    in the combis list.

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
from pygtftk.stats.intersect.modl.dict_learning import Modl
from pygtftk import utils
from pygtftk.stats.intersect.read_bed import read_bed_as_list as read_bed
from pygtftk.utils import chrom_info_as_dict
from pygtftk import arg_formatter

np.random.seed(42)  # Random seed for reproducibility




## Prepare files
bedA = pybedtools.BedTool(QUERY_PATH).sort().merge()
bedsB = [pybedtools.BedTool(bedfilepath).sort().merge() for bedfilepath in MORE_PATHS]


# Do the exclusion manually Generate a fake bed for the entire genome, using the chromsizes
bed_incl = pybedtools.BedTool(EXCL_PATH)
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
flags_matrix = flags_matrix[:,1:] # Remove first column that represents the fake, full genome BED.



# Perform subsampling of the majority class ('0', meaning query is absent)
flags_matrix_with_query = flags_matrix[flags_matrix[:, 0] != 0]
draw_this_many = int(np.around(KEEP_N_TIMES_QUERY * len(flags_matrix_with_query)))
random_rows = flags_matrix[np.random.randint(flags_matrix.shape[0], size=draw_this_many),:]
flags_matrix = np.concatenate((flags_matrix_with_query,random_rows))



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