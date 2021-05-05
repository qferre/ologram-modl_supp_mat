"""
Compare MODL to other algorithms using artificial data matrices.
"""

import numpy as np
import pandas as pd

import time
import subprocess
import warnings
import random

from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
from pygtftk.stats.intersect.modl.apriori import Apriori, matrix_to_list_of_transactions, apriori_results_to_matrix
from pygtftk.stats.intersect.modl.subroutines import learn_dictionary_and_encode
from pygtftk import utils

## ---------------------------- Parameters ---------------------------------- ##

np.random.seed(42)
random.seed(42)

utils.VERBOSITY = 3 # Force debug messages to appear
N_THREADS = 4

## Paths
# Hardcoded for now. It was necessary to launch this script in 
# shell and not via snakemake.script to allow proper redirection of the log
OUTPUT_ROOT = "output/benchmark/comparison/" # Output path
EXT_PATH = "ext/"   # Path to external tools

# # For manual launch
# OUTPUT_ROOT = "../output/benchmark/comparison/" # Output path
# EXT_PATH = "../ext/"   # Path to external tools

# Debug : plots are not shown, only printed, so use Agg backend for matplotlib
import matplotlib
matplotlib.use("Agg")


# Data parameters
NB_SETS = 6
NFLAGS = 10000
NOISE = 0.12
NOISE_HIGH = 0.2
CORR_GROUPS = [(0,1),(0,1,2,3),(4,5)]

# MODL parameters
N_QUERIED_ATOMS = 3

# Other algorithms' parameters
MIN_SUPPORT = 1E-10 # This value is an epsilon for 0. Hence, return all itemsets ordered by support



## --------------------------- Found combinations --------------------------- ##

# Generate data with the AB, ABCD, EF combinations, adding uniform noise
names = [str(i) for i in range(NB_SETS)]
x = test_data_for_modl(nflags = NFLAGS, number_of_sets = NB_SETS,
    noise = NOISE, cor_groups = CORR_GROUPS)

# Convert to list of transactions
transactions = matrix_to_list_of_transactions(x, names) 




## Run MODL

combi_miner = Modl(x, 
    multiple_overlap_target_combi_size = -1,                          # Optional: limit the size of the combinations
    multiple_overlap_max_number_of_combinations = N_QUERIED_ATOMS,    # How many words to find ?
    nb_threads = N_THREADS,                                           # Full multithreading
    smother = True,                                                   # Reduce each row's abundance to its square root. Helps find rarer combinations but magnifies the noise.
    step_1_factor_allowance = 2,                                      # Optional: how many words to ask for in each step 1 rebuilding
    normalize_words = True,                                           # Normalize word sum of square in step 2
    step_2_alpha = None)                                              # Optional: override the sparsity control in step 2
modl_interesting_combis = combi_miner.find_interesting_combinations()


# Run MODL without smothering
x = test_data_for_modl(noise = NOISE_HIGH)
combi_miner = Modl(x, 
    multiple_overlap_max_number_of_combinations=N_QUERIED_ATOMS, 
    smother = False)
modl_interesting_combis_no_smother = combi_miner.find_interesting_combinations()


# Run MODL without word normalization
x = test_data_for_modl(noise = NOISE)
combi_miner = Modl(x, 
    multiple_overlap_max_number_of_combinations=N_QUERIED_ATOMS, 
    normalize_words = False)
modl_interesting_combis_not_normalized = combi_miner.find_interesting_combinations()


print("-- MODL INTERESTING COMBINATIONS --")
print("Standard =", modl_interesting_combis)
print("No smother, more noise =", modl_interesting_combis_no_smother)
print("No word normalization =", modl_interesting_combis_not_normalized)
print("-------------------------------")


## Run Apriori and FP-growth

# Apriori
myminer = Apriori(min_support = MIN_SUPPORT)
myminer.run_apriori(transactions)
results = myminer.produce_results()
apriori_results_df = apriori_results_to_matrix(results, names)

print("-- APRIORI ASSOCIATION RULES --")
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(apriori_results_df)
print("-------------------------------")



# FP-Growth
from mlxtend.frequent_patterns import fpgrowth
x_as_dataframe = pd.DataFrame(x)
fpres = fpgrowth(x_as_dataframe, min_support=MIN_SUPPORT)



pd.set_option('display.max_rows', 1000)
print(fpres)


# Here, output the matrices to BED files
# from matrix_to_bed import intersection_matrix_to_bedfiles
# bedpaths = intersection_matrix_to_bedfiles(x, output_root)

# # Also print it to a tsv, and to a list of transactions
# np.savetxt(OUTPUT_ROOT+'artifmat.tsv', delimiter='\t')

# def format_list_brackets(alist): return '{'+','.join(alist)+'}'
# with open(OUTPUT_ROOT+'artiftransact_brackets.tsv', 'w+') as f:
#     for item in transactions:
#         f.write(format_list_brackets(item)+'\n')


# And in the specific SPMF format
with open(OUTPUT_ROOT+'artiftransact.txt', 'w+') as f:
    for item in transactions: f.write(' '.join(item)+'\n')



## ---------- Manual launch of other itemset miner(s)

## Run LCM
command_line = ["java", "-jar", EXT_PATH + "spmf.jar", # Java SPMF toolset
    "run", "LCM",                        # Run this algorithm
    OUTPUT_ROOT + "artiftransact.txt",   # Query file
    OUTPUT_ROOT + "output_lcm.txt",      # Output
    str(round(MIN_SUPPORT*100))+"%"      # Min support (parameter)
    ] 


# NOTE Remember that in subprocess.run, the command line must be a list of 
# arguments and not a string of the command
stdout_captured = subprocess.run(command_line, 
    stdout=subprocess.PIPE, universal_newlines=True).stdout





## Run CL-Max (approximate itemset miner)

warnings.filterwarnings('ignore') # Silence NumPy warnings for display

# Load it by running the code (kinda hacky but works for now)
exec(open(EXT_PATH + "cl-max.py").read())



# Run a variety of parameters
# Esp. rounding threshold. And random seeds seemed to change results a lot !
RANDOM_SEED = 1234
min_support = 0.1
for num_cluster in [3,6]:
    for rounding_threshold in [0.5,0.9]:

        print("-----------------")
        print("num_cluster =",num_cluster)
        print("rounding_threshold =",rounding_threshold)
            
        np.random.seed(RANDOM_SEED)
        random.seed(RANDOM_SEED)

        cl_max = CL_MAX(min_support, num_cluster, rounding_threshold)

        # Load dataset
        cl_max.dataset = pd.DataFrame(x)
        cl_max.sorted_items = names
        cl_max.transactions = [[int(i) for i in t] for t in transactions if t] # Remove empties
        cl_max.maximum = NB_SETS
        cl_max.remove_non_frequent_single_items()

        # Run the algorithm
        MFIs = cl_max._CL_MAX()