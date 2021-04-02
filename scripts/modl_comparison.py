"""
Compare MODL to other algorithms using artificial data matrices.
"""

import numpy as np
np.random.seed(42)

import pandas as pd

import time
import subprocess

from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
from pygtftk.stats.intersect.modl.apriori import Apriori, matrix_to_list_of_transactions, apriori_results_to_matrix
from pygtftk.stats.intersect.modl.subroutines import learn_dictionary_and_encode
from pygtftk import utils

OUTPUT_ROOT = "output/benchmark/comparison" # Hardcoded for now. It was necessary to launch this script in shell and not via snakemake.script to allow proper redirection of the log

# Debug : plots are not shown, only printed, so use Agg backend for matplotlib
import matplotlib
matplotlib.use("Agg")


## --------------------------- Found combinations --------------------------- ##

# Generate data with the AB, ABCD, EF combinations, adding uniform noise
NB_SETS = 6
names = [str(i) for i in range(NB_SETS)]
x = test_data_for_modl(nflags = 10000, number_of_sets = NB_SETS,
    noise = 0.12, cor_groups = [(0,1),(0,1,2,3),(4,5)])

transactions = matrix_to_list_of_transactions(x, names) # Convert to list of transactions





# Run MODL
utils.VERBOSITY = 3 # Force debug messages to appear

combi_miner = Modl(x, 
    multiple_overlap_target_combi_size = -1,            # Optional : Limit the size of the combinations
    multiple_overlap_max_number_of_combinations = 3,    # How many words to find ?
    nb_threads = 8,                                     # Full multithreading
    smother = True,                                     # Reduce each row's abundance to its square root. Helps find rarer combinations but magnifies the noise.
    step_1_factor_allowance = 2,                        # Optional : How many words to ask for in each step 1 rebuilding
    normalize_words = True,                            # Normalize word sum of square in step 2
    step_2_alpha = None)                                # Override the sparsity control in step 2
modl_interesting_combis = combi_miner.find_interesting_combinations()


# Run MODL without smothering
x = test_data_for_modl(noise = 0.2)
combi_miner = Modl(x, multiple_overlap_max_number_of_combinations=3, smother = False)
modl_interesting_combis_no_smother = combi_miner.find_interesting_combinations()


# Run MODL without word normalization
x = test_data_for_modl(noise = 0.12)
combi_miner = Modl(x, multiple_overlap_max_number_of_combinations=3, normalize_words = False)
modl_interesting_combis_not_normalized = combi_miner.find_interesting_combinations()


print("-- MODL INTERESTING COMBINATIONS --")
print("Standard =", modl_interesting_combis)
print("No smother, more noise =", modl_interesting_combis_no_smother)
print("No word normalization =", modl_interesting_combis_not_normalized)
print("-------------------------------")


## Run Apriori and FP-growth
MIN_SUPPORT = 1E-10 # Epsilon for 0. Return all itemsets, ordered by support

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
fpgrowth(x_as_dataframe, min_support=MIN_SUPPORT)



# Here, output the matrices to BED files
# from matrix_to_bed import intersection_matrix_to_bedfiles
# bedpaths = intersection_matrix_to_bedfiles(x, output_root)

# # Also print it to a tsv, and to a list of transactions
# np.savetxt(OUTPUT_ROOT+'/artifmat.tsv', delimiter='\t')

# def format_list_brackets(alist): return '{'+','.join(alist)+'}'
# with open(OUTPUT_ROOT+'/artiftransact_brackets.tsv', 'w+') as f:
#     for item in transactions:
#         f.write(format_list_brackets(item)+'\n')




# And in the specific SPMF format
# TODO UNHARDCODE THIS
OUTPUT_ROOT = "../output/benchmark/comparison"

with open(OUTPUT_ROOT+'/artiftransact.txt', 'w+') as f:
    for item in transactions: f.write(' '.join(item)+'\n')



## --- Manual itemsets miners

"""
The miners will be in the /ext directory
"""

# TODO add SPMF source as citation in paper for those implementations

# TODO Run LCM
command_line = ["java -jar","./ext/spmf.jar", # Java SPMF toolset
    "run",
    "LCM", # Algorithm
    OUTPUT_ROOT+"/artiftransact.txt", # Query file
    OUTPUT_ROOT+"/output_lcm.txt",
    str(round(MIN_SUPPORT*100))+"%" # Min support
    ] 


# NOTE Remember that in subprocess.run, the command line must be a list of arguments and not a string of the command
stdout_captured = subprocess.run(command_line, stdout=subprocess.PIPE, universal_newlines=True).stdout
