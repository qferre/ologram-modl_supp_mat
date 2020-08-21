"""
Run a benchmark on MODL, dictionary learning and Apriori using artificial data matrices.
"""

import numpy as np
np.random.seed(42)

import pandas as pd
from plotnine import ggplot, aes, geom_point, geom_line, scale_x_continuous

import time

from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
from pygtftk.stats.intersect.modl.apriori import Apriori, matrix_to_list_of_transactions, apriori_results_to_matrix
from pygtftk.stats.intersect.modl.subroutines import learn_dictionary_and_encode
from pygtftk import utils

OUTPUT_ROOT = "output/benchmark/" # Hardcoded for now. It was necessary to launch this script in shell and not via snakemake.script to allow proper redirection of the log


## --------------------------- Found combinations --------------------------- ##

# Generate data with the AB, ABCD, EF combinations, adding 10% uniform noise
NB_SETS = 6
names = [str(i) for i in range(NB_SETS)]
x = test_data_for_modl(nflags = 10000, number_of_sets = NB_SETS,
    noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])

# Run Apriori
transactions = matrix_to_list_of_transactions(x, names)
myminer = Apriori(min_support = 0)
myminer.run_apriori(transactions)
results = myminer.produce_results()
apriori_results_df = apriori_results_to_matrix(results, names)

print("-- APRIORI ASSOCIATION RULES --")
print(apriori_results_df)
print("-------------------------------")

# Run MODL
utils.VERBOSITY = 3 # Force debug messages to appear

combi_miner = Modl(x, 
    multiple_overlap_target_combi_size = -1,            # Optional : Limit the size of the combinations
    multiple_overlap_max_number_of_combinations = 3,    # How many words to find ?
    nb_threads = 8,                                     # Full multithreading
    step_1_factor_allowance = 2)                        # Optional : How many words to ask for in each step 1 rebuilding
modl_interesting_combis = combi_miner.find_interesting_combinations()


print("-- MODL INTERESTING COMBINATIONS --")
print(modl_interesting_combis)
print("-------------------------------")




## ---------------------------- Time benchmarks ----------------------------- ##
utils.VERBOSITY = 0 # We don't want to record debug messages for these tests

SIZES = [6,9,12,25,50]

## Number of sets
df_bench = pd.DataFrame(columns = ['nb_sets','time'])   # Prepare df

for size in SIZES:
    X = test_data_for_modl(nflags = 20000, number_of_sets = size, noise = 0.5)

    start_time = time.time()
    combi_miner = Modl(X, multiple_overlap_max_number_of_combinations = 8, nb_threads = 8)
    modl_interesting_combis = combi_miner.find_interesting_combinations()
    stop_time = time.time()

    df_bench = df_bench.append({'nb_sets':size, 'time': stop_time-start_time}, ignore_index = True)

df_bench['nb_sets'] = df_bench['nb_sets'].astype(int)
p = (ggplot(df_bench) + aes('nb_sets', 'time', color='algo', group='algo')
 + geom_point() + geom_line() + scale_x_continuous())
p.save(filename = OUTPUT_ROOT + "fig1")


## Number of queried words
df_bench = pd.DataFrame(columns = ['step','algo','time'])

STEPS = [3,5,8,10,15,20,25]

X = test_data_for_modl(nflags = 20000, number_of_sets = 8, noise = 0.5)

for step in STEPS:

    start_time = time.time()
    combi_miner = Modl(x, multiple_overlap_max_number_of_combinations = step, nb_threads = 8)
    modl_interesting_combis = combi_miner.find_interesting_combinations()
    stop_time = time.time()

    df_bench = df_bench.append({'step':step, 'algo':'modl', 'time': stop_time-start_time}, ignore_index = True)

df_bench['step'] = df_bench['step'].astype(int)
p = (ggplot(df_bench) + aes('step', 'time', color='algo', group='algo')
 + geom_point() + geom_line() + scale_x_continuous())
p.save(filename = OUTPUT_ROOT + "fig2")



## -------------- Elementary operation benchmark : apriori vs dict learning

# Support and number of queried words
SCALING_FACTORS = [2,5,10,25,40,75,100,150,200]

df_bench = pd.DataFrame(columns = ['scaling_factor','algo','time'])   # Prepare df

NOISE = 0.5
ALPHA = 0.5

for k in SCALING_FACTORS:

    # NOTE Scaling factor of k means :
    # - apriori will have min support of 1/k
    # - DL will learn 2*k words

    X = test_data_for_modl(nflags = 25000, number_of_sets = 13, noise = NOISE)
    names = [str(i) for i in range(X.shape[1])]
    transactions = matrix_to_list_of_transactions(X, names)

    # Apriori
    start_time = time.time()
    myminer = Apriori(min_support = 1/k)
    myminer.run_apriori(transactions)
    results = myminer.produce_results()
    stop_time = time.time()

    df_bench = df_bench.append({'scaling_factor':k, 'algo':'apriori', 'time': stop_time-start_time}, ignore_index = True)   

    # Dict learning
    start_time = time.time()
    U_df, V_df, error = learn_dictionary_and_encode(X, n_atoms = 2*k, alpha = ALPHA, n_jobs = 1)
    stop_time = time.time()

    df_bench = df_bench.append({'scaling_factor':k, 'algo':'DL', 'time': stop_time - start_time}, ignore_index = True)   


df_bench['scaling_factor'] = df_bench['scaling_factor'].astype(int)
p = (ggplot(df_bench) + aes('scaling_factor', 'time', color='algo', group='algo')
 + geom_point() + geom_line() + scale_x_continuous())
p.save(filename = OUTPUT_ROOT + "fig3")





# Number of sets

SETS_NB = [6,8,10,12,14]

df_bench = pd.DataFrame(columns = ['set_nb','algo','time'])   # Prepare df

for size in SETS_NB:

    # Generate data
    X = test_data_for_modl(nflags = 100000, number_of_sets = size, noise = NOISE)
    names = [str(i) for i in range(X.shape[1])]
    transactions = matrix_to_list_of_transactions(X, names)

    # Apriori
    start_time = time.time()
    myminer = Apriori(min_support = 1/100)
    myminer.run_apriori(transactions)
    results = myminer.produce_results()
    stop_time = time.time()

    df_bench = df_bench.append({'set_nb':size, 'algo':'apriori', 'time': stop_time-start_time}, ignore_index = True)   

    # Dict learning
    start_time = time.time()
    U_df, V_df, error = learn_dictionary_and_encode(X, n_atoms = 60, alpha = ALPHA, n_jobs = 1)
    stop_time = time.time()

    df_bench = df_bench.append({'set_nb':size, 'algo':'DL', 'time': stop_time - start_time}, ignore_index = True)   

df_bench['set_nb'] = df_bench['set_nb'].astype(int)
p = (ggplot(df_bench) + aes('set_nb', 'time', color='algo', group='algo')
 + geom_point() + geom_line() + scale_x_continuous())
p.save(filename = OUTPUT_ROOT + "fig4")