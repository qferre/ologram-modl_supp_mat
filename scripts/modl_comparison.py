#!/usr/bin/python
"""
Compare MODL and Apriori using artificial data
"""

import numpy as np
np.random.seed(42) # Seed !
import time

from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
from pygtftk.stats.intersect.modl.apriori import Apriori, matrix_to_list_of_transactions, apriori_results_to_matrix

# Query Snakemake info
result_file = snakemake.output["result"]
log_file = snakemake.output["log"]

# Redirect output to a log
import sys
sys.stdout = open(log_file, 'w')











## --------------------------- Found combinations --------------------------- ##










# Generate data with the AB, ABCD, EF combinations, adding 10% uniform noise
NB_SETS = 7
names = np.array([str(i) for i in range(NB_SETS)]) # Must be an array
np.random.seed(42)
x = test_data_for_modl(nflags = 10000, number_of_sets = NB_SETS,
    noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])

# Run Apriori
transactions = matrix_to_list_of_transactions(x, names)
myminer = Apriori(min_support = 0.05)
myminer.run_apriori(transactions)
results = myminer.produce_results()
"""
TODO : itemsets of 1 are returned, I don't want that !
"""
apriori_results_df = apriori_results_to_matrix(results, names)

# MODL
from pygtftk import utils
utils.VERBOSITY = 3 # Force debug messages
combi_miner = Modl(x, 
    multiple_overlap_target_combi_size = -1,    # Limit the size of the combinations
    multiple_overlap_max_number_of_combinations = 3,    # How many words to find ?
    nb_threads = 8,
    step_1_factor_allowance = 2)    # How many words to ask for in each step 1 rebuilding
modl_interesting_combis = combi_miner.find_interesting_combinations()







## Write to file
with open(result_file, 'a') as fo:

    # Print Apriori rules ordered by support
    fo.write(apriori_results_df.to_string())

    # Print MODL combinations in order of selection
    fo.write(str(modl_interesting_combis))

















## ---------------------------- Time benchmarks ----------------------------- ##

## Parameters to test
# Nb of sets
SIZES = [6,8,10,15,20,30,40,100, 120, 150, 200, 250, 300, 350, 400]
# Support/words steps
STEPS = [1,2,5,7,10,12,15,20,25,30,40,50,55,60]







# DEBUG TRY
SIZES = [6,20,200]

STEPS = [1,5,20]










utils.VERBOSITY = 0 # We don't want to record debug messages for these tests

import pandas as pd
from plotnine import *

import time


## Number of sets
df_bench = pd.DataFrame(columns = ['nb_sets','algo','time'])   # Prepare df

for size in SIZES:
    X = test_data_for_modl(nflags = 1000, number_of_sets = size, noise = 0.5)
    names = ['X' for _ in range(size)]
    transactions = matrix_to_list_of_transactions(X, names)

    # Apriori
    start_time = time.time()
    myminer = Apriori(min_support = 0.05)
    myminer.run_apriori(transactions)
    results = myminer.produce_results()
    stop_time = time.time()

    df_bench = df_bench.append({'nb_sets':size, 'algo':'apriori', 'time': stop_time-start_time}, ignore_index = True)   # Record the time

    # Modl
    start_time = time.time()
    combi_miner = Modl(x, multiple_overlap_max_number_of_combinations = 8, nb_threads = 4)
    modl_interesting_combis = combi_miner.find_interesting_combinations()
    stop_time = time.time()

    df_bench = df_bench.append({'nb_sets':size, 'algo':'modl', 'time': stop_time-start_time}, ignore_index = True)

# Plot
p = (ggplot(df_bench) + aes('parameter', 'time', color='algo', group='algo')
 + geom_point() + geom_line())
p.save(filename = snakemake.output["fig1"])



## Number of queried words and min support
df_bench = pd.DataFrame(columns = ['step','algo','time'])

X = test_data_for_modl(nflags = 100000, number_of_sets = 10, noise = 0.5)
names = ['X' for _ in range(10)]
transactions = matrix_to_list_of_transactions(X, names)

for step in STEPS:

    # Apriori
    start_time = time.time()
    myminer = Apriori(min_support = (1-step/100))
    myminer.run_apriori(transactions)
    results = myminer.produce_results()
    stop_time = time.time()

    df_bench = df_bench.append({'step':step, 'algo':'apriori', 'time': stop_time-start_time}, ignore_index = True)

    # Modl
    start_time = time.time()
    combi_miner = Modl(x, multiple_overlap_max_number_of_combinations = step, nb_threads = 4)
    modl_interesting_combis = combi_miner.find_interesting_combinations()
    stop_time = time.time()

    df_bench = df_bench.append({'step':step, 'algo':'modl', 'time': stop_time-start_time}, ignore_index = True)

# Plot
p = (ggplot(df_bench) + aes('step', 'time', color='algo', group='algo')
 + geom_point() + geom_line())
p.save(filename = snakemake.output["fig2"])








"""
Seconday X axis are hard, but maybe show SCALING SPEED
    MODL : from 1 queried word to 60
    apriori : from 100% min support to 0.4%
"""