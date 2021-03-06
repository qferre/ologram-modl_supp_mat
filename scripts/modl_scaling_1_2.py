"""
Run a scaling benchmark on MODL, dictionary learning and Apriori using artificial data matrices.
"""

import numpy as np
np.random.seed(42)

import pandas as pd
from plotnine import ggplot, aes, geom_point, geom_line, geom_smooth, scale_x_continuous, scale_y_log10, xlab, ylab

import time

from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
from pygtftk.stats.intersect.modl.apriori import Apriori, matrix_to_list_of_transactions, apriori_results_to_matrix
from pygtftk.stats.intersect.modl.subroutines import learn_dictionary_and_encode
from pygtftk import utils


# Debug : plots are not shown, only printed, so use Agg backend for matplotlib
import matplotlib
matplotlib.use("Agg")

OUTPUT_ROOT = "output/benchmark/scaling/" # Hardcoded for now. It was necessary to launch this script in shell and not via snakemake.script to allow proper redirection of the log


## ---------------------------- Parameters ---------------------------------- ##
utils.VERBOSITY = 0 # We don't want to record DEBUG messages for these tests
REPEATS = range(5) # Repeat all operations N times to get the average

N_THREADS = 4

## MODL
SIZES = [6,9,12,18,24,30,40,50]  # Numbers of sets (columns)
STEPS = [3,5,8,10,12,15,18,20,22,25,30] # Numbers of queried atoms

# When varying only one parameter, which is the default for the others ?
DEFAULT_QUERIED_ATOMS = 8
DEFAULT_N_SETS = 8

# Data parameters
NFLAGS = 10000
NOISE = 0.1     # Previously 0.5


## ---------------------------- Time benchmarks ----------------------------- ##

## Number of sets
df_bench = pd.DataFrame(columns = ['nb_sets','time'])   # Prepare dataframe

for _ in REPEATS:

    for size in SIZES:
        X = test_data_for_modl(nflags = NFLAGS, number_of_sets = size, noise = NOISE)

        start_time = time.time()
        combi_miner = Modl(X, multiple_overlap_max_number_of_combinations = DEFAULT_QUERIED_ATOMS, nb_threads = N_THREADS)
        modl_interesting_combis = combi_miner.find_interesting_combinations()
        stop_time = time.time()

        df_bench = df_bench.append({'nb_sets':size, 'time': stop_time-start_time}, ignore_index = True)

df_bench['nb_sets'] = df_bench['nb_sets'].astype(int)
p = (ggplot(df_bench) + aes('nb_sets', 'time')
 + geom_point() + geom_smooth(span=.3) + scale_x_continuous()
 + xlab("Number of sets") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig1")


## Number of queried words
df_bench = pd.DataFrame(columns = ['step','time'])

X = test_data_for_modl(nflags = NFLAGS, number_of_sets = DEFAULT_N_SETS, noise = NOISE)

for _ in REPEATS:
    for step in STEPS:

        start_time = time.time()
        combi_miner = Modl(X, multiple_overlap_max_number_of_combinations = step, nb_threads = N_THREADS)
        modl_interesting_combis = combi_miner.find_interesting_combinations()
        stop_time = time.time()

        df_bench = df_bench.append({'step':step, 'time': stop_time-start_time}, ignore_index = True)

df_bench['step'] = df_bench['step'].astype(int)
p = (ggplot(df_bench) + aes('step', 'time')
 + geom_point() + geom_smooth(span=.3) + scale_x_continuous()
 + xlab("Queried nb. of atoms") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig2")