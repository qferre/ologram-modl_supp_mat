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

from mlxtend.frequent_patterns import fpgrowth, apriori

# Debug : plots are not shown, only printed, so use Agg backend for matplotlib
import matplotlib
matplotlib.use("Agg")

OUTPUT_ROOT = "output/benchmark/scaling/" # Hardcoded for now. It was necessary to launch this script in shell and not via snakemake.script to allow proper redirection of the log

## ---------------------------- Parameters ---------------------------------- ##
utils.VERBOSITY = 0 # We don't want to record debug messages for these tests
REPEATS = range(5) # Repeat all operations N times to get the average


## Elementary operation (DL) vs other itemset miners

# Number of lines (in thousands)
LINES_NB = [1,2,5,8,10,15,20,25,40,50]
NB_SETS = 15  

# Data generation parameters
NOISE = 0.5

# DL parameters
ALPHA = 0
N_ATOMS = 120

# Other algo parameters
MIN_SUPPORT = 1/100


## -------------- Elementary operation benchmark : apriori vs fpgrowth vs dict learning

## Number of lines
df_bench = pd.DataFrame(columns = ['lines','algo','time'])   # Prepare df

for _ in REPEATS:

    for lines in LINES_NB:

        true_line_nb = lines * 1000

        # Generate data
        # Scaling factor is in the thouands
        X = test_data_for_modl(nflags = true_line_nb, number_of_sets = NB_SETS, noise = NOISE)
        names = [str(i) for i in range(X.shape[1])]
        transactions = matrix_to_list_of_transactions(X, names)
        X_as_dataframe = pd.DataFrame(X)

        # Apriori
        start_time = time.time()
        myminer = Apriori(min_support = MIN_SUPPORT)
        myminer.run_apriori(transactions)
        results = myminer.produce_results()
        stop_time = time.time()

        df_bench = df_bench.append({'lines':lines, 'algo':'apriori_pure_python', 'time': stop_time-start_time}, ignore_index = True)   


        # Apriori
        # NOTE: this implementation has high RAM cost for large data or low min_support, be careful with scaling
        if true_line_nb <= 20000:
            start_time = time.time()
            result = apriori(X_as_dataframe, min_support = MIN_SUPPORT)
            stop_time = time.time()
            df_bench = df_bench.append({'lines':lines, 'algo':'apriori', 'time': stop_time-start_time}, ignore_index = True)  
        

        # FP-Growth
        start_time = time.time()
        result = fpgrowth(X_as_dataframe, min_support = MIN_SUPPORT)
        stop_time = time.time()

        df_bench = df_bench.append({'lines':lines, 'algo':'fpgrowth', 'time': stop_time-start_time}, ignore_index = True)  

        
        # Dict learning
        start_time = time.time()
        U_df, V_df, error = learn_dictionary_and_encode(X, n_atoms = N_ATOMS, alpha = ALPHA, n_jobs = 1)
        stop_time = time.time()

        df_bench = df_bench.append({'lines':lines, 'algo':'DL', 'time': stop_time - start_time}, ignore_index = True)   


df_bench['lines'] = df_bench['lines'].astype(int)
p = (ggplot(df_bench) + aes('lines', 'time', color='algo', group='algo')
 + geom_point() + geom_smooth(method='gpr', span=.3) + scale_x_continuous()
 + xlab("Number of lines") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig5")

p = (ggplot(df_bench) + aes('lines', 'time', color='algo', group='algo')
 + geom_point() + geom_smooth(method='gpr', span=.3) + scale_x_continuous() + scale_y_log10()
 + xlab("Number of lines") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig5_log10")


# Normalized time to minimum line number
min_nb_lines = min(LINES_NB)
minimum_time = df_bench[df_bench['lines'] == min_nb_lines][['algo','time']]
df_bench['time_relative'] = df_bench['time'] # Placeholder
for index, row in df_bench.iterrows():
    my_algo = row['algo']
    my_minimum_time = pd.to_numeric(minimum_time.loc[minimum_time['algo'] == my_algo]['time']).min()
    df_bench.at[index,'time_relative'] = row['time']/my_minimum_time

p = (ggplot(df_bench) + aes('lines', 'time_relative', color='algo', group='algo')
 + geom_point() + geom_smooth(method='gpr', span=.3) + scale_x_continuous()
 + xlab("Number of lines") + ylab("Time (relative)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig5_relative")