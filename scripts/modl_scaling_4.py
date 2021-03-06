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

# Number of sets (columns), minimum of 6
SETS_NB = [6,7,8,9,10,11,12,13,14,16,18,20,22,24,26]    

# Data generation parameters
NOISE = 0.5
LOW_NOISE = 0.01
N_FLAGS = 5000

# DL parameters
ALPHA = 0
N_ATOMS = 120

# Other algorithms' parameters
MIN_SUPPORT = 1/100 

## -------------- Elementary operation benchmark : apriori vs fpgrowth vs dict learning



## Number of sets
df_bench = pd.DataFrame(columns = ['set_nb','algo','time'])   # Prepare df

for _ in REPEATS:
        
    for size in SETS_NB:

        # Generate data
        X = test_data_for_modl(nflags = N_FLAGS, number_of_sets = size, noise = NOISE)
        names = [str(i) for i in range(X.shape[1])]
        transactions = matrix_to_list_of_transactions(X, names)
        X_as_dataframe = pd.DataFrame(X)

        X_low_noise = test_data_for_modl(nflags = N_FLAGS, number_of_sets = size, noise = LOW_NOISE)
        X_noiseless = test_data_for_modl(nflags = N_FLAGS, number_of_sets = size, noise = 0)

        # Apriori
        # Cap it to the max number that is not unreasonable
        if size <= 15:
            start_time = time.time()
            myminer = Apriori(min_support = MIN_SUPPORT)
            myminer.run_apriori(transactions)
            results = myminer.produce_results()
            stop_time = time.time()

            df_bench = df_bench.append({'set_nb':size, 'algo':'apriori_pure_python', 'time': stop_time-start_time}, ignore_index = True)   


        # Apriori
        # NOTE: this implementation has high RAM cost for large data or low min_support, be careful
        if size <= 15:
            start_time = time.time()
            result = apriori(X_as_dataframe, min_support = MIN_SUPPORT)
            stop_time = time.time()
            df_bench = df_bench.append({'set_nb':size, 'algo':'apriori', 'time': stop_time-start_time}, ignore_index = True)  


        # FP-Growth
        start_time = time.time()
        result = fpgrowth(X_as_dataframe, min_support = MIN_SUPPORT)
        stop_time = time.time()

        df_bench = df_bench.append({'set_nb':size, 'algo':'fpgrowth', 'time': stop_time-start_time}, ignore_index = True)  

        
        # Dict learning
        start_time = time.time()
        U_df, V_df, error = learn_dictionary_and_encode(X, n_atoms = N_ATOMS, alpha = ALPHA, n_jobs = 1)
        stop_time = time.time()

        df_bench = df_bench.append({'set_nb':size, 'algo':'DL', 'time': stop_time - start_time}, ignore_index = True)   


        # Dict learning - low noise data
        start_time = time.time()
        U_df, V_df, error = learn_dictionary_and_encode(X_low_noise, n_atoms = N_ATOMS, alpha = ALPHA, n_jobs = 1)
        stop_time = time.time()

        df_bench = df_bench.append({'set_nb':size, 'algo':'DL_low_noise_data', 'time': stop_time - start_time}, ignore_index = True)  


        # Dict learning - noiseless data
        start_time = time.time()
        U_df, V_df, error = learn_dictionary_and_encode(X_noiseless, n_atoms = N_ATOMS, alpha = ALPHA, n_jobs = 1)
        stop_time = time.time()

        df_bench = df_bench.append({'set_nb':size, 'algo':'DL_noiseless_data', 'time': stop_time - start_time}, ignore_index = True)  



df_bench['set_nb'] = df_bench['set_nb'].astype(int)
p = (ggplot(df_bench) + aes('set_nb', 'time', color='algo', group='algo')
 + geom_point() + geom_smooth(method='gpr', span=.3) + scale_x_continuous()
 + xlab("Number of sets") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig4")

p = (ggplot(df_bench) + aes('set_nb', 'time', color='algo', group='algo')
 + geom_point() + geom_smooth(method='gpr', span=.3) + scale_x_continuous() + scale_y_log10()
 + xlab("Number of sets") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig4_log10")


# Normalized time to minimum number of sets
min_nb_sets = min(SETS_NB)
minimum_time = df_bench[df_bench['set_nb'] == min_nb_sets][['algo','time']]
df_bench['time_relative'] = df_bench['time'] # Placeholder
for index, row in df_bench.iterrows():
    my_algo = row['algo']
    my_minimum_time = pd.to_numeric(minimum_time.loc[minimum_time['algo'] == my_algo]['time']).min()
    df_bench.at[index,'time_relative'] = row['time']/my_minimum_time

p = (ggplot(df_bench) + aes('set_nb', 'time_relative', color='algo', group='algo')
 + geom_point() + geom_smooth(method='gpr', span=.3) + scale_x_continuous()
 + xlab("Number of sets") + ylab("Time (relative)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig4_relative")


p = (ggplot(df_bench) + aes('set_nb', 'time_relative', color='algo', group='algo')
 + geom_point() + geom_smooth(method='gpr', span=.3) + scale_x_continuous() + scale_y_log10()
 + xlab("Number of sets") + ylab("Time (relative)"))
p.save(filename = OUTPUT_ROOT + "scaling_fig4_relative_log10")