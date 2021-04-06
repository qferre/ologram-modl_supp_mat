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

# Elementary operation (DL)
# Number of words, and 1/min_support for the comparisons
SCALING_FACTORS =  [1,2,5,10,25,30,40,50,75,100,150,200,300,500]    

NOISE = 0.5
ALPHA = 0.5

## -------------- Elementary operation benchmark : apriori vs fpgrowth vs dict learning

## Support and number of queried words
df_bench = pd.DataFrame(columns = ['scaling_factor','algo','time'])   # Prepare df

for _ in REPEATS:

    for k in SCALING_FACTORS:

        # NOTE Scaling factor of k means:
        # - Apriori and FP-growth will have min support of 1/k
        # - Dictionary learning will learn k atoms
        current_min_support = 1/k

        X = test_data_for_modl(nflags = 25000, number_of_sets = 13, noise = NOISE)
        names = [str(i) for i in range(X.shape[1])]
        transactions = matrix_to_list_of_transactions(X, names)
        X_as_dataframe = pd.DataFrame(X)

        # Apriori - my pure Python implementation
        start_time = time.time()
        myminer = Apriori(min_support = current_min_support)
        myminer.run_apriori(transactions)
        results = myminer.produce_results()
        stop_time = time.time()

        df_bench = df_bench.append({'scaling_factor':k, 'algo':'apriori_pure_python', 'time': stop_time-start_time}, ignore_index = True)  


        # Apriori 
        # It seems this particular apriori implementation reserves a very large amount of memory when support is very low, so cap it.
        if current_min_support > 0.05: 
            start_time = time.time()
            apriori(X_as_dataframe, min_support = current_min_support)
            stop_time = time.time()

            df_bench = df_bench.append({'scaling_factor':k, 'algo':'apriori', 'time': stop_time-start_time}, ignore_index = True)  

        # FP-Growth
        start_time = time.time()
        result = fpgrowth(X_as_dataframe, min_support=current_min_support)
        stop_time = time.time()

        df_bench = df_bench.append({'scaling_factor':k, 'algo':'fpgrowth', 'time': stop_time-start_time}, ignore_index = True)  

        # Dict learning
        start_time = time.time()
        U_df, V_df, error = learn_dictionary_and_encode(X, n_atoms = k, alpha = ALPHA, n_jobs = 1)
        stop_time = time.time()

        df_bench = df_bench.append({'scaling_factor':k, 'algo':'DL', 'time': stop_time - start_time}, ignore_index = True)   


df_bench['scaling_factor'] = df_bench['scaling_factor'].astype(int)
p = (ggplot(df_bench) + aes('scaling_factor', 'time', color='algo', group='algo')
 + geom_point() + geom_smooth() + scale_x_continuous()
 + xlab("Scaling factor (k)") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "fig3")

p = (ggplot(df_bench) + aes('scaling_factor', 'time', color='algo', group='algo')
 + geom_point() + geom_smooth() + scale_x_continuous() + scale_y_log10()
 + xlab("Scaling factor (k)") + ylab("Time (seconds)"))
p.save(filename = OUTPUT_ROOT + "fig3_log10")


# Normalized time to scaling factor of 1
minimum_time = df_bench[df_bench['scaling_factor']==1][['algo','time']]
df_bench['time_relative'] = df_bench['time'] # Placeholder
for index, row in df_bench.iterrows():
    my_algo = row['algo']
    my_minimum_time = pd.to_numeric(minimum_time.loc[minimum_time['algo'] == my_algo]['time']).min()
    df_bench.at[index,'time_relative'] = row['time']/my_minimum_time

p = (ggplot(df_bench) + aes('scaling_factor', 'time_relative', color='algo', group='algo')
 + geom_point() + geom_smooth() + scale_x_continuous()
 + xlab("Scaling factor (k)") + ylab("Time (relative)"))
p.save(filename = OUTPUT_ROOT + "fig3_relative")