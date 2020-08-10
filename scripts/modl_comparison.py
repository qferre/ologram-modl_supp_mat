#!/usr/bin/python
"""
Compare MODL and Apriori using artificial data
"""





import numpy as np

from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
from pygtftk.stats.intersect.modl.apriori import Apriori, matrix_to_list_of_transactions, apriori_results_to_matrix



# Query Snakemake info
result_file = snakemake.output[0]




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
results_df = apriori_results_to_matrix(results, names)

print(results_df)





# MODL
combi_miner = Modl(x, 
    multiple_overlap_target_combi_size = -1,    # Limit the size of the combinations
    multiple_overlap_max_number_of_combinations = 3,    # How many words to find ?
    nb_threads = 8,
    step_1_factor_allowance = 2)    # How many words to ask for in each step 1 rebuilding
interesting_combis = combi_miner.find_interesting_combinations()











# Print MODL combinations in order of selection


print(interesting_combis)

# Print Apriori rules ordered by support


# Compare the two
assert_compare(ap_result, modl_result)













# Write to file
with open(result_file, 'a') as fo:
    fo.write("Placeholder result\n")
