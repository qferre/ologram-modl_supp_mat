"""
Adapted from the `ologram_modl_treeify` source code.
"""
import numpy as np
import pandas as pd
import pandas.api.types as pdtypes
from pygtftk.stats.intersect.modl import tree
from functools import partial
import math
from collections import Counter

from plotnine import ggplot, aes, labs, scale_color_gradient, geom_point, geom_abline, scale_x_log10, scale_y_log10, geom_violin, geom_boxplot, position_dodge, scale_y_log10, xlab, ylab

ROOT_PATH = "./output/ologram_result_scatacseq_pbmc/"
DATA_ROOT_PATH = "./output/sc_atac_seq_pbmc_data/"

# For manual execution
# ROOT_PATH = "../output/ologram_result_scatacseq_pbmc/"
# DATA_ROOT_PATH = "../output/sc_atac_seq_pbmc_data/"


## Superclusters: how to regroup the cell types
# For now, this is hardcoded
SUPERCLUSTERS = {
    "CD14+ Monocytes":  "cd14",
    "CD4 Naive":        "cd48",
    "CD8 Naive":        "cd48",
    "Quer":             "NA",
    "...":              "NA",
    "NA":               "NA",
    "pre-B cell":       "preB"
}


# ---------------------------------------------------------------------------- #

## Read the OLOGRAM-MODL result file
# Read dataframe to create the found_combis dictionary
df_res = pd.read_csv(ROOT_PATH+"merged_batches_result.tsv", sep='\t', header=0, index_col=None)
# Pval set to 0 or -1 are changed to 1e-320 and NaN respectively
df_res.loc[df_res['summed_bp_overlaps_pvalue'] == 0, 'summed_bp_overlaps_pvalue'] = 1e-320
df_res.loc[df_res['summed_bp_overlaps_pvalue'] == -1, 'summed_bp_overlaps_pvalue'] = np.nan

# Create a Library tree with the combinations, like in the MODL algorithm itself  
print("Creating a Library for these words. This can be very long for longer words.")
print("In this particular Snakemake, it should take about 12 minutes.")

L = tree.Library()
L.build_nodes_for_words_from_ologram_result_df(df_res)
L.assign_nodes()


# ---------------------------------------------------------------------------- #


# Remember the features names are in L.features_names
# And the interesting attributes of each node are: 
#   node.word, node.children, node.s pval and fc
def word_to_combi(word): return [L.features_names[i]  for i in range(len(word)) if word[i]]

# Retrieves all combis and enrichment
def work_on_this_node(node, graph):
    return word_to_combi(node.word), node.fc, node.pval, node.s

# Iterate over all nodes
mygraphfunc = partial(work_on_this_node, graph=L)
global_results = {}
tree.apply_recursively_to_all_nodes(L.root_node, mygraphfunc, global_results)

df = pd.DataFrame(global_results.values(), columns = ["combination","fc","pval","s"])

df['combi_length'] = (df['combination'].apply(len)) -2
# Substract 2 to the length to discard the "Query" and "..." terms of the combinations

# Now, get the enrichments of all 2-wise combis, all 3-wise, etc.
def entropy(data):
    """
    Compute Shannon entropy.
    """
    if len(data) <= 1: return 0

    # Data counts
    counts = Counter()
    for d in data: counts[d] += 1

    ent = 0
    probs = [float(c) / len(data) for c in counts.values()]
    for p in probs:
        if p > 0.: ent -= p * math.log(p, 2)

    return ent


# ---------------------------------------------------------------------------- #


# Read the translation table as a dictionary and convert
# If df_map is a pandas dataframe with two columns, 'name' and 'category'

df_map = pd.read_csv(DATA_ROOT_PATH + "cell_to_class.tsv", sep='\t', header=0, index_col=None)
df_map.set_index('id', inplace=True)
df_map.transpose()
dict_translate = df_map['class'].to_dict()

dict_translate["Query"] = "NA" # Add the 'Query' class
dict_translate["..."] = "NA"


def combination_to_combi_of_classes(combination):
    new_combination = []
    for element in combination:
        current_key = dict_translate[element]
        element_to_add = SUPERCLUSTERS[current_key]
        if element_to_add != "NA": 
            new_combination += [element_to_add]
    return new_combination

df['combination_classes'] = df['combination'].apply(combination_to_combi_of_classes)


## Compute entropy
df['entropy'] = df['combination_classes'].apply(entropy)
# Round to 2 decimals and make a category
df['entropy'] = df['entropy'].apply(lambda x: round(x,2)).astype(pdtypes.CategoricalDtype())


# Using only combis of each individual length
    
all_lengths = sorted(set(df['combi_length']))

# No point in going beyond, roughly, 12
all_lengths = [l for l in all_lengths if l <= 12]

for length in all_lengths:

    df_filtered = df.loc[df['combi_length'] == length,:]

    try:
        p = (ggplot(data=df_filtered, mapping = aes(x='entropy', y='fc'))
                + geom_violin(position = position_dodge(1), width = 1)
                + geom_boxplot(position = position_dodge(1), width = 0.25)
                + xlab("Entropy") + ylab("Fold change (log2)"))

        p.save(filename = ROOT_PATH + "entropy_graph/entropy_length_" + str(length) + "_fc.png")



        p = (ggplot(data=df_filtered, mapping = aes(x='entropy', y='s'))
                + geom_violin(position = position_dodge(1), width = 1)
                + geom_boxplot(position = position_dodge(1), width = 0.25)
                + xlab("Entropy") + ylab("True total overlapping bp."))

        p.save(filename = ROOT_PATH + "entropy_graph/entropy_length_" + str(length) + "_s.png")


    except:
        print("Skipping length",length,"which caused an error.")