"""
Compare the p-values estimated by OLOGRAM's Negative Binomial to deeply sampled
empirical p-values.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText

# root_dir = './output/supp_fig4/'
# concat_file = root_dir + 'merged_ologram_test.tsv'
# res_file_1 = root_dir + 'S_mean_boxplot.pdf'
# res_file_2 = root_dir + 'S_var_boxplot.pdf'

concat_file = snakemake.input["merged"]         # Concatenated shuffles
truth_file = snakemake.input["truth"]           # Deep sampled empirical p-value
res_file = snakemake.output[0]                  # Result file


# ---------------------------- Read result files ----------------------------- #
def result_ologram_to_list_of_lists(file, is_truth = False):
   
    # Skip even rows except for 0, they are the headers of other files (non-truth data only)
    def logic(index):
        if (index % 2 == 0) and index != 0 : return True
        return False 
    
    if not is_truth: df = pd.read_csv(file, sep='\t', header = 0, index_col = None, skiprows = lambda x: logic(x))
    else: df = pd.read_csv(file, sep='\t', header = 0, index_col = None)

    # Iterate over each row and build a result object
    result = list()
    for i, row in df.iterrows():
        header = row['feature_type'].split('_')

        if not is_truth:
            r = {
                'nb_minibatches': int(header[2]),
                'try_id': int(header[5]),
                'p_val': row['summed_bp_overlaps_pvalue'],
                'log_10(p_val)': np.log10(row['summed_bp_overlaps_pvalue'] + 1E-320),
                'log_10(empirical_pval)': np.log10(row["summed_bp_overlaps_empirical_pvalue"] + 1E-320)
            }
        else:
            r = {
                'p_val': row['summed_bp_overlaps_pvalue'],
                'log_10(p_val)': np.log10(row['summed_bp_overlaps_pvalue'] + 1E-320),
                'log_10(empirical_pval)': np.log10(row["summed_bp_overlaps_empirical_pvalue"] + 1E-320)
                #'log_10(beta_pval)': np.log10(row["ad_hoc_beta_summed_bp_overlaps_pvalue"] + 1E-320)
            }

        result += [r]

    return result


result = result_ologram_to_list_of_lists(concat_file)
result_truth = result_ologram_to_list_of_lists(truth_file, is_truth = True)

# ------------------------------- Plot figure -------------------------------- #
result_df = pd.DataFrame(result)

MINIBATCH_SIZE = 10 # Hardcoded for now, keep it the same in the Snakefile
result_df['Number of shuffles'] = MINIBATCH_SIZE * result_df['nb_minibatches']

meanplot = result_df.boxplot(column=['log_10(p_val)'], by='Number of shuffles')
result_truth_df = pd.DataFrame(result_truth)

# Add a line for the true p-val (as estimated by very deep sampling)
true_pval = list(result_truth_df['log_10(empirical_pval)'])[0] # Hardcoded to read the first line, since there should be only one
plt.axhline(y=true_pval, color='g', linestyle='-')

#beta_pval = list(result_truth_df['log_10(beta_pval)'])[0] # Hardcoded to read the first line, since there should be only one
#plt.axhline(y=beta_pval, color='r', linestyle='-')

at = AnchoredText("True p-val = "+'{0:.4g}'.format(10**true_pval),
                  prop=dict(size=8), frameon=True,
                  loc='upper right',
                  )
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
meanplot.add_artist(at)

plt.savefig(res_file)

