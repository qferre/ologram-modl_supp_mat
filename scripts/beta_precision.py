"""
Evaluate the precision of method-of-moments fitting of a beta distribution.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

from scipy.stats import beta
from pygtftk.stats.beta import fit_beta

# Graph paths
res_file_a = snakemake.output["a"]   
res_file_b = snakemake.output["b"] 
res_file_min = snakemake.output["mini"] 
res_file_max = snakemake.output["maxi"] 
res_file_v = snakemake.output["v"] 

# Create directories if needed
for filepath in [res_file_a, res_file_b, res_file_v]:
    path = Path(filepath).parent # Get parent directory
    Path(path).mkdir(parents=True, exist_ok=True)


# Fix the parameters to be fiitted
a, b = 1., 2.

def beta_var(a,b) : return (a*b)/((a+b)**2 * (a+b+1)) 
true_variance = beta_var(a,b)

np.random.seed(seed=42)

# ---------------------------- Experiment ----------------------------- #
# Make many different fittings with many different samples sizes
K = 200
SAMPLE_SIZES = [100]*K + [200]*K + [1000]*K + [10000]*K + [20000]*K

result = list()

for size in SAMPLE_SIZES:
    obs = beta.rvs(a, b, size=size)  # Generate data
    ahat, bhat, mhat, chat = fit_beta(obs) # Fit

    # Record result
    result += [{
        'samples':size,
        'ahat':ahat,
        'bhat':bhat,
        'mhat':mhat,
        'chat':chat,
        'var':np.var(obs, ddof=1),
        'rebuilt_var':beta_var(ahat, bhat)
    }]


# ------------------------------- Plot figures ------------------------------- #
result_df = pd.DataFrame(result)
result_df['Number of samples'] = result_df['samples'] # Renaming


### For alpha

meanplot = result_df.boxplot(column=['ahat'], by='Number of samples')

# Add a line for the true parameter
plt.axhline(y=a, color='g', linestyle='-')

at = AnchoredText("True α = "+str(a),
                  prop=dict(size=8), frameon=True,
                  loc='upper right',
                  )
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
meanplot.add_artist(at)

plt.savefig(res_file_a)
plt.close()

### Same for beta

meanplot = result_df.boxplot(column=['bhat'], by='Number of samples')

plt.axhline(y=b, color='b', linestyle='-')
at = AnchoredText("True β = "+str(b), prop=dict(size=8), frameon=True, loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
meanplot.add_artist(at)

plt.savefig(res_file_b)
plt.close()


### And min and max

meanplot = result_df.boxplot(column=['mhat'], by='Number of samples')
plt.axhline(y=0, color='r', linestyle='-')
at = AnchoredText("True minimum = "+str(0), prop=dict(size=8), frameon=True, loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
meanplot.add_artist(at)
plt.savefig(res_file_min)
plt.close()

meanplot = result_df.boxplot(column=['chat'], by='Number of samples')
plt.axhline(y=1, color='r', linestyle='-')
at = AnchoredText("True maximum = "+str(1), prop=dict(size=8), frameon=True, loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
meanplot.add_artist(at)
plt.savefig(res_file_max)
plt.close()


### And for variance, to show Neg Binom fitting would thus be better

# Renaming
result_df['Empirical'] = result_df['var']
result_df['Fitted'] = result_df['rebuilt_var']

sns.set_theme(style="whitegrid")
dd = pd.melt(result_df,id_vars=['Number of samples'],value_vars=['Empirical','Fitted'],var_name='Variance')
vplot = sns.boxplot(x='Number of samples', y='value', data=dd, hue='Variance',
    linewidth=1)

vplot.get_figure().savefig(res_file_v)