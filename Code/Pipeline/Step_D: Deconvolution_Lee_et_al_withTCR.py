import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from tqdm import tqdm



"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-18 by Saishi Cui
Program: Step_D: Deconvolution_Lee_et_al_withTCR.py
Purpose: Deconvolution of Lee et al., 2024 using built PBMC atlas with Cancer Reactive T cells.
--------------------------------------------------------------------------------
Data Inputs:
- PBMC atlas with Cancer Reactive T cells and with Granu obtained from Step_C: Build_PBMC_with_CancerReactive_Atlas.py under the folder "Pipeline".
- Lee et al., 2024 data from GSE243529_matrix-processed.tsv under the data folder "Lee_et_al".

Data Outputs:
- Deconvolution results (with TCR) under the data folder "Lee_et_al".

<<<<<<< HEAD
=======
Visualization analysis:
- Boxplot of deconvolution results (with TCR) under the data folder "Lee_et_al".

>>>>>>> bb3e26c (Initial commit)
Functions:
- constrained_ls(): Deconvolution of Lee et al., 2024 using built PBMC atlas with Cancer Reactive T cells.

Note:
- We use the PBMC atlas with Granu to deconvolute the Lee et al., 2024 data.
--------------------------------------------------------------------------------
"""

<<<<<<< HEAD
=======

>>>>>>> bb3e26c (Initial commit)
Atlas_final_df_withGranu = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_final_df_withGranu.csv')
Lee = pd.read_csv('/Users/scui2/DNAmethylation/Lee_et_al/GSE243529_matrix-processed.tsv', sep = '\t', header = 0)


# Preprocess the Lee et al., 2024 data
Lee = Lee.iloc[:, [0] + list(range(1,1049,2))]
Lee_BC = Lee.iloc[:,[0] + list(range(1,257))]
Lee_normal = Lee.iloc[:,[0] + list(range(257,525))]
Lee_clean = pd.concat([Lee_BC, Lee_normal.iloc[:,1:]], axis = 1)
<<<<<<< HEAD
=======

>>>>>>> bb3e26c (Initial commit)
# Keep the first instance of duplicate columns
Lee_clean = Lee_clean.loc[:, ~Lee_clean.columns.duplicated()]
Lee_clean.dropna(inplace = True)
Lee_clean.reset_index(drop = True, inplace = True)
intersect_probe = list(set(Atlas_final_df_withGranu.ProbeID).intersection(set(Lee_clean.ID_REF)))
Lee_atlas = Lee_clean[Lee_clean.ID_REF.isin(intersect_probe)]
Lee_atlas.sort_values(by = "ID_REF", ascending = True, inplace = True)
Lee_atlas.drop_duplicates(subset = "ID_REF", keep = "first", inplace = True)
Lee_atlas.reset_index(drop = True, inplace = True)

Atlas_final_df_withGranu = Atlas_final_df_withGranu[Atlas_final_df_withGranu.ProbeID.isin(intersect_probe)]
Atlas_final_df_withGranu.drop_duplicates(subset = "ProbeID", keep = "first", inplace = True)
Atlas_final_df_withGranu.sort_values(by = "ProbeID", ascending = True, inplace = True)
Atlas_final_df_withGranu.reset_index(drop = True, inplace = True)






X = Atlas_final_df_withGranu.iloc[:, 3:].values  # Target matrix
Y = Lee_atlas.iloc[:, 1:].values.T      # Mixed signal matrix


def constrained_ls(A, b):
    n = A.shape[1]
    
    # Objective function: Least squares
    def objective(x):
        return np.sum((A @ x - b) ** 2)
    
    # Constraints: Sum to 1
    constraints = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
    ]
    
    # Boundary conditions: Non-negative
    bounds = [(0, None) for _ in range(n)]
    
    x0 = (np.ones(n) / n)
    
    # Solve
    result = minimize(objective, x0, 
                     method='SLSQP',
                     bounds=bounds,
                     constraints=constraints)
    
    return result.x

# Initialize result list
nnls_results_Lee = []

# Deconvolution for each sample (add progress bar)
for i in tqdm(range(Y.shape[0]), desc="Processing samples"):
    y = Y[i, :]
    # Deconvolution using constrained least squares
    coeffs = constrained_ls(X, y)
    nnls_results_Lee.append(coeffs)

# Convert results to DataFrame
nnls_df_Lee = pd.DataFrame(nnls_results_Lee, columns=["CD4+T", "CD4+TCR", "CD8+T", "CD8+TCR", "B", "Mono", "NK", "Granu"])
nnls_df_Lee *= 100  # Convert to percentage


# Drop 'Granu' column
nnls_df_Lee_dropped = nnls_df_Lee.drop(columns=['Granu'])

# Re-normalize cell type proportions
nnls_df_Lee_dropped = nnls_df_Lee_dropped.div(
    nnls_df_Lee_dropped.sum(axis=1), axis=0
)

nnls_df_Lee_dropped *= 100







# Add hierarchical labels

labels = ['Breast Cancer'] * 256 + ['Normal'] * 268
nnls_df_Lee_dropped['Group'] = labels

# Convert data to long format for plotting
nnls_long_Lee = pd.melt(nnls_df_Lee_dropped, id_vars='Group', var_name='Cell Type', value_name='Proportion')



# Set light color palette for boxplot
boxplot_palette = ["lightcoral", "lightgreen"]
# Set dark color palette for stripplot
stripplot_palette = ["red", "green"]

# Create figure
plt.figure(figsize=(12, 6))

# Create boxplot, set showfliers=False to hide outliers
sns.boxplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Lee, 
            palette=boxplot_palette, width=0.7, showfliers=False)

# Add individual sample points, using dark color
sns.stripplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Lee,
              palette=stripplot_palette, dodge=True, size=4, alpha=0.6)

# Set labels
plt.xlabel('')
plt.ylabel('Proportion (%)', fontweight='bold')
plt.xticks(rotation=45, ha='right', fontweight='bold')
plt.title('Lee et al., 2024', fontweight='bold')
# Set legend
handles, _ = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[:2], labels=['Breast Cancer (N = 256)', 'Normal (N = 268)'], 
          prop={'weight': 'bold'}, loc = 'upper left')

# Remove top and right borders
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Remove grid lines
plt.grid(False)

# Adjust layout
plt.tight_layout()

# Display figure
plt.show()

nnls_df_Lee_TCR = nnls_df_Lee_dropped.copy()
nnls_df_Lee_TCR.to_csv('/Users/scui2/DNAmethylation/Lee_et_al/nnls_df_Lee_TCR.csv', index = False)
