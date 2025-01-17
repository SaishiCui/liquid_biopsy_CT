import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from tqdm import tqdm



"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-17 by Saishi Cui
Program: Step_D: Deconvolution_Wang_et_al_withTCR.py
Purpose: Deconvolution of Wang et al., 2023 using built PBMC atlas with Cancer Reactive T cells.
--------------------------------------------------------------------------------
Data Inputs:
- PBMC atlas with Cancer Reactive T cells obtained from Step_C: Build_PBMC_with_CancerReactive_Atlas.py under the folder "Pipeline".
- Wang et al., 2023 data from Wang_Meth_data_850K.csv under the data folder "Wang_et_al".

Data Outputs:
- Deconvolution results (with TCR) under the data folder "Wang_et_al".

Functions:
- constrained_ls(): Deconvolution of Wang et al., 2023 using built PBMC atlas with Cancer Reactive T cells.
--------------------------------------------------------------------------------
"""


Atlas_final_df = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_final.csv')
Wang = pd.read_csv('/Users/scui2/DNAmethylation/Wang_et_al/Wang_Meth_data_850K.csv')

# Find intersected probes
intersect_probe = list(set(Atlas_final_df.ProbeID).intersection(set(Wang.cg_ID)))
Wang_atlas = Wang[Wang.cg_ID.isin(intersect_probe)]
Wang_atlas.sort_values(by = "cg_ID", ascending = True, inplace = True)
Wang_atlas.drop_duplicates(subset = "cg_ID", keep = "first", inplace = True)
Wang_atlas.reset_index(drop = True, inplace = True)

Atlas_final_df = Atlas_final_df[Atlas_final_df.ProbeID.isin(intersect_probe)]
Atlas_final_df.drop_duplicates(subset = "ProbeID", keep = "first", inplace = True)
Atlas_final_df.sort_values(by = "ProbeID", ascending = True, inplace = True)
Atlas_final_df.reset_index(drop = True, inplace = True)


X = Atlas_final_df.iloc[:, 3:].iloc[:,].values  # Target matrix
Y = Wang_atlas.iloc[:, 4:].values.T      # Mixed signal matrix



def constrained_ls(A, b):
    n = A.shape[1]
    
    # Objective function: least squares
    def objective(x):
        return np.sum((A @ x - b) ** 2)
    
    # Constraints: sum to 1
    constraints = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
    ]
    
    # Boundaries: non-negative
    bounds = [(0, None) for _ in range(n)]
    
    x0 = (np.ones(n) / n)
    
    # Solve
    result = minimize(objective, x0, 
                     method='SLSQP',
                     bounds=bounds,
                     constraints=constraints)
    
    return result.x

# Initialize result list
nnls_results_Wang = []

# Deconvolution for each sample (with progress bar)
for i in tqdm(range(Y.shape[0]), desc="Processing samples"):
    y = Y[i, :]
    # Deconvolution using constrained least squares
    coeffs = constrained_ls(X, y)
    nnls_results_Wang.append(coeffs)

# Convert results to DataFrame
nnls_df_Wang = pd.DataFrame(nnls_results_Wang, columns=["CD4+T", "CD4+TCR", "CD8+T", "CD8+TCR", "B", "Mono", "NK"])
nnls_df_Wang *= 100  # Convert to percentage



# Add hierarchical labels
labels_Wang = ['Breast Cancer'] * 50 + ['Normal'] * 30
nnls_df_Wang['Group'] = labels_Wang

# Convert data to long format for plotting
nnls_long_Wang = pd.melt(nnls_df_Wang, id_vars='Group', var_name='Cell Type', value_name='Proportion')

# Set light color palette for boxplot
boxplot_palette = ["lightcoral", "lightgreen"]
# Set dark color palette for stripplot
stripplot_palette = ["red", "green"]

# Create figure
plt.figure(figsize=(12, 6))

# Create boxplot, set showfliers=False to hide outliers
sns.boxplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Wang, 
            palette=boxplot_palette, width=0.7, showfliers=False)

# Add individual sample points, using dark color palette
sns.stripplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Wang,
              palette=stripplot_palette, dodge=True, size=4, alpha=0.6)

# Set labels
plt.xlabel('')
plt.ylabel('Proportion (%)', fontweight='bold')
plt.xticks(rotation=45, ha='right', fontweight='bold')
plt.title('Wang et al., 2023', fontweight='bold')
# Set legend
handles, _ = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[:2], labels=['Breast Cancer (N = 50)', 'Normal (N = 30)'], 
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

nnls_df_Wang_TCR = nnls_df_Wang.copy()
nnls_df_Wang_TCR.to_csv('/Users/scui2/DNAmethylation/Wang_et_al/nnls_df_Wang_TCR.csv', index = False)


