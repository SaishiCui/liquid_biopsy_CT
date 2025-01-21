import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from tqdm import tqdm



"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-17 by Saishi Cui
Program: Step_D: Deconvolution_Li_et_al_withoutTCR.py
Purpose: Deconvolution of Li et al., 2024 using built PBMC atlas without Cancer Reactive T cells.
--------------------------------------------------------------------------------
Data Inputs:
- PBMC atlas without Cancer Reactive T cells obtained from Step_C: Build_PBMC_with_CancerReactive_Atlas.py under the folder "Pipeline".
- Li et al., 2024 data from Li_et_al.txt under the data folder "Li_et_al".

Data Outputs:
- Deconvolution results (without TCR) under the data folder "Li_et_al".

<<<<<<< HEAD
=======
Visualization analysis:
- Boxplot of deconvolution results (without TCR) under the data folder "Li_et_al".

>>>>>>> bb3e26c (Initial commit)
Functions:
- constrained_ls(): Deconvolution of Li et al., 2024 using built PBMC atlas without Cancer Reactive T cells.
--------------------------------------------------------------------------------
"""


Atlas_final_df = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_final_df.csv')
Li = pd.read_csv('/Users/scui2/DNAmethylation/Li_et_al/Li_et_al.txt',sep = ",")  

intersect_probe = list(set(Atlas_final_df.ProbeID).intersection(set(Li.probeID)))

Li_atlas = Li[Li.probeID.isin(intersect_probe)]
Li_atlas.sort_values(by = "probeID", ascending = True, inplace = True)
Li_atlas.drop_duplicates(subset = "probeID", keep = "first", inplace = True)
Li_atlas.reset_index(drop = True, inplace = True)

Atlas_final_df = Atlas_final_df[Atlas_final_df.ProbeID.isin(intersect_probe)]
Atlas_final_df.drop_duplicates(subset = "ProbeID", keep = "first", inplace = True)
Atlas_final_df.sort_values(by = "ProbeID", ascending = True, inplace = True)
Atlas_final_df.reset_index(drop = True, inplace = True)



X = Atlas_final_df.iloc[:, 3:].iloc[:,[0,2,4,5,6]].values  # Target matrix
Y = Li_atlas.iloc[:, 1:].iloc[:, :].values.T      # Mixed signal matrix

def constrained_ls(A, b):
    n = A.shape[1]
    
    # Objective function: least squares
    def objective(x):
        return np.sum((A @ x - b) ** 2)
    
    # Constraint: sum to 1
    constraints = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
    ]
    
    # Boundary condition: non-negative
    bounds = [(0, None) for _ in range(n)]
    
    # Initial guess: uniform distribution
    x0 = np.ones(n) / n
    
    # Solve
    result = minimize(objective, x0, 
                     method='SLSQP',
                     bounds=bounds,
                     constraints=constraints)
    
    return result.x

# Initialize result list
nnls_results_Li = []

# Deconvolution for each sample (add progress bar)
for i in tqdm(range(Y.shape[0]), desc="Processing samples"):
    y = Y[i, :]
    # Deconvolution using constrained least squares
    coeffs = constrained_ls(X, y)
    nnls_results_Li.append(coeffs)

# Convert results to DataFrame
nnls_df_Li = pd.DataFrame(nnls_results_Li, columns=["CD4+T", "CD8+T", "B", "Mono", "NK"])
nnls_df_Li *= 100  # Convert to percentage


# Add stratified labels
labels_Li = ['Lung Cancer'] * 35 + ['Normal'] * 50
nnls_df_Li['Group'] = labels_Li

# Convert data to long format for plotting
nnls_long_Li = pd.melt(nnls_df_Li, id_vars='Group', var_name='Cell Type', value_name='Proportion')

# Set light color palette for boxplot
boxplot_palette = ["lightcoral", "lightgreen"]
# Set dark color palette for stripplot
stripplot_palette = ["red", "green"]

# Create figure
plt.figure(figsize=(12, 6))

# Create boxplot, set showfliers=False to hide outliers
sns.boxplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Li, 
            palette=boxplot_palette, width=0.7, showfliers=False)

# Add individual sample points, using dark color palette
sns.stripplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Li,
              palette=stripplot_palette, dodge=True, size=4, alpha=0.6)

# Set labels
plt.xlabel('')
plt.ylabel('Proportion (%)', fontweight='bold')
plt.xticks(rotation=45, ha='right', fontweight='bold')
plt.title('Li et al., 2024', fontweight='bold')
# Set legend
handles, _ = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[:2], labels=['Lung Cancer (N = 35)', 'Normal (N = 50)'], 
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


nnls_df_Li_withoutTCR = nnls_df_Li.copy()
nnls_df_Li_withoutTCR.to_csv('/Users/scui2/DNAmethylation/Li_et_al/nnls_df_Li_withoutTCR.csv', index = False)
