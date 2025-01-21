import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from tqdm import tqdm

"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-17 by Saishi Cui
Program: Step_D: Deconvolution_Xie_et_al_withTCR.py
Purpose: Deconvolution of Xie et al., 2023 using built PBMC atlas with Cancer Reactive T cells.
--------------------------------------------------------------------------------
Data Inputs:
- PBMC atlas with Cancer Reactive T cells obtained from Step_C: Build_PBMC_with_CancerReactive_Atlas.py under the folder "Pipeline".
- Xie et al., 2023 data from Xie_et_al.csv under the data folder "Xie_et_al".

Data Outputs:
- Deconvolution results (with TCR) under the data folder "Xie_et_al".

<<<<<<< HEAD
=======
Visualization analysis:
- Boxplot of deconvolution results (with TCR) under the data folder "Xie_et_al".

>>>>>>> bb3e26c (Initial commit)
Functions:
- constrained_ls(): Deconvolution of Xie et al., 2023 using built PBMC atlas with Cancer Reactive T cells.
--------------------------------------------------------------------------------
"""

# Read data
<<<<<<< HEAD
Atlas_final_df = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_final.csv')
=======
Atlas_final_df = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_final_df.csv')
>>>>>>> bb3e26c (Initial commit)
Xie = pd.read_csv('/Users/scui2/DNAmethylation/Xie_et_al/Xie_et_al.csv')  ### 728,310

# Find intersected probes
intersect_probe = list(set(Atlas_final_df.ProbeID).intersection(set(Xie.probeID)))
Xie_atlas = Xie[Xie.probeID.isin(intersect_probe)]
Xie_atlas.drop_duplicates(subset = "probeID", keep = "first", inplace = True)
Xie_atlas.sort_values(by = "probeID", ascending = True, inplace = True)
Xie_atlas.reset_index(drop = True, inplace = True)

Atlas_final_df = Atlas_final_df[Atlas_final_df.ProbeID.isin(intersect_probe)]
Atlas_final_df.drop_duplicates(subset = "ProbeID", keep = "first", inplace = True)
Atlas_final_df.sort_values(by = "ProbeID", ascending = True, inplace = True)
Atlas_final_df.reset_index(drop = True, inplace = True)





X = Atlas_final_df.iloc[:, 3:].iloc[:,:].values  # Target matrix
Y = Xie_atlas.iloc[:, 1:101].iloc[:, :].values.T      # Mixed signal matrix (need to be deconvoluted)

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

# Initialize the result list
nnls_results_Xie = []

# Process each sample (add progress bar)
for i in tqdm(range(Y.shape[0]), desc="Processing samples"):
    y = Y[i, :]
    # Use constrained least squares to deconvolute
    coeffs = constrained_ls(X, y)
    nnls_results_Xie.append(coeffs)

# Convert the result to DataFrame
nnls_df_Xie = pd.DataFrame(nnls_results_Xie, columns=["CD4+T", "CD4+TCR", "CD8+T", "CD8+TCR", "B", "Mono", "NK"])
nnls_df_Xie *= 100  # Convert to percentage


<<<<<<< HEAD
=======

### Visualization analysis (boxplot)
>>>>>>> bb3e26c (Initial commit)
# Add stratified labels
labels_Xie = ['Colorectal Cancer'] * 50 + ['Normal'] * 50
nnls_df_Xie['Group'] = labels_Xie

# Convert the data to long format for plotting
nnls_long_Xie = pd.melt(nnls_df_Xie, id_vars='Group', var_name='Cell Type', value_name='Proportion')

# Set the palette for boxplot
boxplot_palette = ["lightcoral", "lightgreen"]
# Set the palette for stripplot
stripplot_palette = ["red", "green"]

# Create a figure
plt.figure(figsize=(12, 6))

# Create a boxplot, set showfliers=False to hide outliers
sns.boxplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Xie, 
            palette=boxplot_palette, width=0.7, showfliers=False)

# Add individual sample points, using dark color
sns.stripplot(x='Cell Type', y='Proportion', hue='Group', data=nnls_long_Xie,
              palette=stripplot_palette, dodge=True, size=4, alpha=0.6)

# Set labels
plt.xlabel('')
plt.ylabel('Proportion (%)', fontweight='bold')
plt.xticks(rotation=45, ha='right', fontweight='bold')
plt.title('Xie et al., 2023', fontweight='bold')

# Set legend
handles, _ = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[:2], labels=['Colorectal Cancer (N = 50)', 'Normal (N = 50)'], 
          prop={'weight': 'bold'}, loc = 'upper left')

# Remove top and right borders
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Remove grid lines
plt.grid(False)

# Adjust layout
plt.tight_layout()

# Show the figure
plt.show()

<<<<<<< HEAD

=======
# Save the result
>>>>>>> bb3e26c (Initial commit)
nnls_df_Xie_TCR = nnls_df_Xie.copy()
nnls_df_Xie_TCR.to_csv('/Users/scui2/DNAmethylation/Xie_et_al/nnls_df_Xie_TCR.csv', index = False)


