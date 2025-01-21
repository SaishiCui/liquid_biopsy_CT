import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from tqdm import tqdm

"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-17 by Saishi Cui
Program: Step_D: Deconvolution_Xie_et_al_withoutTCR.py
Purpose: Deconvolution of Xie et al., 2023 without using Cancer Reactive T cells.
--------------------------------------------------------------------------------
Data Inputs:
- PBMC atlas without Cancer Reactive T cells obtained from Step_C: Build_PBMC_with_CancerReactive_Atlas.py under the folder "Pipeline".
- Xie et al., 2023 data from Xie_et_al.csv under the data folder "Xie_et_al".

Data Outputs:
- Deconvolution results (without TCR) under the data folder "Xie_et_al".

<<<<<<< HEAD
=======
Visualization analysis:
- Boxplot of deconvolution results (without TCR) under the data folder "Xie_et_al".

>>>>>>> bb3e26c (Initial commit)
Functions:
- constrained_ls(): Deconvolution of Xie et al., 2023 without using Cancer Reactive T cells.
--------------------------------------------------------------------------------
"""


<<<<<<< HEAD
Atlas_final_df = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_final.csv')
=======
Atlas_final_df = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_final_df.csv')
>>>>>>> bb3e26c (Initial commit)
Xie = pd.read_csv('/Users/scui2/DNAmethylation/Xie_et_al/Xie_et_al.csv')  ### 728,310


intersect_probe = list(set(Atlas_final_df.ProbeID).intersection(set(Xie.probeID)))
Xie_atlas = Xie[Xie.probeID.isin(intersect_probe)]
Xie_atlas.drop_duplicates(subset = "probeID", keep = "first", inplace = True)
Xie_atlas.sort_values(by = "probeID", ascending = True, inplace = True)
Xie_atlas.reset_index(drop = True, inplace = True)

Atlas_final_df = Atlas_final_df[Atlas_final_df.ProbeID.isin(intersect_probe)]
Atlas_final_df.drop_duplicates(subset = "ProbeID", keep = "first", inplace = True)
Atlas_final_df.sort_values(by = "ProbeID", ascending = True, inplace = True)
Atlas_final_df.reset_index(drop = True, inplace = True)





X = Atlas_final_df.iloc[:, 3:].iloc[:,[0,2,4,5,6]].values  # Target matrix
Y = Xie_atlas.iloc[:, 1:101].iloc[:, :].values.T      # Mixed signal matrix

def constrained_ls(A, b):
    n = A.shape[1]
    
    # Objective function: Least squares
    def objective(x):
        return np.sum((A @ x - b) ** 2)
    
    # Constraint: Sum to 1
    constraints = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
    ]
    
    # Boundary condition: Non-negative
    bounds = [(0, None) for _ in range(n)]
    
    # Initial guess: Uniform distribution
    x0 = np.ones(n) / n
    
    # 求解
    result = minimize(objective, x0, 
                     method='SLSQP',
                     bounds=bounds,
                     constraints=constraints)
    
    return result.x

# Initialize result list
nnls_results_Xie = []

# Deconvolution for each sample (with progress bar)
for i in tqdm(range(Y.shape[0]), desc="Processing samples"):
    y = Y[i, :]
    # Deconvolution using constrained least squares
    coeffs = constrained_ls(X, y)
    nnls_results_Xie.append(coeffs)

# Convert results to DataFrame
nnls_df_Xie = pd.DataFrame(nnls_results_Xie, columns=["CD4+T", "CD8+T", "B", "Mono", "NK"])
nnls_df_Xie *= 100  # Convert to percentage



labels_Xie = ['Colorectal Cancer'] * 50 + ['Normal'] * 50
nnls_df_Xie['Group'] = labels_Xie

<<<<<<< HEAD
=======

### Visualization analysis
>>>>>>> bb3e26c (Initial commit)
# Convert data to long format for plotting
nnls_long_Xie = pd.melt(nnls_df_Xie, id_vars='Group', var_name='Cell Type', value_name='Proportion')

# Set light color palette for boxplot
boxplot_palette = ["lightcoral", "lightgreen"]
# Set dark color palette for stripplot
stripplot_palette = ["red", "green"]

# Create figure
plt.figure(figsize=(12, 6))

# Create boxplot, set showfliers=False to hide outliers
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

# Display figure
plt.show()

<<<<<<< HEAD

=======
# Save the result
>>>>>>> bb3e26c (Initial commit)
nnls_df_Xie_withoutTCR = nnls_df_Xie.copy()
nnls_df_Xie_withoutTCR.to_csv('/Users/scui2/DNAmethylation/Xie_et_al/nnls_df_Xie_withoutTCR.csv', index = False)


