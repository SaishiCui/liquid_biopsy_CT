import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import json


"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: Step_B: Calulate_Spearman_correlation_CD8T_testing .py
Purpose: Calculate Spearman correlation between DNA methylation and RNA-seq in the most extreme negative and positive correlated windows for each gene in testing samples.
--------------------------------------------------------------------------------
Data Inputs:

- Matched samples of CD8+T RNA-seq and DNAmethylation data. (Most extreme negative and positive correlated windows)

Data Outputs:
- Spearman correlation between DNA methylation and RNA-seq in the most extreme negative and positive correlated windows for each gene in testing samples.
--------------------------------------------------------------------------------
"""





# Read CD8+T RNA-seq testing data (for most extreme negative and positive correlated windows)  (3 studies combined)
CD8T_RNA_testing_pos_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_testing_pos_window_sorted.csv')
CD8T_RNA_testing_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_testing_neg_window_sorted.csv')
CD8T_Meth_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_testing.csv')
CD8T_RNA_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_testing.csv')



# Negative correlated 
CD8T_Meth_neg_window_testing = pd.DataFrame()
for i in range(CD8T_RNA_testing_neg_window_sorted.shape[0]):

    chr, start, end = CD8T_RNA_testing_neg_window_sorted["Promoter_window"][i].split("_")
    start = int(start)
    end = int(end)
    
    condition = (CD8T_Meth_testing["Chr_hg38"] == chr) & \
               (CD8T_Meth_testing["start_hg38"] >= start) & \
               (CD8T_Meth_testing["end_hg38"] <= end)
    
    filtered_data = CD8T_Meth_testing.loc[condition]
    window_mean_meth = filtered_data.iloc[:,4:].mean(axis=0).to_frame().T
    CD8T_Meth_neg_window_testing = pd.concat([CD8T_Meth_neg_window_testing, window_mean_meth], ignore_index=True)
    print(i)



# Positive correlated 
CD8T_Meth_pos_window_testing = pd.DataFrame()
for i in range(CD8T_RNA_testing_pos_window_sorted.shape[0]):

    chr, start, end = CD8T_RNA_testing_pos_window_sorted["Promoter_window"][i].split("_")
    start = int(start)
    end = int(end)
    
    condition = (CD8T_Meth_testing["Chr_hg38"] == chr) & \
               (CD8T_Meth_testing["start_hg38"] >= start) & \
               (CD8T_Meth_testing["end_hg38"] <= end)
    
    filtered_data = CD8T_Meth_testing.loc[condition]
    window_mean_meth = filtered_data.iloc[:,4:].mean(axis=0).to_frame().T
    CD8T_Meth_pos_window_testing = pd.concat([CD8T_Meth_pos_window_testing, window_mean_meth], ignore_index=True)
    print(i)





# Calculate Spearman correlation coefficient matrix
spearman_corr_matrix_neg = np.zeros((30, 30))
unajusted_pvalue_matrix_neg = np.zeros((30, 30))

for i in range(30):
    for j in range(30):
        spearman_corr_matrix_neg[i, j], unajusted_pvalue_matrix_neg[i, j] = spearmanr(CD8T_Meth_neg_window_testing.iloc[i,:], CD8T_RNA_testing_neg_window_sorted.iloc[j,4:15])
        print(i, j)

# Convert the correlation coefficient matrix to a DataFrame
spearman_corr_df_neg = pd.DataFrame(spearman_corr_matrix_neg, index=CD8T_RNA_testing_neg_window_sorted["Gene_Name"][:30], columns=CD8T_RNA_testing_neg_window_sorted["Gene_Name"][:30])
unajusted_pvalue_df_neg = pd.DataFrame(-np.log10(unajusted_pvalue_matrix_neg), index=CD8T_RNA_testing_neg_window_sorted["Gene_Name"][:30], columns=CD8T_RNA_testing_neg_window_sorted["Gene_Name"][:30])


spearman_corr_matrix_pos = np.zeros((30, 30))
unajusted_pvalue_matrix_pos = np.zeros((30, 30))

for i in range(30):
    for j in range(30):
        spearman_corr_matrix_pos[i, j], unajusted_pvalue_matrix_pos[i, j] = spearmanr(CD8T_Meth_pos_window_testing.iloc[i,:], CD8T_RNA_testing_pos_window_sorted.iloc[j,4:15])
        print(i, j)

# 将相关系数矩阵转换为 DataFrame
spearman_corr_df_pos = pd.DataFrame(spearman_corr_matrix_pos, index=CD8T_RNA_testing_pos_window_sorted["Gene_Name"][:30], columns=CD8T_RNA_testing_pos_window_sorted["Gene_Name"][:30])
unajusted_pvalue_df_pos = pd.DataFrame(-np.log10(unajusted_pvalue_matrix_pos), index=CD8T_RNA_testing_pos_window_sorted["Gene_Name"][:30], columns=CD8T_RNA_testing_pos_window_sorted["Gene_Name"][:30])



# Calculate Spearman Correlation between DNA methylation and RNA-seq for each gene in testing samples


CD8T_spearman_corr_neg_list_testing = []
for i in range(CD8T_RNA_testing_neg_window_sorted.shape[0]):
    spearman_corr, _ = spearmanr(CD8T_RNA_testing_neg_window_sorted.iloc[i,4:15], CD8T_Meth_neg_window_testing.iloc[i,:])
    CD8T_spearman_corr_neg_list_testing.append(spearman_corr)
    print(i)
    

CD8T_spearman_corr_pos_list_testing = []
for i in range(CD8T_RNA_testing_pos_window_sorted.shape[0]):
    spearman_corr, _ = spearmanr(CD8T_RNA_testing_pos_window_sorted.iloc[i,4:15], CD8T_Meth_pos_window_testing.iloc[i,:])
    CD8T_spearman_corr_pos_list_testing.append(spearman_corr)
    print(i)

# Define save path
save_path_pos = '/Users/scui2/DNAmethylation/Corr/CD8T_spearman_corr_pos_list_testing.json'
save_path_neg = '/Users/scui2/DNAmethylation/Corr/CD8T_spearman_corr_neg_list_testing.json'

# Save as JSON file
with open(save_path_pos, 'w') as f:
    json.dump(CD8T_spearman_corr_pos_list_testing, f, indent=4)
with open(save_path_neg, 'w') as f:
    json.dump(CD8T_spearman_corr_neg_list_testing, f, indent=4)




# heatmap
colors = [(0, 'blue'), (0.5, 'white'), (1, 'red')]  # Color from blue to white to red
n_bins = 100  # Number of color steps
cmap_name = 'custom_cmap'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

# Plot heatmap
plt.figure(figsize=(10, 8))  # Increase image size
ax = sns.heatmap(spearman_corr_df_neg.iloc[:30, :30], cmap=custom_cmap, center=0, vmin=-1, vmax=1, cbar_kws={'label': 'Spearman Correlation', "shrink": 0.8, "aspect": 30})
plt.title('CD8+ T cell testing samples (Top 30 negative correlated genes)', fontweight='bold', fontsize= 9)
ax.set_xticklabels(ax.get_xticklabels(), rotation=50, fontsize=6, fontweight='bold')
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6, fontweight='bold')
plt.xlabel('DNA methylation (most extreme negative correlated promoter window)', fontweight='bold', fontsize= 9)
plt.ylabel('RNA-seq (Microarray)', fontweight='bold', fontsize= 9)
plt.show()



