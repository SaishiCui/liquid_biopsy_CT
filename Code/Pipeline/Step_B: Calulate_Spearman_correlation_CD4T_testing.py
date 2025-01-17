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
Program: Step_B: Calulate_Spearman_correlation_CD4T_testing .py
Purpose: Calculate Spearman correlation between DNA methylation and RNA-seq in the most extreme negative and positive correlated windows for each gene in testing samples.
--------------------------------------------------------------------------------
Data Inputs:

- Matched samples of CD4+T RNA-seq and DNAmethylation data. (Most extreme negative and positive correlated windows)

Data Outputs:
- Spearman correlation between DNA methylation and RNA-seq in the most extreme negative and positive correlated windows for each gene in testing samples.
--------------------------------------------------------------------------------
"""



# Read CD4+T RNA-seq testing data (for most extreme negative and positive correlated windows)

CD4T_RNA_testing_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD4T_RNA_testing_neg_window_sorted.csv')
CD4T_RNA_testing_pos_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD4T_RNA_testing_pos_window_sorted.csv')

# Read CD4+T DNA methylation testing data
CD4T_Meth_testing = pd.read_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_Meth_testing.csv', index_col=None)




### Negative correlated 
CD4T_Meth_neg_window_testing = pd.DataFrame()
for i in range(CD4T_RNA_testing_neg_window_sorted.shape[0]):

    chr, start, end = CD4T_RNA_testing_neg_window_sorted["Promoter_window"][i].split("_")
    start = int(start)
    end = int(end)
    
    condition = (CD4T_Meth_testing["Chr_hg38"] == chr) & \
               (CD4T_Meth_testing["start_hg38"] >= start) & \
               (CD4T_Meth_testing["end_hg38"] <= end)
    
    filtered_data = CD4T_Meth_testing.loc[condition]
    window_mean_meth = filtered_data.iloc[:,4:].mean(axis=0).to_frame().T
    CD4T_Meth_neg_window_testing = pd.concat([CD4T_Meth_neg_window_testing, window_mean_meth], ignore_index=True)
    print(i)



### pos correlated 
CD4T_Meth_pos_window_testing = pd.DataFrame()
for i in range(CD4T_RNA_testing_pos_window_sorted.shape[0]):

    chr, start, end = CD4T_RNA_testing_pos_window_sorted["Promoter_window"][i].split("_")
    start = int(start)
    end = int(end)
    
    condition = (CD4T_Meth_testing["Chr_hg38"] == chr) & \
               (CD4T_Meth_testing["start_hg38"] >= start) & \
               (CD4T_Meth_testing["end_hg38"] <= end)
    
    filtered_data = CD4T_Meth_testing.loc[condition]
    window_mean_meth = filtered_data.iloc[:,4:].mean(axis=0).to_frame().T
    CD4T_Meth_pos_window_testing = pd.concat([CD4T_Meth_pos_window_testing, window_mean_meth], ignore_index=True)
    print(i)






# Calculate Spearman correlation coefficient matrix
spearman_corr_matrix_neg = np.zeros((30, 30))
unajusted_pvalue_matrix_neg = np.zeros((30, 30))

for i in range(30):
    for j in range(30):
        spearman_corr_matrix_neg[i, j], unajusted_pvalue_matrix_neg[i, j] = spearmanr(CD4T_Meth_neg_window_testing.iloc[i,:], CD4T_RNA_testing_neg_window_sorted.iloc[j,5:43])
        print(i, j)

# Convert the correlation coefficient matrix to DataFrame
spearman_corr_df_neg = pd.DataFrame(spearman_corr_matrix_neg, index= CD4T_RNA_testing_neg_window_sorted["Gene_Name"][:30], columns=CD4T_RNA_testing_neg_window_sorted["Gene_Name"][:30])
unajusted_pvalue_df_neg = pd.DataFrame(-np.log10(unajusted_pvalue_matrix_neg), index=CD4T_RNA_testing_neg_window_sorted["Gene_Name"][:30], columns=CD4T_RNA_testing_neg_window_sorted["Gene_Name"][:30])


spearman_corr_matrix_pos = np.zeros((30, 30))
unajusted_pvalue_matrix_pos = np.zeros((30, 30))

for i in range(30):
    for j in range(30):
        spearman_corr_matrix_pos[i, j], unajusted_pvalue_matrix_pos[i, j] = spearmanr(CD4T_Meth_pos_window_testing.iloc[i,:], CD4T_RNA_testing_pos_window_sorted.iloc[j,5:43])
        print(i, j)

# Convert the correlation coefficient matrix to DataFrame
spearman_corr_df_pos = pd.DataFrame(spearman_corr_matrix_pos, index=CD4T_RNA_testing_pos_window_sorted["Gene_Name"][:30], columns=CD4T_RNA_testing_pos_window_sorted["Gene_Name"][:30])
unajusted_pvalue_df_pos = pd.DataFrame(-np.log10(unajusted_pvalue_matrix_pos+0.000001), index=CD4T_RNA_testing_pos_window_sorted["Gene_Name"][:30], columns=CD4T_RNA_testing_pos_window_sorted["Gene_Name"][:30])









# Calculate Spearman Correlation for each gene training samples



CD4T_spearman_corr_neg_list_testing = []
for i in range(CD4T_RNA_testing_neg_window_sorted.shape[0]):
    spearman_corr, _ = spearmanr(CD4T_RNA_testing_neg_window_sorted.iloc[i,5:43], CD4T_Meth_neg_window_testing.iloc[i,:])
    CD4T_spearman_corr_neg_list_testing.append(spearman_corr)
    print(i)
    

CD4T_spearman_corr_pos_list_testing = []
for i in range(CD4T_RNA_testing_pos_window_sorted.shape[0]):
    spearman_corr, _ = spearmanr(CD4T_RNA_testing_pos_window_sorted.iloc[i,5:43], CD4T_Meth_pos_window_testing.iloc[i,:])
    CD4T_spearman_corr_pos_list_testing.append(spearman_corr)
    print(i)


# Define the save path for Spearman correlation
save_path_spearman_pos = '/Users/scui2/DNAmethylation/Corr/CD4T_spearman_corr_pos_list_testing.json'
save_path_spearman_neg = '/Users/scui2/DNAmethylation/Corr/CD4T_spearman_corr_neg_list_testing.json'


with open(save_path_spearman_pos, 'w') as f:
    json.dump(CD4T_spearman_corr_pos_list_testing, f, indent=4)
with open(save_path_spearman_neg, 'w') as f:
    json.dump(CD4T_spearman_corr_neg_list_testing, f, indent=4)



# Heatmap for testing samples for visualizing the Spearman correlation

# heatmap spearman correlation
colors = [(0, 'blue'), (0.5, 'white'), (1, 'red')]  # 颜色从蓝色到白色再到红色
n_bins = 100  # 颜色渐变的步数
cmap_name = 'custom_cmap'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

# 绘制热图
plt.figure(figsize=(10, 8))  # 增加图像大小
ax = sns.heatmap(spearman_corr_df_neg.iloc[:30, :30], cmap=custom_cmap, center=0, vmin=-1, vmax=1, cbar_kws={'label': 'Spearman Correlation', "shrink": 0.8, "aspect": 30})
plt.title('CD4+ T cell testing samples (Top 30 negative correlated genes)', fontweight='bold', fontsize= 9)
ax.set_xticklabels(ax.get_xticklabels(), rotation=55, fontsize=8, fontweight='bold')
ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, fontweight='bold')
plt.xlabel('DNA methylation (most extreme negative correlated promoter window)', fontweight='bold', fontsize= 9)
plt.ylabel('RNA-seq (Normalized)', fontweight='bold', fontsize= 9)
plt.show()





# heatmap p-value
colors = [(0, 'blue'), (0.5, 'white'), (1, 'red')]  # 颜色从蓝色到白色再到红色
n_bins = 50  # 颜色渐变的步数
cmap_name = 'custom_cmap'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

# 绘制热图
plt.figure(figsize=(10, 10))  # 增加图像大小
ax = sns.heatmap(unajusted_pvalue_df_neg.iloc[:30, :30], cmap=custom_cmap, center=3, vmin=0, vmax=10, cbar_kws={'label': '-log10(p-value)'})
plt.title('CD4+ T cell testing samples (Top 30 negative correlated genes)', fontweight='bold', fontsize= 10)

ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, fontweight='bold')
ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, fontweight='bold')
plt.xlabel('DNA methylation (promoter window)', fontweight='bold', fontsize= 10)
plt.ylabel('RNA-seq', fontweight='bold', fontsize= 10)
plt.show()