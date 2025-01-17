import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from statsmodels.stats.multitest import multipletests
from collections import Counter
import json
from matplotlib.font_manager import FontProperties


"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: Step_B: Select_Gene_list_B_(CD8T).py
Purpose: Select gene list B (genes whose gene expression is negatively correlated with DNA methylation) after sliding window. (for CD8+T cells)
--------------------------------------------------------------------------------
Data Inputs:

- Spearman correlation between DNA methylation and RNA-seq in the most extreme negative correlated windows for each gene in both training and testing CD4+T cells samples.

Data Outputs:
- Gene list B (genes whose gene expression is negatively correlated with DNA methylation) after sliding window.
--------------------------------------------------------------------------------
"""

# Write paths to Spearman correlation lists

CD8T_path_pos_training = '/Users/scui2/DNAmethylation/Corr/CD8T_spearman_corr_pos_list_training.json'
CD8T_path_neg_training = '/Users/scui2/DNAmethylation/Corr/CD8T_spearman_corr_neg_list_training.json'
CD8T_path_pos_testing = '/Users/scui2/DNAmethylation/Corr/CD8T_spearman_corr_pos_list_testing.json'
CD8T_path_neg_testing = '/Users/scui2/DNAmethylation/Corr/CD8T_spearman_corr_neg_list_testing.json'


# Read Spearman correlation lists
CD8T_RNA_training_pos_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_pos_window.csv')
CD8T_RNA_training_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_neg_window.csv')


with open(CD8T_path_pos_training, 'r') as f:
    CD8T_spearman_corr_pos_list_training = json.load(f)

with open(CD8T_path_neg_training, 'r') as f:
    CD8T_spearman_corr_neg_list_training = json.load(f)

with open(CD8T_path_pos_testing, 'r') as f:
    CD8T_spearman_corr_pos_list_testing = json.load(f)

with open(CD8T_path_neg_testing, 'r') as f:
    CD8T_spearman_corr_neg_list_testing = json.load(f)





CD8T_corr_df = pd.DataFrame({"Pos_training_gene": CD8T_RNA_training_pos_window_sorted["Ensembl_Gene_ID"],
                            "Pos_training_gene_name": CD8T_RNA_training_pos_window_sorted["Gene_Name"],
                             "Pos_training":CD8T_spearman_corr_pos_list_training,
                             "Pos_testing":CD8T_spearman_corr_pos_list_testing,
                             "Neg_training_gene": CD8T_RNA_training_neg_window_sorted["Ensembl_Gene_ID"],
                             "Neg_training_gene_name": CD8T_RNA_training_neg_window_sorted["Gene_Name"],
                             "Neg_training":CD8T_spearman_corr_neg_list_training,
                             "Neg_testing":CD8T_spearman_corr_neg_list_testing})


CD8T_corr_df.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_correlated_genes.csv', index=False)


# Visualize the Spearman correlation between training and testing
x = CD8T_corr_df['Pos_training']
y = CD8T_corr_df['Pos_testing']


# Create scatter plot
plt.figure(figsize=(10, 8))
font_properties = FontProperties(weight='bold')
# Plot all points, default blue
plt.scatter(x, y, color='blue', label='Transcriptional Promoter Windows', alpha=0.8, s=7)


red_points = (x > 0.2) & (y > 0.2)
plt.scatter(x[red_points], y[red_points], color='red', label='Positive correlated ($\\rho > 0.2$) both training and testing', alpha=0.8, s=7)

# Plot top 30 points, green
plt.scatter(x[:30], y[:30], color='green', label='Top 30 positive correlated in training', alpha=0.8, s=7)


# Plot horizontal and vertical lines
plt.axhline(y= 0.2, color='black', linestyle='--')
plt.axvline(x= 0.2, color='black', linestyle='--')

# Add labels and title
plt.xlabel('Training Spearman Correlation ($\\rho$)', fontsize=8, fontweight='bold')
plt.ylabel('Testing Spearman Correlation ($\\rho$)', fontsize=8, fontweight='bold')
plt.xticks(fontsize=8, fontweight='bold')
plt.yticks(fontsize=8, fontweight='bold')
plt.title('Scatter Plot of Training vs Testing Spearman Correlation, CD8+T cells', fontsize=8, fontweight='bold')
plt.legend(fontsize=8, prop=font_properties, loc='upper left')
plt.show()


# Filter genes

CD8T_filtered_pos = CD8T_corr_df[(CD8T_corr_df['Pos_training'] > 0.2) & (CD8T_corr_df['Pos_testing'] > 0.2)]
CD8T_filtered_neg = CD8T_corr_df[(CD8T_corr_df['Neg_training'] < -0.2) & (CD8T_corr_df['Neg_testing'] < -0.2)]

CD8T_filtered_pos.reset_index(drop=True, inplace=True)
CD8T_filtered_neg.reset_index(drop=True, inplace=True)  



CD8T_filtered_pos.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_filtered_pos.csv', index=False)
CD8T_filtered_neg.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_filtered_neg.csv', index=False)


# Each gene's correlation for samples (most extreme window)

CD8T_Meth_window_training = pd.read_csv('/Users/scui2/DNAmethylation/CD4T_CD8T_processed/CD8T_Meth_window_training.csv')
CD8T_Meth_window_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD4T_CD8T_processed/CD8T_Meth_window_testing.csv')

CD8T_RNA_training = pd.read_csv('/Users/scui2/DNAmethylation/CD4T_CD8T_processed/CD8T_RNA_training.csv')
CD8T_RNA_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD4T_CD8T_processed/CD8T_RNA_testing.csv')
CD8T_RNA_training_pos_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_pos_window.csv')
CD8T_RNA_training_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_neg_window.csv')


top1_pos_window_CD8T = CD8T_corr_df["Pos_training_gene"][0]
condition_top1_pos_window_CD8T = (CD8T_RNA_training["gene_id"] == top1_pos_window_CD8T)
top1_pos_window_CD8T_location = CD8T_RNA_training_pos_window_sorted["Promoter_window"][0].split("_")

condition_top1_pos_window_CD8T_meth = (CD8T_Meth_window_training["chr"] == top1_pos_window_CD8T_location[0]) & (CD8T_Meth_window_training["start"] >= int(top1_pos_window_CD8T_location[1])) & (CD8T_Meth_window_training["start"] <= int(top1_pos_window_CD8T_location[2]))
top1_pos_window_Meth_CD8T = CD8T_Meth_window_training[condition_top1_pos_window_CD8T_meth].iloc[:,4:].mean(axis=0).tolist()+CD8T_Meth_window_testing[condition_top1_pos_window_CD8T_meth].iloc[:,4:].mean(axis=0).tolist()
top1_pos_window_RNA_CD8T = CD8T_RNA_training[condition_top1_pos_window_CD8T].iloc[:,2:].values[0].tolist()+CD8T_RNA_testing[condition_top1_pos_window_CD8T].iloc[:,2:].values[0].tolist()


# 创建散点图
plt.figure(figsize=(10, 8))
corr, pval = spearmanr(top1_pos_window_Meth_CD8T, top1_pos_window_RNA_CD8T)
plt.scatter(top1_pos_window_Meth_CD8T, top1_pos_window_RNA_CD8T, color='blue', label='CD8+T cells sample', alpha=1)
coeffs_cubic = np.polyfit(top1_pos_window_Meth_CD8T, top1_pos_window_RNA_CD8T, 3)
# 生成拟合曲线
x_fit = np.linspace(min(top1_pos_window_Meth_CD8T), max(top1_pos_window_Meth_CD8T), 100)
y_fit_cubic = np.polyval(coeffs_cubic, x_fit)

plt.plot(x_fit, y_fit_cubic, color='red', linestyle='-', label='Cubic Fit')
plt.text(0.03, 0.95, f'Spearman $\\rho$ = {corr:.3f}, p-value < 0.001', transform=plt.gca().transAxes, ha='left', va='top', fontsize=12, color='red', fontweight='bold')
# 添加标签和标题
plt.xlabel('Methylation')
plt.ylabel('Gene Expression')
title = f"{top1_pos_window_CD8T} (Most positive correlated window) ({top1_pos_window_CD8T_location[0]}: {top1_pos_window_CD8T_location[1]}-{top1_pos_window_CD8T_location[2]})"
plt.title(title)
plt.legend(loc = "lower right")
plt.grid(True)
plt.show()





####### Plot top 20

fig, axes = plt.subplots(4, 5, figsize=(20, 16))
axes = axes.flatten()  


for i in range(20):
    neg_window_gene = CD8T_corr_df["Neg_training_gene"][i]
    condition_neg_window_gene = (CD8T_RNA_training["gene_id"] == neg_window_gene)
    neg_window_location = CD8T_RNA_training_neg_window_sorted["Promoter_window"][i].split("_")
    condition_neg_window_meth = (
        (CD8T_Meth_window_training["chr"] == neg_window_location[0]) &
        (CD8T_Meth_window_training["start"] >= int(neg_window_location[1])) &
        (CD8T_Meth_window_training["start"] <= int(neg_window_location[2]))
    )
    neg_window_Meth_CD8T = (
        CD8T_Meth_window_training[condition_neg_window_meth].iloc[:, 4:].mean(axis=0).tolist() +
        CD8T_Meth_window_testing[condition_neg_window_meth].iloc[:, 4:].mean(axis=0).tolist()
    )
    neg_window_RNA_CD8T = (
        CD8T_RNA_training[condition_neg_window_gene].iloc[:, 2:].values[0].tolist() +
        CD8T_RNA_testing[condition_neg_window_gene].iloc[:, 2:].values[0].tolist()
    )


    ax = axes[i]

    corr, pval = spearmanr(neg_window_Meth_CD8T, neg_window_RNA_CD8T)
    ax.scatter(neg_window_Meth_CD8T, neg_window_RNA_CD8T, color='green', alpha=1, s=10)
    coeffs_cubic = np.polyfit(neg_window_Meth_CD8T, neg_window_RNA_CD8T, 3)

    # 生成拟合曲线
    x_fit = np.linspace(min(neg_window_Meth_CD8T), max(neg_window_Meth_CD8T), 100)
    y_fit_cubic = np.polyval(coeffs_cubic, x_fit)

    ax.plot(x_fit, y_fit_cubic, color='red', linestyle='-')

    ax.text(0.03, 0.95, f'Spearman $\\rho$ = {corr:.3f}', transform=ax.transAxes, ha='left', va='top', fontsize=5, color='red', fontweight='bold')

    title = f"{neg_window_gene} (Promoter, {neg_window_location[0]}: {neg_window_location[1]}-{neg_window_location[2]})"
    ax.set_title(title, fontsize=6, fontweight='bold', pad = 1)
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=4)
    ax.set_xticklabels([f'{tick:.2f}' for tick in ax.get_xticks()], fontweight='bold', rotation=45)
    ax.set_yticklabels(ax.get_yticks(), fontweight='bold')

plt.tight_layout()
plt.show()




fig, axes = plt.subplots(4, 5, figsize=(20, 16))
axes = axes.flatten()  

for i in range(20):
    pos_window_gene = CD8T_corr_df["Pos_training_gene"][i]
    condition_pos_window_gene = (CD8T_RNA_training["gene_id"] == pos_window_gene)
    pos_window_location = CD8T_RNA_training_pos_window_sorted["Promoter_window"][i].split("_")
    condition_pos_window_meth = (
        (CD8T_Meth_window_training["chr"] == pos_window_location[0]) &
        (CD8T_Meth_window_training["start"] >= int(pos_window_location[1])) &
        (CD8T_Meth_window_training["start"] <= int(pos_window_location[2]))
    )
    pos_window_Meth_CD8T = (
        CD8T_Meth_window_training[condition_pos_window_meth].iloc[:, 4:].mean(axis=0).tolist() +
        CD8T_Meth_window_testing[condition_pos_window_meth].iloc[:, 4:].mean(axis=0).tolist()
    )
    pos_window_RNA_CD8T = (
        CD8T_RNA_training[condition_pos_window_gene].iloc[:, 2:].values[0].tolist() +
        CD8T_RNA_testing[condition_pos_window_gene].iloc[:, 2:].values[0].tolist()
    )

   
    ax = axes[i]

    corr, pval = spearmanr(pos_window_Meth_CD8T, pos_window_RNA_CD8T)
    ax.scatter(pos_window_Meth_CD8T, pos_window_RNA_CD8T, color='blue', alpha=1, s=10)
    coeffs_cubic = np.polyfit(pos_window_Meth_CD8T, pos_window_RNA_CD8T, 3)

    x_fit = np.linspace(min(pos_window_Meth_CD8T), max(pos_window_Meth_CD8T), 100)
    y_fit_cubic = np.polyval(coeffs_cubic, x_fit)

    ax.plot(x_fit, y_fit_cubic, color='red', linestyle='-')

    ax.text(0.03, 0.95, f'Spearman $\\rho$ = {corr:.3f}', transform=ax.transAxes, ha='left', va='top', fontsize=5, color='red', fontweight='bold')

    title = f"{pos_window_gene} (Promoter, {pos_window_location[0]}: {pos_window_location[1]}-{pos_window_location[2]})"
    ax.set_title(title, fontsize=6, fontweight='bold', pad = 1)
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=4)
    ax.set_xticklabels([f'{tick:.2f}' for tick in ax.get_xticks()], fontweight='bold', rotation=45)
    ax.set_yticklabels(ax.get_yticks(), fontweight='bold')


plt.tight_layout()
plt.show()


