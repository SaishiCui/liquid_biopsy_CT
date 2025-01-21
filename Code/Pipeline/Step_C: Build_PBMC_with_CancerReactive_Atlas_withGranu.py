import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy.polynomial.polynomial as poly




"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-17 by Saishi Cui
Program: Step_C: Build_PBMC_with_CancerReactive_Atlas_withGranu.py
Purpose: Build PBMC atlas with Cancer Reactive T cells.
--------------------------------------------------------------------------------
Data Inputs:
- WGBS PBMC DNAmethylation data obtained from WGBS_process.py under the folder "Data_prepare".
- Differential Methylated Probes (DMPs) results of 850K EPIC data obtained from Step_A: Find_DMP.py under the folder "Pipeline".
- CD4+T/CD8+T Cancer Reactive T cells signatures obtained from Step_B: LM_filtering_and_Extrapolation_CD8TCR.py and Step_B: LM_filtering_and_Extrapolation_CD4TCR.py under the folder "Pipeline".

Data Outputs:
- PBMC atlas with Cancer Reactive T cells.

Functions:
- get_850K_ProbeID_mean_beta(): Calculate the mean beta value for each 850K DMP in the Atlas.
- match_850K_ProbeID_TCR(): Match the 850K DMPs with the TCR signatures.
- get_TCRWindow_mean_beta():  Calculate the mean beta value for each matched probe in the Cancer Reactive TCR signatures based on extrapolation.

Analysis:
- Visualize some DMPs in the Atlas (Barplot).
- Check the correlation between the built Atlas.
--------------------------------------------------------------------------------
"""



# Read DMP results
DMP_850K_final_withGranu = pd.read_csv('/Users/scui2/DNAmethylation/DMP/DMP_850K_results_final_withGranu.csv', index_col=None)

# DMP_850K_results_final_withGranu = pd.read_csv('/Users/scui2/DNAmethylation/DMP/DMP_850K_results_final_withGranu.csv', index_col=None) # for DMPs with Granu

DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] = np.abs(DMP_850K_final_withGranu["CD4T_Others_mean_diff"])
DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] = np.abs(DMP_850K_final_withGranu["CD8T_Others_mean_diff"])
DMP_850K_final_withGranu["B_Others_mean_diff_abs"] = np.abs(DMP_850K_final_withGranu["B_Others_mean_diff"])
DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] = np.abs(DMP_850K_final_withGranu["NK_Others_mean_diff"])
DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] = np.abs(DMP_850K_final_withGranu["Mono_Others_mean_diff"])
DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] = np.abs(DMP_850K_final_withGranu["Granu_Others_mean_diff"])
DMP_850K_final_withGranu["CD4T_CD8T_mean_diff_abs"] = np.abs(DMP_850K_final_withGranu["CD4T_CD8T_mean_diff"])


# Filter DMPs based on the criteria (See step A for details)
threshold1= 0.4
threshold2= 0.3

DMP_filter_CD4T_850K_withGranu = DMP_850K_final_withGranu[
    (DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] >= threshold1) &
    (DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["B_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] <= threshold2)].sort_values(
    by=['CD4T_Others_p_value', 'CD4T_Others_mean_diff'],
    ascending=[True, False],
    inplace=False
)[:]


DMP_filter_CD8T_850K_withGranu = DMP_850K_final_withGranu[
    (DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] >= threshold1) &
    (DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["B_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] <= threshold2)
].sort_values(
    by=['CD8T_Others_p_value', 'CD8T_Others_mean_diff'],
    ascending=[True, False],
    inplace=False
)[:]


DMP_filter_B_850K_withGranu = DMP_850K_final_withGranu[
    (DMP_850K_final_withGranu["B_Others_mean_diff_abs"] >= threshold1) &
    (DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] <= threshold2)
].sort_values(
    by=['B_Others_p_value', 'B_Others_mean_diff'],
    ascending=[True, False],
    inplace=False
)[:]


DMP_filter_NK_850K_withGranu = DMP_850K_final_withGranu[
    (DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] >= threshold1) &
    (DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["B_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] <= threshold2)
].sort_values(
    by=['NK_Others_p_value', 'NK_Others_mean_diff'],
    ascending=[True, False],
    inplace=False
)[:]


DMP_filter_Mono_850K_withGranu = DMP_850K_final_withGranu[
    (DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] >= threshold1) &
    (DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["B_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] <= threshold2)
].sort_values(
    by=['Mono_Others_p_value', 'Mono_Others_mean_diff'],
    ascending=[True, False],
    inplace=False
)[:]


DMP_filter_Granu_850K_withGranu = DMP_850K_final_withGranu[
    (DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] >= threshold1) &
    (DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["B_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] <= threshold2)
].sort_values(
    by=['Granu_Others_p_value', 'Granu_Others_mean_diff'],
    ascending=[True, False],
    inplace=False
)[:]




DMP_filter_CD4T_CD8T_850K_withGranu = DMP_850K_final_withGranu[
    (DMP_850K_final_withGranu["CD4T_CD8T_mean_diff_abs"] >= threshold1) &
    (DMP_850K_final_withGranu["CD4T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["CD8T_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["B_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["NK_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Mono_Others_mean_diff_abs"] <= threshold2) &
    (DMP_850K_final_withGranu["Granu_Others_mean_diff_abs"] <= threshold2)
].sort_values(
    by=['CD4T_CD8T_p_value', 'CD4T_CD8T_mean_diff'],
    ascending=[True, False],
    inplace=False
)[:]





DMP_filter_CD4T_850K_withGranu.reset_index(drop=True, inplace=True)
DMP_filter_CD8T_850K_withGranu.reset_index(drop=True, inplace=True)
DMP_filter_B_850K_withGranu.reset_index(drop=True, inplace=True)
DMP_filter_NK_850K_withGranu.reset_index(drop=True, inplace=True)
DMP_filter_Mono_850K_withGranu.reset_index(drop=True, inplace=True)
DMP_filter_Granu_850K_withGranu.reset_index(drop=True, inplace=True)
DMP_filter_CD4T_CD8T_850K_withGranu.reset_index(drop=True, inplace=True)


# Take top 1000 DMPs for each cell type
DMP_filter_CD4T_850K_withGranu   = DMP_filter_CD4T_850K_withGranu.iloc[:1000,:].copy()
DMP_filter_CD8T_850K_withGranu   = DMP_filter_CD8T_850K_withGranu.iloc[:1000,:].copy()
DMP_filter_B_850K_withGranu      = DMP_filter_B_850K_withGranu.iloc[:1000,:].copy()
DMP_filter_NK_850K_withGranu     = DMP_filter_NK_850K_withGranu.iloc[:1000,:].copy()
DMP_filter_Mono_850K_withGranu   = DMP_filter_Mono_850K_withGranu.iloc[:1000,:].copy()
DMP_filter_Granu_850K_withGranu  = DMP_filter_Granu_850K_withGranu.iloc[:1000,:].copy()
DMP_filter_CD4T_CD8T_850K_withGranu  = DMP_filter_CD4T_CD8T_850K_withGranu.iloc[:1000,:].copy()



A = set(DMP_filter_CD4T_850K_withGranu["ProbeID"])
B = set(DMP_filter_CD8T_850K_withGranu["ProbeID"])
C = set(DMP_filter_B_850K_withGranu["ProbeID"])
D = set(DMP_filter_NK_850K_withGranu["ProbeID"])
E = set(DMP_filter_Mono_850K_withGranu["ProbeID"])
F = set(DMP_filter_Granu_850K_withGranu["ProbeID"])
G = set(DMP_filter_CD4T_CD8T_850K_withGranu["ProbeID"])


A_exclusive_list = list(A - (B | C | D | E | F | G))  # exclude the probes that are in other cell types but in CD4+T
B_exclusive_list = list(B - (A | C | D | E | F | G))  # exclude the probes that are in other cell types but in CD8+T
C_exclusive_list = list(C - (A | B | D | E | F | G))  # exclude the probes that are in other cell types but in B
D_exclusive_list = list(D - (A | B | C | E | F | G))  # exclude the probes that are in other cell types but in NK
E_exclusive_list = list(E - (A | B | C | D | F | G))  # exclude the probes that are in other cell types but in Mono
F_exclusive_list = list(F - (A | B | C | D | E | G))  # exclude the probes that are in other cell types but in CD4+T/CD8+T
G_exclusive_list = list(G - (A | B | C | D | E | F))  # exclude the probes that are in other cell types but in Granu


# Visualize some DMPs in the Atlas (Barplot) #

WGBS_allsamples_combined_filtered_withGranu = pd.read_csv('/Users/scui2/DNAmethylation/Loyfer_data/WGBS_allsamples_combined_filtered_withGranu.csv')

fig, axes = plt.subplots(3, 5, figsize=(20, 12))

CpG_ID_CD4T = A_exclusive_list[0:3]
for idx in range(len(CpG_ID_CD4T)):
    chr, start, end = DMP_filter_CD4T_850K_withGranu[DMP_filter_CD4T_850K_withGranu["ProbeID"] == CpG_ID_CD4T[idx]]["Window_Location"].tolist()[0].split("_")
    start = int(start)
    end = int(end)
    condition = (WGBS_allsamples_combined_filtered_withGranu['Chromosome'] == chr) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] >= start) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] <= end)
    CpG_df = WGBS_allsamples_combined_filtered_withGranu[condition]

    
    cell_types = ['CD4+T', 'CD8+T', 'B', 'NK', 'Mono']
    means = []
    sems = []  
    individual_values = []


    colors = ['lightcoral'] + ['lightsteelblue'] * 4

    x = np.arange(len(cell_types)) * 0.7


    for i, cell_type in enumerate(cell_types):
        old_type = cell_type.replace('+', '')  
        values = [CpG_df[f'{old_type}_{j}_Beta'] for j in range(1, 4)]
        means.append(np.mean(values))
        sems.append(np.std(values) / np.sqrt(len(values)))
        individual_values.append(values)


    bars = axes[idx, 0].bar(x, means, yerr=sems, capsize=5, alpha=0.7, color=colors, width=0.4)

 
    baseline = np.mean(means[1:])  
    
  
    axes[idx, 0].axhline(y=baseline, color='gray', linestyle='--', alpha=0.7)

    for i, values in enumerate(individual_values):
        axes[idx, 0].scatter(np.repeat(x[i], len(values)) + np.random.normal(0, 0.05, len(values)), 
                         values, color='black', alpha=0.6)

   
    axes[idx, 0].set_xticks(x, cell_types, rotation=45, fontweight='bold', fontsize=8)
    axes[idx, 0].set_yticks(np.arange(0, 1.1, 0.1))
    axes[idx, 0].set_yticklabels([f'{x:.1f}' for x in np.arange(0, 1.1, 0.1)], fontweight='bold', fontsize=8)
    axes[idx, 0].tick_params(axis='both', which='major', labelsize=8, width=2)
    axes[idx, 0].set_ylabel('Beta Value', fontweight='bold', fontsize=8)
    

    title_bbox = dict(
        facecolor='#FFF7D4',
        edgecolor='black',
        pad=5,
        boxstyle='square,pad=0.5'
    )
    axes[idx, 0].set_title(f'{CpG_ID_CD4T[idx]}', fontweight='bold', fontsize=8,
                          bbox=title_bbox,
                          x=0.5,
                          y=1.02)


    axes[idx, 0].spines['right'].set_visible(False)
    axes[idx, 0].spines['top'].set_visible(False)
    axes[idx, 0].spines['left'].set_linewidth(2)
    axes[idx, 0].spines['bottom'].set_linewidth(2)


    if idx == 2:
        axes[idx, 0].set_xlabel('CD4+T', fontweight='bold')


CpG_ID_CD8T = B_exclusive_list[0:3]
for idx in range(len(CpG_ID_CD8T)):
    chr, start, end = DMP_filter_CD8T_850K_withGranu[DMP_filter_CD8T_850K_withGranu["ProbeID"] == CpG_ID_CD8T[idx]]["Window_Location"].tolist()[0].split("_")
    start = int(start)
    end = int(end)
    condition = (WGBS_allsamples_combined_filtered_withGranu['Chromosome'] == chr) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] >= start) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] <= end)
    CpG_df = WGBS_allsamples_combined_filtered_withGranu[condition]


    cell_types = ['CD4+T', 'CD8+T', 'B', 'NK', 'Mono']
    means = []
    sems = []  
    individual_values = []

    colors = ['lightsteelblue'] + ['lightcoral'] + ['lightsteelblue'] * 3

    x = np.arange(len(cell_types)) * 0.7

    for i, cell_type in enumerate(cell_types):
        old_type = cell_type.replace('+', '')  
        values = [CpG_df[f'{old_type}_{j}_Beta'] for j in range(1, 4)]
        means.append(np.mean(values))
        sems.append(np.std(values) / np.sqrt(len(values)))
        individual_values.append(values)

    bars = axes[idx, 1].bar(x, means, yerr=sems, capsize=5, alpha=0.7, color=colors, width=0.4)

    baseline = np.mean(means[:1]+means[2:])  
    
    axes[idx, 1].axhline(y=baseline, color='gray', linestyle='--', alpha=0.7)

    for i, values in enumerate(individual_values):
        axes[idx, 1].scatter(np.repeat(x[i], len(values)) + np.random.normal(0, 0.05, len(values)), 
                         values, color='black', alpha=0.6)

    axes[idx, 1].set_xticks(x, cell_types, rotation=45, fontweight='bold', fontsize=8)
    axes[idx, 1].set_yticks(np.arange(0, 1.1, 0.1))
    axes[idx, 1].set_yticklabels([f'{x:.1f}' for x in np.arange(0, 1.1, 0.1)], fontweight='bold', fontsize=8)
    axes[idx, 1].tick_params(axis='both', which='major', labelsize=8, width=2)
    axes[idx, 1].set_ylabel('Beta Value', fontweight='bold', fontsize=8)
    
    title_bbox = dict(
        facecolor='#DFF4E0',
        edgecolor='black',
        pad=5,
        boxstyle='square,pad=0.5'
    )
    axes[idx, 1].set_title(f'{CpG_ID_CD8T[idx]}', fontweight='bold', fontsize=8,
                          bbox=title_bbox,
                          x=0.5,
                          y=1.02)

    axes[idx, 1].spines['right'].set_visible(False)
    axes[idx, 1].spines['top'].set_visible(False)
    axes[idx, 1].spines['left'].set_linewidth(2)
    axes[idx, 1].spines['bottom'].set_linewidth(2)

    if idx == 2:
        axes[idx, 1].set_xlabel('CD8+T', fontweight='bold')


CpG_ID_B = C_exclusive_list[0:3]
for idx in range(len(CpG_ID_B)):
    chr, start, end = DMP_filter_B_850K_withGranu[DMP_filter_B_850K_withGranu["ProbeID"] == CpG_ID_B[idx]]["Window_Location"].tolist()[0].split("_")
    start = int(start)
    end = int(end)
    condition = (WGBS_allsamples_combined_filtered_withGranu['Chromosome'] == chr) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] >= start) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] <= end)
    CpG_df = WGBS_allsamples_combined_filtered_withGranu[condition]

    cell_types = ['CD4+T', 'CD8+T', 'B', 'NK', 'Mono']
    means = []
    sems = []  
    individual_values = []

    colors = ['lightsteelblue'] * 2 + ['lightcoral']  + ['lightsteelblue'] * 2

    x = np.arange(len(cell_types)) * 0.7

    for i, cell_type in enumerate(cell_types):
        old_type = cell_type.replace('+', '')  
        values = [CpG_df[f'{old_type}_{j}_Beta'] for j in range(1, 4)]
        means.append(np.mean(values))
        sems.append(np.std(values) / np.sqrt(len(values)))
        individual_values.append(values)

    bars = axes[idx, 2].bar(x, means, yerr=sems, capsize=5, alpha=0.7, color=colors, width=0.4)

    baseline = np.mean(means[:2]+means[3:])  
    
    axes[idx, 2].axhline(y=baseline, color='gray', linestyle='--', alpha=0.7)

    for i, values in enumerate(individual_values):
        axes[idx, 2].scatter(np.repeat(x[i], len(values)) + np.random.normal(0, 0.05, len(values)), 
                         values, color='black', alpha=0.6)

    axes[idx, 2].set_xticks(x, cell_types, rotation=45, fontweight='bold', fontsize=8)
    axes[idx, 2].set_yticks(np.arange(0, 1.1, 0.1))
    axes[idx, 2].set_yticklabels([f'{x:.1f}' for x in np.arange(0, 1.1, 0.1)], fontweight='bold', fontsize=8)
    axes[idx, 2].tick_params(axis='both', which='major', labelsize=8, width=2)
    axes[idx, 2].set_ylabel('Beta Value', fontweight='bold', fontsize=8)
    
    title_bbox = dict(
        facecolor='#E0F4FF',
        edgecolor='black',
        pad=5,
        boxstyle='square,pad=0.5'
    )
    axes[idx, 2].set_title(f'{CpG_ID_B[idx]}', fontweight='bold', fontsize=8,
                          bbox=title_bbox,
                          x=0.5,
                          y=1.02)

    axes[idx, 2].spines['right'].set_visible(False)
    axes[idx, 2].spines['top'].set_visible(False)
    axes[idx, 2].spines['left'].set_linewidth(2)
    axes[idx, 2].spines['bottom'].set_linewidth(2)

    if idx == 2:
        axes[idx, 2].set_xlabel('B', fontweight='bold')


CpG_ID_NK = D_exclusive_list[0:3]
for idx in range(len(CpG_ID_NK)):
    chr, start, end = DMP_filter_NK_850K_withGranu[DMP_filter_NK_850K_withGranu["ProbeID"] == CpG_ID_NK[idx]]["Window_Location"].tolist()[0].split("_")
    start = int(start)
    end = int(end)
    condition = (WGBS_allsamples_combined_filtered_withGranu['Chromosome'] == chr) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] >= start) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] <= end)
    CpG_df = WGBS_allsamples_combined_filtered_withGranu[condition]

    cell_types = ['CD4+T', 'CD8+T', 'B', 'NK', 'Mono']
    means = []
    sems = []  
    individual_values = []

    colors = ['lightsteelblue'] * 3 + ['lightcoral']  + ['lightsteelblue'] * 1

    x = np.arange(len(cell_types)) * 0.7

    for i, cell_type in enumerate(cell_types):
        old_type = cell_type.replace('+', '')  
        values = [CpG_df[f'{old_type}_{j}_Beta'] for j in range(1, 4)]
        means.append(np.mean(values))
        sems.append(np.std(values) / np.sqrt(len(values)))
        individual_values.append(values)

    bars = axes[idx, 3].bar(x, means, yerr=sems, capsize=5, alpha=0.7, color=colors, width=0.4)

    baseline = np.mean(means[:3]+means[4:])  
    
    axes[idx, 3].axhline(y=baseline, color='gray', linestyle='--', alpha=0.7)

    for i, values in enumerate(individual_values):
        axes[idx, 3].scatter(np.repeat(x[i], len(values)) + np.random.normal(0, 0.05, len(values)), 
                         values, color='black', alpha=0.6)

    axes[idx, 3].set_xticks(x, cell_types, rotation=45, fontweight='bold', fontsize=8)
    axes[idx, 3].set_yticks(np.arange(0, 1.1, 0.1))
    axes[idx, 3].set_yticklabels([f'{x:.1f}' for x in np.arange(0, 1.1, 0.1)], fontweight='bold', fontsize=8)
    axes[idx, 3].tick_params(axis='both', which='major', labelsize=8, width=2)
    axes[idx, 3].set_ylabel('Beta Value', fontweight='bold', fontsize=8)
    
    title_bbox = dict(
        facecolor='#F0E6FF',
        edgecolor='black',
        pad=5,
        boxstyle='square,pad=0.5'
    )
    axes[idx, 3].set_title(f'{CpG_ID_NK[idx]}', fontweight='bold', fontsize=8,
                          bbox=title_bbox,
                          x=0.5,
                          y=1.02)

    axes[idx, 3].spines['right'].set_visible(False)
    axes[idx, 3].spines['top'].set_visible(False)
    axes[idx, 3].spines['left'].set_linewidth(2)
    axes[idx, 3].spines['bottom'].set_linewidth(2)

    if idx == 2:
        axes[idx, 3].set_xlabel('NK', fontweight='bold')

CpG_ID_Mono = E_exclusive_list[0:3]
for idx in range(len(CpG_ID_Mono)):
    chr, start, end = DMP_filter_Mono_850K_withGranu[DMP_filter_Mono_850K_withGranu["ProbeID"] == CpG_ID_Mono[idx]]["Window_Location"].tolist()[0].split("_")
    start = int(start)
    end = int(end)
    condition = (WGBS_allsamples_combined_filtered_withGranu['Chromosome'] == chr) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] >= start) & \
                (WGBS_allsamples_combined_filtered_withGranu['Position'] <= end)
    CpG_df = WGBS_allsamples_combined_filtered_withGranu[condition]

    cell_types = ['CD4+T', 'CD8+T', 'B', 'NK', 'Mono']
    means = []
    sems = []  
    individual_values = []

    colors = ['lightsteelblue'] * 4 + ['lightcoral']  

    x = np.arange(len(cell_types)) * 0.7

    for i, cell_type in enumerate(cell_types):
        old_type = cell_type.replace('+', '')  
        values = [CpG_df[f'{old_type}_{j}_Beta'] for j in range(1, 4)]
        means.append(np.mean(values))
        sems.append(np.std(values) / np.sqrt(len(values)))
        individual_values.append(values)

    bars = axes[idx, 4].bar(x, means, yerr=sems, capsize=5, alpha=0.7, color=colors, width=0.4)

    baseline = np.mean(means[:4])  
    
    axes[idx, 4].axhline(y=baseline, color='gray', linestyle='--', alpha=0.7)

    for i, values in enumerate(individual_values):
        axes[idx, 4].scatter(np.repeat(x[i], len(values)) + np.random.normal(0, 0.05, len(values)), 
                         values, color='black', alpha=0.6)

    axes[idx, 4].set_xticks(x, cell_types, rotation=45, fontweight='bold', fontsize=8)
    axes[idx, 4].set_yticks(np.arange(0, 1.1, 0.1))
    axes[idx, 4].set_yticklabels([f'{x:.1f}' for x in np.arange(0, 1.1, 0.1)], fontweight='bold', fontsize=8)
    axes[idx, 4].tick_params(axis='both', which='major', labelsize=8, width=2)
    axes[idx, 4].set_ylabel('Beta Value', fontweight='bold', fontsize=8)
    
    title_bbox = dict(
        facecolor='#FFE6E6',
        edgecolor='black',
        pad=5,
        boxstyle='square,pad=0.5'
    )
    axes[idx, 4].set_title(f'{CpG_ID_Mono[idx]}', fontweight='bold', fontsize=8,
                          bbox=title_bbox,
                          x=0.5,
                          y=1.02)

    axes[idx, 4].spines['right'].set_visible(False)
    axes[idx, 4].spines['top'].set_visible(False)
    axes[idx, 4].spines['left'].set_linewidth(2)
    axes[idx, 4].spines['bottom'].set_linewidth(2)

    if idx == 2:
        axes[idx, 4].set_xlabel('Mono', fontweight='bold')


plt.tight_layout()
plt.show()



# Build the Atlas for traditional PBMC cell types (CD4+T, CD8+T, B, NK, Mono) (for version 2, we added Granu)
Final_850K_probeID_withGranu = A_exclusive_list + B_exclusive_list + C_exclusive_list + D_exclusive_list + E_exclusive_list + F_exclusive_list + G_exclusive_list


classifications_withGranu = ['CD4+T diff marker'] * len(A_exclusive_list) + \
             ['CD8+T diff marker'] * len(B_exclusive_list) + \
             ['B diff marker'] * len(C_exclusive_list) + \
             ['NK diff marker'] * len(D_exclusive_list) + \
             ['Mono diff marker'] * len(E_exclusive_list) + \
             ['Granu diff marker'] * len(F_exclusive_list) + \
             ['CD4+T vs CD8+T diff marker'] * len(G_exclusive_list)

# Create the Atlas for traditional PBMC cell types (CD4+T, CD8+T, B, NK, Mono) (for version 2, we added Granu)
Atlas_850K_withGranu = pd.DataFrame({
    'ProbeID': Final_850K_probeID_withGranu,
    'Class': classifications_withGranu
})


Atlas_850K_withGranu = pd.merge(Atlas_850K_withGranu, DMP_850K_final_withGranu, on='ProbeID', how='left').iloc[:,[0,2,1]]
Atlas_850K_withGranu.to_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_850K_withGranu.csv', index=False)


# Read the WGBS data and the Atlas
WGBS_allsamples_combined_filtered_withGranu = pd.read_csv('/Users/scui2/DNAmethylation/Loyfer_data/WGBS_allsamples_combined_filtered_withGranu.csv')
Atlas_850K_withGranu = pd.read_csv('/Users/scui2/DNAmethylation/Atlas/Atlas_850K_withGranu.csv', index_col=None)



# Calculate the mean beta value for each DMR in the traditional PBMC Atlas
def get_850K_ProbeID_mean_beta(WGBS_df):
    
    output_mean_beta_list_CD4T = []
    output_mean_beta_list_CD8T = []
    output_mean_beta_list_B = []
    output_mean_beta_list_NK = []
    output_mean_beta_list_Mono = []
    output_mean_beta_list_Granu = []
    output_DMR_list = []
    output_class_list = []
    output_ProbeID_list = []
    for i in range(len(Atlas_850K_withGranu)):
        chr, start, end = Atlas_850K_withGranu.iloc[i,1].split("_")
        start = int(start)
        end = int(end)
        condition = (WGBS_df['Chromosome'] == chr) & (WGBS_df['Position'] >= start) & (WGBS_df['Position'] <= end)

        DMR_mean_beta_CD4T = WGBS_df[condition].iloc[:,2:5].mean().mean().round(3)
        DMR_mean_beta_CD8T = WGBS_df[condition].iloc[:,5:8].mean().mean().round(3)
        DMR_mean_beta_B = WGBS_df[condition].iloc[:,8:11].mean().mean().round(3)
        DMR_mean_beta_NK = WGBS_df[condition].iloc[:,11:14].mean().mean().round(3)
        DMR_mean_beta_Mono = WGBS_df[condition].iloc[:,14:17].mean().mean().round(3)
        DMR_mean_beta_Granu = WGBS_df[condition].iloc[:,17:20].mean().mean().round(3)
        output_mean_beta_list_CD4T += [DMR_mean_beta_CD4T]
        output_mean_beta_list_CD8T += [DMR_mean_beta_CD8T]
        output_mean_beta_list_B += [DMR_mean_beta_B]
        output_mean_beta_list_NK += [DMR_mean_beta_NK]
        output_mean_beta_list_Mono += [DMR_mean_beta_Mono]
        output_mean_beta_list_Granu += [DMR_mean_beta_Granu]
        output_ProbeID_list += [Atlas_850K_withGranu.iloc[i,0]]
        output_DMR_list += [Atlas_850K_withGranu.iloc[i,1]]
        output_class_list += [Atlas_850K_withGranu.iloc[i,2]]

        print(i)

    output_df = pd.DataFrame({"ProbeID":output_ProbeID_list,
    "Window_Location":output_DMR_list,
    "Class":output_class_list,
    "DMR_mean_beta_CD4T":output_mean_beta_list_CD4T, 
    "DMR_mean_beta_CD4T_CR":output_mean_beta_list_CD4T,
    "DMR_mean_beta_CD8T":output_mean_beta_list_CD8T,
    "DMR_mean_beta_CD8T_CR":output_mean_beta_list_CD8T,
    "DMR_mean_beta_B":output_mean_beta_list_B, 
    "DMR_mean_beta_NK":output_mean_beta_list_NK, 
    "DMR_mean_beta_Mono":output_mean_beta_list_Mono,
    "DMR_mean_beta_Granu":output_mean_beta_list_Granu})
    return output_df



Atlas_850K_withGranu_df = get_850K_ProbeID_mean_beta(WGBS_allsamples_combined_filtered_withGranu)
Atlas_850K_withGranu_df.to_csv("/Users/scui2/DNAmethylation/Atlas/Atlas_850K_withGranu_df.csv", index=None)






CD4T_850K_TCR_df = pd.read_csv("/Users/scui2/DNAmethylation/Atlas/CD4T_850K_TCR_df.csv", index_col=None)
CD8T_850K_TCR_df = pd.read_csv("/Users/scui2/DNAmethylation/Atlas/CD8T_850K_TCR_df.csv", index_col=None)


# Calculate the mean beta value for each matched probe in the Cancer Reactive TCR signatures
def get_TCRWindow_mean_beta(WGBS_df, TCR_df, cell_type):
    
    output_mean_beta_list_CD4T = []
    output_mean_beta_list_CD4T_CR = []
    output_mean_beta_list_CD8T = []
    output_mean_beta_list_CD8T_CR = []
    output_mean_beta_list_B = []
    output_mean_beta_list_NK = []
    output_mean_beta_list_Mono = []
    output_mean_beta_list_Granu = []
    output_850K_probeID_list = []
    output_850K_window_list = []
    output_class_list = []
    for i in range(TCR_df.shape[0]):
        chr, start, end = TCR_df.iloc[i,1].split("_")
        start = int(start)
        end = int(end)
        if cell_type == "CD4+T":
            percentage_change_CD4T = TCR_df.iloc[i,5]
            percentage_change_CD8T = 1
        elif cell_type == "CD8+T":
            percentage_change_CD8T = TCR_df.iloc[i,5]
            percentage_change_CD4T = 1
        condition = (WGBS_df['Chromosome'] == chr) & (WGBS_df['Position'] >= start) & (WGBS_df['Position'] <= end)
        if condition.sum() > 0:
            TCR_mean_beta_CD4T = WGBS_df[condition].iloc[:,3:6].mean().mean().round(3)
            TCR_mean_beta_CD4T_CR = np.clip((TCR_mean_beta_CD4T * percentage_change_CD4T).round(3), 0, 1)
            TCR_mean_beta_CD8T = WGBS_df[condition].iloc[:,6:9].mean().mean().round(3)
            TCR_mean_beta_CD8T_CR = np.clip((TCR_mean_beta_CD8T * percentage_change_CD8T).round(3), 0, 1)
            TCR_mean_beta_B = WGBS_df[condition].iloc[:,9:12].mean().mean().round(3)
            TCR_mean_beta_NK = WGBS_df[condition].iloc[:,12:15].mean().mean().round(3)
            TCR_mean_beta_Mono = WGBS_df[condition].iloc[:,15:18].mean().mean().round(3)
            TCR_mean_beta_Granu = WGBS_df[condition].iloc[:,18:21].mean().mean().round(3)
            output_mean_beta_list_CD4T += [TCR_mean_beta_CD4T]
            output_mean_beta_list_CD4T_CR += [TCR_mean_beta_CD4T_CR]
            output_mean_beta_list_CD8T += [TCR_mean_beta_CD8T]
            output_mean_beta_list_CD8T_CR += [TCR_mean_beta_CD8T_CR]
            output_mean_beta_list_B += [TCR_mean_beta_B]
            output_mean_beta_list_NK += [TCR_mean_beta_NK]
            output_mean_beta_list_Mono += [TCR_mean_beta_Mono]
            output_mean_beta_list_Granu += [TCR_mean_beta_Granu]
            output_850K_probeID_list += [TCR_df.iloc[i,0]]
            output_850K_window_list += [TCR_df.iloc[i,1]]
            output_class_list += [cell_type + " Cancer Reactive"]
        else:
            output_mean_beta_list_CD4T += [np.nan]
            output_mean_beta_list_CD4T_CR += [np.nan]
            output_mean_beta_list_CD8T += [np.nan]
            output_mean_beta_list_CD8T_CR += [np.nan]
            output_mean_beta_list_B += [np.nan]
            output_mean_beta_list_NK += [np.nan]
            output_mean_beta_list_Mono += [np.nan]
            output_mean_beta_list_Granu += [np.nan]
            output_850K_probeID_list += [TCR_df.iloc[i,0]]
            output_850K_window_list += [TCR_df.iloc[i,1]]
            output_class_list += [cell_type + " Cancer Reactive"]
        print(i)
    output_df = pd.DataFrame({"ProbeID":output_850K_probeID_list, 
    "Window_Location":output_850K_window_list, 
    "Class":output_class_list,
    "DMR_mean_beta_CD4T":output_mean_beta_list_CD4T,
    "DMR_mean_beta_CD4T_CR":output_mean_beta_list_CD4T_CR,
    "DMR_mean_beta_CD8T":output_mean_beta_list_CD8T,
    "DMR_mean_beta_CD8T_CR":output_mean_beta_list_CD8T_CR,
    "DMR_mean_beta_B":output_mean_beta_list_B, 
    "DMR_mean_beta_NK":output_mean_beta_list_NK, 
    "DMR_mean_beta_Mono":output_mean_beta_list_Mono,
    "DMR_mean_beta_Granu":output_mean_beta_list_Granu})
    return output_df



Atlas_CD4TCR_df_withGranu = get_TCRWindow_mean_beta(WGBS_allsamples_combined_filtered_withGranu, CD4T_850K_TCR_df, "CD4+T")
Atlas_CD4TCR_df_withGranu.dropna(inplace=True)
Atlas_CD4TCR_df_withGranu.reset_index(drop=True, inplace=True)
Atlas_CD4TCR_df_withGranu.to_csv("/Users/scui2/DNAmethylation/Atlas/Atlas_CD4TCR_withGranu.csv", index=None)


Atlas_CD8TCR_df_withGranu = get_TCRWindow_mean_beta(WGBS_allsamples_combined_filtered_withGranu, CD8T_850K_TCR_df, "CD8+T")
Atlas_CD8TCR_df_withGranu.dropna(inplace=True)
Atlas_CD8TCR_df_withGranu.reset_index(drop=True, inplace=True)
Atlas_CD8TCR_df_withGranu.to_csv("/Users/scui2/DNAmethylation/Atlas/Atlas_CD8TCR_withGranu.csv", index=None)




# Remove the probes that are already in the traditional PBMC Atlas
intersected_list = list(set(Atlas_850K_withGranu_df["ProbeID"]).intersection(set(Atlas_CD4TCR_df_withGranu["ProbeID"]))) + list(set(Atlas_850K_withGranu_df["ProbeID"]).intersection(set(Atlas_CD8TCR_df_withGranu["ProbeID"])))
need_to_drop_index = Atlas_850K_withGranu_df[Atlas_850K_withGranu_df["ProbeID"].isin(intersected_list)].index
Atlas_850K_withGranu_df.drop(index=need_to_drop_index, inplace=True)
Atlas_850K_withGranu_df.reset_index(drop=True, inplace=True)




Atlas_final_df_withGranu = pd.concat([Atlas_850K_withGranu_df, Atlas_CD4TCR_df_withGranu, Atlas_CD8TCR_df_withGranu], axis=0, ignore_index=True)
Atlas_final_df_withGranu.to_csv("/Users/scui2/DNAmethylation/Atlas/Atlas_final_df_withGranu.csv", index=None)


# Check the built Atlas' correlation


Atlas_final_df = pd.read_csv("/Users/scui2/DNAmethylation/Atlas/Atlas_final_df.csv", index_col=None)

correlation_columns = ['DMR_mean_beta_CD4T', 'DMR_mean_beta_CD4T_CR', 
                      'DMR_mean_beta_CD8T', 'DMR_mean_beta_CD8T_CR',
                      'DMR_mean_beta_B', 'DMR_mean_beta_NK', 
                      'DMR_mean_beta_Mono']

corr_matrix = Atlas_final_df[correlation_columns].corr(method='spearman')


plt.figure(figsize=(10, 8))


mask = np.triu(np.ones_like(corr_matrix), k=1)  
sns.heatmap(corr_matrix, 
            annot=True,  
            fmt='.2f',   
            cmap='RdBu_r', 
            center=0,     
            square=True,  
            mask=mask,    
            cbar_kws={'label': 'Spearman Correlation'},
            xticklabels=[col.replace('DMR_mean_beta_', '') for col in correlation_columns],
            yticklabels=[col.replace('DMR_mean_beta_', '') for col in correlation_columns])

plt.xticks(rotation=45, ha='right', weight='bold')
plt.yticks(rotation=0, weight='bold')

plt.tight_layout()

plt.title('Correlation Matrix of DNA Methylation Levels\nacross Different Cell Types', 
          pad=20, fontsize=12, fontweight='bold')

plt.show()