import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from statsmodels.stats.multitest import multipletests
from matplotlib.font_manager import FontProperties



"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: Step_B: Sliding_Window_CD8T.py
Purpose: Using matched samples of CD8+T RNA-seq and DNAmethylation data, to find the best positive and negative correlated windows.
--------------------------------------------------------------------------------
Data Inputs:

- Matched samples of CD8+T RNA-seq and DNAmethylation data.

Data Outputs:
- The best positive and negative correlated windows with their information (gene location, window location, spearman correlation, spearman p-value).

Notes:
- Three kinds of output, one for CD8+T 3 studies combined training positive and negative windows, one for CD8+T study2 only positive and negative windows.
--------------------------------------------------------------------------------
"""



## Read external CD8T RNA-seq data and DNAmethylation data for both training and testing data
CD8T_RNA_training = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_training.csv', index_col=None)
CD8T_Meth_training = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_training.csv', index_col=None)
CD8T_RNA_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_testing.csv', index_col=None)
CD8T_Meth_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_testing.csv', index_col=None)



## Read external CD8T RNA-seq data for study2
CD8T_RNA_study2 = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_study2.csv', index_col=None)



stride = 500 
Spearman_Corr_matrix_study2 = pd.DataFrame(np.full((19, CD8T_RNA_study2.shape[0]), np.nan))
Spearman_Corr_Pvalue_matrix_study2 = pd.DataFrame(np.full((19, CD8T_RNA_study2.shape[0]), np.nan))

# Sliding window approach
for j in range(CD8T_RNA_study2.shape[0]):

    chr = "chr" + str(CD8T_RNA_study2["Chromosome"][j])
    TSS = CD8T_RNA_study2["Transcription_Start_Site"][j]

    for i in range(19):

        window_start = TSS-5000 + stride*i
        window_end = window_start + 1000
 
        condition = (CD8T_Meth_training["Chr_hg38"] == chr) & (CD8T_Meth_training["start_hg38"] >= window_start) & (CD8T_Meth_training["end_hg38"] <= window_end)
        
        
        if condition.sum() == 0:
            print(j,i)
            continue

        filtered_data = CD8T_Meth_training.loc[condition]

        window_mean_meth = filtered_data.iloc[:, 4:].mean(axis=0).to_list()
    
        
        
        corr, p_value = spearmanr(CD8T_RNA_study2.iloc[j, 4:], window_mean_meth)
        Spearman_Corr_matrix_study2.iloc[i,j] = corr
        Spearman_Corr_Pvalue_matrix_study2.iloc[i,j] = p_value
        print(j,i, "Yes")





bestpos_corr_window_index_study2 = Spearman_Corr_matrix_study2.idxmax(axis=0).to_list()
bestneg_corr_window_index_study2 = Spearman_Corr_matrix_study2.idxmin(axis=0).to_list()


# Positive correlated windows

window_info_list_pos_study2 = []
Spearman_Corr_list_pos_study2 = []
Spearman_Corr_Pvalue_list_pos_study2 = []
index_list_pos_study2 = []
for i in range(len(bestpos_corr_window_index_study2)):

    index = bestpos_corr_window_index_study2[i]
    chr = "chr" + str(CD8T_RNA_study2["Chromosome"][i])
    TSS = CD8T_RNA_study2["Transcription_Start_Site"][i]

    if np.isnan(index):
        window_info_list_pos_study2.append(np.nan)
        Spearman_Corr_list_pos_study2.append(np.nan)
        Spearman_Corr_Pvalue_list_pos_study2.append(np.nan)
        index_list_pos_study2.append(np.nan)
        continue
    else:
        index = int(index)
        index_list_pos_study2.append(index)

    window_start = TSS-5000 + stride*index
    window_end = window_start + 1000

    window_info = chr + "_" + str(int(window_start)) + "_" + str(int(window_end))
    window_info_list_pos_study2.append(window_info)

    Spearman_Corr = Spearman_Corr_matrix_study2.iloc[index, i] 
    Spearman_Corr_list_pos_study2.append(Spearman_Corr)
    Spearman_Corr_Pvalue = Spearman_Corr_Pvalue_matrix_study2.iloc[index, i]
    Spearman_Corr_Pvalue_list_pos_study2.append(Spearman_Corr_Pvalue)

    print(i)



CD8T_RNA_training_pos_window_study2 = CD8T_RNA_study2.copy()
CD8T_RNA_training_neg_window_study2 = CD8T_RNA_study2.copy()


CD8T_RNA_training_pos_window_study2["Gene_Name"] = CD8T_RNA_study2["Gene_Name"]
CD8T_RNA_training_pos_window_study2["Promoter_window"] = window_info_list_pos_study2
CD8T_RNA_training_pos_window_study2["Spearman_Correlation"] = Spearman_Corr_list_pos_study2
CD8T_RNA_training_pos_window_study2["Spearman_Correlation_pvalue"] = Spearman_Corr_Pvalue_list_pos_study2
CD8T_RNA_training_pos_window_study2["index"] = index_list_pos_study2
CD8T_RNA_training_pos_window_study2_sorted = CD8T_RNA_training_pos_window_study2.sort_values(by="Spearman_Correlation", ascending=False)
CD8T_RNA_training_pos_window_study2_sorted = CD8T_RNA_training_pos_window_study2_sorted.dropna()
CD8T_RNA_training_pos_window_study2_sorted = CD8T_RNA_training_pos_window_study2_sorted.reset_index(drop=True)




# Negative correlated windows

window_info_list_neg_study2 = []
Spearman_Corr_list_neg_study2 = []
Spearman_Corr_Pvalue_list_neg_study2 = []
index_list_neg_study2 = []
for i in range(len(bestneg_corr_window_index_study2)):

    index = bestneg_corr_window_index_study2[i]
    chr = "chr" + str(CD8T_RNA_study2["Chromosome"][i])
    TSS = CD8T_RNA_study2["Transcription_Start_Site"][i]

    if np.isnan(index):
        window_info_list_neg_study2.append(np.nan)
        Spearman_Corr_list_neg_study2.append(np.nan)
        Spearman_Corr_Pvalue_list_neg_study2.append(np.nan)
        index_list_neg_study2.append(np.nan)
        continue
    else:
        index = int(index)
        index_list_neg_study2.append(index)
    window_start = TSS-5000 + stride*index    
    window_end = window_start + 1000

    window_info = chr+ "_" + str(int(window_start)) + "_" + str(int(window_end))
    window_info_list_neg_study2.append(window_info)

    Spearman_Corr = Spearman_Corr_matrix_study2.iloc[index, i] 
    Spearman_Corr_list_neg_study2.append(Spearman_Corr)
    Spearman_Corr_Pvalue = Spearman_Corr_Pvalue_matrix_study2.iloc[index, i]
    Spearman_Corr_Pvalue_list_neg_study2.append(Spearman_Corr_Pvalue)

    print(i)



CD8T_RNA_training_neg_window_study2["Gene_Name"] = CD8T_RNA_study2["Gene_Name"]
CD8T_RNA_training_neg_window_study2["Promoter_window"] = window_info_list_neg_study2 
CD8T_RNA_training_neg_window_study2["Spearman_Correlation"] = Spearman_Corr_list_neg_study2
CD8T_RNA_training_neg_window_study2["Spearman_Correlation_pvalue"] = Spearman_Corr_Pvalue_list_neg_study2
CD8T_RNA_training_neg_window_study2["index"] = index_list_neg_study2
CD8T_RNA_training_neg_window_study2_sorted = CD8T_RNA_training_neg_window_study2.sort_values(by="Spearman_Correlation", ascending=True)
CD8T_RNA_training_neg_window_study2_sorted = CD8T_RNA_training_neg_window_study2_sorted.dropna()
CD8T_RNA_training_neg_window_study2_sorted = CD8T_RNA_training_neg_window_study2_sorted.reset_index(drop=True)


CD8T_RNA_training_neg_window_study2_sorted.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_neg_window_study2_sorted.csv', index=False)
CD8T_RNA_training_pos_window_study2_sorted.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_pos_window_study2_sorted.csv', index=False)







### For all three studies combined


CD8T_RNA_training_neg_window = CD8T_RNA_training.copy()
CD8T_RNA_training_pos_window = CD8T_RNA_training.copy()
CD8T_RNA_testing_neg_window = CD8T_RNA_testing.copy()
CD8T_RNA_testing_pos_window = CD8T_RNA_testing.copy()


Spearman_Corr_matrix_training = pd.DataFrame(np.full((19, CD8T_RNA_training.shape[0]), np.nan))
Spearman_Corr_Pvalue_matrix_training = pd.DataFrame(np.full((19, CD8T_RNA_training.shape[0]), np.nan))


# Sliding window approach
for j in range(CD8T_RNA_training.shape[0]):

    chr = "chr" + str(CD8T_RNA_training["Chromosome"][j])
    TSS = CD8T_RNA_training["Transcription_Start_Site"][j]

    for i in range(19):

        window_start = TSS-5000 + stride*i
        window_end = window_start + 1000
 
        condition = (CD8T_Meth_training["Chr_hg38"] == chr) & (CD8T_Meth_training["start_hg38"] >= window_start) & (CD8T_Meth_training["end_hg38"] <= window_end)
        
        
        if condition.sum() == 0:
            print(j,i)
            continue

        filtered_data = CD8T_Meth_training.loc[condition]

        window_mean_meth = filtered_data.iloc[:, 4:].mean(axis=0).to_list()
    
        
        
        corr, p_value = spearmanr(CD8T_RNA_training.iloc[j, 4:], window_mean_meth)
        Spearman_Corr_matrix_training.iloc[i,j] = corr
        Spearman_Corr_Pvalue_matrix_training.iloc[i,j] = p_value
        print(j,i, "Yes")



bestpos_corr_window_index_training = Spearman_Corr_matrix_training.idxmax(axis=0).to_list()
bestneg_corr_window_index_training = Spearman_Corr_matrix_training.idxmin(axis=0).to_list()



window_info_list_pos_training = []
Spearman_Corr_list_pos_training = []
Spearman_Corr_Pvalue_list_pos_training = []
index_list_pos_training = []
for i in range(len(bestpos_corr_window_index_training)):

    index = bestpos_corr_window_index_training[i]
    chr = "chr" + str(CD8T_RNA_training["Chromosome"][i])
    TSS = CD8T_RNA_training["Transcription_Start_Site"][i]

    if np.isnan(index):
        window_info_list_pos_training.append(np.nan)
        Spearman_Corr_list_pos_training.append(np.nan)
        Spearman_Corr_Pvalue_list_pos_training.append(np.nan)
        index_list_pos_training.append(np.nan)
        continue
    else:
        index = int(index)
        index_list_pos_training.append(index)

    window_start = TSS-5000 + stride*index
    window_end = window_start + 1000

    window_info = chr + "_" + str(int(window_start)) + "_" + str(int(window_end))
    window_info_list_pos_training.append(window_info)

    Spearman_Corr = Spearman_Corr_matrix_training.iloc[index, i] 
    Spearman_Corr_list_pos_training.append(Spearman_Corr)
    Spearman_Corr_Pvalue = Spearman_Corr_Pvalue_matrix_training.iloc[index, i]
    Spearman_Corr_Pvalue_list_pos_training.append(Spearman_Corr_Pvalue)

    print(i)



CD8T_RNA_training_pos_window["Gene_Name"] = CD8T_RNA_training["Gene_Name"]
CD8T_RNA_training_pos_window["Promoter_window"] = window_info_list_pos_training
CD8T_RNA_training_pos_window["Spearman_Correlation"] = Spearman_Corr_list_pos_training
CD8T_RNA_training_pos_window["Spearman_Correlation_pvalue"] = Spearman_Corr_Pvalue_list_pos_training
CD8T_RNA_training_pos_window["index"] = index_list_pos_training
CD8T_RNA_training_pos_window_sorted = CD8T_RNA_training_pos_window.sort_values(by="Spearman_Correlation", ascending=False)
CD8T_RNA_training_pos_window_sorted = CD8T_RNA_training_pos_window_sorted.dropna()
CD8T_RNA_training_pos_window_sorted = CD8T_RNA_training_pos_window_sorted.reset_index(drop=True)



CD8T_RNA_testing_pos_window["Gene_Name"] = CD8T_RNA_testing["Gene_Name"]
CD8T_RNA_testing_pos_window["Promoter_window"] = window_info_list_pos_training
CD8T_RNA_testing_pos_window["Spearman_Correlation"] = Spearman_Corr_list_pos_training
CD8T_RNA_testing_pos_window["Spearman_Correlation_pvalue"] = Spearman_Corr_Pvalue_list_pos_training
CD8T_RNA_testing_pos_window["index"] = index_list_pos_training
CD8T_RNA_testing_pos_window_sorted = CD8T_RNA_testing_pos_window.sort_values(by="Spearman_Correlation", ascending=False)
CD8T_RNA_testing_pos_window_sorted = CD8T_RNA_testing_pos_window_sorted.dropna()
CD8T_RNA_testing_pos_window_sorted = CD8T_RNA_testing_pos_window_sorted.reset_index(drop=True)




window_info_list_neg_training = []
Spearman_Corr_list_neg_training = []
Spearman_Corr_Pvalue_list_neg_training = []
index_list_neg_training = []
for i in range(len(bestneg_corr_window_index_training)):

    index = bestneg_corr_window_index_training[i]
    chr = "chr" + str(CD8T_RNA_training["Chromosome"][i])
    TSS = CD8T_RNA_training["Transcription_Start_Site"][i]

    if np.isnan(index):
        window_info_list_neg_training.append(np.nan)
        Spearman_Corr_list_neg_training.append(np.nan)
        Spearman_Corr_Pvalue_list_neg_training.append(np.nan)
        index_list_neg_training.append(np.nan)
        continue
    else:
        index = int(index)
        index_list_neg_training.append(index)

    window_start = TSS-5000 + stride*index
    window_end = window_start + 1000

    window_info = chr + "_" + str(int(window_start)) + "_" + str(int(window_end))
    window_info_list_neg_training.append(window_info)

    Spearman_Corr = Spearman_Corr_matrix_training.iloc[index, i] 
    Spearman_Corr_list_neg_training.append(Spearman_Corr)
    Spearman_Corr_Pvalue = Spearman_Corr_Pvalue_matrix_training.iloc[index, i]
    Spearman_Corr_Pvalue_list_neg_training.append(Spearman_Corr_Pvalue)

    print(i)



CD8T_RNA_training_neg_window["Gene_Name"] = CD8T_RNA_training["Gene_Name"]
CD8T_RNA_training_neg_window["Promoter_window"] = window_info_list_neg_training
CD8T_RNA_training_neg_window["Spearman_Correlation"] = Spearman_Corr_list_neg_training
CD8T_RNA_training_neg_window["Spearman_Correlation_pvalue"] = Spearman_Corr_Pvalue_list_neg_training
CD8T_RNA_training_neg_window["index"] = index_list_neg_training
CD8T_RNA_training_neg_window_sorted = CD8T_RNA_training_neg_window.sort_values(by="Spearman_Correlation", ascending=True)
CD8T_RNA_training_neg_window_sorted = CD8T_RNA_training_neg_window_sorted.dropna()
CD8T_RNA_training_neg_window_sorted = CD8T_RNA_training_neg_window_sorted.reset_index(drop=True)



CD8T_RNA_testing_neg_window["Gene_Name"] = CD8T_RNA_testing["Gene_Name"]
CD8T_RNA_testing_neg_window["Promoter_window"] = window_info_list_neg_training
CD8T_RNA_testing_neg_window["Spearman_Correlation"] = Spearman_Corr_list_neg_training
CD8T_RNA_testing_neg_window["Spearman_Correlation_pvalue"] = Spearman_Corr_Pvalue_list_neg_training
CD8T_RNA_testing_neg_window["index"] = index_list_neg_training
CD8T_RNA_testing_neg_window_sorted = CD8T_RNA_testing_neg_window.sort_values(by="Spearman_Correlation", ascending=True)
CD8T_RNA_testing_neg_window_sorted = CD8T_RNA_testing_neg_window_sorted.dropna()
CD8T_RNA_testing_neg_window_sorted = CD8T_RNA_testing_neg_window_sorted.reset_index(drop=True)



CD8T_RNA_training_neg_window_sorted.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_neg_window_sorted.csv', index=False)
CD8T_RNA_training_pos_window_sorted.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_pos_window_sorted.csv', index=False)
CD8T_RNA_testing_neg_window_sorted.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_testing_neg_window_sorted.csv', index=False)
CD8T_RNA_testing_pos_window_sorted.to_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_testing_pos_window_sorted.csv', index=False)




### Summarize the location of extreme positive and negative correlated windows

CD8T_RNA_training_pos_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_pos_window_sorted.csv')
CD8T_RNA_training_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_neg_window_sorted.csv')


index_list_neg_clean = [i for i in CD8T_RNA_training_neg_window_sorted['index'].to_list() if not np.isnan(i)]
Spearman_Corr_list_neg_clean = [i for i in CD8T_RNA_training_neg_window_sorted['Spearman_Correlation'].to_list() if not np.isnan(i)]
index_list_pos_clean = [i for i in CD8T_RNA_training_pos_window_sorted['index'].to_list() if not np.isnan(i)]
Spearman_Corr_list_pos_clean = [i for i in CD8T_RNA_training_pos_window_sorted['Spearman_Correlation'].to_list() if not np.isnan(i)]

result_neg = np.array(index_list_neg_clean) + 1 - 10
result_neg = result_neg.astype(int)

result_pos = np.array(index_list_pos_clean) + 1 - 10
result_pos = result_pos.astype(int)

x_labels = np.arange(-9, 10)

counts_neg = np.bincount(result_neg+9, minlength=19)  
counts_pos = np.bincount(result_pos+9, minlength=19)  


bar_width = 0.35

total_neg = np.sum(counts_neg)
total_pos = np.sum(counts_pos)
percent_neg = counts_neg / total_neg * 100
percent_pos = counts_pos / total_pos * 100

plt.figure(figsize=(15, 8))

font_properties = FontProperties(weight='bold')
bars_neg = plt.bar(x_labels - bar_width/2, counts_neg, bar_width, label='Negative Correlation', 
                  color='#1f77b4', edgecolor='black')  
bars_pos = plt.bar(x_labels + bar_width/2, counts_pos, bar_width, label='Positive Correlation', 
                  color='#d62728', edgecolor='black')  

plt.title('Distribution of Extreme Correlated Windows\' Locations (CD8+T training, n = 5,948)', fontsize=12, fontweight='bold')
plt.xlabel('Distance from window center to TSS ($\\times$500bp)\n$\\leftarrow$ upstream $\\leftarrow$ TSS $\\rightarrow$ downstream $\\rightarrow$', fontsize=10, labelpad=15, fontweight='bold')
plt.ylabel('Frequency', fontsize=10, fontweight='bold')
plt.xticks(x_labels, fontsize=10, fontweight='bold')   
plt.yticks(fontsize=10, fontweight='bold')
plt.legend(fontsize=10, prop=font_properties) 
plt.grid(axis='y', alpha=0.75)


plt.plot(x_labels - (bar_width/2), counts_neg, marker='o', 
         color='#1f77b4', linestyle='-', linewidth=1, alpha=0.7, markersize=4)
plt.plot(x_labels + (bar_width/2), counts_pos, marker='o', 
         color='#d62728', linestyle='-', linewidth=1, alpha=0.7, markersize=4)



def add_percentage_labels(bars, percentages, x_offset):
    for bar, percentage in zip(bars, percentages):
        height = bar.get_height()
        x_position = bar.get_x() + bar.get_width()/2 + x_offset
        plt.text(x_position, height + 1,
                 f'{percentage:.1f}%',
                 ha='center', va='bottom', fontsize=7, fontweight='bold')


add_percentage_labels(bars_neg, percent_neg, 0)  
add_percentage_labels(bars_pos, percent_pos, 0)  


plt.ylim(0, plt.ylim()[1] * 1.15)  

plt.tight_layout()
plt.show()
plt.show()




####### Split visualization 
np.quantile(Spearman_Corr_list_neg_clean, [0.33, 0.66])
np.quantile(Spearman_Corr_list_pos_clean, [0.33, 0.66])

index_list_neg_clean_high = [idx for idx, corr in zip(index_list_neg_clean, Spearman_Corr_list_neg_clean) if corr >= -0.34]
index_list_neg_clean_mid = [idx for idx, corr in zip(index_list_neg_clean, Spearman_Corr_list_neg_clean) if -0.54 <= corr < -0.34]
index_list_neg_clean_low = [idx for idx, corr in zip(index_list_neg_clean, Spearman_Corr_list_neg_clean) if corr < -0.54]


index_list_pos_clean_high = [idx for idx, corr in zip(index_list_pos_clean, Spearman_Corr_list_pos_clean) if corr > 0.47]
index_list_pos_clean_mid = [idx for idx, corr in zip(index_list_pos_clean, Spearman_Corr_list_pos_clean) if 0.25 < corr <= 0.47]
index_list_pos_clean_low = [idx for idx, corr in zip(index_list_pos_clean, Spearman_Corr_list_pos_clean) if corr <= 0.25]


result_neg_high = np.array(index_list_neg_clean_high) + 1 - 10
result_neg_high = result_neg_high.astype(int)

result_neg_mid = np.array(index_list_neg_clean_mid) + 1 - 10
result_neg_mid = result_neg_mid.astype(int)

result_neg_low = np.array(index_list_neg_clean_low) + 1 - 10
result_neg_low = result_neg_low.astype(int)

result_pos_high = np.array(index_list_pos_clean_high) + 1 - 10
result_pos_high = result_pos_high.astype(int)

result_pos_mid = np.array(index_list_pos_clean_mid) + 1 - 10
result_pos_mid = result_pos_mid.astype(int)

result_pos_low = np.array(index_list_pos_clean_low) + 1 - 10
result_pos_low = result_pos_low.astype(int)



counts_neg_high = np.bincount(result_neg_high+9, minlength=19)  
counts_neg_mid = np.bincount(result_neg_mid+9, minlength=19)  
counts_neg_low = np.bincount(result_neg_low+9, minlength=19)  
counts_pos_high = np.bincount(result_pos_high+9, minlength=19)  
counts_pos_mid = np.bincount(result_pos_mid+9, minlength=19)  
counts_pos_low = np.bincount(result_pos_low+9, minlength=19)  



bar_width = 0.35
x_labels = np.arange(-9, 10)
x_positions = x_labels * 1.5  
bar_width = 0.4 

plt.figure(figsize=(20, 8))  

font_properties = FontProperties(weight='bold')

color_blue = ['#aec7e8', '#6baed6', '#2171b5']
color_red = ['#fcae91', '#fb6a4a', '#cb181d'] 

bars_pos_low = plt.bar(x_positions - bar_width, counts_pos_low, bar_width, 
                       label='Spearman $\\rho$ < 0.25', color=color_red[0], edgecolor='black')  
bars_pos_mid = plt.bar(x_positions, counts_pos_mid, bar_width, 
                       label='0.25 $\\leq$ Spearman $\\rho$ < 0.47', color=color_red[1], edgecolor='black')  
bars_pos_high = plt.bar(x_positions + bar_width, counts_pos_high, bar_width, 
                        label='Spearman $\\rho$ $\\geq$ 0.47', color=color_red[2], edgecolor='black')  

plt.title('Distribution of Extreme Positive Correlated Windows\' Locations (CD8+T training, n = 5,948)', fontsize=12, fontweight='bold')
plt.xlabel('Distance from window center to TSS ($\\times$500bp)\n$\\leftarrow$ upstream $\\leftarrow$ TSS $\\rightarrow$ downstream $\\rightarrow$', fontsize=10, labelpad=15, fontweight='bold')
plt.ylabel('Frequency', fontsize=10, fontweight='bold')
plt.xticks(x_positions, x_labels, fontsize=10, fontweight='bold')  
plt.yticks(fontsize=10, fontweight='bold')
plt.legend(fontsize=10, loc='upper left', prop=font_properties)
plt.grid(axis='y', alpha=0.75)



def add_percentage_labels(bars, percentages, x_offset):
    for bar, percentage in zip(bars, percentages):
        height = bar.get_height()
        x_position = bar.get_x() + bar.get_width()/2 + x_offset
        plt.text(x_position, height + 1,
                 f'{percentage:.1f}%',
                 ha='center', va='bottom', fontsize=7, fontweight='bold')


total = np.sum(counts_pos_low + counts_pos_mid + counts_pos_high)
percent_pos_low = counts_pos_low / total * 100
percent_pos_mid = counts_pos_mid / total * 100
percent_pos_high = counts_pos_high / total * 100

add_percentage_labels(bars_pos_low, percent_pos_low, -0.1)
add_percentage_labels(bars_pos_mid, percent_pos_mid, 0.0)
add_percentage_labels(bars_pos_high, percent_pos_high, 0.1)

plt.ylim(0, plt.ylim()[1] * 1.15)  

plt.tight_layout()
plt.show()