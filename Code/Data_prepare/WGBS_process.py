import pandas as pd
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from scipy import stats

"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: WGBS_process.py
Purpose: Process WGBS data from Loyfer et al., 2023.
--------------------------------------------------------------------------------
Data Inputs:

- Path to the bigwig files of the WGBS PBMC DNAmethylation data. (3 CD4+T samples, 3 CD8+T samples, 3 B cells samples, 3 NK cells samples, 3 Monocytes samples, 3 Granulocytes samples)
- Path to the beta files of the WGBS PBMC DNAmethylation data. (3 CD4+T samples, 3 CD8+T samples, 3 B cells samples, 3 NK cells samples, 3 Monocytes samples, 3 Granulocytes samples)

Data Outputs:
- Filtered and combined averaged beta values of WGBS PBMC DNAmethylation data for each cell type.
- Filtered and combined separate samples beta values of WGBS PBMC DNAmethylation data for each cell type.

--------------------------------------------------------------------------------
Functions:
- get_beta_value: Extract beta values from bigwig files and beta files.
- categorize_coverage: Categorize coverage to remove low coverage CpG sites.
--------------------------------------------------------------------------------
Notes:
- There are two versions, one for the averaged beta values, and one for the 3 samples beta values.
- There are two versions of combined WGBS PBMC DNAmethylation data: one with Granu, one without Granu.
-------------------------------------------------------------------------------
"""



### Function of getting beta values 

def get_beta_value(bw_file_path, beta_file_path):
    bw = pyBigWig.open(bw_file_path)
    chromosomes = bw.chroms()
    data = []

    non_missing_count_by_chrom = []
  # 遍历每个染色体
    for chrom in chromosomes:
        # 获取染色体的长度
        chrom_length = chromosomes[chrom]
        
        # 获取该染色体上的所有数据点
        values = np.array(bw.values(chrom, 0, chrom_length))
        
        # 非nan且不为-1的位置
        non_nan_positions = np.arange(chrom_length)[(~np.isnan(values))][::2]
        non_missing_count_by_chrom.append(len(non_nan_positions))

        values = np.where(values == -1, np.nan, values)

        # 使用 zip 来同时迭代 coverage 和 non_nan_positions
        data.extend([(chrom, pos+1, round(values[pos],3)) for pos in non_nan_positions])

        print("finished chrom: ", chrom)
    # 关闭 bigWig 文件
    bw.close()
    beta_file = np.fromfile(beta_file_path, dtype=np.uint8).reshape((-1, 2))
    df = pd.DataFrame(data, columns=['Chromosome', 'Position', 'Beta_Value'])
    df["Coverage"] = beta_file[:,1]
    return df



### Paths to the bigwig files and beta files
bw_file_path_list_CD4T = ['/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_CD4T_1_hg38.bigwig',
                         '/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_CD4T_2_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_CD4T_3_hg38.bigwig']


bw_file_path_list_CD8T = ['/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_CD8T_1_hg38.bigwig',
                            '/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_CD8T_2_hg38.bigwig',
                            '/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_CD8T_3_hg38.bigwig']


bw_file_path_list_B = ['/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_B_1_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_B_2_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_B_3_hg38.bigwig']


bw_file_path_list_NK = ['/Users/scui2/DNAmethylation/Loyfer_data/NK_data/Loyfer_NK_1_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/NK_data/Loyfer_NK_2_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/NK_data/Loyfer_NK_3_hg38.bigwig']


bw_file_path_list_Mono = ['/Users/scui2/DNAmethylation/Loyfer_data/Mono_data/Loyfer_Mono_1_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/Mono_data/Loyfer_Mono_2_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/Mono_data/Loyfer_Mono_3_hg38.bigwig']


bw_file_path_list_Granu = ['/Users/scui2/DNAmethylation/Loyfer_data/Granu_data/Loyfer_Granu_1_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/Granu_data/Loyfer_Granu_2_hg38.bigwig',
                        '/Users/scui2/DNAmethylation/Loyfer_data/Granu_data/Loyfer_Granu_3_hg38.bigwig']


beta_file_path_list_CD4T = ["/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_CD4T_1_hg38.beta",
                            "/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_CD4T_2_hg38.beta",
                            "/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_CD4T_3_hg38.beta"]


beta_file_path_list_CD8T = ['/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_CD8T_1_hg38.beta',
                            '/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_CD8T_2_hg38.beta',
                            '/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_CD8T_3_hg38.beta']

beta_file_path_list_B = ['/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_B_1_hg38.beta',
                        '/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_B_2_hg38.beta',
                        '/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_B_3_hg38.beta']

beta_file_path_list_NK = ['/Users/scui2/DNAmethylation/Loyfer_data/NK_data/Loyfer_NK_1_hg38.beta',
                        '/Users/scui2/DNAmethylation/Loyfer_data/NK_data/Loyfer_NK_2_hg38.beta',
                        '/Users/scui2/DNAmethylation/Loyfer_data/NK_data/Loyfer_NK_3_hg38.beta']

beta_file_path_list_Mono = ['/Users/scui2/DNAmethylation/Loyfer_data/Mono_data/Loyfer_Mono_1_hg38.beta',
                        '/Users/scui2/DNAmethylation/Loyfer_data/Mono_data/Loyfer_Mono_2_hg38.beta',
                        '/Users/scui2/DNAmethylation/Loyfer_data/Mono_data/Loyfer_Mono_3_hg38.beta']


beta_file_path_list_Granu = ['/Users/scui2/DNAmethylation/Loyfer_data/Granu_data/Loyfer_Granu_1_hg38.beta',
                             '/Users/scui2/DNAmethylation/Loyfer_data/Granu_data/Loyfer_Granu_2_hg38.beta',
                            '/Users/scui2/DNAmethylation/Loyfer_data/Granu_data/Loyfer_Granu_3_hg38.beta'] 


### Get beta values for each cell type and each sample

CD4T_1_df = get_beta_value(bw_file_path = bw_file_path_list_CD4T[0], beta_file_path = beta_file_path_list_CD4T[0])
CD4T_2_df = get_beta_value(bw_file_path = bw_file_path_list_CD4T[1], beta_file_path = beta_file_path_list_CD4T[1])
CD4T_3_df = get_beta_value(bw_file_path = bw_file_path_list_CD4T[2], beta_file_path = beta_file_path_list_CD4T[2])

CD8T_1_df = get_beta_value(bw_file_path = bw_file_path_list_CD8T[0], beta_file_path = beta_file_path_list_CD8T[0])
CD8T_2_df = get_beta_value(bw_file_path = bw_file_path_list_CD8T[1], beta_file_path = beta_file_path_list_CD8T[1])
CD8T_3_df = get_beta_value(bw_file_path = bw_file_path_list_CD8T[2], beta_file_path = beta_file_path_list_CD8T[2])

B_1_df = get_beta_value(bw_file_path = bw_file_path_list_B[0], beta_file_path = beta_file_path_list_B[0])
B_2_df = get_beta_value(bw_file_path = bw_file_path_list_B[1], beta_file_path = beta_file_path_list_B[1])
B_3_df = get_beta_value(bw_file_path = bw_file_path_list_B[2], beta_file_path = beta_file_path_list_B[2])

NK_1_df = get_beta_value(bw_file_path = bw_file_path_list_NK[0], beta_file_path = beta_file_path_list_NK[0])
NK_2_df = get_beta_value(bw_file_path = bw_file_path_list_NK[1], beta_file_path = beta_file_path_list_NK[1])
NK_3_df = get_beta_value(bw_file_path = bw_file_path_list_NK[2], beta_file_path = beta_file_path_list_NK[2])

Mono_1_df = get_beta_value(bw_file_path = bw_file_path_list_Mono[0], beta_file_path = beta_file_path_list_Mono[0])
Mono_2_df = get_beta_value(bw_file_path = bw_file_path_list_Mono[1], beta_file_path = beta_file_path_list_Mono[1])
Mono_3_df = get_beta_value(bw_file_path = bw_file_path_list_Mono[2], beta_file_path = beta_file_path_list_Mono[2])

Granu_1_df = get_beta_value(bw_file_path = bw_file_path_list_Granu[0], beta_file_path = beta_file_path_list_Granu[0])
Granu_2_df = get_beta_value(bw_file_path = bw_file_path_list_Granu[1], beta_file_path = beta_file_path_list_Granu[1])
Granu_3_df = get_beta_value(bw_file_path = bw_file_path_list_Granu[2], beta_file_path = beta_file_path_list_Granu[2])



### Get averaged beta values for each cell type ---- XXX_avg_df means averaged beta values for XXX cell type

CD4T_df = pd.concat([CD4T_1_df.iloc[:,:2],
    CD4T_1_df.iloc[:,2],
    CD4T_2_df.iloc[:,2],
    CD4T_3_df.iloc[:,2],
    CD4T_1_df.iloc[:,3],
    CD4T_2_df.iloc[:,3],
    CD4T_3_df.iloc[:,3]], axis=1)

sample_coverage_mean = CD4T_df.iloc[:,5:].mean(axis=1)
sample_beta_mean = CD4T_df.iloc[:,2:5].mean(axis=1)
CD4T_avg_df = pd.concat([CD4T_df.iloc[:,:2], sample_coverage_mean, sample_beta_mean], axis=1)
CD4T_avg_df.columns = ['Chromosome', 'Position', 'Coverage_Mean', 'Beta_Mean']
sample_coverage_mean.describe()


CD8T_df = pd.concat([CD8T_1_df.iloc[:,:2],
    CD8T_1_df.iloc[:,2],
    CD8T_2_df.iloc[:,2],
    CD8T_3_df.iloc[:,2],
    CD8T_1_df.iloc[:,3],
    CD8T_2_df.iloc[:,3],
    CD8T_3_df.iloc[:,3]], axis=1)
sample_coverage_mean = CD8T_df.iloc[:,5:].mean(axis=1)
sample_beta_mean = CD8T_df.iloc[:,2:5].mean(axis=1)
CD8T_avg_df = pd.concat([CD8T_df.iloc[:,:2], sample_coverage_mean, sample_beta_mean], axis=1)
CD8T_avg_df.columns = ['Chromosome', 'Position', 'Coverage_Mean', 'Beta_Mean']
sample_coverage_mean.describe()




B_df = pd.concat([B_1_df.iloc[:,:2],
    B_1_df.iloc[:,2],
    B_2_df.iloc[:,2],
    B_3_df.iloc[:,2],
    B_1_df.iloc[:,3],
    B_2_df.iloc[:,3],
    B_3_df.iloc[:,3]], axis=1)
sample_coverage_mean = B_df.iloc[:,5:].mean(axis=1)
sample_beta_mean = B_df.iloc[:,2:5].mean(axis=1)
B_avg_df = pd.concat([B_df.iloc[:,:2], sample_coverage_mean, sample_beta_mean], axis=1)
B_avg_df.columns = ['Chromosome', 'Position', 'Coverage_Mean', 'Beta_Mean']
sample_coverage_mean.describe()





NK_df = pd.concat([NK_1_df.iloc[:,:2],
    NK_1_df.iloc[:,2],
    NK_2_df.iloc[:,2],
    NK_3_df.iloc[:,2],
    NK_1_df.iloc[:,3],
    NK_2_df.iloc[:,3],
    NK_3_df.iloc[:,3]], axis=1)
sample_coverage_mean = NK_df.iloc[:,5:].mean(axis=1)
sample_beta_mean = NK_df.iloc[:,2:5].mean(axis=1)
NK_avg_df = pd.concat([NK_df.iloc[:,:2], sample_coverage_mean, sample_beta_mean], axis=1)
NK_avg_df.columns = ['Chromosome', 'Position', 'Coverage_Mean', 'Beta_Mean']
sample_coverage_mean.describe()





Mono_df = pd.concat([Mono_1_df.iloc[:,:2],
    Mono_1_df.iloc[:,2],
    Mono_2_df.iloc[:,2],
    Mono_3_df.iloc[:,2],
    Mono_1_df.iloc[:,3],
    Mono_2_df.iloc[:,3],
    Mono_3_df.iloc[:,3]], axis=1)
sample_coverage_mean = Mono_df.iloc[:,5:].mean(axis=1)
sample_beta_mean = Mono_df.iloc[:,2:5].mean(axis=1)
Mono_avg_df = pd.concat([Mono_df.iloc[:,:2], sample_coverage_mean, sample_beta_mean], axis=1)
Mono_avg_df.columns = ['Chromosome', 'Position', 'Coverage_Mean', 'Beta_Mean']
sample_coverage_mean.describe()





Granu_df = pd.concat([Granu_1_df.iloc[:,:2],
    Granu_1_df.iloc[:,2],
    Granu_2_df.iloc[:,2],
    Granu_3_df.iloc[:,2],
    Granu_1_df.iloc[:,3],
    Granu_2_df.iloc[:,3],
    Granu_3_df.iloc[:,3]], axis=1)
sample_coverage_mean = Granu_df.iloc[:,5:].mean(axis=1)
sample_beta_mean = Granu_df.iloc[:,2:5].mean(axis=1)
Granu_avg_df = pd.concat([Granu_df.iloc[:,:2], sample_coverage_mean, sample_beta_mean], axis=1)
Granu_avg_df.columns = ['Chromosome', 'Position', 'Coverage_Mean', 'Beta_Mean']
sample_coverage_mean.describe()




### Get 3 samples beta values separately for each cell type, not averaged


CD4T_df_3samples = pd.concat([CD4T_1_df.iloc[:,:2],
                              CD4T_1_df.iloc[:,2],
                              CD4T_1_df.iloc[:,3],
                              CD4T_2_df.iloc[:,2],
                              CD4T_2_df.iloc[:,3],
                              CD4T_3_df.iloc[:,2],
                              CD4T_3_df.iloc[:,3]], axis=1)

CD4T_df_3samples.columns = ['Chromosome', 'Position', "CD4T_1_Beta", "CD4T_1_Coverage", "CD4T_2_Beta", "CD4T_2_Coverage", "CD4T_3_Beta", "CD4T_3_Coverage"]



CD8T_df_3samples = pd.concat([CD8T_1_df.iloc[:,:2],
    CD8T_1_df.iloc[:,2],
    CD8T_1_df.iloc[:,3],
    CD8T_2_df.iloc[:,2],
    CD8T_2_df.iloc[:,3],
    CD8T_3_df.iloc[:,2],
    CD8T_3_df.iloc[:,3]], axis=1)

CD8T_df_3samples.columns = ['Chromosome', 'Position', "CD8T_1_Beta", "CD8T_1_Coverage", "CD8T_2_Beta", "CD8T_2_Coverage", "CD8T_3_Beta", "CD8T_3_Coverage"]



B_df_3samples = pd.concat([B_1_df.iloc[:,:2],
    B_1_df.iloc[:,2],
    B_1_df.iloc[:,3],
    B_2_df.iloc[:,2],
    B_2_df.iloc[:,3],
    B_3_df.iloc[:,2],
    B_3_df.iloc[:,3]], axis=1)

B_df_3samples.columns = ['Chromosome', 'Position', "B_1_Beta", "B_1_Coverage", "B_2_Beta", "B_2_Coverage", "B_3_Beta", "B_3_Coverage"]




NK_df_3samples = pd.concat([NK_1_df.iloc[:,:2],
    NK_1_df.iloc[:,2],
    NK_1_df.iloc[:,3],
    NK_2_df.iloc[:,2],
    NK_2_df.iloc[:,3],
    NK_3_df.iloc[:,2],
    NK_3_df.iloc[:,3]], axis=1)

NK_df_3samples.columns = ['Chromosome', 'Position', "NK_1_Beta", "NK_1_Coverage", "NK_2_Beta", "NK_2_Coverage", "NK_3_Beta", "NK_3_Coverage"]


Mono_df_3samples = pd.concat([Mono_1_df.iloc[:,:2],
    Mono_1_df.iloc[:,2],
    Mono_1_df.iloc[:,3],
    Mono_2_df.iloc[:,2],
    Mono_2_df.iloc[:,3],
    Mono_3_df.iloc[:,2],
    Mono_3_df.iloc[:,3]], axis=1)

Mono_df_3samples.columns = ['Chromosome', 'Position', "Mono_1_Beta", "Mono_1_Coverage", "Mono_2_Beta", "Mono_2_Coverage", "Mono_3_Beta", "Mono_3_Coverage"]



Granu_df_3samples = pd.concat([Granu_1_df.iloc[:,:2],
    Granu_1_df.iloc[:,2],
    Granu_1_df.iloc[:,3],
    Granu_2_df.iloc[:,2],
    Granu_2_df.iloc[:,3],
    Granu_3_df.iloc[:,2],
    Granu_3_df.iloc[:,3]], axis=1)

Granu_df_3samples.columns = ['Chromosome', 'Position', "Granu_1_Beta", "Granu_1_Coverage", "Granu_2_Beta", "Granu_2_Coverage", "Granu_3_Beta", "Granu_3_Coverage"]





### Categorize coverage to remove low coverage CpG sites

def categorize_coverage(coverage):
    if 0 <= coverage < 10:
        return '<10'
    elif 10 <= coverage < 20:
        return '10-20'
    elif 20 <= coverage < 30:
        return '20-30'
    elif 30 <= coverage <= 40:
        return '30-40'
    else:
        return '>40'


### Apply the function to categorize coverage

CD4T_avg_df['Coverage_Category'] = CD4T_avg_df['Coverage_Mean'].apply(categorize_coverage)
CD4T_avg_df.iloc[:,[2,3]] = CD4T_avg_df.iloc[:,[2,3]].round(3)
CD4T_avg_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_WGBS_mean_CD4T_df_all.csv', index=False)

CD8T_avg_df['Coverage_Category'] = CD8T_avg_df['Coverage_Mean'].apply(categorize_coverage)
CD8T_avg_df.iloc[:,[2,3]] = CD8T_avg_df.iloc[:,[2,3]].round(3)
CD8T_avg_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_WGBS_mean_CD8T_df_all.csv', index=False)

B_avg_df['Coverage_Category'] = B_avg_df['Coverage_Mean'].apply(categorize_coverage)
B_avg_df.iloc[:,[2,3]] = B_avg_df.iloc[:,[2,3]].round(3)
B_avg_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_WGBS_mean_B_df_all.csv', index=False)

NK_avg_df['Coverage_Category'] = NK_avg_df['Coverage_Mean'].apply(categorize_coverage)
NK_avg_df.iloc[:,[2,3]] = NK_avg_df.iloc[:,[2,3]].round(3)

Mono_avg_df['Coverage_Category'] = Mono_avg_df['Coverage_Mean'].apply(categorize_coverage)
Mono_avg_df.iloc[:,[2,3]] = Mono_avg_df.iloc[:,[2,3]].round(3)


Granu_avg_df['Coverage_Category'] = Granu_avg_df['Coverage_Mean'].apply(categorize_coverage)
Granu_avg_df.iloc[:,[2,3]] = Granu_avg_df.iloc[:,[2,3]].round(3)


### Visualize the beta value distribution by coverage category

sns.set_style("whitegrid")
colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8']
order = ['<10', '10-20', '20-30', '30-40', '>40']

category_counts = CD4T_avg_df['Coverage_Category'].value_counts()
total_count = category_counts.sum()
category_percentages = (category_counts / total_count * 100).round(1)


new_labels = [f"{cat} ({category_percentages[cat]}%)" for cat in order]


plt.figure(figsize=(12, 6))
sns.violinplot(x='Coverage_Category', y='Beta_Mean', data=CD4T_avg_df,
               order=order,
               hue='Coverage_Category',  
               palette=colors,  
               inner='box', 
               cut=0,  
               legend=False) 


plt.xticks(range(5), new_labels, rotation=45, ha='right')
plt.title('CD4+ T Cells Beta Value Distribution by Coverage Category (WGBS)', fontsize=12, fontweight='bold')
plt.xlabel('Coverage Category', fontsize=12, fontweight='bold')
plt.ylabel('Beta Value', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.show()

### Remove low coverage CpG sites (coverage < 10)
Loyfer_WGBS_mean_CD4T_df = CD4T_avg_df[CD4T_avg_df['Coverage_Mean'] >= 10]
Loyfer_WGBS_mean_CD8T_df = CD8T_avg_df[CD8T_avg_df['Coverage_Mean'] >= 10]
Loyfer_WGBS_mean_B_df = B_avg_df[B_avg_df['Coverage_Mean'] >= 10]
Loyfer_WGBS_mean_NK_df = NK_avg_df[NK_avg_df['Coverage_Mean'] >= 10]
Loyfer_WGBS_mean_Mono_df = Mono_avg_df[Mono_avg_df['Coverage_Mean'] >= 10]
Loyfer_WGBS_mean_Granu_df = Granu_avg_df[Granu_avg_df['Coverage_Mean'] >= 10]

Loyfer_WGBS_mean_CD4T_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/CD4T_data/Loyfer_WGBS_mean_CD4T_df_10.csv', index=False)
Loyfer_WGBS_mean_CD8T_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/CD8T_data/Loyfer_WGBS_mean_CD8T_df_10.csv', index=False)
Loyfer_WGBS_mean_B_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/B_data/Loyfer_WGBS_mean_B_df_10.csv', index=False)
Loyfer_WGBS_mean_NK_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/NK_data/Loyfer_WGBS_mean_NK_df_10.csv', index=False)
Loyfer_WGBS_mean_Mono_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/Mono_data/Loyfer_WGBS_mean_Mono_df_10.csv', index=False)
Loyfer_WGBS_mean_Granu_df.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/Granu_data/Loyfer_WGBS_mean_Granu_df_10.csv', index=False)


######## Combine all separate samples, not averaged


WGBS_allsamples_combined = pd.concat([CD4T_df_3samples, CD8T_df_3samples.iloc[:,2:], B_df_3samples.iloc[:,2:],
 NK_df_3samples.iloc[:,2:], Mono_df_3samples.iloc[:,2:], Granu_df_3samples.iloc[:,2:]], axis=1)



#### remove all samples with coverage less than 10

condition = (WGBS_allsamples_combined["CD4T_1_Coverage"] >= 10) & (WGBS_allsamples_combined["CD4T_2_Coverage"] >= 10) & \
(WGBS_allsamples_combined["CD4T_3_Coverage"] >= 10) & (WGBS_allsamples_combined["CD8T_1_Coverage"] >= 10) & \
(WGBS_allsamples_combined["CD8T_2_Coverage"] >= 10) & (WGBS_allsamples_combined["CD8T_3_Coverage"] >= 10) & \
(WGBS_allsamples_combined["B_1_Coverage"] >= 10) & (WGBS_allsamples_combined["B_2_Coverage"] >= 10) & \
(WGBS_allsamples_combined["B_3_Coverage"] >= 10) & (WGBS_allsamples_combined["NK_1_Coverage"] >= 10) & \
(WGBS_allsamples_combined["NK_2_Coverage"] >= 10) & (WGBS_allsamples_combined["NK_3_Coverage"] >= 10) & \
(WGBS_allsamples_combined["Mono_1_Coverage"] >= 10) & (WGBS_allsamples_combined["Mono_2_Coverage"] >= 10) & \
(WGBS_allsamples_combined["Mono_3_Coverage"] >= 10) & (WGBS_allsamples_combined["Granu_1_Coverage"] >= 10) & \
(WGBS_allsamples_combined["Granu_2_Coverage"] >= 10) & (WGBS_allsamples_combined["Granu_3_Coverage"] >= 10)


### Two versions of the filtered data: one with Granu, one without Granu
WGBS_allsamples_combined_filtered = WGBS_allsamples_combined[condition][["Chromosome", "Position", "CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta"]]
WGBS_allsamples_combined_filtered_withGranu = WGBS_allsamples_combined[condition][["Chromosome", "Position", "CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]]


WGBS_allsamples_combined_filtered.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/WGBS_allsamples_combined_filtered.csv', index=False)
WGBS_allsamples_combined_filtered_withGranu.to_csv('/Users/scui2/DNAmethylation/Loyfer_data/WGBS_allsamples_combined_filtered_withGranu.csv', index=False)
