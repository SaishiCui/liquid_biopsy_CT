import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: External_CD8T_DNAMeth_data_prepare.py
Purpose: Prepare external CD8+T healthy samples data for DNAmethylation data, and split into training and testing data.
--------------------------------------------------------------------------------
Data Inputs:

- Three studies' CD8+T healthy samples data, including samples matched DNAmethylation data.
- Study 1: Rodriguez et al., 2017 ( N = 6)
- Study 2: Ventham et al., 2016 ( N = 12)
- Study 3: Mamrut et al., 2015 ( N = 5)

Data Outputs:
- CD8+T healthy samples data for DNAmethylation data, and split into training and testing data.

Notes:
- Training data is using the 12 samples from Ventham et al., 2016.
- Testing data is using the 11 samples from the other two studies.
--------------------------------------------------------------------------------
"""


### Study 1 (N = 6) （log2 transformed and quantile normalized）

GSE83156_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE83156_Meth.txt', sep='\t', header=0)


### Study 2 (N = 12) （quantile normalized）

GSM2336845_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336845.txt', sep='\t', header=0)
GSM2336851_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336851.txt', sep='\t', header=0)
GSM2336874_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336874.txt', sep='\t', header=0)
GSM2336911_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336911.txt', sep='\t', header=0)


GSM2336922_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336922.txt', sep='\t', header=0)
GSM2336927_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336927.txt', sep='\t', header=0)
GSM2336935_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336935.txt', sep='\t', header=0)
GSM2336941_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336941.txt', sep='\t', header=0)
GSM2336952_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336952.txt', sep='\t', header=0)
GSM2336953_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2336953.txt', sep='\t', header=0)
GSM2337002_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2337002.txt', sep='\t', header=0)
GSM2337046_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE87640/GSM2337046.txt', sep='\t', header=0)




### Study 3 (N = 5) （normalized Average Beta）

GSM1831330_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE71245/GSM1831330.txt', sep='\t', header=0)
GSM1831331_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE71245/GSM1831331.txt', sep='\t', header=0)
GSM1831332_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE71245/GSM1831332.txt', sep='\t', header=0)
GSM1831333_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE71245/GSM1831333.txt', sep='\t', header=0)
GSM1831334_Meth_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/GSE71245/GSM1831334.txt', sep='\t', header=0)



####### Merge dna methylation data #####

intersect_ID_REF = set(GSM2336845_Meth_df["ID_REF"]).intersection(set(GSE83156_Meth_df["ID_REF"])).intersection(set(GSM1831330_Meth_df["ID_REF"]))

filtered_intersect_ID_REF = [id for id in intersect_ID_REF if id.startswith('cg')]
filtered_intersect_ID_REF = pd.DataFrame(filtered_intersect_ID_REF, columns=["cg_ID"])



CD8T_Meth_info = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_info.csv', index_col=None)




merge6 = pd.merge(CD8T_Meth_info, GSE83156_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge6 = merge6[["cg_ID", "Chr_hg38", "start_hg38", "end_hg38", "CD8_NAIVE_1", "CD8_NAIVE_2",
 "CD8_EM_1", "CD8_EM_2", "CD8_TEMRA_1", "CD8_TEMRA_2"]]

merge7 = pd.merge(merge6, GSM2336845_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge7 = merge7.drop(merge7.columns[-2], axis=1)
merge7.rename(columns={"VALUE": "GSM2336845"}, inplace=True)

merge8 = pd.merge(merge7, GSM2336851_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge8 = merge8.drop(merge8.columns[-2], axis=1)
merge8.rename(columns={"VALUE": "GSM2336851"}, inplace=True)

merge9 = pd.merge(merge8, GSM2336874_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge9 = merge9.drop(merge9.columns[-2], axis=1)
merge9.rename(columns={"VALUE": "GSM2336874"}, inplace=True)

merge10 = pd.merge(merge9, GSM2336911_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge10 = merge10.drop(merge10.columns[-2], axis=1)
merge10.rename(columns={"VALUE": "GSM2336911"}, inplace=True)


merge11 = pd.merge(merge10, GSM2336922_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge11 = merge11.drop(merge11.columns[-2], axis=1)
merge11.rename(columns={"VALUE": "GSM2336922"}, inplace=True)


merge12 = pd.merge(merge11, GSM2336927_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge12 = merge12.drop(merge12.columns[-2], axis=1)
merge12.rename(columns={"VALUE": "GSM2336927"}, inplace=True)

merge13 = pd.merge(merge12, GSM2336935_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge13 = merge13.drop(merge13.columns[-2], axis=1)
merge13.rename(columns={"VALUE": "GSM2336935"}, inplace=True)


merge14 = pd.merge(merge13, GSM2336941_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge14 = merge14.drop(merge14.columns[-2], axis=1)
merge14.rename(columns={"VALUE": "GSM2336941"}, inplace=True)

merge15 = pd.merge(merge14, GSM2336952_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge15 = merge15.drop(merge15.columns[-2], axis=1)
merge15.rename(columns={"VALUE": "GSM2336952"}, inplace=True)

merge16 = pd.merge(merge15, GSM2336953_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge16 = merge16.drop(merge16.columns[-2], axis=1)
merge16.rename(columns={"VALUE": "GSM2336953"}, inplace=True)

merge17 = pd.merge(merge16, GSM2337002_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge17 = merge17.drop(merge17.columns[-2], axis=1)
merge17.rename(columns={"VALUE": "GSM2337002"}, inplace=True)

merge18 = pd.merge(merge17, GSM2337046_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge18 = merge18.drop(merge18.columns[-2], axis=1)
merge18.rename(columns={"VALUE": "GSM2337046"}, inplace=True)

merge19 = pd.merge(merge18, GSM1831330_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge19 = merge19.drop(merge19.columns[[-1,-3]], axis=1)
merge19.rename(columns={"VALUE": "GSM1831330"}, inplace=True)

merge20 = pd.merge(merge19, GSM1831331_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge20 = merge20.drop(merge20.columns[[-1,-3]], axis=1)
merge20.rename(columns={"VALUE": "GSM1831331"}, inplace=True)

merge21 = pd.merge(merge20, GSM1831332_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge21 = merge21.drop(merge21.columns[[-1,-3]], axis=1)
merge21.rename(columns={"VALUE": "GSM1831332"}, inplace=True)

merge22 = pd.merge(merge21, GSM1831333_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge22 = merge22.drop(merge22.columns[[-1,-3]], axis=1)
merge22.rename(columns={"VALUE": "GSM1831333"}, inplace=True)

merge23 = pd.merge(merge22, GSM1831334_Meth_df, left_on="cg_ID", right_on="ID_REF", how="left")
merge23 = merge23.drop(merge23.columns[[-1,-3]], axis=1)
merge23.rename(columns={"VALUE": "GSM1831334"}, inplace=True)


merge23.iloc[:, 4:] = merge23.iloc[:, 4:].round(4)

CD8T_Meth_allsamples = merge23.copy()
CD8T_Meth_training = merge23.iloc[:, [0,1,2,3,10,11,12,13,14,15,16,17,18,19,20,21]].copy()
CD8T_Meth_testing = merge23.iloc[:, [0,1,2,3,4,5,6,7,8,9,22,23,24,25,26]].copy()

CD8T_Meth_allsamples.to_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_allsamples.csv', index=False)
CD8T_Meth_training.to_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_training.csv', index=False)
CD8T_Meth_testing.to_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_testing.csv', index=False)




###### Check the sample mean and sd for three studies ######

sample_mean_beta_data1 = merge23.iloc[:,4:10].mean(axis=1)
sample_sd_beta_data1 = merge23.iloc[:,4:10].std(axis=1)

sample_mean_beta_data2 = merge23.iloc[:,10:22].mean(axis=1)
sample_sd_beta_data2 = merge23.iloc[:,10:22].std(axis=1)

sample_mean_beta_data3 = merge23.iloc[:,22:].mean(axis=1)
sample_sd_beta_data3 = merge23.iloc[:,22:].std(axis=1)

sample_mean_beta_all = merge23.iloc[:,4:].mean(axis=1)
sample_sd_beta_all = merge23.iloc[:,4:].std(axis=1)


sample_mean_beta_data1.describe()
sample_mean_beta_data2.describe()
sample_mean_beta_data3.describe()



df_sorted = np.sort(merge23.iloc[:,3:21].values, axis=0)
df_mean = np.mean(df_sorted, axis=1)
df_ranked = merge23.iloc[:,21:].rank(method="average").fillna(0).astype(int) - 1
merge_data3_quantile = pd.DataFrame(df_mean[df_ranked], index=merge23.iloc[:,21:].index, columns=merge23.iloc[:,21:].columns)


sample_mean_beta_data3_quantile = merge_data3_quantile.mean(axis=1)
sample_sd_beta_data3_quantile = merge_data3_quantile.std(axis=1)


plt.figure(figsize=(10, 8))

ax = sns.histplot(sample_mean_beta_data1, bins=50, kde=True)
plt.title('CD8+T Samples Mean Beta Value of all 450k Probes (N=23)')
plt.xlabel('Sample Mean Beta Value')
plt.ylabel('Frequency')

# 禁用自动缩放
ax.autoscale(enable=False, axis='y')

plt.tight_layout()
plt.show()




plt.figure(figsize=(10, 8))




sns.kdeplot(data=sample_sd_beta_data1, color='blue', label='GSE83156 (N=6)', shade=False)

sns.kdeplot(data=sample_sd_beta_data2, color='red', label='GSE87640 (N=12)', shade=False)

sns.kdeplot(data=sample_sd_beta_data3_quantile, color='green', label='GSE71244 (N=5)', shade=False)

plt.title('After quantile normalization of GSE71244')
plt.xlabel('Samples\' Standard Deviation of Beta Value')
plt.ylabel('Density')
plt.legend()
plt.tight_layout()
plt.show()





rho, p_value = stats.spearmanr(sample_mean_beta_data1, sample_mean_beta_data3, nan_policy='omit')

print(f"Spearman correlation coefficient (rho): {rho:.4f}")
print(f"P-value: {p_value:.4e}")



plt.figure(figsize=(10, 10))
scatter = plt.scatter(sample_mean_beta_data1, sample_mean_beta_data3, alpha=0.5, label='450K Probes')
plt.plot([0, 1], [0, 1], color='blue', linestyle='--', linewidth=2, label='Identity Line')
plt.text(0.15, 0.95, f'Spearman rho: {rho:.4f}', transform=plt.gca().transAxes, fontsize=12, ha='center', fontweight='bold', color='red')
plt.xlabel('GSE83156 (N=6) (log2 transformed and quantile normalized)')
plt.ylabel('GSE71244 (N=5) (normalized Average Beta)')
plt.title('Comparison of mean beta value between GSE83156 (N=6) and GSE71244 (N=5)')
plt.legend(loc = 'lower right')
plt.grid(True)
plt.show()


sample_mean_beta_data2.describe()
sample_sd_beta_data1.describe()




