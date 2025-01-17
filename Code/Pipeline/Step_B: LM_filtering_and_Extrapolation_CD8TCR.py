import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from statsmodels.stats.multitest import multipletests
import requests
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy.polynomial.polynomial as poly
import statsmodels.api as sm
from numpy.polynomial import Polynomial
from sklearn.metrics import r2_score
from matplotlib.font_manager import FontProperties


"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-17 by Saishi Cui
Program: Step_B: LM_filtering_and_Extrapolation_CD8TCR.py
Purpose: Filter genes by linear model (R squared >= 0.1) and extrapolate DNA methylation for CD8+TCR signatures.
--------------------------------------------------------------------------------
Data Inputs:

- Potential CD8+TCR signatures.
- Negative correlated genes for CD8+T samples.
- Matched samples of CD8+T RNA-seq and DNAmethylation data. (Most extreme negative correlated windows)

Data Outputs:
- CD8+TCR signatures filtered by linear model (R squared >= 0.1).
- Extrapolated DNA methylation for those filtered CD8+TCR signatures.
--------------------------------------------------------------------------------
"""


# Read potential CD8+TCR signatures
method_CD8TCR_sigs = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/method_CD8TCR_sigs.csv", index_col=None)


# Previously identified CD8+TCR signatures
CD8_sigs = pd.read_csv("/Users/scui2/DNAmethylation/Cancer_reactive_T/CD8_sigs.txt", sep="\t", index_col=None)
CD8_sigs.columns = ["Gene_Name", "effect"]

a = set(method_CD8TCR_sigs["Gene_Name"]) & set(CD8_sigs["Gene_Name"])
len(a)

b = set(method_CD8TCR_sigs["Gene_Name"]) - a

method_CD8TCR_sigs_final = pd.concat([method_CD8TCR_sigs[method_CD8TCR_sigs["Gene_Name"].isin(b)], CD8_sigs], axis=0)


# Read CD8T negative correlation genes

def query_geneSymbol(Symbol):
    base_url = "http://www.ensembl.org/biomart/martservice"
    Symbol_string = ','.join(Symbol)
    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "external_gene_name" value = "{Symbol_string}"/>
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "transcript_is_canonical" />
            <Attribute name = "gene_biotype" />
            <Attribute name = "transcript_biotype" />
            <Attribute name = "external_gene_name" />
            <Attribute name = "transcription_start_site" />
            <Attribute name = "chromosome_name" />
        </Dataset>
    </Query>
    """
    response = requests.get(base_url, params={'query': xml_query})
    return response.text.strip()



Symbol_CD8T1 = method_CD8TCR_sigs_final["Gene_Name"].to_list()[:500]
Symbol_CD8T2 = method_CD8TCR_sigs_final["Gene_Name"].to_list()[500:1000]
Symbol_CD8T3 = method_CD8TCR_sigs_final["Gene_Name"].to_list()[1000:1500]
Symbol_CD8T4 = method_CD8TCR_sigs_final["Gene_Name"].to_list()[1500:2000]
Symbol_CD8T5 = method_CD8TCR_sigs_final["Gene_Name"].to_list()[2000:]


a= query_geneSymbol(Symbol_CD8T1)
b= query_geneSymbol(Symbol_CD8T2)
c= query_geneSymbol(Symbol_CD8T3)
d= query_geneSymbol(Symbol_CD8T4)
e= query_geneSymbol(Symbol_CD8T5)


results_CD8T1 = a.split('\n')
results_CD8T2 = b.split('\n')
results_CD8T3 = c.split('\n')
results_CD8T4 = d.split('\n')
results_CD8T5 = e.split('\n')

results_CD8T = results_CD8T1 + results_CD8T2 + results_CD8T3 + results_CD8T4 + results_CD8T5


columns = ['Ensembl_Gene_ID', 'Ensembl_Transcript_ID',  "Is_Canonical", "Is_protein_coding_gene",
"Is_protein_coding_transcript", "Gene_Name", "Transcription_start_site", "Chromosome"]

df_CD8T = pd.DataFrame([line.split('\t') for line in results_CD8T],columns=columns)
# Filter genes based on Ensembl attributes
df_CD8T = df_CD8T[(df_CD8T["Is_Canonical"] == "1") & (df_CD8T["Is_protein_coding_gene"] == "protein_coding") & (df_CD8T["Is_protein_coding_transcript"] == "protein_coding")].copy()

# Select relevant columns
df_CD8T = df_CD8T[["Ensembl_Gene_ID", "Gene_Name","Chromosome", "Transcription_start_site"]]

# Filter genes on chromosomes 1-22, X, and Y
df_CD8T = df_CD8T[df_CD8T["Chromosome"].isin(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"])]
df_CD8T.drop_duplicates(subset="Gene_Name", keep="first", inplace=True)
df_CD8T = pd.merge(df_CD8T, method_CD8TCR_sigs_final, left_on="Gene_Name", right_on="Gene_Name", how="left")
df_CD8T.dropna(inplace=True)



# Read CD8T filtered negative correlation genes

CD8T_filtered_neg = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_filtered_neg.csv')


condition1_CD8T = df_CD8T["Ensembl_Gene_ID"].isin(set(CD8T_filtered_neg["Neg_training_gene"]).intersection(set(df_CD8T["Ensembl_Gene_ID"])))
condition2_CD8T = CD8T_filtered_neg["Neg_training_gene"].isin(set(CD8T_filtered_neg["Neg_training_gene"]).intersection(set(df_CD8T["Ensembl_Gene_ID"])))

CD8T_signatures = pd.merge(df_CD8T[condition1_CD8T], CD8T_filtered_neg[condition2_CD8T], left_on="Ensembl_Gene_ID", right_on="Neg_training_gene", how="left").sort_values(by="Neg_training", ascending=True)[["Ensembl_Gene_ID", "Gene_Name", "Chromosome", "Transcription_start_site", "effect",  "Neg_training", "Neg_testing"]].reset_index(drop=True)



# Linear model for extrapolation

CD8T_RNA_testing_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_testing_neg_window_sorted.csv')
CD8T_Meth_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_testing.csv')
CD8T_RNA_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_testing.csv')

CD8T_RNA_training_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_neg_window_sorted.csv')
CD8T_Meth_training = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_training.csv')
CD8T_RNA_training = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_training.csv')


CD8T_Meth_neg_window_testing = pd.DataFrame()
CD8T_Meth_neg_window_training = pd.DataFrame()

for i in range(CD8T_RNA_testing_neg_window_sorted.shape[0]):

    chr, start, end = CD8T_RNA_testing_neg_window_sorted["Promoter_window"][i].split("_")
    start = int(start)
    end = int(end)
    
    condition = (CD8T_Meth_testing["Chr_hg38"] == chr) & \
               (CD8T_Meth_testing["start_hg38"] >= start) & \
               (CD8T_Meth_testing["end_hg38"] <= end)
    
    filtered_data_testing = CD8T_Meth_testing.loc[condition]
    filtered_data_training = CD8T_Meth_training.loc[condition]
    window_mean_meth_testing = filtered_data_testing.iloc[:,4:].mean(axis=0).to_frame().T
    window_mean_meth_training = filtered_data_training.iloc[:,4:].mean(axis=0).to_frame().T
    CD8T_Meth_neg_window_testing = pd.concat([CD8T_Meth_neg_window_testing, window_mean_meth_testing], ignore_index=True)
    CD8T_Meth_neg_window_training = pd.concat([CD8T_Meth_neg_window_training, window_mean_meth_training], ignore_index=True)
    print(i)


CD8T_predicting_RNA_testing = CD8T_RNA_testing_neg_window_sorted[CD8T_RNA_testing_neg_window_sorted["Ensembl_Gene_ID"].isin(CD8T_signatures["Ensembl_Gene_ID"])][["Ensembl_Gene_ID", "Gene_Name", "Promoter_window", 'N1', 'N2', 'TEMRA1', 'TEMRA2', 'EM1', 'EM2', 'GSM1831330', 'GSM1831331', 'GSM1831332', 'GSM1831333', 'GSM1831334']]
CD8T_predicting_RNA_testing.reset_index(drop=True, inplace=True)
CD8T_predicting_RNA_testing["effect"] = CD8T_signatures["effect"]

CD8T_predicting_RNA_testing = pd.concat([CD8T_predicting_RNA_testing.iloc[:,:3], CD8T_predicting_RNA_testing.iloc[:,[-1]], CD8T_predicting_RNA_testing.iloc[:,3:-1]], axis=1)


CD8T_predicting_RNA_training = CD8T_RNA_training_neg_window_sorted[CD8T_RNA_training_neg_window_sorted["Ensembl_Gene_ID"].isin(CD8T_signatures["Ensembl_Gene_ID"])][["Ensembl_Gene_ID", "Gene_Name", "Promoter_window", 'GSM2336845', 'GSM2336851', 'GSM2336874', 'GSM2336911', 'GSM2336922', 'GSM2336927', 'GSM2336935', 'GSM2336941', 'GSM2336952', 'GSM2336953', 'GSM2337002', 'GSM2337046']]
CD8T_predicting_RNA_training.reset_index(drop=True, inplace=True)
CD8T_predicting_RNA_training["effect"] = CD8T_signatures["effect"]
CD8T_predicting_RNA_training = pd.concat([CD8T_predicting_RNA_training.iloc[:,:3], CD8T_predicting_RNA_training.iloc[:,[-1]], CD8T_predicting_RNA_training.iloc[:,3:-1]], axis=1)

CD8T_predicting_RNA = pd.concat([CD8T_predicting_RNA_training,CD8T_predicting_RNA_testing.iloc[:,4:]], axis=1)


CD8T_predicting_Meth_testing = CD8T_Meth_neg_window_testing[CD8T_RNA_testing_neg_window_sorted["Ensembl_Gene_ID"].isin(CD8T_signatures["Ensembl_Gene_ID"])]
CD8T_predicting_Meth_training = CD8T_Meth_neg_window_training[CD8T_RNA_training_neg_window_sorted["Ensembl_Gene_ID"].isin(CD8T_signatures["Ensembl_Gene_ID"])]

CD8T_predicting_Meth = pd.concat([CD8T_predicting_Meth_training, CD8T_predicting_Meth_testing], axis=1)




# Filter genes by linear model (R squared >= 0.1)


need_to_keep_genes = []
need_to_keep_genes_index = []
for i in range(608):
    Ensembl_Gene_ID = CD8T_predicting_RNA_training.iloc[i, 0]
    gene_name = CD8T_predicting_RNA_training.iloc[i, 1]
    window = CD8T_predicting_RNA_training.iloc[i, 2]
    chrom, start, end = window.split('_')
    x_training = CD8T_predicting_RNA_training.iloc[i, 4:].astype(float)
    y_training = CD8T_predicting_Meth_training.iloc[i, :].astype(float)
    x_testing = CD8T_predicting_RNA_testing.iloc[i, 4:].astype(float)
    y_testing = CD8T_predicting_Meth_testing.iloc[i, :].astype(float)


    # 多项式拟合
    p1_training = Polynomial.fit(x_training, y_training, 1)
    p1_testing = Polynomial.fit(x_testing, y_testing, 1)
    r_squared_training = r2_score(y_training, p1_training(x_training))
    r_squared_testing = r2_score(y_testing, p1_testing(x_testing))
    
    if r_squared_training >= 0.1 and r_squared_testing >= 0.1:
        need_to_keep_genes.append(gene_name)
        need_to_keep_genes_index.append(i)


len(need_to_keep_genes)
CD8T_predicting_RNA = CD8T_predicting_RNA.iloc[need_to_keep_genes_index, :].reset_index(drop=True)
CD8T_predicting_Meth = CD8T_predicting_Meth.iloc[need_to_keep_genes_index, :].reset_index(drop=True)
negative_list_CD8TCR = CD8T_predicting_RNA[CD8T_predicting_RNA["effect"] == "negative"]["Gene_Name"].tolist()


# Plotting
fig, axes = plt.subplots(3, 5, figsize=(15, 9))
fig.subplots_adjust(hspace=0.3, wspace=0.3)

axes_flat = axes.flatten()
font_properties = FontProperties(weight='bold', size=6)

for i in range(15):
    Ensembl_Gene_ID = CD8T_predicting_RNA.iloc[i, 0]
    gene_name = CD8T_predicting_RNA.iloc[i, 1]
    if gene_name in negative_list_CD8TCR:
        gene_name = gene_name + "*"
    window = CD8T_predicting_RNA.iloc[i, 2]
    chrom, start, end = window.split('_')
    x_training = CD8T_predicting_RNA.iloc[i, 4:16].astype(float)
    y_training = CD8T_predicting_Meth.iloc[i, :12].astype(float)
    x_testing = CD8T_predicting_RNA.iloc[i, 16:].astype(float)
    y_testing = CD8T_predicting_Meth.iloc[i, 12:].astype(float)

    # 计算 Spearman 相关系数
    rho_training, p_value_training = spearmanr(x_training, y_training)
    rho_testing, p_value_testing = spearmanr(x_testing, y_testing)

    # 多项式拟合
    p1_training = Polynomial.fit(x_training, y_training, 1)
    p1_testing = Polynomial.fit(x_testing, y_testing, 1)
    x_min_training, x_max_training = np.min(x_training), np.max(x_training)
    x_min_testing, x_max_testing = np.min(x_testing), np.max(x_testing)
    x_range_training = x_max_training - x_min_training
    x_range_testing = x_max_testing - x_min_testing
    x_fit_training = np.linspace(x_min_training - 0.2 * x_range_training, x_max_training + 0.2 * x_range_training, 200)
    x_fit_testing = np.linspace(x_min_testing - 0.2 * x_range_testing, x_max_testing + 0.2 * x_range_testing, 200)
    y_fit_training = p1_training(x_fit_training)
    y_fit_testing = p1_testing(x_fit_testing)
    r_squared_training = r2_score(y_training, p1_training(x_training))
    r_squared_testing = r2_score(y_testing, p1_testing(x_testing))

    # 计算95% CI
    n_training = len(x_training)
    n_testing = len(x_testing)
    x_mean_training = np.mean(x_training)
    x_mean_testing = np.mean(x_testing)
    
    # 计算预测值的标准误差
    y_pred_training = p1_training(x_training)
    y_pred_testing = p1_testing(x_testing)
    mse_training = np.sum((y_training.values - y_pred_training.values) ** 2) / (n_training - 2)
    mse_testing = np.sum((y_testing.values - y_pred_testing.values) ** 2) / (n_testing - 2)
    std_err_training = np.sqrt(mse_training * (1/n_training + (x_fit_training - x_mean_training)**2 / np.sum((x_training.values - x_mean_training)**2)))
    std_err_testing = np.sqrt(mse_testing * (1/n_testing + (x_fit_testing - x_mean_testing)**2 / np.sum((x_testing.values - x_mean_testing)**2)))
    
    # 计算95% CI
    ci_training = 1.96 * std_err_training
    ci_testing = 1.96 * std_err_testing

    # 在对应的子图中绘制
    ax = axes_flat[i]

    # 设置边框，只显示左边和下边
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # 绘制散点图和拟合线
    ax.scatter(x_training, y_training, color='green', label='CD8+T training', alpha=0.8, s=8)
    ax.scatter(x_testing, y_testing, color='red', label='CD8+T testing', alpha=0.8, s=8)
    ax.plot(x_fit_training, y_fit_training, color='blue', label='Linear fit (training)')
    ax.plot(x_fit_testing, y_fit_testing, color='red', label='Linear fit (testing)')

    # 添加95% CI灰色区域
    ax.fill_between(x_fit_training, y_fit_training - ci_training, y_fit_training + ci_training, color='gray', alpha=0.6, label='95% CI (training)')
    ax.fill_between(x_fit_testing, y_fit_testing - ci_testing, y_fit_testing + ci_testing, color='gray', alpha=0.6, label='95% CI (testing)')

    # 添加注释
    ax.annotate(f'$\\rho$ train = {rho_training:.3f}, $R^2$ train = {r_squared_training:.3f}\n$\\rho$ test = {rho_testing:.3f}, $R^2$ test = {r_squared_testing:.3f}',
                xy=(0.6, 0.95), xycoords='axes fraction', color='blue',
                fontsize=6, fontweight='bold', ha='left', va='top')

    ax.tick_params(axis='both', which='major', labelsize=6)
    # 设置标题
    ax.set_title(f'{gene_name}', 
                 fontsize=10, fontweight='bold')

    # 只为左下角子图(index=6)添加标签
    if i != 10:
        ax.set_xlabel('')
        ax.set_ylabel('')
    else:
        ax.set_xlabel('Gene expression (RNA-seq, Microarray)', fontsize=8, fontweight='bold')
        ax.set_ylabel('Beta value (Methylation)', fontsize=8, fontweight='bold')

    # 设置图例
    if i == 10:
        ax.legend(loc='lower left', prop=font_properties)  



plt.tight_layout()
plt.show()


axes_flat[-1].set_visible(False)







percentage_change_list_CD8T_10pct = []
percentage_change_list_CD8T_20pct = []
for i in range(CD8T_predicting_RNA.shape[0]):
    x_training = CD8T_predicting_RNA.iloc[i, 4:16].astype(float)
    x_testing = CD8T_predicting_RNA.iloc[i, 16:].astype(float)
    y_training = CD8T_predicting_Meth.iloc[i, :12].astype(float)
    y_testing = CD8T_predicting_Meth.iloc[i, 12:].astype(float)

    # 多项式拟合
    p1_training = Polynomial.fit(x_training, y_training, 1)
    p1_testing = Polynomial.fit(x_testing, y_testing, 1)
    x_min_training, x_max_training = np.min(x_training), np.max(x_training)
    x_min_testing, x_max_testing = np.min(x_testing), np.max(x_testing)
    x_range_training = x_max_training - x_min_training
    x_range_testing = x_max_testing - x_min_testing
    x_fit_training_10pct = np.linspace(x_min_training - 0.1 * x_range_training, x_max_training + 0.1 * x_range_training, 200)
    y_fit_training_10pct = p1_training(x_fit_training_10pct)
    x_fit_testing_10pct = np.linspace(x_min_testing - 0.1 * x_range_testing, x_max_testing + 0.1 * x_range_testing, 200)
    y_fit_testing_10pct = p1_testing(x_fit_testing_10pct)

    x_fit_training_20pct = np.linspace(x_min_training - 0.2 * x_range_training, x_max_training + 0.2 * x_range_training, 200)
    y_fit_training_20pct = p1_training(x_fit_training_20pct)
    x_fit_testing_20pct = np.linspace(x_min_testing - 0.2 * x_range_testing, x_max_testing + 0.2 * x_range_testing, 200)
    y_fit_testing_20pct = p1_testing(x_fit_testing_20pct)


    y_fit_max_training_10pct = y_fit_training_10pct.max()
    y_fit_min_training_10pct = y_fit_training_10pct.min()
    y_fit_max_training_20pct = y_fit_training_20pct.max()
    y_fit_min_training_20pct = y_fit_training_20pct.min()
    y_fit_max_testing_10pct = y_fit_testing_10pct.max()
    y_fit_min_testing_10pct = y_fit_testing_10pct.min()
    y_fit_max_testing_20pct = y_fit_testing_20pct.max()
    y_fit_min_testing_20pct = y_fit_testing_20pct.min()

    if y_fit_max_training_10pct < 0:
        y_fit_max_training_10pct = 0
    if y_fit_min_training_10pct > 1:
        y_fit_min_training_10pct = 1
    if y_fit_min_training_20pct > 1:
        y_fit_min_training_20pct = 1
    if y_fit_max_testing_10pct < 0:
        y_fit_max_testing_10pct = 0
    if y_fit_min_testing_10pct > 1:
        y_fit_min_testing_10pct = 1
    if y_fit_max_testing_20pct < 0:
        y_fit_max_testing_20pct = 0
    if y_fit_min_testing_20pct > 1:
        y_fit_min_testing_20pct = 1

    if CD8T_predicting_RNA["effect"][i] == "positive":
        percentage_change_training_10pct =  1-(np.mean(y_training) - y_fit_min_training_10pct) / np.mean(y_training)
        percentage_change_testing_10pct =  1-(np.mean(y_testing) - y_fit_min_testing_10pct) / np.mean(y_testing)
        percentage_change_training_20pct =  1-(np.mean(y_training) - y_fit_min_training_20pct) / np.mean(y_training)
        percentage_change_testing_20pct =  1-(np.mean(y_testing) - y_fit_min_testing_20pct) / np.mean(y_testing)
    else:
        percentage_change_training_10pct =  1+(y_fit_max_training_10pct - np.mean(y_training)) / np.mean(y_training)
        percentage_change_testing_10pct =  1+(y_fit_max_testing_10pct - np.mean(y_testing)) / np.mean(y_testing)
        percentage_change_training_20pct =  1+(y_fit_max_training_20pct - np.mean(y_training)) / np.mean(y_training)
        percentage_change_testing_20pct =  1+(y_fit_max_testing_20pct - np.mean(y_testing)) / np.mean(y_testing)
    
    percentage_change_10pct = np.mean([percentage_change_training_10pct, percentage_change_testing_10pct])
    percentage_change_20pct = np.mean([percentage_change_training_20pct, percentage_change_testing_20pct])
    percentage_change_list_CD8T_10pct.append(percentage_change_10pct)
    percentage_change_list_CD8T_20pct.append(percentage_change_20pct)



CD8T_predicting_RNA["Percentage_change_10pct"] = percentage_change_list_CD8T_10pct
CD8T_predicting_RNA["Percentage_change_20pct"] = percentage_change_list_CD8T_20pct

CD8T_sigs_change = CD8T_predicting_RNA.iloc[:,[0,1,2,3,27,28]]


CD8T_sigs_change.to_csv('/Users/scui2/DNAmethylation/Cancer_Reactive_T/CD8T_sigs_change.csv', index=False)





##### CD8T study 2

CD8T_RNA_training_neg_window_study2_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD8T_RNA_training_neg_window_study2_sorted.csv')

CD8T_sigs_change = pd.read_csv('/Users/scui2/DNAmethylation/Cancer_Reactive_T/CD8T_sigs_change.csv')
CD8T_RNA_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_testing.csv')
CD8T_Meth_training = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_training.csv')
CD8T_Meth_testing = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_DNAMeth_External_data_3studies/CD8T_Meth_testing.csv')
        
Study2_CD8T_neg_genes = CD8T_RNA_training_neg_window_study2_sorted[CD8T_RNA_training_neg_window_study2_sorted["Spearman_Correlation"] < -0.8]["Gene_Name"].to_list()
All_CD8T_CR_sig_genes = method_CD8TCR_sigs_final["Gene_Name"].to_list()
Old_CD8T_CR_sig_genes = CD8T_sigs_change["Gene_Name"].to_list()
CD8T_genes_in_testing = CD8T_RNA_testing["Gene_Name"].to_list()

# Find new signatures, criteria: 
# 1. genes are negative correlated with CD8+T in study 2 (< -0.8)
# 2. genes are in the list of all potential CD8+T CR signatures
# 3. genes are not in the list of old CR signatures
# 4. genes are not in the list of genes in testing set

new_sigs_need_to_add =set(Study2_CD8T_neg_genes).intersection(set(All_CD8T_CR_sig_genes)).difference(set(Old_CD8T_CR_sig_genes)).difference(set(CD8T_genes_in_testing))


new_sigs_df = CD8T_RNA_training_neg_window_study2_sorted[CD8T_RNA_training_neg_window_study2_sorted["Gene_Name"].isin(new_sigs_need_to_add)].reset_index(drop=True).iloc[:,[0,1,16]]

new_sigs_effect = method_CD8TCR_sigs_final["effect"][method_CD8TCR_sigs_final["Gene_Name"].isin(new_sigs_need_to_add)].reset_index(drop=True)

new_sigs_df = pd.concat([new_sigs_df, new_sigs_effect], axis=1)

new_sigs_meth = pd.DataFrame()
for i in range(new_sigs_df.shape[0]):
    chr, start, end = new_sigs_df["Promoter_window"][i].split("_")
    start = int(start)
    end = int(end)

    condition = (CD8T_Meth_training["Chr_hg38"] == chr) & \
               (CD8T_Meth_training["start_hg38"] >= start) & \
               (CD8T_Meth_training["end_hg38"] <= end)
    
    filtered_data = CD8T_Meth_training.loc[condition]
    window_mean_meth = filtered_data.iloc[:,4:].mean(axis=0).to_frame().T
    new_sigs_meth = pd.concat([new_sigs_meth, window_mean_meth], ignore_index=True)
    print(i)


CD8T_RNA_new = CD8T_RNA_training_neg_window_study2_sorted[CD8T_RNA_training_neg_window_study2_sorted["Gene_Name"].isin(new_sigs_need_to_add)].reset_index(drop=True).iloc[:,[0,1,16] + list(range(4,16))]
CD8T_RNA_new = pd.concat([CD8T_RNA_new.iloc[:,:3], new_sigs_df["effect"], CD8T_RNA_new.iloc[:,3:]], axis=1)



# Plotting
fig, axes = plt.subplots(3, 5, figsize=(15, 10))
fig.subplots_adjust(hspace=0.3, wspace=0.3)

axes_flat = axes.flatten()
font_properties = FontProperties(weight='bold', size=6)

for i in range(15):
    Ensembl_Gene_ID = CD8T_RNA_new["Ensembl_Gene_ID"][i]
    gene_name = CD8T_RNA_new["Gene_Name"][i]
    window = CD8T_RNA_new["Promoter_window"][i]
    chrom, start, end = window.split('_')
    x_training = CD8T_RNA_new.iloc[i, 5:17].astype(float)
    y_training = new_sigs_meth.iloc[i, :].astype(float)

    # 计算 Spearman 相关系数
    rho_training, p_value_training = spearmanr(x_training, y_training)

    # 多项式拟合
    p1_training = Polynomial.fit(x_training, y_training, 1)
    x_min_training, x_max_training = np.min(x_training), np.max(x_training)
    x_range_training = x_max_training - x_min_training
    x_fit_training = np.linspace(x_min_training - 0.2 * x_range_training, x_max_training + 0.2 * x_range_training, 200)
    y_fit_training = p1_training(x_fit_training)
    r_squared_training = r2_score(y_training, p1_training(x_training))

    # 计算95% CI
    n_training = len(x_training)
    x_mean_training = np.mean(x_training)
    
    # 计算预测值的标准误差
    y_pred_training = p1_training(x_training)
    mse_training = np.sum((y_training.values - y_pred_training.values) ** 2) / (n_training - 2)
    std_err_training = np.sqrt(mse_training * (1/n_training + (x_fit_training - x_mean_training)**2 / np.sum((x_training.values - x_mean_training)**2)))
    
    # 计算95% CI
    ci_training = 1.96 * std_err_training

    # 在对应的子图中绘制
    ax = axes_flat[i]

    # 设置边框，只显示左边和下边
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # 绘制散点图和拟合线
    ax.scatter(x_training, y_training, color='green', label='CD8+T training', alpha=0.8, s=8)
    ax.plot(x_fit_training, y_fit_training, color='blue', label='Linear fit (training)')

    # 添加95% CI灰色区域
    ax.fill_between(x_fit_training, y_fit_training - ci_training, y_fit_training + ci_training, color='gray', alpha=0.6, label='95% CI (training)')

    # 添加注释
    ax.annotate(f'$\\rho$ train = {rho_training:.3f}\n$R^2$ train = {r_squared_training:.3f}',
                xy=(0.6, 0.95), xycoords='axes fraction', color='blue',
                fontsize=6, fontweight='bold', ha='left', va='top')

    ax.tick_params(axis='both', which='major', labelsize=6)
    # 设置标题
    ax.set_title(f'{gene_name}', 
                 fontsize=10, fontweight='bold')

    # 只为左下角子图(index=6)添加标签
    if i != 10:
        ax.set_xlabel('')
        ax.set_ylabel('')
    else:
        ax.set_xlabel('Gene expression (RNA-seq, Microarray)', fontsize=8, fontweight='bold')
        ax.set_ylabel('Beta value (Methylation)', fontsize=8, fontweight='bold')

    # 设置图例
    if i == 10:
        ax.legend(loc='lower left', prop=font_properties)


plt.tight_layout()
plt.show()




# Extrapolation

percentage_change_list_CD8T_new_sigs_10pct = []
percentage_change_list_CD8T_new_sigs_20pct = []
for i in range(CD8T_RNA_new.shape[0]):
    x_new_sigs = CD8T_RNA_new.iloc[i, 4:16].astype(float)
    y_new_sigs = new_sigs_meth.iloc[i, :].astype(float)

    # 多项式拟合
    p1_new_sigs = Polynomial.fit(x_new_sigs, y_new_sigs, 1)
    x_min_new_sigs, x_max_new_sigs = np.min(x_new_sigs), np.max(x_new_sigs)
    x_range_new_sigs = x_max_new_sigs - x_min_new_sigs
    x_fit_new_sigs_20pct = np.linspace(x_min_new_sigs - 0.2 * x_range_new_sigs, x_max_new_sigs + 0.2 * x_range_new_sigs, 200)
    y_fit_new_sigs_20pct = p1_new_sigs(x_fit_new_sigs_20pct)

    x_fit_new_sigs_10pct = np.linspace(x_min_new_sigs - 0.1 * x_range_new_sigs, x_max_new_sigs + 0.1 * x_range_new_sigs, 200)
    y_fit_new_sigs_10pct = p1_new_sigs(x_fit_new_sigs_10pct)



    y_fit_max_new_sigs_20pct = y_fit_new_sigs_20pct.max()
    y_fit_min_new_sigs_20pct = y_fit_new_sigs_20pct.min()

    y_fit_max_new_sigs_10pct = y_fit_new_sigs_10pct.max()
    y_fit_min_new_sigs_10pct = y_fit_new_sigs_10pct.min()

    if y_fit_max_new_sigs_20pct < 0:
        y_fit_max_new_sigs_20pct = 0
    if y_fit_min_new_sigs_20pct > 1:
        y_fit_min_new_sigs_20pct = 1
    if y_fit_max_new_sigs_10pct < 0:
        y_fit_max_new_sigs_10pct = 0
    if y_fit_min_new_sigs_10pct > 1:
        y_fit_min_new_sigs_10pct = 1

    if CD8T_RNA_new["effect"][i] == "positive":
        percentage_change_new_sigs_10pct =  1-(np.mean(y_new_sigs) - y_fit_min_new_sigs_10pct) / np.mean(y_new_sigs)
        percentage_change_new_sigs_20pct =  1-(np.mean(y_new_sigs) - y_fit_min_new_sigs_20pct) / np.mean(y_new_sigs)
    else:
        percentage_change_new_sigs_10pct =  1+(y_fit_max_new_sigs_10pct - np.mean(y_new_sigs)) / np.mean(y_new_sigs)
        percentage_change_new_sigs_20pct =  1+(y_fit_max_new_sigs_20pct - np.mean(y_new_sigs)) / np.mean(y_new_sigs)


    percentage_change_new_sigs_10pct = np.mean([percentage_change_new_sigs_10pct, percentage_change_new_sigs_20pct])
    percentage_change_new_sigs_20pct = np.mean([percentage_change_new_sigs_10pct, percentage_change_new_sigs_20pct])
    percentage_change_list_CD8T_new_sigs_10pct.append(percentage_change_new_sigs_10pct)
    percentage_change_list_CD8T_new_sigs_20pct.append(percentage_change_new_sigs_20pct)



CD8T_RNA_new["Percentage_change_new_sigs_10pct"] = percentage_change_list_CD8T_new_sigs_10pct
CD8T_RNA_new["Percentage_change_new_sigs_20pct"] = percentage_change_list_CD8T_new_sigs_20pct

CD8T_sigs_change_new = CD8T_RNA_new.iloc[:,[0,1,2,3,16,17]]


CD8T_sigs_change_new.rename(columns={"Percentage_change_new_sigs_10pct":"Percentage_change_10pct", \
                                    "Percentage_change_new_sigs_20pct":"Percentage_change_20pct"}, inplace=True)

CD8T_sigs_change_combined = pd.concat([CD8T_sigs_change, CD8T_sigs_change_new], axis=0, ignore_index=True)


CD8T_sigs_change_combined.to_csv('/Users/scui2/DNAmethylation/Cancer_Reactive_T/CD8T_sigs_change_combined.csv', index=False)
