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
Program: Step_B: LM_filtering_and_Extrapolation_CD4TCR.py
Purpose: Filter genes by linear model (R squared >= 0.1) and extrapolate DNA methylation for CD4+TCR signatures.
--------------------------------------------------------------------------------
Data Inputs:

- Potential CD4+TCR signatures.
- Negative correlated genes for CD4+T samples.
- Matched samples of CD4+T RNA-seq and DNAmethylation data. (Most extreme negative correlated windows)

Data Outputs:
- CD4+TCR signatures filtered by linear model (R squared >= 0.1).
- Extrapolated DNA methylation for those filtered CD4+TCR signatures.
--------------------------------------------------------------------------------
"""


# Read potential CD4+TCR signatures
method_CD4TCR_sigs = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/method_CD4TCR_sigs.csv", index_col=None)


# Previously identified CD4+TCR signatures
CD4_sigs = pd.read_csv("/Users/scui2/DNAmethylation/Cancer_reactive_T/CD4_sigs.txt", sep="\t", index_col=None)
CD4_sigs.columns = ["Gene_Name", "effect"]

# Find intersection between potential CD4+TCR signatures and previously identified CD4+TCR signatures
a = set(method_CD4TCR_sigs["Gene_Name"]) & set(CD4_sigs["Gene_Name"])
len(a)

# Find genes in potential CD4+TCR signatures that are not in previously identified CD4+TCR signatures
b = set(method_CD4TCR_sigs["Gene_Name"]) - a

# Combine potential CD4+TCR signatures and previously identified CD4+TCR signatures
method_CD4TCR_sigs_final = pd.concat([method_CD4TCR_sigs[method_CD4TCR_sigs["Gene_Name"].isin(b)], CD4_sigs], axis=0)
method_CD4TCR_sigs_final["effect"].value_counts()
method_CD4TCR_sigs_final.reset_index(drop=True, inplace=True)



# Query gene information from Ensembl
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

Symbol_CD4T = method_CD4TCR_sigs_final["Gene_Name"].to_list()
a1=query_geneSymbol(Symbol_CD4T[:500])
a2=query_geneSymbol(Symbol_CD4T[500:1000])
a3=query_geneSymbol(Symbol_CD4T[1000:1500])
a4=query_geneSymbol(Symbol_CD4T[1500:])


results_CD4T = a1.split('\n') + a2.split('\n') + a3.split('\n') + a4.split('\n')
columns = ['Ensembl_Gene_ID', 'Ensembl_Transcript_ID',  "Is_Canonical", "Is_protein_coding_gene",
"Is_protein_coding_transcript", "Symbol", "Transcription_start_site", "Chromosome"]
df_CD4T = pd.DataFrame([line.split('\t') for line in results_CD4T],columns=columns)

# Filter genes based on Ensembl attributes
df_CD4T = df_CD4T[(df_CD4T["Is_Canonical"] == "1") & (df_CD4T["Is_protein_coding_gene"] == "protein_coding") & (df_CD4T["Is_protein_coding_transcript"] == "protein_coding")].copy()

# Select relevant columns
df_CD4T = df_CD4T[["Ensembl_Gene_ID", "Symbol","Chromosome", "Transcription_start_site"]]

# Filter genes on chromosomes 1-22, X, and Y
df_CD4T = df_CD4T[df_CD4T["Chromosome"].isin(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"])]

df_CD4T.drop_duplicates(subset="Symbol", keep="first", inplace=True)


df_CD4T = pd.merge(df_CD4T, method_CD4TCR_sigs, left_on="Symbol", right_on="Gene_Name", how="left").iloc[:,[0,1,2,3,5]]

df_CD4T.dropna(inplace=True)



# Read CD4T filtered negative correlation genes

CD4T_filtered_neg = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD4T_filtered_neg.csv')

condition1_CD4T = df_CD4T["Ensembl_Gene_ID"].isin(set(CD4T_filtered_neg["Neg_training_gene"]).intersection(set(df_CD4T["Ensembl_Gene_ID"])))
condition2_CD4T = CD4T_filtered_neg["Neg_training_gene"].isin(set(CD4T_filtered_neg["Neg_training_gene"]).intersection(set(df_CD4T["Ensembl_Gene_ID"])))

CD4T_signatures = pd.merge(df_CD4T[condition1_CD4T], CD4T_filtered_neg[condition2_CD4T], left_on="Ensembl_Gene_ID", right_on="Neg_training_gene", how="left").sort_values(by="Neg_training_Spearman", ascending=True)[["Ensembl_Gene_ID", "Symbol", "Chromosome", "Transcription_start_site", "effect","Neg_training_Spearman", "Neg_testing_Spearman"]].reset_index(drop=True)



# Linear model for extrapolation


CD4T_RNA_training_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD4T_RNA_training_neg_window_sorted.csv')
CD4T_Meth_training = pd.read_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_Meth_training.csv', index_col=None)
CD4T_RNA_testing_neg_window_sorted = pd.read_csv('/Users/scui2/DNAmethylation/Corr/CD4T_RNA_testing_neg_window_sorted.csv')
CD4T_Meth_testing = pd.read_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_Meth_testing.csv', index_col=None)


CD4T_Meth_neg_window_testing = pd.DataFrame()
CD4T_Meth_neg_window_training = pd.DataFrame()

for i in range(CD4T_RNA_testing_neg_window_sorted.shape[0]):

    chr, start, end = CD4T_RNA_testing_neg_window_sorted["Promoter_window"][i].split("_")
    start = int(start)
    end = int(end)
    
    condition = (CD4T_Meth_testing["Chr_hg38"] == chr) & \
               (CD4T_Meth_testing["start_hg38"] >= start) & \
               (CD4T_Meth_testing["end_hg38"] <= end)
    
    filtered_data_testing = CD4T_Meth_testing.loc[condition]
    filtered_data_training = CD4T_Meth_training.loc[condition]
    window_mean_meth_testing = filtered_data_testing.iloc[:,4:].mean(axis=0).to_frame().T
    window_mean_meth_training = filtered_data_training.iloc[:,4:].mean(axis=0).to_frame().T
    CD4T_Meth_neg_window_testing = pd.concat([CD4T_Meth_neg_window_testing, window_mean_meth_testing], ignore_index=True)
    CD4T_Meth_neg_window_training = pd.concat([CD4T_Meth_neg_window_training, window_mean_meth_training], ignore_index=True)
    print(i)


CD4T_predicting_RNA_testing = CD4T_RNA_testing_neg_window_sorted[CD4T_RNA_testing_neg_window_sorted["Ensembl_Gene_ID"].isin(CD4T_signatures["Ensembl_Gene_ID"])].iloc[:, [0,2,43] + list(range(5, 43))]
CD4T_predicting_RNA_testing.reset_index(drop=True, inplace=True)
CD4T_predicting_RNA_testing["effect"] = CD4T_signatures["effect"]

CD4T_predicting_RNA_testing = pd.concat([CD4T_predicting_RNA_testing.iloc[:,:3], CD4T_predicting_RNA_testing.iloc[:,-1], CD4T_predicting_RNA_testing.iloc[:,3:-1]], axis=1)


CD4T_predicting_RNA_training = CD4T_RNA_training_neg_window_sorted[CD4T_RNA_training_neg_window_sorted["Ensembl_Gene_ID"].isin(CD4T_signatures["Ensembl_Gene_ID"])].iloc[:, [0,2,93] + list(range(5, 93))]
CD4T_predicting_RNA_training.reset_index(drop=True, inplace=True)
CD4T_predicting_RNA_training["effect"] = CD4T_signatures["effect"]

CD4T_predicting_RNA_training = pd.concat([CD4T_predicting_RNA_training.iloc[:,:3], CD4T_predicting_RNA_training.iloc[:,-1], CD4T_predicting_RNA_training.iloc[:,3:-1]], axis=1)

CD4T_predicting_RNA = pd.concat([CD4T_predicting_RNA_training,CD4T_predicting_RNA_testing.iloc[:,4:]], axis=1)


CD4T_predicting_Meth_testing = CD4T_Meth_neg_window_testing[CD4T_RNA_testing_neg_window_sorted["Ensembl_Gene_ID"].isin(CD4T_signatures["Ensembl_Gene_ID"])]
CD4T_predicting_Meth_training = CD4T_Meth_neg_window_training[CD4T_RNA_training_neg_window_sorted["Ensembl_Gene_ID"].isin(CD4T_signatures["Ensembl_Gene_ID"])]

CD4T_predicting_Meth = pd.concat([CD4T_predicting_Meth_training, CD4T_predicting_Meth_testing], axis=1)

CD4T_predicting_RNA.reset_index(drop=True, inplace=True)
CD4T_predicting_Meth.reset_index(drop=True, inplace=True)


# Filter genes by linear model (R squared >= 0.1)

need_to_keep_genes = []
need_to_keep_genes_index = []
for i in range(257):
    Ensembl_Gene_ID = CD4T_predicting_RNA.iloc[i, 0]
    gene_name = CD4T_predicting_RNA.iloc[i, 1]
    window = CD4T_predicting_RNA.iloc[i, 2]
    chrom, start, end = window.split('_')
    x = CD4T_predicting_RNA.iloc[i, 4:].astype(float)
    y = CD4T_predicting_Meth.iloc[i, :].astype(float)

    # Calculate Spearman correlation coefficient
    rho, p_value = spearmanr(x, y)

    # Polynomial fitting
    p1 = Polynomial.fit(x, y, 1)
    x_min, x_max = np.min(x), np.max(x)
    x_range = x_max - x_min
    x_fit = np.linspace(x_min - 0.2 * x_range, x_max + 0.2 * x_range, 200)
    y_fit = p1(x_fit)
    r_squared = r2_score(y, p1(x))
    
    if r_squared >= 0.1:
        need_to_keep_genes.append(gene_name)
        need_to_keep_genes_index.append(i)


CD4T_predicting_RNA = CD4T_predicting_RNA.iloc[need_to_keep_genes_index, :].reset_index(drop=True)
CD4T_predicting_Meth = CD4T_predicting_Meth.iloc[need_to_keep_genes_index, :].reset_index(drop=True)
negative_list_CD4TCR = CD4T_predicting_RNA[CD4T_predicting_RNA["effect"] == "negative"]["Gene_Name"].tolist()


# Plotting
fig, axes = plt.subplots(3, 5, figsize=(15, 9))
fig.subplots_adjust(hspace=0.3, wspace=0.3)

axes_flat = axes.flatten()
font_properties = FontProperties(weight='bold', size=6)

for i in range(15):
    Ensembl_Gene_ID = CD4T_predicting_RNA.iloc[i, 0]
    gene_name = CD4T_predicting_RNA.iloc[i, 1]
    if gene_name in negative_list_CD4TCR:
        gene_name = gene_name + "*"
    window = CD4T_predicting_RNA.iloc[i, 2]
    chrom, start, end = window.split('_')
    x = CD4T_predicting_RNA.iloc[i, 4:].astype(float)
    y = CD4T_predicting_Meth.iloc[i, :].astype(float)

    # 计算 Spearman 相关系数
    rho, p_value = spearmanr(x, y)

    # 多项式拟合
    p1 = Polynomial.fit(x, y, 1)
    x_min, x_max = np.min(x), np.max(x)
    x_range = x_max - x_min
    x_fit = np.linspace(x_min - 0.2 * x_range, x_max + 0.2 * x_range, 200)
    y_fit = p1(x_fit)
    r_squared = r2_score(y, p1(x))

    # 计算95% CI
    n = len(x)
    x_mean = np.mean(x)
    
    # 计算预测值的标准误差
    y_pred = p1(x)
    mse = np.sum((y.values - y_pred.values) ** 2) / (n - 2)
    std_err = np.sqrt(mse * (1/n + (x_fit - x_mean)**2 / np.sum((x.values - x_mean)**2)))
    
    # 计算95% CI
    ci = 1.96 * std_err

    # 在对应的子图中绘制
    ax = axes_flat[i]

    # 设置边框，只显示左边和下边
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # 绘制散点图和拟合线
    ax.scatter(x[:88], y[:88], color='green', label='CD4+T training', alpha=0.8, s=8)
    ax.scatter(x[88:], y[88:], color='red', label='CD4+T testing', alpha=0.8, s=8)
    ax.plot(x_fit, y_fit, color='blue', label='Linear fit')

    # 添加95% CI灰色区域
    ax.fill_between(x_fit, y_fit - ci, y_fit + ci, color='gray', alpha=0.6, label='95% CI')

    # 添加注释
    ax.annotate(f'$\\rho$ = {rho:.3f}\n$R^2$ = {r_squared:.3f}',
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
        ax.set_xlabel('Gene expression (RNA-seq, Normalized)', fontsize=8, fontweight='bold')
        ax.set_ylabel('Beta value (Methylation)', fontsize=8, fontweight='bold')

    # 设置图例
    if i == 10:
        ax.legend(loc='lower left', prop=font_properties)  

plt.tight_layout()
plt.show()



axes_flat[-1].set_visible(False)



CD4TCR_genelist = CD4T_predicting_RNA.iloc[:,:4].copy()
CD4TCR_genelist.to_csv("/Users/scui2/DNAmethylation/Cancer_Reactive_T/CD4TCR_genelist.csv", index=False)


# Extrapolation

percentage_change_list_CD4T_20pct = []
percentage_change_list_CD4T_30pct = []
percentage_change_list_CD4T_40pct = []
percentage_change_list_CD4T_50pct = []
for i in range(CD4T_predicting_RNA.shape[0]):
    x = CD4T_predicting_RNA.iloc[i, 4:].astype(float)
    y = CD4T_predicting_Meth.iloc[i, :].astype(float)

    # 多项式拟合
    p1 = Polynomial.fit(x, y, 1)
    x_min, x_max = np.min(x), np.max(x)
    x_range = x_max - x_min
    x_fit_20pct = np.linspace(x_min - 0.2 * x_range, x_max + 0.2 * x_range, 200)
    y_fit_20pct = p1(x_fit_20pct)

    x_fit_30pct = np.linspace(x_min - 0.3 * x_range, x_max + 0.3 * x_range, 200)
    y_fit_30pct = p1(x_fit_30pct)

    x_fit_40pct = np.linspace(x_min - 0.4 * x_range, x_max + 0.4 * x_range, 200)
    y_fit_40pct = p1(x_fit_40pct)

    x_fit_50pct = np.linspace(x_min - 0.5 * x_range, x_max + 0.5 * x_range, 200)
    y_fit_50pct = p1(x_fit_50pct)

    y_fit_max_20pct = y_fit_20pct.max()
    y_fit_min_20pct = y_fit_20pct.min()
    y_fit_max_30pct = y_fit_30pct.max()
    y_fit_min_30pct = y_fit_30pct.min()
    y_fit_max_40pct = y_fit_40pct.max()
    y_fit_min_40pct = y_fit_40pct.min()
    y_fit_max_50pct = y_fit_50pct.max()
    y_fit_min_50pct = y_fit_50pct.min()


    if y_fit_max_20pct > 1:
        y_fit_max_20pct = 1
    if y_fit_min_20pct < 0:
        y_fit_min_20pct = 0
    if y_fit_max_30pct > 1:
        y_fit_max_30pct = 1
    if y_fit_min_30pct < 0:
        y_fit_min_30pct = 0
    if y_fit_max_40pct > 1:
        y_fit_max_40pct = 1
    if y_fit_min_40pct < 0:
        y_fit_min_40pct = 0
    if y_fit_max_50pct > 1:
        y_fit_max_50pct = 1
    if y_fit_min_50pct < 0:
        y_fit_min_50pct = 0

    if CD4T_predicting_RNA["effect"][i] == "positive":
        percentage_change_20pct =  1-(np.mean(y) - y_fit_min_20pct) / np.mean(y)
        percentage_change_30pct =  1-(np.mean(y) - y_fit_min_30pct) / np.mean(y)
        percentage_change_40pct =  1-(np.mean(y) - y_fit_min_40pct) / np.mean(y)
        percentage_change_50pct =  1-(np.mean(y) - y_fit_min_50pct) / np.mean(y)
    else:
        percentage_change_20pct =  1+(y_fit_max_20pct - np.mean(y)) / np.mean(y)
        percentage_change_30pct =  1+(y_fit_max_30pct - np.mean(y)) / np.mean(y)
        percentage_change_40pct =  1+(y_fit_max_40pct - np.mean(y)) / np.mean(y)
        percentage_change_50pct =  1+(y_fit_max_50pct - np.mean(y)) / np.mean(y)
    
    percentage_change_list_CD4T_20pct.append(percentage_change_20pct)
    percentage_change_list_CD4T_30pct.append(percentage_change_30pct)
    percentage_change_list_CD4T_40pct.append(percentage_change_40pct)
    percentage_change_list_CD4T_50pct.append(percentage_change_50pct)


CD4T_predicting_RNA["Percentage_change_20pct"] = percentage_change_list_CD4T_20pct  
CD4T_predicting_RNA["Percentage_change_30pct"] = percentage_change_list_CD4T_30pct
CD4T_predicting_RNA["Percentage_change_40pct"] = percentage_change_list_CD4T_40pct
CD4T_predicting_RNA["Percentage_change_50pct"] = percentage_change_list_CD4T_50pct
CD4T_sigs_change = CD4T_predicting_RNA.iloc[:,[0,1,2,3,130,131,132,133]]


CD4T_sigs_change.to_csv('/Users/scui2/DNAmethylation/Cancer_Reactive_T/CD4T_sigs_change.csv', index=False)





