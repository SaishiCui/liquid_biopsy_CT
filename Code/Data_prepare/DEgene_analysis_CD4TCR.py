import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from collections import Counter
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
from scipy.stats import linregress, spearmanr

"""
--------------------------------------------------------------------------------
<<<<<<< HEAD
Date last modified: 2025-01-17 by Saishi Cui
=======
Date last modified: 2025-01-16 by Saishi Cui
>>>>>>> bb3e26c (Initial commit)
Program: DEgene_analysis_CD4TCR.py
Purpose: Using external single-cell RNA-seq data of CD4+T and CD4+T cancer reactive samples,
analyze the DE genes between CD4+T cancer reactive and non-cancer reactive CD4+T samples.
--------------------------------------------------------------------------------
Data Inputs:

- Single-cell RNA-seq data of CD4+T and CD4+T cancer reactive samples, and their label information.

Data Outputs:

- DE genes between CD4+T cancer reactive and non-cancer reactive CD4+T samples.
--------------------------------------------------------------------------------
"""


# Read the label information of CD4+T cancer reactive and non-cancer reactive samples

CD4T_cell_status = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/CD4_cancer_reactive_status.tsv", sep="\t")


# We only focus on samples with more than 50 normal CD4+T and 50 CD4+T cancer reactive samples.

L_sample1_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'L-P20181123-T']
L_sample2_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'L-P20190404-T']
BC_sample1_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'BC-P20190403-T']
BC_sample2_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'BC-P20190123-T']
ESCA_sample1_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'ESCA-P20190411-T']
ESCA_sample2_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'ESCA-P20190613-T']
ESCA_sample3_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'ESCA-P20181114-T']
ESCA_sample4_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'ESCA-P20181123-T']
ESCA_sample5_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'ESCA-P20190410-T']
ESCA_sample6_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'ESCA-P20190404-T']
PACA_sample1_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'PACA-P20190306-T']
PACA_sample2_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'PACA-P20181121-T']
THCA_sample1_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'THCA-P20190730-T']
THCA_sample2_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'THCA-P20190816-T']
THCA_sample3_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'THCA-P20190703-T']
THCA_sample4_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'THCA-P20190122-T']
THCA_sample5_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'THCA-P20181226-T']
UCEC_sample1_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'UCEC-P20190910-T']
UCEC_sample2_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'UCEC-P20190312-T']
UCEC_sample3_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'UCEC-P20190305-T']
UCEC_sample4_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'UCEC-P20190213-T']
UCEC_sample5_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'UCEC-P20190625-T']
UCEC_sample6_df = CD4T_cell_status[CD4T_cell_status['library.id'] == 'UCEC-P20181122-T']




# Read the single-cell RNA-seq data for each cancer type


BCL_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/BCL.txt", sep="\t")
BC_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/BC.txt", sep="\t")
ESCA_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/ESCA.txt", sep="\t")
MM_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/MM.txt", sep="\t")
PACA_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/PACA.txt", sep="\t")
RC_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/RC.txt", sep="\t")
THCA_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/THCA.txt", sep="\t")
UCEC_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/UCEC.txt", sep="\t")






# Create the pseudobulk data for each sample.
# BCL
cells_label = BCL_df.T.index.tolist()
genes_label = BCL_df["Unnamed: 0"].tolist()
BCL_df_transposed = BCL_df.T
BCL_df_transposed = BCL_df_transposed.drop(BCL_df_transposed.index[0], axis=0)
BCL_df_transposed = BCL_df_transposed.drop(BCL_df_transposed.columns[0], axis=1)
BCL_df_transposed.columns = genes_label[1:]  
BCL_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  
BCL_df_transposed.reset_index(drop=True, inplace=True)

BCL_scRNA_sample1_df = pd.merge(L_sample1_df, BCL_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
BCL_scRNA_sample1_df = BCL_scRNA_sample1_df.iloc[:, [4,3]+list(range(5,len(BCL_scRNA_sample1_df.columns)))]
BCL_scRNA_sample1_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
BCL_scRNA_sample1_df.reset_index(drop=True, inplace=True)

pseudobulk_sample1_CD4TCR = BCL_scRNA_sample1_df[BCL_scRNA_sample1_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample1_CD4T = BCL_scRNA_sample1_df[BCL_scRNA_sample1_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


BCL_scRNA_sample2_df = pd.merge(L_sample2_df, BCL_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
BCL_scRNA_sample2_df = BCL_scRNA_sample2_df.iloc[:, [4,3]+list(range(5,len(BCL_scRNA_sample2_df.columns)))]
BCL_scRNA_sample2_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
BCL_scRNA_sample2_df.reset_index(drop=True, inplace=True)

pseudobulk_sample2_CD4TCR = BCL_scRNA_sample2_df[BCL_scRNA_sample2_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample2_CD4T = BCL_scRNA_sample2_df[BCL_scRNA_sample2_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)




# BC

cells_label = BC_df.T.index.tolist()
genes_label = BC_df["Unnamed: 0"].tolist()
BC_df_transposed = BC_df.T
BC_df_transposed = BC_df_transposed.drop(BC_df_transposed.index[0], axis=0)
BC_df_transposed = BC_df_transposed.drop(BC_df_transposed.columns[0], axis=1)
<<<<<<< HEAD
BC_df_transposed.columns = genes_label[1:]  # 去掉第一个元素，因为它是列名
BC_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  # 去掉第一个元素，因为它是行名
=======
BC_df_transposed.columns = genes_label[1:]  
BC_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  
>>>>>>> bb3e26c (Initial commit)
BC_df_transposed.reset_index(drop=True, inplace=True)


BC_scRNA_sample1_df = pd.merge(BC_sample1_df, BC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
BC_scRNA_sample1_df = BC_scRNA_sample1_df.iloc[:, [4,3]+list(range(5,len(BC_scRNA_sample1_df.columns)))]
BC_scRNA_sample1_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
BC_scRNA_sample1_df.reset_index(drop=True, inplace=True)

pseudobulk_sample3_CD4TCR = BC_scRNA_sample1_df[BC_scRNA_sample1_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample3_CD4T = BC_scRNA_sample1_df[BC_scRNA_sample1_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


BC_scRNA_sample2_df = pd.merge(BC_sample2_df, BC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
BC_scRNA_sample2_df = BC_scRNA_sample2_df.iloc[:, [4,3]+list(range(5,len(BC_scRNA_sample2_df.columns)))]
BC_scRNA_sample2_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
BC_scRNA_sample2_df.reset_index(drop=True, inplace=True)

pseudobulk_sample4_CD4TCR = BC_scRNA_sample2_df[BC_scRNA_sample2_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample4_CD4T = BC_scRNA_sample2_df[BC_scRNA_sample2_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



# ESCA    


cells_label = ESCA_df.T.index.tolist()
genes_label = ESCA_df["Unnamed: 0"].tolist()
ESCA_df_transposed = ESCA_df.T
ESCA_df_transposed = ESCA_df_transposed.drop(ESCA_df_transposed.index[0], axis=0)
ESCA_df_transposed = ESCA_df_transposed.drop(ESCA_df_transposed.columns[0], axis=1)
<<<<<<< HEAD
ESCA_df_transposed.columns = genes_label[1:]  # 去掉第一个元素，因为它是列名
ESCA_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  # 去掉第一个元素，因为它是行名
=======
ESCA_df_transposed.columns = genes_label[1:]  
ESCA_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  
>>>>>>> bb3e26c (Initial commit)
ESCA_df_transposed.reset_index(drop=True, inplace=True)



ESCA_scRNA_sample1_df = pd.merge(ESCA_sample1_df, ESCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample1_df = ESCA_scRNA_sample1_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample1_df.columns)))]
ESCA_scRNA_sample1_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample1_df.reset_index(drop=True, inplace=True)

pseudobulk_sample5_CD4TCR = ESCA_scRNA_sample1_df[ESCA_scRNA_sample1_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample5_CD4T = ESCA_scRNA_sample1_df[ESCA_scRNA_sample1_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


ESCA_scRNA_sample2_df = pd.merge(ESCA_sample2_df, ESCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample2_df = ESCA_scRNA_sample2_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample2_df.columns)))]
ESCA_scRNA_sample2_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample2_df.reset_index(drop=True, inplace=True)

pseudobulk_sample6_CD4TCR = ESCA_scRNA_sample2_df[ESCA_scRNA_sample2_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample6_CD4T = ESCA_scRNA_sample2_df[ESCA_scRNA_sample2_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



ESCA_scRNA_sample3_df = pd.merge(ESCA_sample3_df, ESCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample3_df = ESCA_scRNA_sample3_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample3_df.columns)))]
ESCA_scRNA_sample3_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample3_df.reset_index(drop=True, inplace=True)

pseudobulk_sample7_CD4TCR = ESCA_scRNA_sample3_df[ESCA_scRNA_sample3_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample7_CD4T = ESCA_scRNA_sample3_df[ESCA_scRNA_sample3_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



ESCA_scRNA_sample4_df = pd.merge(ESCA_sample4_df, ESCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample4_df = ESCA_scRNA_sample4_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample4_df.columns)))]
ESCA_scRNA_sample4_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample4_df.reset_index(drop=True, inplace=True)


pseudobulk_sample8_CD4TCR = ESCA_scRNA_sample4_df[ESCA_scRNA_sample4_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample8_CD4T = ESCA_scRNA_sample4_df[ESCA_scRNA_sample4_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



ESCA_scRNA_sample5_df = pd.merge(ESCA_sample5_df, ESCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
ESCA_scRNA_sample5_df = ESCA_scRNA_sample5_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample5_df.columns)))]
ESCA_scRNA_sample5_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample5_df.reset_index(drop=True, inplace=True)

pseudobulk_sample9_CD4TCR = ESCA_scRNA_sample5_df[ESCA_scRNA_sample5_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample9_CD4T = ESCA_scRNA_sample5_df[ESCA_scRNA_sample5_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


ESCA_scRNA_sample6_df = pd.merge(ESCA_sample6_df, ESCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
ESCA_scRNA_sample6_df = ESCA_scRNA_sample6_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample6_df.columns)))]
ESCA_scRNA_sample6_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample6_df.reset_index(drop=True, inplace=True)

pseudobulk_sample10_CD4TCR = ESCA_scRNA_sample6_df[ESCA_scRNA_sample6_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample10_CD4T = ESCA_scRNA_sample6_df[ESCA_scRNA_sample6_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


# PACA

cells_label = PACA_df.T.index.tolist()
genes_label = PACA_df["Unnamed: 0"].tolist()
PACA_df_transposed = PACA_df.T
PACA_df_transposed = PACA_df_transposed.drop(PACA_df_transposed.index[0], axis=0)
PACA_df_transposed = PACA_df_transposed.drop(PACA_df_transposed.columns[0], axis=1)
<<<<<<< HEAD
PACA_df_transposed.columns = genes_label[1:]  # 去掉第一个元素，因为它是列名
PACA_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  # 去掉第一个元素，因为它是行名
=======
PACA_df_transposed.columns = genes_label[1:]  
PACA_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  
>>>>>>> bb3e26c (Initial commit)
PACA_df_transposed.reset_index(drop=True, inplace=True)


PACA_scRNA_sample1_df = pd.merge(PACA_sample1_df, PACA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
PACA_scRNA_sample1_df = PACA_scRNA_sample1_df.iloc[:, [4,3]+list(range(5,len(PACA_scRNA_sample1_df.columns)))]
PACA_scRNA_sample1_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
PACA_scRNA_sample1_df.reset_index(drop=True, inplace=True)

pseudobulk_sample11_CD4TCR = PACA_scRNA_sample1_df[PACA_scRNA_sample1_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample11_CD4T = PACA_scRNA_sample1_df[PACA_scRNA_sample1_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


PACA_scRNA_sample2_df = pd.merge(PACA_sample2_df, PACA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
PACA_scRNA_sample2_df = PACA_scRNA_sample2_df.iloc[:, [4,3]+list(range(5,len(PACA_scRNA_sample2_df.columns)))]
PACA_scRNA_sample2_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
PACA_scRNA_sample2_df.reset_index(drop=True, inplace=True)

pseudobulk_sample12_CD4TCR = PACA_scRNA_sample2_df[PACA_scRNA_sample2_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample12_CD4T = PACA_scRNA_sample2_df[PACA_scRNA_sample2_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


# THCA

cells_label = THCA_df.T.index.tolist()
genes_label = THCA_df["Unnamed: 0"].tolist()
THCA_df_transposed = THCA_df.T
THCA_df_transposed = THCA_df_transposed.drop(THCA_df_transposed.index[0], axis=0)
THCA_df_transposed = THCA_df_transposed.drop(THCA_df_transposed.columns[0], axis=1)
<<<<<<< HEAD
THCA_df_transposed.columns = genes_label[1:]  # 去掉第一个元素，因为它是列名
THCA_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  # 去掉第一个元素，因为它是行名
=======
THCA_df_transposed.columns = genes_label[1:]  
THCA_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  
>>>>>>> bb3e26c (Initial commit)
THCA_df_transposed.reset_index(drop=True, inplace=True)


THCA_scRNA_sample1_df = pd.merge(THCA_sample1_df, THCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample1_df = THCA_scRNA_sample1_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample1_df.columns)))]
THCA_scRNA_sample1_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample1_df.reset_index(drop=True, inplace=True)

pseudobulk_sample13_CD4TCR = THCA_scRNA_sample1_df[THCA_scRNA_sample1_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample13_CD4T = THCA_scRNA_sample1_df[THCA_scRNA_sample1_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


THCA_scRNA_sample2_df = pd.merge(THCA_sample2_df, THCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample2_df = THCA_scRNA_sample2_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample2_df.columns)))]
THCA_scRNA_sample2_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample2_df.reset_index(drop=True, inplace=True)

pseudobulk_sample14_CD4TCR = THCA_scRNA_sample2_df[THCA_scRNA_sample2_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample14_CD4T = THCA_scRNA_sample2_df[THCA_scRNA_sample2_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


THCA_scRNA_sample3_df = pd.merge(THCA_sample3_df, THCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample3_df = THCA_scRNA_sample3_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample3_df.columns)))]
THCA_scRNA_sample3_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample3_df.reset_index(drop=True, inplace=True)

pseudobulk_sample15_CD4TCR = THCA_scRNA_sample3_df[THCA_scRNA_sample3_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample15_CD4T = THCA_scRNA_sample3_df[THCA_scRNA_sample3_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


THCA_scRNA_sample4_df = pd.merge(THCA_sample4_df, THCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample4_df = THCA_scRNA_sample4_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample4_df.columns)))]
THCA_scRNA_sample4_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample4_df.reset_index(drop=True, inplace=True)

pseudobulk_sample16_CD4TCR = THCA_scRNA_sample4_df[THCA_scRNA_sample4_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample16_CD4T = THCA_scRNA_sample4_df[THCA_scRNA_sample4_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


THCA_scRNA_sample5_df = pd.merge(THCA_sample5_df, THCA_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample5_df = THCA_scRNA_sample5_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample5_df.columns)))]
THCA_scRNA_sample5_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample5_df.reset_index(drop=True, inplace=True)

pseudobulk_sample17_CD4TCR = THCA_scRNA_sample5_df[THCA_scRNA_sample5_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample17_CD4T = THCA_scRNA_sample5_df[THCA_scRNA_sample5_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


# UCEC


cells_label = UCEC_df.T.index.tolist()
genes_label = UCEC_df["Unnamed: 0"].tolist()
UCEC_df_transposed = UCEC_df.T
UCEC_df_transposed = UCEC_df_transposed.drop(UCEC_df_transposed.index[0], axis=0)
UCEC_df_transposed = UCEC_df_transposed.drop(UCEC_df_transposed.columns[0], axis=1)
<<<<<<< HEAD
UCEC_df_transposed.columns = genes_label[1:]  # 去掉第一个元素，因为它是列名
UCEC_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  # 去掉第一个元素，因为它是行名
=======
UCEC_df_transposed.columns = genes_label[1:]  
UCEC_df_transposed.insert(0, 'Cell_Label', cells_label[1:])  
>>>>>>> bb3e26c (Initial commit)
UCEC_df_transposed.reset_index(drop=True, inplace=True)


UCEC_scRNA_sample1_df = pd.merge(UCEC_sample1_df, UCEC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample1_df = UCEC_scRNA_sample1_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample1_df.columns)))]
UCEC_scRNA_sample1_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample1_df.reset_index(drop=True, inplace=True)

pseudobulk_sample18_CD4TCR = UCEC_scRNA_sample1_df[UCEC_scRNA_sample1_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample18_CD4T = UCEC_scRNA_sample1_df[UCEC_scRNA_sample1_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample2_df = pd.merge(UCEC_sample2_df, UCEC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample2_df = UCEC_scRNA_sample2_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample2_df.columns)))]
UCEC_scRNA_sample2_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample2_df.reset_index(drop=True, inplace=True)

pseudobulk_sample19_CD4TCR = UCEC_scRNA_sample2_df[UCEC_scRNA_sample2_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample19_CD4T = UCEC_scRNA_sample2_df[UCEC_scRNA_sample2_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample3_df = pd.merge(UCEC_sample3_df, UCEC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample3_df = UCEC_scRNA_sample3_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample3_df.columns)))]
UCEC_scRNA_sample3_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample3_df.reset_index(drop=True, inplace=True)

pseudobulk_sample20_CD4TCR = UCEC_scRNA_sample3_df[UCEC_scRNA_sample3_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample20_CD4T = UCEC_scRNA_sample3_df[UCEC_scRNA_sample3_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample4_df = pd.merge(UCEC_sample4_df, UCEC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample4_df = UCEC_scRNA_sample4_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample4_df.columns)))]
UCEC_scRNA_sample4_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample4_df.reset_index(drop=True, inplace=True)

pseudobulk_sample21_CD4TCR = UCEC_scRNA_sample4_df[UCEC_scRNA_sample4_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample21_CD4T = UCEC_scRNA_sample4_df[UCEC_scRNA_sample4_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample5_df = pd.merge(UCEC_sample5_df, UCEC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
UCEC_scRNA_sample5_df = UCEC_scRNA_sample5_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample5_df.columns)))]
UCEC_scRNA_sample5_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample5_df.reset_index(drop=True, inplace=True)

pseudobulk_sample22_CD4TCR = UCEC_scRNA_sample5_df[UCEC_scRNA_sample5_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample22_CD4T = UCEC_scRNA_sample5_df[UCEC_scRNA_sample5_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample6_df = pd.merge(UCEC_sample6_df, UCEC_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
UCEC_scRNA_sample6_df = UCEC_scRNA_sample6_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample6_df.columns)))]
UCEC_scRNA_sample6_df.sort_values(by='CD4_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample6_df.reset_index(drop=True, inplace=True)

pseudobulk_sample23_CD4TCR = UCEC_scRNA_sample6_df[UCEC_scRNA_sample6_df['CD4_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample23_CD4T = UCEC_scRNA_sample6_df[UCEC_scRNA_sample6_df['CD4_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



# Combine the pseudobulk data across all cancer types.

intersected_genes = list(set(pseudobulk_sample1_CD4TCR.index) & set(pseudobulk_sample2_CD4TCR.index) & set(pseudobulk_sample3_CD4TCR.index) & set(pseudobulk_sample4_CD4TCR.index) & set(pseudobulk_sample5_CD4TCR.index) & set(pseudobulk_sample6_CD4TCR.index) & set(pseudobulk_sample7_CD4TCR.index) & set(pseudobulk_sample8_CD4TCR.index) & set(pseudobulk_sample9_CD4TCR.index) & set(pseudobulk_sample10_CD4TCR.index) & set(pseudobulk_sample11_CD4TCR.index) & set(pseudobulk_sample12_CD4TCR.index) & set(pseudobulk_sample13_CD4TCR.index) & set(pseudobulk_sample14_CD4TCR.index) & set(pseudobulk_sample15_CD4TCR.index) & set(pseudobulk_sample16_CD4TCR.index) & set(pseudobulk_sample17_CD4TCR.index) & set(pseudobulk_sample18_CD4TCR.index) & set(pseudobulk_sample19_CD4TCR.index) & set(pseudobulk_sample20_CD4TCR.index) & set(pseudobulk_sample21_CD4TCR.index) & set(pseudobulk_sample22_CD4TCR.index) & set(pseudobulk_sample23_CD4TCR.index))
<<<<<<< HEAD
pseudobulk_CD4TCR = pd.concat([pseudobulk_sample1_CD4TCR[intersected_genes], pseudobulk_sample2_CD4TCR[intersected_genes], pseudobulk_sample3_CD4TCR[intersected_genes], pseudobulk_sample4_CD4TCR[intersected_genes], pseudobulk_sample5_CD4TCR[intersected_genes], pseudobulk_sample6_CD4TCR[intersected_genes], pseudobulk_sample7_CD4TCR[intersected_genes], pseudobulk_sample8_CD4TCR[intersected_genes], pseudobulk_sample9_CD4TCR[intersected_genes], pseudobulk_sample10_CD4TCR[intersected_genes], pseudobulk_sample11_CD4TCR[intersected_genes], pseudobulk_sample12_CD4TCR[intersected_genes], pseudobulk_sample13_CD4TCR[intersected_genes], pseudobulk_sample14_CD4TCR[intersected_genes], pseudobulk_sample15_CD4TCR[intersected_genes], pseudobulk_sample16_CD4TCR[intersected_genes], pseudobulk_sample17_CD4TCR[intersected_genes], pseudobulk_sample18_CD4TCR[intersected_genes], pseudobulk_sample19_CD4TCR[intersected_genes], pseudobulk_sample20_CD4TCR[intersected_genes], pseudobulk_sample21_CD4TCR[intersected_genes], pseudobulk_sample22_CD4TCR[intersected_genes], pseudobulk_sample23_CD4TCR[intersected_genes]], axis=1)
=======

pseudobulk_CD4TCR = pd.concat([pseudobulk_sample1_CD4TCR[intersected_genes], pseudobulk_sample2_CD4TCR[intersected_genes], pseudobulk_sample3_CD4TCR[intersected_genes], pseudobulk_sample4_CD4TCR[intersected_genes], pseudobulk_sample5_CD4TCR[intersected_genes], pseudobulk_sample6_CD4TCR[intersected_genes], pseudobulk_sample7_CD4TCR[intersected_genes], pseudobulk_sample8_CD4TCR[intersected_genes], pseudobulk_sample9_CD4TCR[intersected_genes], pseudobulk_sample10_CD4TCR[intersected_genes], pseudobulk_sample11_CD4TCR[intersected_genes], pseudobulk_sample12_CD4TCR[intersected_genes], pseudobulk_sample13_CD4TCR[intersected_genes], pseudobulk_sample14_CD4TCR[intersected_genes], pseudobulk_sample15_CD4TCR[intersected_genes], pseudobulk_sample16_CD4TCR[intersected_genes], pseudobulk_sample17_CD4TCR[intersected_genes], pseudobulk_sample18_CD4TCR[intersected_genes], pseudobulk_sample19_CD4TCR[intersected_genes], pseudobulk_sample20_CD4TCR[intersected_genes], pseudobulk_sample21_CD4TCR[intersected_genes], pseudobulk_sample22_CD4TCR[intersected_genes], pseudobulk_sample23_CD4TCR[intersected_genes]], axis=1)

>>>>>>> bb3e26c (Initial commit)
pseudobulk_CD4T = pd.concat([pseudobulk_sample1_CD4T[intersected_genes], pseudobulk_sample2_CD4T[intersected_genes], pseudobulk_sample3_CD4T[intersected_genes], pseudobulk_sample4_CD4T[intersected_genes], pseudobulk_sample5_CD4T[intersected_genes], pseudobulk_sample6_CD4T[intersected_genes], pseudobulk_sample7_CD4T[intersected_genes], pseudobulk_sample8_CD4T[intersected_genes], pseudobulk_sample9_CD4T[intersected_genes], pseudobulk_sample10_CD4T[intersected_genes], pseudobulk_sample11_CD4T[intersected_genes], pseudobulk_sample12_CD4T[intersected_genes], pseudobulk_sample13_CD4T[intersected_genes], pseudobulk_sample14_CD4T[intersected_genes], pseudobulk_sample15_CD4T[intersected_genes], pseudobulk_sample16_CD4T[intersected_genes], pseudobulk_sample17_CD4T[intersected_genes], pseudobulk_sample18_CD4T[intersected_genes], pseudobulk_sample19_CD4T[intersected_genes], pseudobulk_sample20_CD4T[intersected_genes], pseudobulk_sample21_CD4T[intersected_genes], pseudobulk_sample22_CD4T[intersected_genes], pseudobulk_sample23_CD4T[intersected_genes]], axis=1)


# Check 75% percentile

log2_read_depth_list = []
log2_read_depth_75percentile_list = []
for i in range(23):
    log2_read_depth_list.append(np.log2(pseudobulk_CD4TCR.iloc[:,i].sum()))
    log2_read_depth_list.append(np.log2(pseudobulk_CD4T.iloc[:,i].sum()))
    log2_read_depth_75percentile_list.append(np.log2(pseudobulk_CD4TCR.iloc[:,i].cumsum().quantile(0.75)))
    log2_read_depth_75percentile_list.append(np.log2(pseudobulk_CD4T.iloc[:,i].cumsum().quantile(0.75)))
    print(i)

# Linear regression
slope, intercept, r_value, p_value, std_err = linregress(log2_read_depth_list, log2_read_depth_75percentile_list)

# Calculate Spearman correlation coefficient
rho, _ = spearmanr(log2_read_depth_list, log2_read_depth_75percentile_list)

# Create figure
plt.figure(figsize=(10, 8))

# Plot scatter plot
for i in range(46):
    if i % 2 == 0:
        plt.scatter(log2_read_depth_list[i], log2_read_depth_75percentile_list[i], color='red', label='Pseudobulk cancer reactive CD4+T sample' if i == 0 else "")
    else:
        plt.scatter(log2_read_depth_list[i], log2_read_depth_75percentile_list[i], color='blue', label='Pseudobulk normal CD4+T sample' if i == 1 else "")

# Plot linear regression line
x_vals = np.array(plt.gca().get_xlim())
y_vals = intercept + slope * x_vals
plt.plot(x_vals, y_vals, color='black', label='Linear fit')

# Set border
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Set axis label and tick
plt.xlabel('Log2 Read Depth', fontsize=10, fontweight='bold')
plt.ylabel('Log2 75th Percentile Read Depth', fontsize=10, fontweight='bold')
plt.xticks(fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')

# Add annotation
plt.annotate(f'Spearman $\\rho$: {rho:.2f}\n$R^2$: {r_value**2:.2f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, fontweight='bold', ha='left', va='top', color='blue')

# Show legend
plt.legend(loc='lower right')

# Show figure
plt.show()


# Quality control

pseudobulk_CD4T_CD4TCR = pd.concat([pseudobulk_CD4T, pseudobulk_CD4TCR], axis=1)

non_zero_counts = (pseudobulk_CD4T_CD4TCR != 0).sum(axis=1)
total_samples = pseudobulk_CD4T_CD4TCR.shape[1]
pseudobulk_CD4T_CD4TCR_qc = pseudobulk_CD4T_CD4TCR[non_zero_counts > 0.1 * total_samples]




# Create sample metadata, ensure sample names are unique
sample_names = [f"normal_{i}" for i in range(1, 24)] + [f"cancer_reactive_{i}" for i in range(1, 24)]
sample_metadata = pd.DataFrame({
    'Sample': sample_names,
    'Condition': ['normal'] * 23 + ['cancer_reactive'] * 23,
    'Pair': list(range(1, 24)) * 2,
    'Log2_75th_Quantile': log2_read_depth_75percentile_list    
})

# Set index
sample_metadata.set_index('Sample', inplace=True)

# Convert Pair to string type factor
sample_metadata['Pair'] = sample_metadata['Pair'].astype(str).astype('category')

# Center and standardize Log2_75th_Quantile
sample_metadata['Log2_75th_Quantile'] = (sample_metadata['Log2_75th_Quantile'] - sample_metadata['Log2_75th_Quantile'].mean()) / sample_metadata['Log2_75th_Quantile'].std()

# Ensure gene_expression columns are in the same order as sample_metadata rows
pseudobulk_CD4T_CD4TCR_qc.columns = sample_metadata.index

<<<<<<< HEAD
=======

pseudobulk_CD4T_CD4TCR_qc.loc["GSTM3"][:23].mean()
pseudobulk_CD4T_CD4TCR_qc.loc["GSTM3"][23:].mean()


>>>>>>> bb3e26c (Initial commit)
# Activate pandas and R data frame automatic conversion
pandas2ri.activate()

# Import DESeq2
deseq2 = importr('DESeq2')

# Convert data to R data frame
r_count_data = pandas2ri.py2rpy(pseudobulk_CD4T_CD4TCR_qc)
r_col_data = pandas2ri.py2rpy(sample_metadata)

# Create DESeq2 data set, add pair information and covariate
dds_method = deseq2.DESeqDataSetFromMatrix(countData=r_count_data,
                                    colData=r_col_data,
                                    design=ro.Formula('~ Pair + Condition + offset(Log2_75th_Quantile)'))




# Run DESeq2
dds_method = deseq2.DESeq(dds_method)

# Obtain results
res_method = deseq2.results(dds_method)

# Convert RS4 object to R DataFrame
r_dataframe_method = ro.r['as.data.frame'](res_method)

# Convert R DataFrame to Pandas DataFrame
res_df_method = pandas2ri.rpy2py(r_dataframe_method)

res_df_method.dropna(inplace=True)
res_df_method.sort_values(by='pvalue', ascending=True, inplace=True)

# Bonferroni correction
method_CD4T_CD4TCR_DEgenes = res_df_method[res_df_method["pvalue"] < 0.05/12691]

method_CD4T_CD4TCR_DEgenes["Gene_Name"] = method_CD4T_CD4TCR_DEgenes.index
method_CD4T_CD4TCR_DEgenes.reset_index(drop=True, inplace=True)


method_CD4TCR_sigs = pd.DataFrame({"Gene_Name": method_CD4T_CD4TCR_DEgenes["Gene_Name"], "effect": ["positive"]*len(method_CD4T_CD4TCR_DEgenes)})
method_CD4TCR_sigs["effect"][method_CD4T_CD4TCR_DEgenes["log2FoldChange"] > 0] = "negative"

method_CD4TCR_sigs["effect"][method_CD4T_CD4TCR_DEgenes["log2FoldChange"] > 0] = "negative"


method_CD4TCR_sigs["effect"].value_counts()



method_CD4TCR_sigs.to_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD4TCR_Data/method_CD4TCR_sigs.csv", index=None)

