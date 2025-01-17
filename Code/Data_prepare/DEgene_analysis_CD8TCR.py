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
Date last modified: 2025-01-16 by Saishi Cui
Program: DEgene_analysis_CD8TCR.py
Purpose: Using external single-cell RNA-seq data of CD8+T and CD8+T cancer reactive samples,
analyze the DE genes between CD8+T cancer reactive and non-cancer reactive CD8+T samples.
--------------------------------------------------------------------------------
Data Inputs:

- Single-cell RNA-seq data of CD8+T and CD8+T cancer reactive samples, and their label information.

Data Outputs:

- DE genes between CD8+T cancer reactive and non-cancer reactive CD8+T samples.
--------------------------------------------------------------------------------
"""


# Read the label information of CD8+T cancer reactive and non-cancer reactive samples

CD8T_cell_status = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/CD8_cancer_reactive_status.tsv", sep="\t")


# We only focus on samples with more than 50 normal CD8+T and 50 CD8+T cancer reactive samples.

L_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'L-P20181123-T']
L_sample2_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'L-P20190404-T']
BC_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'BC-P20190403-T']
ESCA_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'ESCA-P20190411-T']
ESCA_sample2_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'ESCA-P20190613-T']
ESCA_sample3_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'ESCA-P20181114-T']
ESCA_sample4_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'ESCA-P20181123-T']
ESCA_sample5_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'ESCA-P20190410-T']
ESCA_sample6_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'ESCA-P20190404-T']
MM_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'MM-P20190122-T']
PACA_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'PACA-P20181128-T']
PACA_sample2_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'PACA-P20190306-T']
PACA_sample3_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'PACA-P20181121-T']
RC_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'RC-P20190925-T']
RC_sample2_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'RC-P20190923-T']
THCA_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'THCA-P20190730-T']
THCA_sample2_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'THCA-P20190816-T']
THCA_sample3_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'THCA-P20190703-T']
UCEC_sample1_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'UCEC-P20190910-T']
UCEC_sample2_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'UCEC-P20190312-T']
UCEC_sample3_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'UCEC-P20190305-T']
UCEC_sample4_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'UCEC-P20190213-T']
UCEC_sample5_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'UCEC-P20190911-T']
UCEC_sample6_df = CD8T_cell_status[CD8T_cell_status['library.id'] == 'UCEC-P20190625-T']




# Read the single-cell RNA-seq data for each cancer type



BCL_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/BCL.txt", sep="\t")
BC_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/BC.txt", sep="\t")
ESCA_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/ESCA.txt", sep="\t")
MM_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/MM.txt", sep="\t")
PACA_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/PACA.txt", sep="\t")
RC_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/RC.txt", sep="\t")
THCA_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/THCA.txt", sep="\t")
UCEC_CD8T_df = pd.read_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/UCEC.txt", sep="\t")






# Create the pseudobulk data for each sample.
# BCL
cells_CD8T_label = BCL_CD8T_df.T.index.tolist()
genes_CD8T_label = BCL_CD8T_df["Unnamed: 0"].tolist()
BCL_CD8T_df_transposed = BCL_CD8T_df.T
BCL_CD8T_df_transposed = BCL_CD8T_df_transposed.drop(BCL_CD8T_df_transposed.index[0], axis=0)
BCL_CD8T_df_transposed = BCL_CD8T_df_transposed.drop(BCL_CD8T_df_transposed.columns[0], axis=1)
BCL_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
BCL_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
BCL_CD8T_df_transposed.reset_index(drop=True, inplace=True)

BCL_scRNA_sample1_CD8T_df = pd.merge(L_sample1_df, BCL_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
BCL_scRNA_sample1_CD8T_df = BCL_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(BCL_scRNA_sample1_CD8T_df.columns)))]
BCL_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
BCL_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample1_CD8TCR = BCL_scRNA_sample1_CD8T_df[BCL_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample1_CD8T = BCL_scRNA_sample1_CD8T_df[BCL_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


BCL_scRNA_sample2_CD8T_df = pd.merge(L_sample2_df, BCL_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
BCL_scRNA_sample2_CD8T_df = BCL_scRNA_sample2_CD8T_df.iloc[:, [4,3]+list(range(5,len(BCL_scRNA_sample2_CD8T_df.columns)))]
BCL_scRNA_sample2_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
BCL_scRNA_sample2_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample2_CD8TCR = BCL_scRNA_sample2_CD8T_df[BCL_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample2_CD8T = BCL_scRNA_sample2_CD8T_df[BCL_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)





# BC


cells_CD8T_label = BC_CD8T_df.T.index.tolist()
genes_CD8T_label = BC_CD8T_df["Unnamed: 0"].tolist()
BC_CD8T_df_transposed = BC_CD8T_df.T
BC_CD8T_df_transposed = BC_CD8T_df_transposed.drop(BC_CD8T_df_transposed.index[0], axis=0)
BC_CD8T_df_transposed = BC_CD8T_df_transposed.drop(BC_CD8T_df_transposed.columns[0], axis=1)
BC_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
BC_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
BC_CD8T_df_transposed.reset_index(drop=True, inplace=True)


BC_scRNA_sample1_CD8T_df = pd.merge(BC_sample1_df, BC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
BC_scRNA_sample1_CD8T_df = BC_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(BC_scRNA_sample1_CD8T_df.columns)))]
BC_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
BC_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample3_CD8TCR = BC_scRNA_sample1_CD8T_df[BC_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample3_CD8T = BC_scRNA_sample1_CD8T_df[BC_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



# ESCA    


cells_CD8T_label = ESCA_CD8T_df.T.index.tolist()
genes_CD8T_label = ESCA_CD8T_df["Unnamed: 0"].tolist()
ESCA_CD8T_df_transposed = ESCA_CD8T_df.T
ESCA_CD8T_df_transposed = ESCA_CD8T_df_transposed.drop(ESCA_CD8T_df_transposed.index[0], axis=0)
ESCA_CD8T_df_transposed = ESCA_CD8T_df_transposed.drop(ESCA_CD8T_df_transposed.columns[0], axis=1)
ESCA_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
ESCA_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
ESCA_CD8T_df_transposed.reset_index(drop=True, inplace=True)


ESCA_scRNA_sample1_CD8T_df = pd.merge(ESCA_sample1_df, ESCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample1_CD8T_df = ESCA_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample1_CD8T_df.columns)))]
ESCA_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample4_CD8TCR = ESCA_scRNA_sample1_CD8T_df[ESCA_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample4_CD8T = ESCA_scRNA_sample1_CD8T_df[ESCA_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)

ESCA_scRNA_sample2_CD8T_df = pd.merge(ESCA_sample2_df, ESCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample2_CD8T_df = ESCA_scRNA_sample2_CD8T_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample2_CD8T_df.columns)))]
ESCA_scRNA_sample2_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample2_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample5_CD8TCR = ESCA_scRNA_sample2_CD8T_df[ESCA_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample5_CD8T = ESCA_scRNA_sample2_CD8T_df[ESCA_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



ESCA_scRNA_sample3_CD8T_df = pd.merge(ESCA_sample3_df, ESCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample3_CD8T_df = ESCA_scRNA_sample3_CD8T_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample3_CD8T_df.columns)))]
ESCA_scRNA_sample3_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample3_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample6_CD8TCR = ESCA_scRNA_sample3_CD8T_df[ESCA_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample6_CD8T = ESCA_scRNA_sample3_CD8T_df[ESCA_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


ESCA_scRNA_sample4_CD8T_df = pd.merge(ESCA_sample4_df, ESCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
ESCA_scRNA_sample4_CD8T_df = ESCA_scRNA_sample4_CD8T_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample4_CD8T_df.columns)))]
ESCA_scRNA_sample4_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample4_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample7_CD8TCR = ESCA_scRNA_sample4_CD8T_df[ESCA_scRNA_sample4_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample7_CD8T = ESCA_scRNA_sample4_CD8T_df[ESCA_scRNA_sample4_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


ESCA_scRNA_sample5_CD8T_df = pd.merge(ESCA_sample5_df, ESCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
ESCA_scRNA_sample5_CD8T_df = ESCA_scRNA_sample5_CD8T_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample5_CD8T_df.columns)))]
ESCA_scRNA_sample5_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample5_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample8_CD8TCR = ESCA_scRNA_sample5_CD8T_df[ESCA_scRNA_sample5_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample8_CD8T = ESCA_scRNA_sample5_CD8T_df[ESCA_scRNA_sample5_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



ESCA_scRNA_sample6_CD8T_df = pd.merge(ESCA_sample6_df, ESCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
ESCA_scRNA_sample6_CD8T_df = ESCA_scRNA_sample6_CD8T_df.iloc[:, [4,3]+list(range(5,len(ESCA_scRNA_sample6_CD8T_df.columns)))]
ESCA_scRNA_sample6_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
ESCA_scRNA_sample6_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample9_CD8TCR = ESCA_scRNA_sample6_CD8T_df[ESCA_scRNA_sample6_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample9_CD8T = ESCA_scRNA_sample6_CD8T_df[ESCA_scRNA_sample6_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



# MM 


cells_CD8T_label = MM_CD8T_df.T.index.tolist()
genes_CD8T_label = MM_CD8T_df["Unnamed: 0"].tolist()
MM_CD8T_df_transposed = MM_CD8T_df.T
MM_CD8T_df_transposed = MM_CD8T_df_transposed.drop(MM_CD8T_df_transposed.index[0], axis=0)
MM_CD8T_df_transposed = MM_CD8T_df_transposed.drop(MM_CD8T_df_transposed.columns[0], axis=1)
MM_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
MM_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
MM_CD8T_df_transposed.reset_index(drop=True, inplace=True)

MM_scRNA_sample1_CD8T_df = pd.merge(MM_sample1_df, MM_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
MM_scRNA_sample1_CD8T_df = MM_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(MM_scRNA_sample1_CD8T_df.columns)))]
MM_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
MM_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample10_CD8TCR = MM_scRNA_sample1_CD8T_df[MM_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample10_CD8T = MM_scRNA_sample1_CD8T_df[MM_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



# PACA


cells_CD8T_label = PACA_CD8T_df.T.index.tolist()
genes_CD8T_label = PACA_CD8T_df["Unnamed: 0"].tolist()
PACA_CD8T_df_transposed = PACA_CD8T_df.T
PACA_CD8T_df_transposed = PACA_CD8T_df_transposed.drop(PACA_CD8T_df_transposed.index[0], axis=0)
PACA_CD8T_df_transposed = PACA_CD8T_df_transposed.drop(PACA_CD8T_df_transposed.columns[0], axis=1)
PACA_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
PACA_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
PACA_CD8T_df_transposed.reset_index(drop=True, inplace=True)


PACA_scRNA_sample1_CD8T_df = pd.merge(PACA_sample1_df, PACA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
PACA_scRNA_sample1_CD8T_df = PACA_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(PACA_scRNA_sample1_CD8T_df.columns)))]
PACA_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
PACA_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample11_CD8TCR = PACA_scRNA_sample1_CD8T_df[PACA_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample11_CD8T = PACA_scRNA_sample1_CD8T_df[PACA_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



PACA_scRNA_sample2_CD8T_df = pd.merge(PACA_sample2_df, PACA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
PACA_scRNA_sample2_CD8T_df = PACA_scRNA_sample2_CD8T_df.iloc[:, [4,3]+list(range(5,len(PACA_scRNA_sample2_CD8T_df.columns)))]
PACA_scRNA_sample2_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
PACA_scRNA_sample2_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample12_CD8TCR = PACA_scRNA_sample2_CD8T_df[PACA_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample12_CD8T = PACA_scRNA_sample2_CD8T_df[PACA_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



PACA_scRNA_sample3_CD8T_df = pd.merge(PACA_sample3_df, PACA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
PACA_scRNA_sample3_CD8T_df = PACA_scRNA_sample3_CD8T_df.iloc[:, [4,3]+list(range(5,len(PACA_scRNA_sample3_CD8T_df.columns)))]
PACA_scRNA_sample3_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
PACA_scRNA_sample3_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample13_CD8TCR = PACA_scRNA_sample3_CD8T_df[PACA_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample13_CD8T = PACA_scRNA_sample3_CD8T_df[PACA_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


# RC


cells_CD8T_label = RC_CD8T_df.T.index.tolist()
genes_CD8T_label = RC_CD8T_df["Unnamed: 0"].tolist()
RC_CD8T_df_transposed = RC_CD8T_df.T
RC_CD8T_df_transposed = RC_CD8T_df_transposed.drop(RC_CD8T_df_transposed.index[0], axis=0)
RC_CD8T_df_transposed = RC_CD8T_df_transposed.drop(RC_CD8T_df_transposed.columns[0], axis=1)
RC_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
RC_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
RC_CD8T_df_transposed.reset_index(drop=True, inplace=True)


RC_scRNA_sample1_CD8T_df = pd.merge(RC_sample1_df, RC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
RC_scRNA_sample1_CD8T_df = RC_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(RC_scRNA_sample1_CD8T_df.columns)))]
RC_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
RC_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample14_CD8TCR = RC_scRNA_sample1_CD8T_df[RC_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample14_CD8T = RC_scRNA_sample1_CD8T_df[RC_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



RC_scRNA_sample2_CD8T_df = pd.merge(RC_sample2_df, RC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
RC_scRNA_sample2_CD8T_df = RC_scRNA_sample2_CD8T_df.iloc[:, [4,3]+list(range(5,len(RC_scRNA_sample2_CD8T_df.columns)))]
RC_scRNA_sample2_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
RC_scRNA_sample2_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample15_CD8TCR = RC_scRNA_sample2_CD8T_df[RC_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample15_CD8T = RC_scRNA_sample2_CD8T_df[RC_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)



# THCA


cells_CD8T_label = THCA_CD8T_df.T.index.tolist()
genes_CD8T_label = THCA_CD8T_df["Unnamed: 0"].tolist()
THCA_CD8T_df_transposed = THCA_CD8T_df.T
THCA_CD8T_df_transposed = THCA_CD8T_df_transposed.drop(THCA_CD8T_df_transposed.index[0], axis=0)
THCA_CD8T_df_transposed = THCA_CD8T_df_transposed.drop(THCA_CD8T_df_transposed.columns[0], axis=1)
THCA_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
THCA_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
THCA_CD8T_df_transposed.reset_index(drop=True, inplace=True)


THCA_scRNA_sample1_CD8T_df = pd.merge(THCA_sample1_df, THCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample1_CD8T_df = THCA_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample1_CD8T_df.columns)))]
THCA_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample16_CD8TCR = THCA_scRNA_sample1_CD8T_df[THCA_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample16_CD8T = THCA_scRNA_sample1_CD8T_df[THCA_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


THCA_scRNA_sample2_CD8T_df = pd.merge(THCA_sample2_df, THCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample2_CD8T_df = THCA_scRNA_sample2_CD8T_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample2_CD8T_df.columns)))]
THCA_scRNA_sample2_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample2_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample17_CD8TCR = THCA_scRNA_sample2_CD8T_df[THCA_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample17_CD8T = THCA_scRNA_sample2_CD8T_df[THCA_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


THCA_scRNA_sample3_CD8T_df = pd.merge(THCA_sample3_df, THCA_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
THCA_scRNA_sample3_CD8T_df = THCA_scRNA_sample3_CD8T_df.iloc[:, [4,3]+list(range(5,len(THCA_scRNA_sample3_CD8T_df.columns)))]
THCA_scRNA_sample3_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
THCA_scRNA_sample3_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample18_CD8TCR = THCA_scRNA_sample3_CD8T_df[THCA_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample18_CD8T = THCA_scRNA_sample3_CD8T_df[THCA_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


# UCEC


cells_CD8T_label = UCEC_CD8T_df.T.index.tolist()
genes_CD8T_label = UCEC_CD8T_df["Unnamed: 0"].tolist()
UCEC_CD8T_df_transposed = UCEC_CD8T_df.T
UCEC_CD8T_df_transposed = UCEC_CD8T_df_transposed.drop(UCEC_CD8T_df_transposed.index[0], axis=0)
UCEC_CD8T_df_transposed = UCEC_CD8T_df_transposed.drop(UCEC_CD8T_df_transposed.columns[0], axis=1)
UCEC_CD8T_df_transposed.columns = genes_CD8T_label[1:]  
UCEC_CD8T_df_transposed.insert(0, 'Cell_Label', cells_CD8T_label[1:])  
UCEC_CD8T_df_transposed.reset_index(drop=True, inplace=True)


UCEC_scRNA_sample1_CD8T_df = pd.merge(UCEC_sample1_df, UCEC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample1_CD8T_df = UCEC_scRNA_sample1_CD8T_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample1_CD8T_df.columns)))]
UCEC_scRNA_sample1_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample1_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample19_CD8TCR = UCEC_scRNA_sample1_CD8T_df[UCEC_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample19_CD8T = UCEC_scRNA_sample1_CD8T_df[UCEC_scRNA_sample1_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample2_CD8T_df = pd.merge(UCEC_sample2_df, UCEC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample2_CD8T_df = UCEC_scRNA_sample2_CD8T_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample2_CD8T_df.columns)))]
UCEC_scRNA_sample2_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample2_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample20_CD8TCR = UCEC_scRNA_sample2_CD8T_df[UCEC_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample20_CD8T = UCEC_scRNA_sample2_CD8T_df[UCEC_scRNA_sample2_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample3_CD8T_df = pd.merge(UCEC_sample3_df, UCEC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample3_CD8T_df = UCEC_scRNA_sample3_CD8T_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample3_CD8T_df.columns)))]
UCEC_scRNA_sample3_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample3_CD8T_df.reset_index(drop=True, inplace=True)


pseudobulk_sample21_CD8TCR = UCEC_scRNA_sample3_CD8T_df[UCEC_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample21_CD8T = UCEC_scRNA_sample3_CD8T_df[UCEC_scRNA_sample3_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample4_CD8T_df = pd.merge(UCEC_sample4_df, UCEC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left")
UCEC_scRNA_sample4_CD8T_df = UCEC_scRNA_sample4_CD8T_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample4_CD8T_df.columns)))]
UCEC_scRNA_sample4_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample4_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample22_CD8TCR = UCEC_scRNA_sample4_CD8T_df[UCEC_scRNA_sample4_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample22_CD8T = UCEC_scRNA_sample4_CD8T_df[UCEC_scRNA_sample4_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample5_CD8T_df = pd.merge(UCEC_sample5_df, UCEC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
UCEC_scRNA_sample5_CD8T_df = UCEC_scRNA_sample5_CD8T_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample5_CD8T_df.columns)))]
UCEC_scRNA_sample5_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample5_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample23_CD8TCR = UCEC_scRNA_sample5_CD8T_df[UCEC_scRNA_sample5_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample23_CD8T = UCEC_scRNA_sample5_CD8T_df[UCEC_scRNA_sample5_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)


UCEC_scRNA_sample6_CD8T_df = pd.merge(UCEC_sample6_df, UCEC_CD8T_df_transposed, left_on="barcode", right_on="Cell_Label", how="left") 
UCEC_scRNA_sample6_CD8T_df = UCEC_scRNA_sample6_CD8T_df.iloc[:, [4,3]+list(range(5,len(UCEC_scRNA_sample6_CD8T_df.columns)))]
UCEC_scRNA_sample6_CD8T_df.sort_values(by='CD8_cancer_reactive', ascending=False, inplace=True)    
UCEC_scRNA_sample6_CD8T_df.reset_index(drop=True, inplace=True)

pseudobulk_sample24_CD8TCR = UCEC_scRNA_sample6_CD8T_df[UCEC_scRNA_sample6_CD8T_df['CD8_cancer_reactive'] == "yes"].iloc[:, 2:].sum(axis=0)
pseudobulk_sample24_CD8T = UCEC_scRNA_sample6_CD8T_df[UCEC_scRNA_sample6_CD8T_df['CD8_cancer_reactive'] == "no"].iloc[:, 2:].sum(axis=0)






# Combine the pseudobulk data across all samples.

intersected_genes = list(set(pseudobulk_sample1_CD8TCR.index) & set(pseudobulk_sample2_CD8TCR.index) & set(pseudobulk_sample3_CD8TCR.index) & set(pseudobulk_sample4_CD8TCR.index) & set(pseudobulk_sample5_CD8TCR.index) & set(pseudobulk_sample6_CD8TCR.index) & set(pseudobulk_sample7_CD8TCR.index) & set(pseudobulk_sample8_CD8TCR.index) & set(pseudobulk_sample9_CD8TCR.index) & set(pseudobulk_sample10_CD8TCR.index) & set(pseudobulk_sample11_CD8TCR.index) & set(pseudobulk_sample12_CD8TCR.index) & set(pseudobulk_sample13_CD8TCR.index) & set(pseudobulk_sample14_CD8TCR.index) & set(pseudobulk_sample15_CD8TCR.index) & set(pseudobulk_sample16_CD8TCR.index) & set(pseudobulk_sample17_CD8TCR.index) & set(pseudobulk_sample18_CD8TCR.index) & set(pseudobulk_sample19_CD8TCR.index) & set(pseudobulk_sample20_CD8TCR.index) & set(pseudobulk_sample21_CD8TCR.index) & set(pseudobulk_sample22_CD8TCR.index) & set(pseudobulk_sample23_CD8TCR.index) & set(pseudobulk_sample24_CD8TCR.index))
pseudobulk_CD8TCR = pd.concat([pseudobulk_sample1_CD8TCR[intersected_genes], pseudobulk_sample2_CD8TCR[intersected_genes], pseudobulk_sample3_CD8TCR[intersected_genes], pseudobulk_sample4_CD8TCR[intersected_genes], pseudobulk_sample5_CD8TCR[intersected_genes], pseudobulk_sample6_CD8TCR[intersected_genes], pseudobulk_sample7_CD8TCR[intersected_genes], pseudobulk_sample8_CD8TCR[intersected_genes], pseudobulk_sample9_CD8TCR[intersected_genes], pseudobulk_sample10_CD8TCR[intersected_genes], pseudobulk_sample11_CD8TCR[intersected_genes], pseudobulk_sample12_CD8TCR[intersected_genes], pseudobulk_sample13_CD8TCR[intersected_genes], pseudobulk_sample14_CD8TCR[intersected_genes], pseudobulk_sample15_CD8TCR[intersected_genes], pseudobulk_sample16_CD8TCR[intersected_genes], pseudobulk_sample17_CD8TCR[intersected_genes], pseudobulk_sample18_CD8TCR[intersected_genes], pseudobulk_sample19_CD8TCR[intersected_genes], pseudobulk_sample20_CD8TCR[intersected_genes], pseudobulk_sample21_CD8TCR[intersected_genes], pseudobulk_sample22_CD8TCR[intersected_genes], pseudobulk_sample23_CD8TCR[intersected_genes],pseudobulk_sample24_CD8TCR[intersected_genes]], axis=1)
pseudobulk_CD8T = pd.concat([pseudobulk_sample1_CD8T[intersected_genes], pseudobulk_sample2_CD8T[intersected_genes], pseudobulk_sample3_CD8T[intersected_genes], pseudobulk_sample4_CD8T[intersected_genes], pseudobulk_sample5_CD8T[intersected_genes], pseudobulk_sample6_CD8T[intersected_genes], pseudobulk_sample7_CD8T[intersected_genes], pseudobulk_sample8_CD8T[intersected_genes], pseudobulk_sample9_CD8T[intersected_genes], pseudobulk_sample10_CD8T[intersected_genes], pseudobulk_sample11_CD8T[intersected_genes], pseudobulk_sample12_CD8T[intersected_genes], pseudobulk_sample13_CD8T[intersected_genes], pseudobulk_sample14_CD8T[intersected_genes], pseudobulk_sample15_CD8T[intersected_genes], pseudobulk_sample16_CD8T[intersected_genes], pseudobulk_sample17_CD8T[intersected_genes], pseudobulk_sample18_CD8T[intersected_genes], pseudobulk_sample19_CD8T[intersected_genes], pseudobulk_sample20_CD8T[intersected_genes], pseudobulk_sample21_CD8T[intersected_genes], pseudobulk_sample22_CD8T[intersected_genes], pseudobulk_sample23_CD8T[intersected_genes], pseudobulk_sample24_CD8T[intersected_genes]], axis=1)



##### Check 75% percentile

log2_read_depth_list = []
log2_read_depth_75percentile_list = []
for i in range(24):
    log2_read_depth_list.append(np.log2(pseudobulk_CD8TCR.iloc[:,i].sum()))
    log2_read_depth_list.append(np.log2(pseudobulk_CD8T.iloc[:,i].sum()))
    log2_read_depth_75percentile_list.append(np.log2(pseudobulk_CD8TCR.iloc[:,i].cumsum().quantile(0.75)))
    log2_read_depth_75percentile_list.append(np.log2(pseudobulk_CD8T.iloc[:,i].cumsum().quantile(0.75)))
    print(i)

#  Linear regression
slope, intercept, r_value, p_value, std_err = linregress(log2_read_depth_list, log2_read_depth_75percentile_list)

#  Calculate Spearman correlation coefficient
rho, _ = spearmanr(log2_read_depth_list, log2_read_depth_75percentile_list)

#  Create a graph
plt.figure(figsize=(10, 8))

#  Plot scatter plot
for i in range(48):
    if i % 2 == 0:
        plt.scatter(log2_read_depth_list[i], log2_read_depth_75percentile_list[i], color='red', label='Pseudobulk cancer reactive CD4+T sample' if i == 0 else "")
    else:
        plt.scatter(log2_read_depth_list[i], log2_read_depth_75percentile_list[i], color='blue', label='Pseudobulk normal CD4+T sample' if i == 1 else "")

#  Plot linear regression line
x_vals = np.array(plt.gca().get_xlim())
y_vals = intercept + slope * x_vals
plt.plot(x_vals, y_vals, color='black', label='Linear fit')

#  Set the border
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Set the axis label and tick
plt.xlabel('Log2 Read Depth', fontsize=10, fontweight='bold')
plt.ylabel('Log2 75th Percentile Read Depth', fontsize=10, fontweight='bold')
plt.xticks(fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')

#  Add annotation
plt.annotate(f'Spearman $\\rho$: {rho:.2f}\n$R^2$: {r_value**2:.2f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, fontweight='bold', ha='left', va='top', color='blue')

#  Display legend
plt.legend(loc='lower right')

#  Display graph
plt.show()


####### Quality control

pseudobulk_CD8T_CD8TCR = pd.concat([pseudobulk_CD8T, pseudobulk_CD8TCR], axis=1)
non_zero_counts = (pseudobulk_CD8T_CD8TCR != 0).sum(axis=1)
total_samples = pseudobulk_CD8T_CD8TCR.shape[1]
pseudobulk_CD8T_CD8TCR_qc = pseudobulk_CD8T_CD8TCR[non_zero_counts > 0.1 * total_samples]




#  Create sample metadata, ensure sample name is unique
sample_names_CD8T = [f"normal_{i}" for i in range(1, 25)] + [f"cancer_reactive_{i}" for i in range(1, 25)]
sample_metadata_CD8T = pd.DataFrame({
    'Sample': sample_names_CD8T,
    'Condition': ['normal'] * 24 + ['cancer_reactive'] * 24,
    'Pair': list(range(1, 25)) * 2,
    'Log2_75th_Quantile': log2_read_depth_75percentile_list    
})

#  Set the index
sample_metadata_CD8T.set_index('Sample', inplace=True)

#  Convert Pair to string type factor
sample_metadata_CD8T['Pair'] = sample_metadata_CD8T['Pair'].astype(str).astype('category')

#  Center and standardize Log2_75th_Quantile
sample_metadata_CD8T['Log2_75th_Quantile'] = (sample_metadata_CD8T['Log2_75th_Quantile'] - sample_metadata_CD8T['Log2_75th_Quantile'].mean()) / sample_metadata_CD8T['Log2_75th_Quantile'].std()

#  Ensure gene_expression column order is consistent with sample_metadata row order
pseudobulk_CD8T_CD8TCR_qc.columns = sample_metadata_CD8T.index

#  Activate pandas and R data frame automatic conversion
pandas2ri.activate()

#  Import DESeq2
deseq2 = importr('DESeq2')

#  Convert data to R data frame
r_count_data_CD8T = pandas2ri.py2rpy(pseudobulk_CD8T_CD8TCR_qc)
r_col_data_CD8T = pandas2ri.py2rpy(sample_metadata_CD8T)

#  Create DESeq2 data set, add pair information and covariate
dds_method_CD8T = deseq2.DESeqDataSetFromMatrix(countData=r_count_data_CD8T,
                                    colData=r_col_data_CD8T,
                                    design=ro.Formula('~ Pair + Condition + offset(Log2_75th_Quantile)'))




#  Run DESeq2
dds_method_CD8T = deseq2.DESeq(dds_method_CD8T)

#  Get results
res_method_CD8T = deseq2.results(dds_method_CD8T)

#  Convert RS4 object to R data frame
r_dataframe_method_CD8T = ro.r['as.data.frame'](res_method_CD8T)

#  Convert R data frame to Pandas DataFrame
res_df_method_CD8T = pandas2ri.rpy2py(r_dataframe_method_CD8T)

res_df_method_CD8T.dropna(inplace=True)
res_df_method_CD8T.sort_values(by='pvalue', ascending=True, inplace=True)

method_CD8T_CD8TCR_DEgenes = res_df_method_CD8T[res_df_method_CD8T["pvalue"] < 0.05/13403]



method_CD8T_CD8TCR_DEgenes["Gene_Name"] = method_CD8T_CD8TCR_DEgenes.index
method_CD8T_CD8TCR_DEgenes.reset_index(drop=True, inplace=True)


method_CD8TCR_sigs = pd.DataFrame({"Gene_Name": method_CD8T_CD8TCR_DEgenes["Gene_Name"], "effect": ["positive"]*len(method_CD8T_CD8TCR_DEgenes)})
method_CD8TCR_sigs["effect"][method_CD8T_CD8TCR_DEgenes["log2FoldChange"] > 0] = "negative"




method_CD8TCR_sigs.to_csv("/Users/scui2/DNAmethylation/External_scRNAseq_CD8TCR_Data/method_CD8TCR_sigs.csv", index=None)






