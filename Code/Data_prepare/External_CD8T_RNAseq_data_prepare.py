import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
import requests
import xml.etree.ElementTree as ET
from tqdm import tqdm


"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: External_CD8T_RNAseq_data_prepare.py
Purpose: Prepare external CD8+T healthy samples data for RNA-seq data, and split into training and testing data.
--------------------------------------------------------------------------------
Data Inputs:

- Three studies' CD8+T healthy samples data, including samples matched RNA-seq data.
- Study 1: Rodriguez et al., 2017 ( N = 6)
- Study 2: Ventham et al., 2016 ( N = 12)
- Study 3: Mamrut et al., 2015 ( N = 5)

Data Outputs:
- CD8+T healthy samples data for RNA-seq data, and split into training and testing data.
--------------------------------------------------------------------------------
Functions:
-  chunk_list: Split the Ensembl_ID into chunks, to avoid network issue.
-  query_ensembl_with_retry: Query gene information from Ensembl.
--------------------------------------------------------------------------------
Notes:
- There are two versions for saving the data:
    - CD8T_RNA_study2.csv: only include the 12 samples from Ventham et al., 2016.
    - CD8T_RNA_allsamples.csv: include all the 23 samples from three studies.
<<<<<<< HEAD
=======
- Training data is using the 12 samples from Ventham et al., 2016.
- Testing data is using the 11 samples from the other two studies.
>>>>>>> bb3e26c (Initial commit)
-------------------------------------------------------------------------------
"""


### Study 1 (N = 6) （log2 transformed and quantile normalized）

GSE83156_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE83156_RNA.txt', sep='\t', header=0)

GSE83156_RNAseq_df.iloc[:,1:].sum(axis=1)
GSE83156_RNAseq_df.iloc[:,1:].sum(axis=0)

#### Study 2 (N = 12) （normalized signal）

## Age= 41.5, gender = M
GSM2336845_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336845.txt', sep='\t', header=0)

## Age = 30.8, gender = M
GSM2336851_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336851.txt', sep='\t', header=0)

## Age= 33.7, gender = M
GSM2336874_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336874.txt', sep='\t', header=0)

## Age = 32.6, gender = M 
GSM2336911_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336911.txt', sep='\t', header=0)

## Age = 43.9, gender = F
GSM2336922_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336922.txt', sep='\t', header=0)

## Age = 31.1, gender = M
GSM2336927_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336927.txt', sep='\t', header=0)

## Age = 37.4, gender = M
GSM2336935_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336935.txt', sep='\t', header=0)

## Age = 20.5, gender = M
GSM2336941_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336941.txt', sep='\t', header=0)

## Age = 30.7, gender = M
GSM2336952_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336952.txt', sep='\t', header=0)

## Age = 42.8, gender = M
GSM2336953_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2336953.txt', sep='\t', header=0)


## Age = 30.0, gender = F
GSM2337002_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2337002.txt', sep='\t', header=0)

## Age = 58.5, gender = F
GSM2337046_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE87640/GSM2337046.txt', sep='\t', header=0)


### Study 3 (N = 5) （normalized signal）

GSM1831330_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE71245/GSM1831330.txt', sep='\t', header=0)
GSM1831331_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE71245/GSM1831331.txt', sep='\t', header=0)
GSM1831332_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE71245/GSM1831332.txt', sep='\t', header=0)
GSM1831333_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE71245/GSM1831333.txt', sep='\t', header=0)
GSM1831334_RNAseq_df = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/GSE71245/GSM1831334.txt', sep='\t', header=0)


### Read the Illumina RNA-seq probe information
Illumina_RNA = pd.read_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/Illumina_RNA.txt', sep='\t', header=0)



intersect_ID_3studies = set(Illumina_RNA["ID"]).intersection(set(GSE83156_RNAseq_df["ID_REF"])).intersection(set(GSM2336845_RNAseq_df["ID_REF"])).intersection(set(GSM1831330_RNAseq_df["ID_REF"]))
intersect_df_3studies = pd.DataFrame(list(intersect_ID_3studies), columns=["ID"])
Illumina_RNA_intersect_3studies = pd.merge(intersect_df_3studies, Illumina_RNA, left_on="ID", right_on="ID", how="left")
unique_gene_3studies = Illumina_RNA_intersect_3studies["ILMN_Gene"].unique()


### focus on study 2 ###

intersect_ID_study2 = set(Illumina_RNA["ID"]).intersection(set(GSM2336845_RNAseq_df["ID_REF"]))
intersect_df_study2 = pd.DataFrame(list(intersect_ID_study2), columns=["ID"])

Illumina_RNA_intersect_study2 = pd.merge(intersect_df_study2, Illumina_RNA, left_on="ID", right_on="ID", how="left")
unique_gene_study2 = Illumina_RNA_intersect_study2["ILMN_Gene"].unique()


def chunk_list(lst, chunk_size):
    
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]



def query_ensembl(unique_gene_batch):
    base_url = "http://www.ensembl.org/biomart/martservice"
    unique_gene_string = ','.join(unique_gene_batch)
    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "external_gene_name" value = "{unique_gene_string}"/>
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "external_gene_name" />
            <Attribute name = "gene_biotype" />
            <Attribute name = "transcript_biotype" />
            <Attribute name = "chromosome_name" />
            <Attribute name = "transcription_start_site" />
            <Attribute name = "transcript_is_canonical" />
        </Dataset>
    </Query>
    """
    response = requests.get(base_url, params={'query': xml_query})
    return response.text.strip()

batches_study2 = chunk_list(unique_gene_study2, 200)
all_results_study2 = []


for batch in tqdm(batches_study2, desc="Querying Ensembl"):
    result = query_ensembl(batch)
    if result:
        all_results_study2.extend(result.split('\n'))


columns = ['Ensembl_Gene_ID', 'Ensembl_Transcript_ID', 'Gene_Name', "Gene_Biotype", "Transcript_Biotype", "Chromosome", "Transcription_Start_Site", 'Is_Canonical']
df_study2 = pd.DataFrame([line.split('\t') for line in all_results_study2], columns=columns)

df_study2 = df_study2[(df_study2["Is_Canonical"] == "1") & (df_study2["Gene_Biotype"] == "protein_coding") & (df_study2["Chromosome"].isin(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]))].copy()

df_study2 = df_study2[["Ensembl_Gene_ID", "Gene_Name", "Chromosome", "Transcription_Start_Site"]]



df_study2 = pd.merge(df_study2, Illumina_RNA_intersect_study2, left_on="Gene_Name", right_on="ILMN_Gene", how="left")
df_study2 = df_study2[["Ensembl_Gene_ID", "Gene_Name", "Chromosome_x", "Transcription_Start_Site", "ID"]]
df_study2.rename(columns={"Chromosome_x": "Chromosome", "ID": "PROBEID"}, inplace=True)



merge_study2_1 = pd.merge(df_study2, GSM2336845_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_1 = merge_study2_1.drop(merge_study2_1.columns[-2], axis=1)
merge_study2_1.rename(columns={"VALUE": "GSM2336845"}, inplace=True)


merge_study2_2 = pd.merge(merge_study2_1, GSM2336851_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_2 = merge_study2_2.drop(merge_study2_2.columns[-2], axis=1)
merge_study2_2.rename(columns={"VALUE": "GSM2336851"}, inplace=True)

merge_study2_3 = pd.merge(merge_study2_2, GSM2336874_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_3 = merge_study2_3.drop(merge_study2_3.columns[-2], axis=1)
merge_study2_3.rename(columns={"VALUE": "GSM2336874"}, inplace=True)

merge_study2_4 = pd.merge(merge_study2_3, GSM2336911_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_4 = merge_study2_4.drop(merge_study2_4.columns[-2], axis=1)
merge_study2_4.rename(columns={"VALUE": "GSM2336911"}, inplace=True)

merge_study2_5 = pd.merge(merge_study2_4, GSM2336922_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_5 = merge_study2_5.drop(merge_study2_5.columns[-2], axis=1)
merge_study2_5.rename(columns={"VALUE": "GSM2336922"}, inplace=True)

merge_study2_6 = pd.merge(merge_study2_5, GSM2336927_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_6 = merge_study2_6.drop(merge_study2_6.columns[-2], axis=1)
merge_study2_6.rename(columns={"VALUE": "GSM2336927"}, inplace=True)

merge_study2_7 = pd.merge(merge_study2_6, GSM2336935_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_7 = merge_study2_7.drop(merge_study2_7.columns[-2], axis=1)
merge_study2_7.rename(columns={"VALUE": "GSM2336935"}, inplace=True)

merge_study2_8 = pd.merge(merge_study2_7, GSM2336941_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_8 = merge_study2_8.drop(merge_study2_8.columns[-2], axis=1)
merge_study2_8.rename(columns={"VALUE": "GSM2336941"}, inplace=True)

merge_study2_9 = pd.merge(merge_study2_8, GSM2336952_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_9 = merge_study2_9.drop(merge_study2_9.columns[-2], axis=1)
merge_study2_9.rename(columns={"VALUE": "GSM2336952"}, inplace=True)

merge_study2_10 = pd.merge(merge_study2_9, GSM2336953_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_10 = merge_study2_10.drop(merge_study2_10.columns[-2], axis=1)
merge_study2_10.rename(columns={"VALUE": "GSM2336953"}, inplace=True)


merge_study2_11 = pd.merge(merge_study2_10, GSM2337002_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_11 = merge_study2_11.drop(merge_study2_11.columns[-2], axis=1)
merge_study2_11.rename(columns={"VALUE": "GSM2337002"}, inplace=True)

merge_study2_12 = pd.merge(merge_study2_11, GSM2337046_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_study2_12 = merge_study2_12.drop(merge_study2_12.columns[-2], axis=1)
merge_study2_12.rename(columns={"VALUE": "GSM2337046"}, inplace=True)

numeric_columns = merge_study2_12.columns[5:]
grouped_means = merge_study2_12.groupby('Gene_Name')[numeric_columns].mean()
additional_info = merge_study2_12.groupby('Gene_Name').first().reset_index()[['Ensembl_Gene_ID', 'Gene_Name', 'Chromosome', 'Transcription_Start_Site']]
result_study2 = pd.merge(grouped_means, additional_info, on='Gene_Name')


CD8T_RNA_study2 = pd.concat([result_study2.iloc[:,[13,0,14,15]], result_study2.iloc[:,1:13]], axis=1)

CD8T_RNA_study2.to_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_study2.csv', index=False)


#######

batches_3studies = chunk_list(unique_gene_3studies, 200)
all_results_3studies = []


for batch in tqdm(batches_3studies, desc="Querying Ensembl"):
    result = query_ensembl(batch)
    if result:
        all_results_3studies.extend(result.split('\n'))


columns = ['Ensembl_Gene_ID', 'Ensembl_Transcript_ID', 'Gene_Name', "Gene_Biotype", "Transcript_Biotype", "Chromosome", "Transcription_Start_Site", 'Is_Canonical']
df_3studies = pd.DataFrame([line.split('\t') for line in all_results_3studies], columns=columns)

df_3studies = df_3studies[(df_3studies["Is_Canonical"] == "1") & (df_3studies["Gene_Biotype"] == "protein_coding") & (df_3studies["Chromosome"].isin(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]))].copy()

df_3studies = df_3studies[["Ensembl_Gene_ID", "Gene_Name", "Chromosome", "Transcription_Start_Site"]]



df_3studies = pd.merge(df_3studies, Illumina_RNA_intersect_3studies, left_on="Gene_Name", right_on="ILMN_Gene", how="left")
df_3studies = df_3studies[["Ensembl_Gene_ID", "Gene_Name", "Chromosome_x", "Transcription_Start_Site", "ID"]]
df_3studies.rename(columns={"Chromosome_x": "Chromosome", "ID": "PROBEID"}, inplace=True)


merge_3studies_6 = pd.merge(df_3studies, GSE83156_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_6 = merge_3studies_6[["Ensembl_Gene_ID", "Gene_Name", "Chromosome", "Transcription_Start_Site", "PROBEID", "N1", "N2", "TEMRA1", "TEMRA2", "EM1", "EM2"]]


merge_3studies_7 = pd.merge(merge_3studies_6, GSM2336845_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_7 = merge_3studies_7.drop(merge_3studies_7.columns[-2], axis=1)
merge_3studies_7.rename(columns={"VALUE": "GSM2336845"}, inplace=True)

merge_3studies_8 = pd.merge(merge_3studies_7, GSM2336851_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_8 = merge_3studies_8.drop(merge_3studies_8.columns[-2], axis=1)
merge_3studies_8.rename(columns={"VALUE": "GSM2336851"}, inplace=True)

merge_3studies_9 = pd.merge(merge_3studies_8, GSM2336874_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_9 = merge_3studies_9.drop(merge_3studies_9.columns[-2], axis=1)
merge_3studies_9.rename(columns={"VALUE": "GSM2336874"}, inplace=True)

merge_3studies_10 = pd.merge(merge_3studies_9, GSM2336911_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_10 = merge_3studies_10.drop(merge_3studies_10.columns[-2], axis=1)
merge_3studies_10.rename(columns={"VALUE": "GSM2336911"}, inplace=True)

merge_3studies_11 = pd.merge(merge_3studies_10, GSM2336922_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_11 = merge_3studies_11.drop(merge_3studies_11.columns[-2], axis=1)
merge_3studies_11.rename(columns={"VALUE": "GSM2336922"}, inplace=True)

merge_3studies_12 = pd.merge(merge_3studies_11, GSM2336927_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_12 = merge_3studies_12.drop(merge_3studies_12.columns[-2], axis=1)
merge_3studies_12.rename(columns={"VALUE": "GSM2336927"}, inplace=True)

merge_3studies_13 = pd.merge(merge_3studies_12, GSM2336935_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_13 = merge_3studies_13.drop(merge_3studies_13.columns[-2], axis=1)
merge_3studies_13.rename(columns={"VALUE": "GSM2336935"}, inplace=True)

merge_3studies_14 = pd.merge(merge_3studies_13, GSM2336941_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_14 = merge_3studies_14.drop(merge_3studies_14.columns[-2], axis=1)
merge_3studies_14.rename(columns={"VALUE": "GSM2336941"}, inplace=True)


merge_3studies_15 = pd.merge(merge_3studies_14, GSM2336952_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_15 = merge_3studies_15.drop(merge_3studies_15.columns[-2], axis=1)
merge_3studies_15.rename(columns={"VALUE": "GSM2336952"}, inplace=True)

merge_3studies_16 = pd.merge(merge_3studies_15, GSM2336953_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_16 = merge_3studies_16.drop(merge_3studies_16.columns[-2], axis=1)
merge_3studies_16.rename(columns={"VALUE": "GSM2336953"}, inplace=True)


merge_3studies_17 = pd.merge(merge_3studies_16, GSM2337002_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_17 = merge_3studies_17.drop(merge_3studies_17.columns[-2], axis=1)
merge_3studies_17.rename(columns={"VALUE": "GSM2337002"}, inplace=True)

merge_3studies_18 = pd.merge(merge_3studies_17, GSM2337046_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_18 = merge_3studies_18.drop(merge_3studies_18.columns[-2], axis=1)
merge_3studies_18.rename(columns={"VALUE": "GSM2337046"}, inplace=True)


merge_3studies_19 = pd.merge(merge_3studies_18, GSM1831330_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_19 = merge_3studies_19.drop(merge_3studies_19.columns[-2], axis=1)
merge_3studies_19.rename(columns={"VALUE": "GSM1831330"}, inplace=True)

merge_3studies_20 = pd.merge(merge_3studies_19, GSM1831331_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_20 = merge_3studies_20.drop(merge_3studies_20.columns[-2], axis=1)
merge_3studies_20.rename(columns={"VALUE": "GSM1831331"}, inplace=True)

merge_3studies_21 = pd.merge(merge_3studies_20, GSM1831332_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_21 = merge_3studies_21.drop(merge_3studies_21.columns[-2], axis=1)
merge_3studies_21.rename(columns={"VALUE": "GSM1831332"}, inplace=True)

merge_3studies_22 = pd.merge(merge_3studies_21, GSM1831333_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_22 = merge_3studies_22.drop(merge_3studies_22.columns[-2], axis=1)
merge_3studies_22.rename(columns={"VALUE": "GSM1831333"}, inplace=True)

merge_3studies_23 = pd.merge(merge_3studies_22, GSM1831334_RNAseq_df, left_on="PROBEID", right_on="ID_REF", how="left")
merge_3studies_23 = merge_3studies_23.drop(merge_3studies_23.columns[-2], axis=1)
merge_3studies_23.rename(columns={"VALUE": "GSM1831334"}, inplace=True)

numeric_columns = merge_3studies_23.columns[5:]


grouped_means = merge_3studies_23.groupby('Gene_Name')[numeric_columns].mean()
additional_info = merge_3studies_23.groupby('Gene_Name').first().reset_index()[['Ensembl_Gene_ID', 'Gene_Name', 'Chromosome', 'Transcription_Start_Site']]
result_3studies = pd.merge(grouped_means, additional_info, on='Gene_Name')


merge_3studies_final = pd.concat([result_3studies.iloc[:,[24,0,25,26]], result_3studies.iloc[:,1:24]], axis=1)


CD8T_RNA_allsamples = merge_3studies_final.copy()

CD8T_RNA_training = merge_3studies_final.iloc[:, [0,1,2,3, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]].copy()
CD8T_RNA_testing = merge_3studies_final.iloc[:, [0,1,2,3, 4,5,6,7,8,9, 22, 23, 24, 25, 26]].copy()

CD8T_RNA_allsamples.to_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_allsamples.csv', index=False)
CD8T_RNA_training.to_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_training.csv', index=False)
CD8T_RNA_testing.to_csv('/Users/scui2/DNAmethylation/CD8T_RNAseq_External_data_3studies/CD8T_RNA_testing.csv', index=False)


###### Check the sample mean and sd #####

sample_mean_rna_data1 = merge_3studies_final.iloc[:,4:10].mean(axis=1)
sample_sd_rna_data1 = merge_3studies_final.iloc[:,4:10].std(axis=1)

sample_mean_rna_data2 = merge_3studies_final.iloc[:,10:22].mean(axis=1)
sample_sd_rna_data2 = merge_3studies_final.iloc[:,10:22].std(axis=1)

sample_mean_rna_data3 = merge_3studies_final.iloc[:,22:].mean(axis=1)
sample_sd_rna_data3 = merge_3studies_final.iloc[:,22:].std(axis=1)

sample_mean_rna_all = merge_3studies_final.iloc[:,4:].mean(axis=1)
sample_sd_rna_all = merge_3studies_final.iloc[:,4:].std(axis=1)

df_sorted = np.sort(merge_3studies_final.iloc[:, list(range(4,10)) + list(range(22,27))].values, axis=0)
df_mean = np.mean(df_sorted, axis=1)
df_ranked = merge_3studies_final.iloc[:,10:22].rank(method="average").fillna(0).astype(int) - 1
merge_data2_quantile = pd.DataFrame(df_mean[df_ranked], index=merge_3studies_final.iloc[:,10:22].index, columns=merge_3studies_final.iloc[:,10:22].columns)


sample_mean_rna_data2_quantile = merge_data2_quantile.mean(axis=1)
sample_sd_rna_data2_quantile = merge_data2_quantile.std(axis=1)





# Compare the distribution of RNA-seq data between the three studies

plt.figure(figsize=(10, 10))
sns.kdeplot(data=sample_mean_rna_data1, color='blue', label='GSE83156 (N=6)', shade=False)
sns.kdeplot(data=sample_mean_rna_data2, color='red', label='GSE87640 (N=12)', shade=False)
sns.kdeplot(data=sample_mean_rna_data3, color='green', label='GSE71244 (N=5)', shade=False)

plt.title('After quantile normalization of GSE87640')
plt.xlabel('Samples\' mean of gene expression')
plt.ylabel('Density')
plt.legend()
plt.tight_layout()
plt.show()


sample_mean_rna_data1.describe()
sample_mean_rna_data2.describe()
sample_mean_rna_data3.describe()



rho, p_value = stats.spearmanr(sample_mean_rna_data2_quantile, sample_mean_rna_data3, nan_policy='omit')

print(f"Spearman correlation coefficient (rho): {rho:.4f}")
print(f"P-value: {p_value:.4e}")



plt.figure(figsize=(10, 10))
scatter = plt.scatter(sample_mean_rna_data2, sample_mean_rna_data3, alpha=0.5, label='RNAseq')
plt.plot([6, 15], [6, 15], color='blue', linestyle='--', linewidth=2, label='Identity Line')
plt.text(0.15, 0.95, f'Spearman rho: {rho:.4f}', transform=plt.gca().transAxes, fontsize=12, ha='center', fontweight='bold', color='red')
plt.xlabel('GSE87640 (N=12) (normalized signal)')
plt.ylabel('GSE71244 (N=5) (normalized signal)')
plt.title('Comparison of mean gene expression between GSE87640 (N=12) and GSE71244 (N=5)')
plt.legend(loc = 'lower right')
plt.grid(True)
plt.show()


merge_3studies_final_quantile =pd.concat([merge_3studies_final.iloc[:,list(range(0,10))], merge_data2_quantile, merge_3studies_final.iloc[:,22:]], axis=1)
merge_3studies_final_quantile.iloc[:,4:] = merge_3studies_final_quantile.iloc[:,4:].round(4)

CD8T_RNA_allsamples = merge_3studies_final_quantile.copy()
CD8T_RNA_allsamples = CD8T_RNA_allsamples.iloc[:,list(range(0,6)) + [8,9] + [6,7] +list(range(10,27))].copy()
CD8T_RNA_allsamples.rename(columns={"N1": "CD8_NAIVE_1",
                                    "N2": "CD8_NAIVE_2",
                                    "TEMRA1": "CD8_TEMRA_1",
                                    "TEMRA2": "CD8_TEMRA_2",
                                    "EM1": "CD8_EM_1",
                                    "EM2": "CD8_EM_2"}, inplace=True)




