import pandas as pd
import numpy as np
import random
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from statsmodels.stats.multitest import multipletests
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm
from scipy import stats
from sklearn.metrics import r2_score
import time




"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: External_CD4T_data_prepare.py
Purpose: Prepare external CD4+T healthy samples data for both RNA-seq and DNAmethylation data, and split into training and testing data.
--------------------------------------------------------------------------------
Data Inputs:

- Chen et al.,2016's CD4+T healthy samples data, including samples matched methylation and RNA-seq data.

Data Outputs:
- CD4+T healthy samples data for both RNA-seq and DNAmethylation data, and split into training and testing data.
--------------------------------------------------------------------------------
Functions:
-  chunk_list: Split the Ensembl_ID into chunks, to avoid network issue.
-  query_ensembl_with_retry: Query gene information from Ensembl.
--------------------------------------------------------------------------------
Notes:
- 
-------------------------------------------------------------------------------
"""





## Read External healthy CD4+T samples data, including methylation and RNA-seq data. (Samples matched)

CD4T_Meth_data = pd.read_csv('/Users/scui2/DNAmethylation/Chen_Meth_data/CD4T_Meth_data.csv', index_col=None)
CD4T_RNA_data = pd.read_csv('/Users/scui2/DNAmethylation/Chen_RNA_seq_data/RNA_seq_data.csv')
CD4T_RNA_meta_data = pd.read_csv('/Users/scui2/DNAmethylation/Chen_RNA_seq_data/Meta_RNAseq.csv')


counts = CD4T_RNA_data.iloc[:, 1:]  
Ensembl_IDs = CD4T_RNA_data.iloc[:, 0] 


non_zero_proportion = (counts > 0).sum(axis=1)/ counts.shape[1]

### remove genes with less than 20% non-zero proportion
counts = counts[non_zero_proportion > 0.2].reset_index(drop=True)
Ensembl_IDs = Ensembl_IDs[non_zero_proportion > 0.2].reset_index(drop=True)

### CPM normalization
library_sizes = counts.sum(axis=0)  # 每个样本的总reads数
log_cpm = np.log2( (counts * 1e6 / library_sizes) + 1 )
CD4T_RNA_data = pd.concat([Ensembl_IDs, log_cpm], axis=1)



##### Find gene information ######
CD4T_RNA_data.rename(columns={"Unnamed: 0":"Ensembl_ID"}, inplace=True)
Ensembl_ID = CD4T_RNA_data["Ensembl_ID"]

### Split the Ensembl_ID into chunks, to avoid network issue.
def chunk_list(lst, chunk_size):
    
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

### Query gene information from Ensembl
def query_ensembl_with_retry(unique_gene_batch, max_retries=100, timeout=10):
    base_url = "http://www.ensembl.org/biomart/martservice"
    unique_gene_string = ','.join(unique_gene_batch)
    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "ensembl_gene_id" value = "{unique_gene_string}"/>
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
    
   
    retry_strategy = Retry(
        total=max_retries,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session = requests.Session()
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    
    for attempt in range(max_retries):
        try:
            response = session.get(base_url, 
                                 params={'query': xml_query}, 
                                 timeout=timeout)
            response.raise_for_status()
            
            # 检查返回的内容
            content = response.text.strip()
            if content.startswith('<html>'):
                print(f"Rate limit hit, waiting for 10 seconds...")
                time.sleep(10)  # 等待60秒
                continue
                
            if not content:  # 如果内容为空
                raise requests.exceptions.RequestException("Empty response")
                
            # 验证第一行是否符合预期格式（包含tab分隔的字段）
            first_line = content.split('\n')[0]
            if len(first_line.split('\t')) < 5:  # 应该有至少5个字段
                raise requests.exceptions.RequestException("Invalid response format")
                
            return content
            
        except requests.exceptions.RequestException as e:
            if attempt == max_retries - 1:
                print(f"Failed after {max_retries} attempts for batch: {unique_gene_string[:50]}...")
                return None
            time.sleep(2 ** attempt)
            continue

batches = chunk_list(Ensembl_ID, 200)
all_results = []


for batch in tqdm(batches, desc="Querying Ensembl"):
    result = query_ensembl_with_retry(batch)
    if result:
        print(result.split('\n')[0])
        all_results.extend(result.split('\n'))
    else:
        print(f"Skipping failed batch")
    
    
    time.sleep(5)  



### Create a dataframe
columns = ['Ensembl_Gene_ID', 'Ensembl_Transcript_ID', 'Gene_Name', "Gene_Biotype", "Transcript_Biotype", "Chromosome", "Transcription_Start_Site", 'Is_Canonical']
df = pd.DataFrame([line.split('\t') for line in all_results], columns=columns)


df = df[(df["Is_Canonical"] == "1") & (df["Gene_Biotype"] == "protein_coding") & (df["Chromosome"].isin(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]))].copy()


df = df[["Ensembl_Gene_ID", "Ensembl_Transcript_ID", "Gene_Name", "Chromosome", "Transcription_Start_Site"]].reset_index(drop=True)


CD4T_RNA_data = CD4T_RNA_data[CD4T_RNA_data["Ensembl_ID"].isin(df["Ensembl_Gene_ID"])].reset_index(drop=True)



# Match the samples


CD4T_Meth_allsubjects_name = CD4T_Meth_data.columns[4:].tolist()


match_df = pd.DataFrame()
for i in range(len(CD4T_Meth_allsubjects_name)):
    condition = (CD4T_RNA_meta_data["individual"] == CD4T_Meth_allsubjects_name[i]) & (CD4T_RNA_meta_data["cell_type"] == "CD4-positive, alpha-beta T cell")
    if CD4T_RNA_meta_data[condition].iloc[:,0].tolist() == []:
        continue
    selected_sample = random.choice(CD4T_RNA_meta_data[condition].iloc[:,0].tolist())
    match_df = pd.concat([match_df, pd.DataFrame([[CD4T_Meth_allsubjects_name[i], selected_sample]])], axis=0)

match_df.rename(columns={0:"methyl_sample", 1:"rna_sample"}, inplace=True)

### 126 samples matched

CD4T_RNA_allsamples_matched = pd.concat([df, CD4T_RNA_data[match_df["rna_sample"].tolist()]], axis=1)
CD4T_Meth_allsamples_matched = CD4T_Meth_data[CD4T_Meth_data.columns[0:4].tolist() + match_df["methyl_sample"].tolist()]



##### Split into training and testing data ######


CD4T_RNA_training = CD4T_RNA_allsamples_matched.iloc[:, [0,1,2,3,4] + list(np.arange(5,93))]
CD4T_RNA_testing = CD4T_RNA_allsamples_matched.iloc[:, [0,1,2,3,4] + list(np.arange(93,131))]

CD4T_Meth_training = CD4T_Meth_allsamples_matched.iloc[:, [0,1,2,3] + list(np.arange(4,92))]
CD4T_Meth_testing = CD4T_Meth_allsamples_matched.iloc[:, [0,1,2,3] + list(np.arange(92,130))]


CD4T_RNA_allsamples_matched.to_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_RNA_allsamples_matched.csv', index=False)
CD4T_Meth_allsamples_matched.to_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_Meth_allsamples_matched.csv', index=False)

CD4T_RNA_training.to_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_RNA_training.csv', index=False)
CD4T_RNA_testing.to_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_RNA_testing.csv', index=False)
CD4T_Meth_training.to_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_Meth_training.csv', index=False)
CD4T_Meth_testing.to_csv('/Users/scui2/DNAmethylation/Chen_CD4T_External_data/CD4T_Meth_testing.csv', index=False)



