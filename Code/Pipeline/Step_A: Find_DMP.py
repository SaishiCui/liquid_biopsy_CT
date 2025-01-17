import multiprocessing
from multiprocessing import Pool
import pandas as pd
from scipy.stats import mannwhitneyu
from tqdm import tqdm

"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-16 by Saishi Cui
Program: Step_A: Find_DMP.py
Purpose: Find differentially methylated Probes (DMPs) between different cell types in WGBS PBMC DNAmethylation data.
--------------------------------------------------------------------------------
Data Inputs:

- Pre-processed WGBS PBMC DNAmethylation data.
- EPIC 850K DNAmethylation data.

Data Outputs:
- Differential Methylated Probes (DMPs) between different cell types in WGBS PBMC DNAmethylation data, not filtered, just raw results with difference and p-value.
--------------------------------------------------------------------------------
Functions:
- process_dmps_chunk: First match WGBS data with EPIC 850K data, then calculate the difference and p-value.
- Parallel_process_DMP_850K: Parallel processing, to avoid memory issue.
--------------------------------------------------------------------------------
Notes:
- There are two versions of output, one for DMPs with Granu, one for DMPs without Granu.
-------------------------------------------------------------------------------
"""

# Set up multiprocessing
multiprocessing.set_start_method('fork', force=True)

# Read pre-processed WGBS PBMC DNAmethylation data
WGBS_allsamples_combined_filtered_withGranu = pd.read_csv('/Users/scui2/DNAmethylation/Loyfer_data/WGBS_allsamples_combined_filtered_withGranu.csv')
WGBS_allsamples_combined_filtered = pd.read_csv('/Users/scui2/DNAmethylation/Loyfer_data/WGBS_allsamples_combined_filtered.csv')


# Read 850K EPIC data

Xie = pd.read_csv('/Users/scui2/DNAmethylation/Xie_et_al/Xie_et_al.csv') 
Epic_850K = Xie[["probeID", "CHR_hg38", "Start_hg38", "End_hg38"]]

Epic_850K_clean = Epic_850K.dropna(inplace=False)
Epic_850K_clean.reset_index(drop=True, inplace=True)


## Match WGBS data with EPIC 850K data, then calculate the difference and p-value.
def process_dmps_chunk_v1(args_chunk, WGBS_data = WGBS_allsamples_combined_filtered_withGranu):

    results = []
    for i in args_chunk:
        chr, start, end = Epic_850K_clean.iloc[i, 1], Epic_850K_clean.iloc[i, 2], Epic_850K_clean.iloc[i, 3]
        start = int(float(start))
        end = int(float(end))
        
        condition = (WGBS_data['Chromosome'] == chr) & \
                   (WGBS_data['Position'] >= start) & \
                   (WGBS_data['Position'] <= end)
        
        if condition.sum() == 0:
            continue

        filtered_df = WGBS_data[condition]

        CD4T_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta"]].mean(axis =0)     
        CD8T_filtered_df = filtered_df[["CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta"]].mean(axis =0)
        B_filtered_df = filtered_df[["B_1_Beta", "B_2_Beta", "B_3_Beta"]].mean(axis =0)
        NK_filtered_df = filtered_df[["NK_1_Beta", "NK_2_Beta", "NK_3_Beta"]].mean(axis =0)
        Mono_filtered_df = filtered_df[["Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta"]].mean(axis =0)
        Granu_filtered_df = filtered_df[["Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis =0)
        Others_except_CD4T_filtered_df = filtered_df[["CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_CD8T_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_B_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_NK_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_Mono_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_Granu_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta"]].mean(axis=0)


        # Perform Wilcoxon signed-rank test
        try:
            # Perform Wilcoxon signed-rank test
            _, p_value_CD4T_Others = mannwhitneyu(CD4T_filtered_df.values, Others_except_CD4T_filtered_df.values, alternative='two-sided')
            _, p_value_CD8T_Others = mannwhitneyu(CD8T_filtered_df.values, Others_except_CD8T_filtered_df.values, alternative='two-sided')
            _, p_value_B_Others = mannwhitneyu(B_filtered_df.values, Others_except_B_filtered_df.values, alternative='two-sided')
            _, p_value_NK_Others = mannwhitneyu(NK_filtered_df.values, Others_except_NK_filtered_df.values, alternative='two-sided')
            _, p_value_Mono_Others = mannwhitneyu(Mono_filtered_df.values, Others_except_Mono_filtered_df.values, alternative='two-sided')
            _, p_value_Granu_Others = mannwhitneyu(Granu_filtered_df.values, Others_except_Granu_filtered_df.values, alternative='two-sided')
            _, p_value_CD4T_CD8T = mannwhitneyu(CD4T_filtered_df.values, CD8T_filtered_df.values, alternative='two-sided')
            
        except ValueError as e:
            print(f"Error in window {chr}_{start}_{end}: {e}")
            continue

        Diff_CD4T_Others = CD4T_filtered_df.mean() - Others_except_CD4T_filtered_df.mean()
        Diff_CD8T_Others = CD8T_filtered_df.mean() - Others_except_CD8T_filtered_df.mean()
        Diff_B_Others = B_filtered_df.mean() - Others_except_B_filtered_df.mean()
        Diff_NK_Others = NK_filtered_df.mean() - Others_except_NK_filtered_df.mean()
        Diff_Mono_Others = Mono_filtered_df.mean() - Others_except_Mono_filtered_df.mean()
        Diff_Granu_Others = Granu_filtered_df.mean() - Others_except_Granu_filtered_df.mean()
        Diff_CD4T_CD8T = CD4T_filtered_df.mean() - CD8T_filtered_df.mean()
        


        results.append({
            "ProbeID": Epic_850K_clean.iloc[i, 0],
            "Window_Location": f"{chr}_{start}_{end}",
            "CD4T_Others_mean_diff": Diff_CD4T_Others,
            "CD4T_Others_p_value": p_value_CD4T_Others,
            "CD8T_Others_mean_diff": Diff_CD8T_Others,
            "CD8T_Others_p_value": p_value_CD8T_Others,
            "B_Others_mean_diff": Diff_B_Others,
            "B_Others_p_value": p_value_B_Others,
            "NK_Others_mean_diff": Diff_NK_Others,
            "NK_Others_p_value": p_value_NK_Others,
            "Mono_Others_mean_diff": Diff_Mono_Others,
            "Mono_Others_p_value": p_value_Mono_Others,
            "Granu_Others_mean_diff": Diff_Granu_Others,
            "Granu_Others_p_value": p_value_Granu_Others,
            "CD4T_CD8T_mean_diff": Diff_CD4T_CD8T,
            "CD4T_CD8T_p_value": p_value_CD4T_CD8T
        })
    return results

def process_dmps_chunk_v2(args_chunk, WGBS_data = WGBS_allsamples_combined_filtered):

    results = []
    for i in args_chunk:
        chr, start, end = Epic_850K_clean.iloc[i, 1], Epic_850K_clean.iloc[i, 2], Epic_850K_clean.iloc[i, 3]
        start = int(float(start))
        end = int(float(end))
        
        condition = (WGBS_data['Chromosome'] == chr) & \
                   (WGBS_data['Position'] >= start) & \
                   (WGBS_data['Position'] <= end)
        
        if condition.sum() == 0:
            continue

        filtered_df = WGBS_data[condition]

        CD4T_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta"]].mean(axis =0)     
        CD8T_filtered_df = filtered_df[["CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta"]].mean(axis =0)
        B_filtered_df = filtered_df[["B_1_Beta", "B_2_Beta", "B_3_Beta"]].mean(axis =0)
        NK_filtered_df = filtered_df[["NK_1_Beta", "NK_2_Beta", "NK_3_Beta"]].mean(axis =0)
        Mono_filtered_df = filtered_df[["Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta"]].mean(axis =0)
        Granu_filtered_df = filtered_df[["Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis =0)
        Others_except_CD4T_filtered_df = filtered_df[["CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_CD8T_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_B_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_NK_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_Mono_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Granu_1_Beta", "Granu_2_Beta", "Granu_3_Beta"]].mean(axis=0)
        Others_except_Granu_filtered_df = filtered_df[["CD4T_1_Beta", "CD4T_2_Beta", "CD4T_3_Beta", "CD8T_1_Beta", "CD8T_2_Beta", "CD8T_3_Beta", "B_1_Beta", "B_2_Beta", "B_3_Beta", "NK_1_Beta", "NK_2_Beta", "NK_3_Beta", "Mono_1_Beta", "Mono_2_Beta", "Mono_3_Beta"]].mean(axis=0)


        # Perform Wilcoxon signed-rank test
        try:
            # Perform Wilcoxon signed-rank test
            _, p_value_CD4T_Others = mannwhitneyu(CD4T_filtered_df.values, Others_except_CD4T_filtered_df.values, alternative='two-sided')
            _, p_value_CD8T_Others = mannwhitneyu(CD8T_filtered_df.values, Others_except_CD8T_filtered_df.values, alternative='two-sided')
            _, p_value_B_Others = mannwhitneyu(B_filtered_df.values, Others_except_B_filtered_df.values, alternative='two-sided')
            _, p_value_NK_Others = mannwhitneyu(NK_filtered_df.values, Others_except_NK_filtered_df.values, alternative='two-sided')
            _, p_value_Mono_Others = mannwhitneyu(Mono_filtered_df.values, Others_except_Mono_filtered_df.values, alternative='two-sided')
            _, p_value_Granu_Others = mannwhitneyu(Granu_filtered_df.values, Others_except_Granu_filtered_df.values, alternative='two-sided')
            _, p_value_CD4T_CD8T = mannwhitneyu(CD4T_filtered_df.values, CD8T_filtered_df.values, alternative='two-sided')
            
        except ValueError as e:
            print(f"Error in window {chr}_{start}_{end}: {e}")
            continue

        Diff_CD4T_Others = CD4T_filtered_df.mean() - Others_except_CD4T_filtered_df.mean()
        Diff_CD8T_Others = CD8T_filtered_df.mean() - Others_except_CD8T_filtered_df.mean()
        Diff_B_Others = B_filtered_df.mean() - Others_except_B_filtered_df.mean()
        Diff_NK_Others = NK_filtered_df.mean() - Others_except_NK_filtered_df.mean()
        Diff_Mono_Others = Mono_filtered_df.mean() - Others_except_Mono_filtered_df.mean()
        Diff_Granu_Others = Granu_filtered_df.mean() - Others_except_Granu_filtered_df.mean()
        Diff_CD4T_CD8T = CD4T_filtered_df.mean() - CD8T_filtered_df.mean()
        


        results.append({
            "ProbeID": Epic_850K_clean.iloc[i, 0],
            "Window_Location": f"{chr}_{start}_{end}",
            "CD4T_Others_mean_diff": Diff_CD4T_Others,
            "CD4T_Others_p_value": p_value_CD4T_Others,
            "CD8T_Others_mean_diff": Diff_CD8T_Others,
            "CD8T_Others_p_value": p_value_CD8T_Others,
            "B_Others_mean_diff": Diff_B_Others,
            "B_Others_p_value": p_value_B_Others,
            "NK_Others_mean_diff": Diff_NK_Others,
            "NK_Others_p_value": p_value_NK_Others,
            "Mono_Others_mean_diff": Diff_Mono_Others,
            "Mono_Others_p_value": p_value_Mono_Others,
            "Granu_Others_mean_diff": Diff_Granu_Others,
            "Granu_Others_p_value": p_value_Granu_Others,
            "CD4T_CD8T_mean_diff": Diff_CD4T_CD8T,
            "CD4T_CD8T_p_value": p_value_CD4T_CD8T
        })
    return results

## Parallel processing, to avoid memory issue.
def Parallel_process_DMP_850K_v1(start_chunk=0):

    total_dmps = len(Epic_850K_clean)
    chunk_size = 2  
    dmp_indices = list(range(total_dmps))
    chunks = [dmp_indices[i:i + chunk_size] for i in range(0, total_dmps, chunk_size)]
    
    remaining_chunks = chunks[start_chunk:]
    all_results = []
    
    print(f"Total {len(chunks)} chunks to process")
    print(f"Starting from chunk {start_chunk}")
    print("Press Ctrl+C to pause processing")
    
    try:
        with Pool(processes=multiprocessing.cpu_count()-1) as pool:
            for i, chunk_result in enumerate(tqdm(
                pool.imap(process_dmps_chunk_v1, remaining_chunks),
                total=len(remaining_chunks),
                desc="Processing DMPs"
            )):
                all_results.extend(chunk_result)
                
                if (i + 1) % 10000 == 0:
                    temp_df = pd.DataFrame(all_results)
                    temp_df.to_csv(f'/Users/scui2/DNAmethylation/DMP/DMP_850K_results_temp_withGranu_{start_chunk + i}.csv', index=False)
                    print(f"\nSaved temporary results up to chunk {start_chunk + i}")
                
    except KeyboardInterrupt:
        print(f"\nProcessing interrupted at chunk {start_chunk + i}")
        print("You can resume processing by calling check_DMP_850K_parallel(start_chunk={})".format(start_chunk + i))
        return all_results, start_chunk + i

    final_df = pd.DataFrame(all_results)
    final_df.to_csv('/Users/scui2/DNAmethylation/DMP/DMP_850K_results_final_withGranu.csv', index=False)
    print("\nProcessing complete, results saved to DMP_results_final.csv")
    
    return all_results, len(chunks)

def Parallel_process_DMP_850K_v2(start_chunk=0):

    total_dmps = len(Epic_850K_clean)
    chunk_size = 2  
    dmp_indices = list(range(total_dmps))
    chunks = [dmp_indices[i:i + chunk_size] for i in range(0, total_dmps, chunk_size)]
    
    remaining_chunks = chunks[start_chunk:]
    all_results = []
    
    print(f"Total {len(chunks)} chunks to process")
    print(f"Starting from chunk {start_chunk}")
    print("Press Ctrl+C to pause processing")
    
    try:
        with Pool(processes=multiprocessing.cpu_count()-1) as pool:
            for i, chunk_result in enumerate(tqdm(
                pool.imap(process_dmps_chunk_v2, remaining_chunks),
                total=len(remaining_chunks),
                desc="Processing DMPs"
            )):
                all_results.extend(chunk_result)
                
                if (i + 1) % 10000 == 0:
                    temp_df = pd.DataFrame(all_results)
                    temp_df.to_csv(f'/Users/scui2/DNAmethylation/DMP/DMP_850K_results_temp_withGranu_{start_chunk + i}.csv', index=False)
                    print(f"\nSaved temporary results up to chunk {start_chunk + i}")
                
    except KeyboardInterrupt:
        print(f"\nProcessing interrupted at chunk {start_chunk + i}")
        print("You can resume processing by calling check_DMP_850K_parallel(start_chunk={})".format(start_chunk + i))
        return all_results, start_chunk + i

    final_df = pd.DataFrame(all_results)
    final_df.to_csv('/Users/scui2/DNAmethylation/DMP/DMP_850K_results_final_withGranu.csv', index=False)
    print("\nProcessing complete, results saved to DMP_results_final.csv")
    
    return all_results, len(chunks)



if __name__ == '__main__': 
    Parallel_process_DMP_850K_v1(start_chunk=0)
    Parallel_process_DMP_850K_v2(start_chunk=0)


