# liquid biopsy by cell type proportions
# Project Name

## Overview
This project is designed to [project purpose]. It involves a series of steps to achieve [project functionality].

![Overview Image](/liquid_biopsy_CT/Images/WorkFlow.png)

## Table of Contents
1. [Data Preparation](#data-preparation)
2. [Step A: Find DMP](#step-a-find-dmp)
3. [Step B: Calculate Correlation and Filtering](#step-b-calculate-correlation-and-filtering)
4. [Step C: Build PBMC with Cancer Reactive Atlas](#step-c-build-pbmc-with-cancer-reactive-atlas)
5. [Step D: Deconvolution and Prediction Models](#step-d-deconvolution-and-prediction-models)

## Data Preparation
Before starting, ensure that your data is prepared. Use the following scripts for data preparation:

- `DEgene_analysis_CD4TCR.py`
- `DEgene_analysis_CD8TCR.py`
- `External_CD4T_data_prepare.py`
- `External_CD8T_DNAMeth_data_prepare.py`
- `External_CD8T_RNAseq_data_prepare.py`
- `WGBS_process.py`

These scripts are used to prepare your data.

## Step A: Find DMP
After data preparation, proceed with Step A using the `Find_DMP.py` script. This step involves identifying differentially methylated positions (DMPs).

![Step A Image](/liquid_biopsy_CT/Images/Step_A.png)

## Step B: Calculate Correlation and Filtering
Step B includes several sub-steps involving Spearman correlation calculation, linear model filtering, and gene list selection. The scripts involved are:

![Step B Image](/liquid_biopsy_CT/Images/Step_B.png)

- `Calulate_Spearman_correlation_CD4T_testing.py`
- `Calulate_Spearman_correlation_CD4T_training.py`
- `Calulate_Spearman_correlation_CD8T_testing.py`
- `Calulate_Spearman_correlation_CD8T_training.py`
- `LM_filtering_and_Extrapolation_CD4TCR.py`
- `LM_filtering_and_Extrapolation_CD8TCR.py`
- `Select_Gene_list_A_(CD4T).py`
- `Select_Gene_list_B_(CD8T).py`
- `Sliding_Window_CD4T.py`
- `Sliding_Window_CD8T.py`

These scripts are used for analysis and filtering.

## Step C: Build PBMC with Cancer Reactive Atlas
In Step C, the `Build_PBMC_with_CancerReactive_Atlas.py` script is used to construct a PBMC with a cancer-reactive atlas.
![Step C Image](/liquid_biopsy_CT/Images/Step_C.png)

## Step D: Deconvolution and Prediction Models
Step D involves deconvolution and prediction models. The scripts available for this step are:

- `Deconvolution_Lee_et_al_withTCR.py`
- `Deconvolution_Lee_et_al_withoutTCR.py`
- `Deconvolution_Li_et_al_withTCR.py`
- `Deconvolution_Li_et_al_withoutTCR.py`
- `Deconvolution_Wang_et_al_withTCR.py`
- `Deconvolution_Wang_et_al_withoutTCR.py`
- `Deconvolution_Xie_et_al_withTCR.py`
- `Deconvolution_Xie_et_al_withoutTCR.py`
- `Predicting_Models_withTCR.py`
- `Prediction_Models_withoutTCR.py`

These scripts are used for deconvolution and making predictions.

## Contribution
If you wish to contribute to this project, please submit a Pull Request or contact the project maintainers.

## License
This project is licensed under the [License Name] License.
