import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import statsmodels.api as sm
import numpy as np
import pandas as pd


"""
--------------------------------------------------------------------------------
Date last modified: 2025-01-17 by Saishi Cui
Program: Step_D: Predicting_Models_withTCR.py
Purpose: Predicting models using deconvolution results with TCR.
--------------------------------------------------------------------------------
Data Inputs:
- Deconvolution results (with TCR) obtained from Step_D: Deconvolution_xxx_et_al_withTCR.py under the folder "Pipeline".

Data Outputs:
- ROC curves of predicting models using deconvolution results with TCR.
- Logistic regression to predict cancer vs normal using deconvolution results with TCR.
- Machine learning models using deconvolution results with TCR.
--------------------------------------------------------------------------------
"""

nnls_df_Lee_TCR = pd.read_csv('/Users/scui2/DNAmethylation/Lee_et_al/nnls_df_Lee_TCR.csv')
nnls_df_Li_TCR = pd.read_csv('/Users/scui2/DNAmethylation/Li_et_al/nnls_df_Li_TCR.csv')
nnls_df_Wang_TCR = pd.read_csv('/Users/scui2/DNAmethylation/Wang_et_al/nnls_df_Wang_TCR.csv')
nnls_df_Xie_TCR = pd.read_csv('/Users/scui2/DNAmethylation/Xie_et_al/nnls_df_Xie_TCR.csv')

nnls_df_combined_TCR = pd.concat([nnls_df_Lee_TCR, nnls_df_Li_TCR, nnls_df_Wang_TCR, nnls_df_Xie_TCR], axis = 0)

nnls_df_combined_TCR["Group"][nnls_df_combined_TCR["Group"] != "Normal"] = "Cancer"
nnls_df_combined_TCR["Group"][nnls_df_combined_TCR["Group"] == "Normal"] = "Normal"


####### Logistic regression
### CD4+T 

X_Xie = nnls_df_Xie_TCR.iloc[:, [0,1,2,3,4,5,6]]
y_Xie = np.array(['Cancer'] * 50 + ['Normal'] * 50)
y_num_Xie = (y_Xie == 'Cancer').astype(int)


X_Wang = nnls_df_Wang_TCR.iloc[:, [0,1,2,3,4,5,6]]
y_Wang = np.array(['Cancer'] * 50 + ['Normal'] * 30)
y_num_Wang = (y_Wang == 'Cancer').astype(int)


X_Li = nnls_df_Li_TCR.iloc[:, [0,1,2,3,4,5,6]]
y_Li = np.array(['Cancer'] * 35 + ['Normal'] * 50)
y_num_Li = (y_Li == 'Cancer').astype(int)

X_Lee = nnls_df_Lee_TCR.iloc[:, [0,1,2,3,4,5,6]]
y_Lee = np.array(['Cancer'] * 256 + ['Normal'] * 268)
y_num_Lee = (y_Lee == 'Cancer').astype(int)



plt.figure(figsize=(10, 8))

X_CD4T_Xie = X_Xie[['CD4+T']].values
X_sm_CD4T_Xie = sm.add_constant(X_CD4T_Xie)
logit_CD4T_Xie = sm.Logit(y_num_Xie, X_sm_CD4T_Xie)
results_CD4T_Xie = logit_CD4T_Xie.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4T_Xie = results_CD4T_Xie.pvalues[1]
y_score_CD4T_Xie = results_CD4T_Xie.predict(X_sm_CD4T_Xie)
fpr_CD4T_Xie, tpr_CD4T_Xie, _ = roc_curve(y_num_Xie, y_score_CD4T_Xie)
roc_auc_CD4T_Xie = auc(fpr_CD4T_Xie, tpr_CD4T_Xie)


X_CD4T_Wang = X_Wang[['CD4+T']].values
X_sm_CD4T_Wang = sm.add_constant(X_CD4T_Wang)
logit_CD4T_Wang = sm.Logit(y_num_Wang, X_sm_CD4T_Wang)
results_CD4T_Wang = logit_CD4T_Wang.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4T_Wang = results_CD4T_Wang.pvalues[1]
y_score_CD4T_Wang = results_CD4T_Wang.predict(X_sm_CD4T_Wang)
fpr_CD4T_Wang, tpr_CD4T_Wang, _ = roc_curve(y_num_Wang, y_score_CD4T_Wang)
roc_auc_CD4T_Wang = auc(fpr_CD4T_Wang, tpr_CD4T_Wang)


X_CD4T_Li = X_Li[['CD4+T']].values
X_sm_CD4T_Li = sm.add_constant(X_CD4T_Li)
logit_CD4T_Li = sm.Logit(y_num_Li, X_sm_CD4T_Li)
results_CD4T_Li = logit_CD4T_Li.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4T_Li = results_CD4T_Li.pvalues[1]
y_score_CD4T_Li = results_CD4T_Li.predict(X_sm_CD4T_Li)
fpr_CD4T_Li, tpr_CD4T_Li, _ = roc_curve(y_num_Li, y_score_CD4T_Li)
roc_auc_CD4T_Li = auc(fpr_CD4T_Li, tpr_CD4T_Li)



X_CD4T_Lee = X_Lee[['CD4+T']].values
X_sm_CD4T_Lee = sm.add_constant(X_CD4T_Lee)
logit_CD4T_Lee = sm.Logit(y_num_Lee, X_sm_CD4T_Lee)
results_CD4T_Lee = logit_CD4T_Lee.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4T_Lee = results_CD4T_Lee.pvalues[1]
y_score_CD4T_Lee = results_CD4T_Lee.predict(X_sm_CD4T_Lee)
fpr_CD4T_Lee, tpr_CD4T_Lee, _ = roc_curve(y_num_Lee, y_score_CD4T_Lee)
roc_auc_CD4T_Lee = auc(fpr_CD4T_Lee, tpr_CD4T_Lee)



# Plot ROC curve
plt.plot(fpr_CD4T_Xie, tpr_CD4T_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_CD4T_Xie:.2f}, p = {p_CD4T_Xie:.3})')
plt.plot(fpr_CD4T_Wang, tpr_CD4T_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_CD4T_Wang:.2f}, p = {p_CD4T_Wang:.3})')
plt.plot(fpr_CD4T_Li, tpr_CD4T_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_CD4T_Li:.2f}, p = {p_CD4T_Li:.3})')
plt.plot(fpr_CD4T_Lee, tpr_CD4T_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_CD4T_Lee:.2f}, p = {p_CD4T_Lee:.3})')

# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('CD4+T cells', fontweight='bold')
plt.legend(loc='lower right')
plt.show()


### CD4+TCR

plt.figure(figsize=(10, 8))

X_CD4TCR_Xie = X_Xie[['CD4+TCR']].values
X_sm_CD4TCR_Xie = sm.add_constant(X_CD4TCR_Xie)
logit_CD4TCR_Xie = sm.Logit(y_num_Xie, X_sm_CD4TCR_Xie)
results_CD4TCR_Xie = logit_CD4TCR_Xie.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4TCR_Xie = results_CD4TCR_Xie.pvalues[1]
y_score_CD4TCR_Xie = results_CD4TCR_Xie.predict(X_sm_CD4TCR_Xie)
fpr_CD4TCR_Xie, tpr_CD4TCR_Xie, _ = roc_curve(y_num_Xie, y_score_CD4TCR_Xie)
roc_auc_CD4TCR_Xie = auc(fpr_CD4TCR_Xie, tpr_CD4TCR_Xie)


X_CD4TCR_Wang = X_Wang[['CD4+TCR']].values
X_sm_CD4TCR_Wang = sm.add_constant(X_CD4TCR_Wang)
logit_CD4TCR_Wang = sm.Logit(y_num_Wang, X_sm_CD4TCR_Wang)
results_CD4TCR_Wang = logit_CD4TCR_Wang.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4TCR_Wang = results_CD4TCR_Wang.pvalues[1]
y_score_CD4TCR_Wang = results_CD4TCR_Wang.predict(X_sm_CD4TCR_Wang)
fpr_CD4TCR_Wang, tpr_CD4TCR_Wang, _ = roc_curve(y_num_Wang, y_score_CD4TCR_Wang)
roc_auc_CD4TCR_Wang = auc(fpr_CD4TCR_Wang, tpr_CD4TCR_Wang)


X_CD4TCR_Li = X_Li[['CD4+TCR']].values
X_sm_CD4TCR_Li = sm.add_constant(X_CD4TCR_Li)
logit_CD4TCR_Li = sm.Logit(y_num_Li, X_sm_CD4TCR_Li)
results_CD4TCR_Li = logit_CD4TCR_Li.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4TCR_Li = results_CD4TCR_Li.pvalues[1]
y_score_CD4TCR_Li = results_CD4TCR_Li.predict(X_sm_CD4TCR_Li)
fpr_CD4TCR_Li, tpr_CD4TCR_Li, _ = roc_curve(y_num_Li, y_score_CD4TCR_Li)
roc_auc_CD4TCR_Li = auc(fpr_CD4TCR_Li, tpr_CD4TCR_Li)



X_CD4TCR_Lee = X_Lee[['CD4+TCR']].values
X_sm_CD4TCR_Lee = sm.add_constant(X_CD4TCR_Lee)
logit_CD4TCR_Lee = sm.Logit(y_num_Lee, X_sm_CD4TCR_Lee)
results_CD4TCR_Lee = logit_CD4TCR_Lee.fit(method='bfgs', maxiter=1000, disp=False)
p_CD4TCR_Lee = results_CD4TCR_Lee.pvalues[1]
y_score_CD4TCR_Lee = results_CD4TCR_Lee.predict(X_sm_CD4TCR_Lee)
fpr_CD4TCR_Lee, tpr_CD4TCR_Lee, _ = roc_curve(y_num_Lee, y_score_CD4TCR_Lee)
roc_auc_CD4TCR_Lee = auc(fpr_CD4TCR_Lee, tpr_CD4TCR_Lee)



# Plot ROC curve
plt.plot(fpr_CD4TCR_Xie, tpr_CD4TCR_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_CD4TCR_Xie:.2f}, p = {p_CD4TCR_Xie:.3})')
plt.plot(fpr_CD4TCR_Wang, tpr_CD4TCR_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_CD4TCR_Wang:.2f}, p = {p_CD4TCR_Wang:.3})')
plt.plot(fpr_CD4TCR_Li, tpr_CD4TCR_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_CD4TCR_Li:.2f}, p = {p_CD4TCR_Li:.3})')
plt.plot(fpr_CD4TCR_Lee, tpr_CD4TCR_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_CD4TCR_Lee:.2f}, p = {p_CD4TCR_Lee:.3})')

# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('CD4+T cancer reactive cells', fontweight='bold')
plt.legend(loc='lower right')
plt.show()





### CD8+T

plt.figure(figsize=(10, 8))

X_CD8T_Xie = X_Xie[['CD8+T']].values
X_sm_CD8T_Xie = sm.add_constant(X_CD8T_Xie)
logit_CD8T_Xie = sm.Logit(y_num_Xie, X_sm_CD8T_Xie)
results_CD8T_Xie = logit_CD8T_Xie.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8T_Xie = results_CD8T_Xie.pvalues[1]
y_score_CD8T_Xie = results_CD8T_Xie.predict(X_sm_CD8T_Xie)
fpr_CD8T_Xie, tpr_CD8T_Xie, _ = roc_curve(y_num_Xie, y_score_CD8T_Xie)
roc_auc_CD8T_Xie = auc(fpr_CD8T_Xie, tpr_CD8T_Xie)


X_CD8T_Wang = X_Wang[['CD8+T']].values
X_sm_CD8T_Wang = sm.add_constant(X_CD8T_Wang)
logit_CD8T_Wang = sm.Logit(y_num_Wang, X_sm_CD8T_Wang)
results_CD8T_Wang = logit_CD8T_Wang.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8T_Wang = results_CD8T_Wang.pvalues[1]
y_score_CD8T_Wang = results_CD8T_Wang.predict(X_sm_CD8T_Wang)
fpr_CD8T_Wang, tpr_CD8T_Wang, _ = roc_curve(y_num_Wang, y_score_CD8T_Wang)
roc_auc_CD8T_Wang = auc(fpr_CD8T_Wang, tpr_CD8T_Wang)


X_CD8T_Li = X_Li[['CD8+T']].values
X_sm_CD8T_Li = sm.add_constant(X_CD8T_Li)
logit_CD8T_Li = sm.Logit(y_num_Li, X_sm_CD8T_Li)
results_CD8T_Li = logit_CD8T_Li.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8T_Li = results_CD8T_Li.pvalues[1]
y_score_CD8T_Li = results_CD8T_Li.predict(X_sm_CD8T_Li)
fpr_CD8T_Li, tpr_CD8T_Li, _ = roc_curve(y_num_Li, y_score_CD8T_Li)
roc_auc_CD8T_Li = auc(fpr_CD8T_Li, tpr_CD8T_Li)



X_CD8T_Lee = X_Lee[['CD8+T']].values
X_sm_CD8T_Lee = sm.add_constant(X_CD8T_Lee)
logit_CD8T_Lee = sm.Logit(y_num_Lee, X_sm_CD8T_Lee)
results_CD8T_Lee = logit_CD8T_Lee.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8T_Lee = results_CD8T_Lee.pvalues[1]
y_score_CD8T_Lee = results_CD8T_Lee.predict(X_sm_CD8T_Lee)
fpr_CD8T_Lee, tpr_CD8T_Lee, _ = roc_curve(y_num_Lee, y_score_CD8T_Lee)
roc_auc_CD8T_Lee = auc(fpr_CD8T_Lee, tpr_CD8T_Lee)



# Plot ROC curve
plt.plot(fpr_CD8T_Xie, tpr_CD8T_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_CD8T_Xie:.2f}, p = {p_CD8T_Xie:.3})')
plt.plot(fpr_CD8T_Wang, tpr_CD8T_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_CD8T_Wang:.2f}, p = {p_CD8T_Wang:.3})')
plt.plot(fpr_CD8T_Li, tpr_CD8T_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_CD8T_Li:.2f}, p = {p_CD8T_Li:.3})')
plt.plot(fpr_CD8T_Lee, tpr_CD8T_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_CD8T_Lee:.2f}, p = {p_CD8T_Lee:.3})')

# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('CD8+T cells', fontweight='bold')
plt.legend(loc='lower right')
plt.show()






### CD8+TCR

plt.figure(figsize=(10, 8))

X_CD8TCR_Xie = X_Xie[['CD8+TCR']].values
X_sm_CD8TCR_Xie = sm.add_constant(X_CD8TCR_Xie)
logit_CD8TCR_Xie = sm.Logit(y_num_Xie, X_sm_CD8TCR_Xie)
results_CD8TCR_Xie = logit_CD8TCR_Xie.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8TCR_Xie = results_CD8TCR_Xie.pvalues[1]
y_score_CD8TCR_Xie = results_CD8TCR_Xie.predict(X_sm_CD8TCR_Xie)
fpr_CD8TCR_Xie, tpr_CD8TCR_Xie, _ = roc_curve(y_num_Xie, y_score_CD8TCR_Xie)
roc_auc_CD8TCR_Xie = auc(fpr_CD8TCR_Xie, tpr_CD8TCR_Xie)


X_CD8TCR_Wang = X_Wang[['CD8+TCR']].values
X_sm_CD8TCR_Wang = sm.add_constant(X_CD8TCR_Wang)
logit_CD8TCR_Wang = sm.Logit(y_num_Wang, X_sm_CD8TCR_Wang)
results_CD8TCR_Wang = logit_CD8TCR_Wang.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8TCR_Wang = results_CD8TCR_Wang.pvalues[1]
y_score_CD8TCR_Wang = results_CD8TCR_Wang.predict(X_sm_CD8TCR_Wang)
fpr_CD8TCR_Wang, tpr_CD8TCR_Wang, _ = roc_curve(y_num_Wang, y_score_CD8TCR_Wang)
roc_auc_CD8TCR_Wang = auc(fpr_CD8TCR_Wang, tpr_CD8TCR_Wang)


X_CD8TCR_Li = X_Li[['CD8+TCR']].values
X_sm_CD8TCR_Li = sm.add_constant(X_CD8TCR_Li)
logit_CD8TCR_Li = sm.Logit(y_num_Li, X_sm_CD8TCR_Li)
results_CD8TCR_Li = logit_CD8TCR_Li.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8TCR_Li = results_CD8TCR_Li.pvalues[1]
y_score_CD8TCR_Li = results_CD8TCR_Li.predict(X_sm_CD8TCR_Li)
fpr_CD8TCR_Li, tpr_CD8TCR_Li, _ = roc_curve(y_num_Li, y_score_CD8TCR_Li)
roc_auc_CD8TCR_Li = auc(fpr_CD8TCR_Li, tpr_CD8TCR_Li)



X_CD8TCR_Lee = X_Lee[['CD8+TCR']].values
X_sm_CD8TCR_Lee = sm.add_constant(X_CD8TCR_Lee)
logit_CD8TCR_Lee = sm.Logit(y_num_Lee, X_sm_CD8TCR_Lee)
results_CD8TCR_Lee = logit_CD8TCR_Lee.fit(method='bfgs', maxiter=1000, disp=False)
p_CD8TCR_Lee = results_CD8TCR_Lee.pvalues[1]
y_score_CD8TCR_Lee = results_CD8TCR_Lee.predict(X_sm_CD8TCR_Lee)
fpr_CD8TCR_Lee, tpr_CD8TCR_Lee, _ = roc_curve(y_num_Lee, y_score_CD8TCR_Lee)
roc_auc_CD8TCR_Lee = auc(fpr_CD8TCR_Lee, tpr_CD8TCR_Lee)



# Plot ROC curve
plt.plot(fpr_CD8TCR_Xie, tpr_CD8TCR_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_CD8TCR_Xie:.2f}, p = {p_CD8TCR_Xie:.3})')
plt.plot(fpr_CD8TCR_Wang, tpr_CD8TCR_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_CD8TCR_Wang:.2f}, p = {p_CD8TCR_Wang:.3})')
plt.plot(fpr_CD8TCR_Li, tpr_CD8TCR_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_CD8TCR_Li:.2f}, p = {p_CD8TCR_Li:.3})')
plt.plot(fpr_CD8TCR_Lee, tpr_CD8TCR_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_CD8TCR_Lee:.2f}, p = {p_CD8TCR_Lee:.3})')

# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('CD8+T cancer reactive cells', fontweight='bold')
plt.legend(loc='lower right')
plt.show()





### B cells

plt.figure(figsize=(10, 8))

X_B_Xie = X_Xie[['B']].values
X_sm_B_Xie = sm.add_constant(X_B_Xie)
logit_B_Xie = sm.Logit(y_num_Xie, X_sm_B_Xie)
results_B_Xie = logit_B_Xie.fit(method='bfgs', maxiter=1000, disp=False)
p_B_Xie = results_B_Xie.pvalues[1]
y_score_B_Xie = results_B_Xie.predict(X_sm_B_Xie)
fpr_B_Xie, tpr_B_Xie, _ = roc_curve(y_num_Xie, y_score_B_Xie)
roc_auc_B_Xie = auc(fpr_B_Xie, tpr_B_Xie)


X_B_Wang = X_Wang[['B']].values
X_sm_B_Wang = sm.add_constant(X_B_Wang)
logit_B_Wang = sm.Logit(y_num_Wang, X_sm_B_Wang)
results_B_Wang = logit_B_Wang.fit(method='bfgs', maxiter=1000, disp=False)
p_B_Wang = results_B_Wang.pvalues[1]
y_score_B_Wang = results_B_Wang.predict(X_sm_B_Wang)
fpr_B_Wang, tpr_B_Wang, _ = roc_curve(y_num_Wang, y_score_B_Wang)
roc_auc_B_Wang = auc(fpr_B_Wang, tpr_B_Wang)


X_B_Li = X_Li[['B']].values
X_sm_B_Li = sm.add_constant(X_B_Li)
logit_B_Li = sm.Logit(y_num_Li, X_sm_B_Li)
results_B_Li = logit_B_Li.fit(method='bfgs', maxiter=1000, disp=False)
p_B_Li = results_B_Li.pvalues[1]
y_score_B_Li = results_B_Li.predict(X_sm_B_Li)
fpr_B_Li, tpr_B_Li, _ = roc_curve(y_num_Li, y_score_B_Li)
roc_auc_B_Li = auc(fpr_B_Li, tpr_B_Li)



X_B_Lee = X_Lee[['B']].values
X_sm_B_Lee = sm.add_constant(X_B_Lee)
logit_B_Lee = sm.Logit(y_num_Lee, X_sm_B_Lee)
results_B_Lee = logit_B_Lee.fit(method='bfgs', maxiter=1000, disp=False)
p_B_Lee = results_B_Lee.pvalues[1]
y_score_B_Lee = results_B_Lee.predict(X_sm_B_Lee)
fpr_B_Lee, tpr_B_Lee, _ = roc_curve(y_num_Lee, y_score_B_Lee)
roc_auc_B_Lee = auc(fpr_B_Lee, tpr_B_Lee)



# Plot ROC curve
plt.plot(fpr_B_Xie, tpr_B_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_B_Xie:.2f}, p = {p_B_Xie:.3})')
plt.plot(fpr_B_Wang, tpr_B_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_B_Wang:.2f}, p = {p_B_Wang:.3})')
plt.plot(fpr_B_Li, tpr_B_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_B_Li:.2f}, p = {p_B_Li:.3})')
plt.plot(fpr_B_Lee, tpr_B_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_B_Lee:.2f}, p = {p_B_Lee:.3})')

# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('B cells', fontweight='bold')
plt.legend(loc='lower right')
plt.show()






### Monocytes

plt.figure(figsize=(10, 8))

X_Mono_Xie = X_Xie[['Mono']].values
X_sm_Mono_Xie = sm.add_constant(X_Mono_Xie)
logit_Mono_Xie = sm.Logit(y_num_Xie, X_sm_Mono_Xie)
results_Mono_Xie = logit_Mono_Xie.fit(method='bfgs', maxiter=1000, disp=False)
p_Mono_Xie = results_Mono_Xie.pvalues[1]
y_score_Mono_Xie = results_Mono_Xie.predict(X_sm_Mono_Xie)
fpr_Mono_Xie, tpr_Mono_Xie, _ = roc_curve(y_num_Xie, y_score_Mono_Xie)
roc_auc_Mono_Xie = auc(fpr_Mono_Xie, tpr_Mono_Xie)


X_Mono_Wang = X_Wang[['Mono']].values
X_sm_Mono_Wang = sm.add_constant(X_Mono_Wang)
logit_Mono_Wang = sm.Logit(y_num_Wang, X_sm_Mono_Wang)
results_Mono_Wang = logit_Mono_Wang.fit(method='bfgs', maxiter=1000, disp=False)
p_Mono_Wang = results_Mono_Wang.pvalues[1]
y_score_Mono_Wang = results_Mono_Wang.predict(X_sm_Mono_Wang)
fpr_Mono_Wang, tpr_Mono_Wang, _ = roc_curve(y_num_Wang, y_score_Mono_Wang)
roc_auc_Mono_Wang = auc(fpr_Mono_Wang, tpr_Mono_Wang)


X_Mono_Li = X_Li[['Mono']].values
X_sm_Mono_Li = sm.add_constant(X_Mono_Li)
logit_Mono_Li = sm.Logit(y_num_Li, X_sm_Mono_Li)
results_Mono_Li = logit_Mono_Li.fit(method='bfgs', maxiter=1000, disp=False)
p_Mono_Li = results_Mono_Li.pvalues[1]
y_score_Mono_Li = results_Mono_Li.predict(X_sm_Mono_Li)
fpr_Mono_Li, tpr_Mono_Li, _ = roc_curve(y_num_Li, y_score_Mono_Li)
roc_auc_Mono_Li = auc(fpr_Mono_Li, tpr_Mono_Li)



X_Mono_Lee = X_Lee[['Mono']].values
X_sm_Mono_Lee = sm.add_constant(X_Mono_Lee)
logit_Mono_Lee = sm.Logit(y_num_Lee, X_sm_Mono_Lee)
results_Mono_Lee = logit_Mono_Lee.fit(method='bfgs', maxiter=1000, disp=False)
p_Mono_Lee = results_Mono_Lee.pvalues[1]
y_score_Mono_Lee = results_Mono_Lee.predict(X_sm_Mono_Lee)
fpr_Mono_Lee, tpr_Mono_Lee, _ = roc_curve(y_num_Lee, y_score_Mono_Lee)
roc_auc_Mono_Lee = auc(fpr_Mono_Lee, tpr_Mono_Lee)



# Plot ROC curve
plt.plot(fpr_Mono_Xie, tpr_Mono_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_Mono_Xie:.2f}, p = {p_Mono_Xie:.3})')
plt.plot(fpr_Mono_Wang, tpr_Mono_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_Mono_Wang:.2f}, p = {p_Mono_Wang:.3})')
plt.plot(fpr_Mono_Li, tpr_Mono_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_Mono_Li:.2f}, p = {p_Mono_Li:.3})')
plt.plot(fpr_Mono_Lee, tpr_Mono_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_Mono_Lee:.2f}, p = {p_Mono_Lee:.3})')

# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('Monocytes', fontweight='bold')
plt.legend(loc='lower right')
plt.show()







### NK cells

plt.figure(figsize=(10, 8))

X_NK_Xie = X_Xie[['NK']].values
X_sm_NK_Xie = sm.add_constant(X_NK_Xie)
logit_NK_Xie = sm.Logit(y_num_Xie, X_sm_NK_Xie)
results_NK_Xie = logit_NK_Xie.fit(method='bfgs', maxiter=1000, disp=False)
p_NK_Xie = results_NK_Xie.pvalues[1]
y_score_NK_Xie = results_NK_Xie.predict(X_sm_NK_Xie)
fpr_NK_Xie, tpr_NK_Xie, _ = roc_curve(y_num_Xie, y_score_NK_Xie)
roc_auc_NK_Xie = auc(fpr_NK_Xie, tpr_NK_Xie)


X_NK_Wang = X_Wang[['NK']].values
X_sm_NK_Wang = sm.add_constant(X_NK_Wang)
logit_NK_Wang = sm.Logit(y_num_Wang, X_sm_NK_Wang)
results_NK_Wang = logit_NK_Wang.fit(method='bfgs', maxiter=1000, disp=False)
p_NK_Wang = results_NK_Wang.pvalues[1]
y_score_NK_Wang = results_NK_Wang.predict(X_sm_NK_Wang)
fpr_NK_Wang, tpr_NK_Wang, _ = roc_curve(y_num_Wang, y_score_NK_Wang)
roc_auc_NK_Wang = auc(fpr_NK_Wang, tpr_NK_Wang)


X_NK_Li = X_Li[['NK']].values
X_sm_NK_Li = sm.add_constant(X_NK_Li)
logit_NK_Li = sm.Logit(y_num_Li, X_sm_NK_Li)
results_NK_Li = logit_NK_Li.fit(method='bfgs', maxiter=1000, disp=False)
p_NK_Li = results_NK_Li.pvalues[1]
y_score_NK_Li = results_NK_Li.predict(X_sm_NK_Li)
fpr_NK_Li, tpr_NK_Li, _ = roc_curve(y_num_Li, y_score_NK_Li)
roc_auc_NK_Li = auc(fpr_NK_Li, tpr_NK_Li)



X_NK_Lee = X_Lee[['NK']].values
X_sm_NK_Lee = sm.add_constant(X_NK_Lee)
logit_NK_Lee = sm.Logit(y_num_Lee, X_sm_NK_Lee)
results_NK_Lee = logit_NK_Lee.fit(method='bfgs', maxiter=1000, disp=False)
p_NK_Lee = results_NK_Lee.pvalues[1]
y_score_NK_Lee = results_NK_Lee.predict(X_sm_NK_Lee)
fpr_NK_Lee, tpr_NK_Lee, _ = roc_curve(y_num_Lee, y_score_NK_Lee)
roc_auc_NK_Lee = auc(fpr_NK_Lee, tpr_NK_Lee)



# Plot ROC curve
plt.plot(fpr_NK_Xie, tpr_NK_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_NK_Xie:.2f}, p = {p_NK_Xie:.3})')
plt.plot(fpr_NK_Wang, tpr_NK_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_NK_Wang:.2f}, p = {p_NK_Wang:.3})')
plt.plot(fpr_NK_Li, tpr_NK_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_NK_Li:.2f}, p = {p_NK_Li:.3})')
plt.plot(fpr_NK_Lee, tpr_NK_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_NK_Lee:.2f}, p = {p_NK_Lee:.3})')

# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('NK cells', fontweight='bold')
plt.legend(loc='lower right')
plt.show()






### Combined without CD4+T and CD8+T


plt.figure(figsize=(10, 8))

X_Combined_Xie = X_Xie[["CD4+T", "CD4+TCR", "CD8+T", "CD8+TCR", "Mono", "NK"]].values
X_sm_Combined_Xie = sm.add_constant(X_Combined_Xie)
logit_Combined_Xie = sm.Logit(y_num_Xie, X_sm_Combined_Xie)
results_Combined_Xie = logit_Combined_Xie.fit(method='bfgs', maxiter=1000, disp=False)
y_score_Combined_Xie = results_Combined_Xie.predict(X_sm_Combined_Xie)
fpr_Combined_Xie, tpr_Combined_Xie, _ = roc_curve(y_num_Xie, y_score_Combined_Xie)
roc_auc_Combined_Xie = auc(fpr_Combined_Xie, tpr_Combined_Xie)


X_Combined_Wang = X_Wang[["CD4+T", "CD4+TCR", "CD8+T", "CD8+TCR", "Mono", "NK"]].values
X_sm_Combined_Wang = sm.add_constant(X_Combined_Wang)
logit_Combined_Wang = sm.Logit(y_num_Wang, X_sm_Combined_Wang)
results_Combined_Wang = logit_Combined_Wang.fit(method='bfgs', maxiter=1000, disp=False)
y_score_Combined_Wang = results_Combined_Wang.predict(X_sm_Combined_Wang)
fpr_Combined_Wang, tpr_Combined_Wang, _ = roc_curve(y_num_Wang, y_score_Combined_Wang)
roc_auc_Combined_Wang = auc(fpr_Combined_Wang, tpr_Combined_Wang)


X_Combined_Li = X_Li[["CD4+T", "CD4+TCR", "CD8+T", "CD8+TCR", "Mono", "NK"]].values
X_sm_Combined_Li = sm.add_constant(X_Combined_Li)
logit_Combined_Li = sm.Logit(y_num_Li, X_sm_Combined_Li)
results_Combined_Li = logit_Combined_Li.fit(method='bfgs', maxiter=1000, disp=False)
y_score_Combined_Li = results_Combined_Li.predict(X_sm_Combined_Li)
fpr_Combined_Li, tpr_Combined_Li, _ = roc_curve(y_num_Li, y_score_Combined_Li)
roc_auc_Combined_Li = auc(fpr_Combined_Li, tpr_Combined_Li)



X_Combined_Lee = X_Lee[["CD4+T", "CD4+TCR", "CD8+T", "CD8+TCR", "Mono", "NK"]].values
X_sm_Combined_Lee = sm.add_constant(X_Combined_Lee)
logit_Combined_Lee = sm.Logit(y_num_Lee, X_sm_Combined_Lee)
results_Combined_Lee = logit_Combined_Lee.fit(method='bfgs', maxiter=1000, disp=False)
y_score_Combined_Lee = results_Combined_Lee.predict(X_sm_Combined_Lee)
fpr_Combined_Lee, tpr_Combined_Lee, _ = roc_curve(y_num_Lee, y_score_Combined_Lee)
roc_auc_Combined_Lee = auc(fpr_Combined_Lee, tpr_Combined_Lee)



# Plot ROC curve
plt.plot(fpr_Combined_Xie, tpr_Combined_Xie, lw=2, color='royalblue', label=f'Xie et al., 2023 (CRC, N = 100) (AUC = {roc_auc_Combined_Xie:.2f})')
plt.plot(fpr_Combined_Wang, tpr_Combined_Wang, lw=2, color='darkorange', label=f'Wang et al., 2023 (BC, N = 80) (AUC = {roc_auc_Combined_Wang:.2f})')
plt.plot(fpr_Combined_Li, tpr_Combined_Li, lw=2, color='forestgreen', label=f'Li et al., 2024 (LUAD, N = 85) (AUC = {roc_auc_Combined_Li:.2f})')
plt.plot(fpr_Combined_Lee, tpr_Combined_Lee, lw=2, color='firebrick', label=f'Lee et al., 2024 (BC, N = 524) (AUC = {roc_auc_Combined_Lee:.2f})')


# Plot the diagonal line
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

# Set plot limits and labels
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontweight='bold')
plt.ylabel('True Positive Rate', fontweight='bold')
plt.title('Joint model without CD4+T and CD8+T', fontweight='bold')
plt.legend(loc='lower right')
plt.show()




########### Combined 4 studies and marginal model 


study_labels = ['Xie'] * X_Xie.shape[0] + ['Wang'] * X_Wang.shape[0] + ['Li'] * X_Li.shape[0] + ['Lee'] * X_Lee.shape[0]
study_df = pd.DataFrame(study_labels, columns=['Study'])
study_dummies = pd.get_dummies(study_df['Study'], prefix='Study')
study_dummies = study_dummies.drop(columns=['Study_Lee'])

X_Xie_marginal = X_Xie["CD4+TCR"]
X_Wang_marginal = X_Wang["CD4+TCR"]
X_Li_marginal = X_Li["CD4+TCR"]
X_Lee_marginal = X_Lee["CD4+TCR"]


# 将 Series 转换为 numpy 数组并确保是二维的
X_Xie_marginal = X_Xie_marginal.to_numpy().reshape(-1, 1) if X_Xie_marginal.ndim == 1 else X_Xie_marginal.to_numpy()
X_Wang_marginal = X_Wang_marginal.to_numpy().reshape(-1, 1) if X_Wang_marginal.ndim == 1 else X_Wang_marginal.to_numpy()
X_Li_marginal = X_Li_marginal.to_numpy().reshape(-1, 1) if X_Li_marginal.ndim == 1 else X_Li_marginal.to_numpy()
X_Lee_marginal = X_Lee_marginal.to_numpy().reshape(-1, 1) if X_Lee_marginal.ndim == 1 else X_Lee_marginal.to_numpy()

# 连接所有部分
X_combined = np.concatenate([X_Xie_marginal, X_Wang_marginal, X_Li_marginal, X_Lee_marginal], axis=0)

# 连接 study_dummies
X_combined = np.concatenate([X_combined, study_dummies.values], axis=1)


y_combined = np.concatenate([y_Xie, y_Wang, y_Li, y_Lee], axis = 0)
y_num_combined = (y_combined == 'Cancer').astype(int)   



X_sm_Combined = sm.add_constant(X_combined)
logit_Combined = sm.Logit(y_num_combined, X_sm_Combined)
results_Combined = logit_Combined.fit(method='bfgs', maxiter=1000, disp=False)
y_score_Combined = results_Combined.predict(X_sm_Combined)
fpr_Combined, tpr_Combined, _ = roc_curve(y_num_combined, y_score_Combined)
roc_auc_Combined = auc(fpr_Combined, tpr_Combined)



# 创建统计结果表格
feature_names_Combined = ['Intercept'] + [ "CD4+TCR", "Study_Xie", 'Study_Wang', 'Study_Li']
stats_df_Combined = pd.DataFrame({
    'Feature': feature_names_Combined,
    'Coefficient': results_Combined.params,
    'Std Error': results_Combined.bse,
    'P-value': results_Combined.pvalues,
    'OR': np.exp(results_Combined.params),
    'CI 2.5%': np.exp(results_Combined.conf_int()[:,0]),
    'CI 97.5%': np.exp(results_Combined.conf_int()[:,1])
})

# 格式化输出
stats_df_Combined['P-value'] = stats_df_Combined['P-value'].apply(lambda x: f"{x:.2e}")
stats_df_Combined['OR'] = stats_df_Combined['OR'].apply(lambda x: f"{x:.2f}")
stats_df_Combined['CI 2.5%'] = stats_df_Combined['CI 2.5%'].apply(lambda x: f"{x:.2f}")
stats_df_Combined['CI 97.5%'] = stats_df_Combined['CI 97.5%'].apply(lambda x: f"{x:.2f}")
stats_df_Combined['OR (95% CI)'] = stats_df_Combined.apply(lambda x: f"{x['OR']} ({x['CI 2.5%']}-{x['CI 97.5%']})", axis=1)


# 打印结果
print("\nLogistic Regression Statistics:")
print(stats_df_Combined[['Feature', 'Coefficient', 'P-value', 'OR (95% CI)']])




############# ML models 
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, f1_score, classification_report



study_labels = ['Xie'] * X_Xie.shape[0] + ['Wang'] * X_Wang.shape[0] + ['Li'] * X_Li.shape[0] 
study_df = pd.DataFrame(study_labels, columns=['Study'])
study_dummies = pd.get_dummies(study_df['Study'], prefix='Study')
study_dummies = study_dummies.drop(columns=['Study_Li'])


X_Xie_joint = X_Xie[["CD4+T", "CD8+T", "CD8+TCR", "Mono", "NK"]]
X_Wang_joint = X_Wang[["CD4+T", "CD8+T", "CD8+TCR", "Mono", "NK"]]
X_Li_joint = X_Li[["CD4+T", "CD8+T", "CD8+TCR", "Mono", "NK"]]
X_Lee_joint = X_Lee[["CD4+T", "CD8+T", "CD8+TCR", "Mono", "NK"]]

# 将 Series 转换为 numpy 数组并确保是二维的
X_Xie_joint = X_Xie_joint.to_numpy().reshape(-1, 1) if X_Xie_joint.ndim == 1 else X_Xie_joint.to_numpy()
X_Wang_joint = X_Wang_joint.to_numpy().reshape(-1, 1) if X_Wang_joint.ndim == 1 else X_Wang_joint.to_numpy()
X_Li_joint = X_Li_joint.to_numpy().reshape(-1, 1) if X_Li_joint.ndim == 1 else X_Li_joint.to_numpy()
X_Lee_joint = X_Lee_joint.to_numpy().reshape(-1, 1) if X_Lee_joint.ndim == 1 else X_Lee_joint.to_numpy()

# 连接所有部分
X_combined = np.concatenate([X_Xie_joint, X_Wang_joint, X_Li_joint], axis=0)

# 连接 study_dummies
X_combined = np.concatenate([X_combined, study_dummies.values], axis=1)


y_combined = np.concatenate([y_Xie, y_Wang, y_Li], axis = 0)
y_num_combined = (y_combined == 'Cancer').astype(int)   





# Combine Lee et al.
X_combined = np.concatenate([X_Lee_joint], axis=0)
y_combined = np.concatenate([y_Lee], axis = 0)
y_num_combined = (y_combined == 'Cancer').astype(int)   



### logistic regression


X_train, X_test, y_train, y_test = train_test_split(X_combined, y_num_combined, test_size=0.25, random_state=45)

logistic_model = LogisticRegression(max_iter=1000)
logistic_model.fit(X_train, y_train)
y_train_pred_logistic = logistic_model.predict(X_train)
y_test_pred_logistic = logistic_model.predict(X_test)


print("Logistic Regression Training Accuracy:", accuracy_score(y_train, y_train_pred_logistic))
print("Logistic Regression Training F1 Score:", f1_score(y_train, y_train_pred_logistic))
print("Logistic Regression Testing Accuracy:", accuracy_score(y_test, y_test_pred_logistic))
print("Logistic Regression Testing F1 Score:", f1_score(y_test, y_test_pred_logistic))



### Random Forest


rf_params = {'n_estimators': [50, 100, 200], 'max_depth': [None, 10, 20]}
rf_model = RandomForestClassifier(random_state=42)
rf_cv = GridSearchCV(rf_model, rf_params, cv=5)
rf_cv.fit(X_train, y_train)
y_train_pred_rf = rf_cv.predict(X_train)
y_test_pred_rf = rf_cv.predict(X_test)

# Random Forest Metrics
print("Random Forest Best Parameters:", rf_cv.best_params_)
print("Random Forest Training Accuracy:", accuracy_score(y_train, y_train_pred_rf))
print("Random Forest Training F1 Score:", f1_score(y_train, y_train_pred_rf))
print("Random Forest Testing Accuracy:", accuracy_score(y_test, y_test_pred_rf))
print("Random Forest Testing F1 Score:", f1_score(y_test, y_test_pred_rf))



### SVM


svm_params = {'C': [0.1, 1, 10], 'gamma': ['scale', 'auto']}
svm_model = SVC(kernel='rbf', random_state=42)
svm_cv = GridSearchCV(svm_model, svm_params, cv=5)
svm_cv.fit(X_train, y_train)
y_train_pred_svm = svm_cv.predict(X_train)
y_test_pred_svm = svm_cv.predict(X_test)

# SVM Metrics
print("SVM Best Parameters:", svm_cv.best_params_)
print("SVM Training Accuracy:", accuracy_score(y_train, y_train_pred_svm))
print("SVM Training F1 Score:", f1_score(y_train, y_train_pred_svm))
print("SVM Testing Accuracy:", accuracy_score(y_test, y_test_pred_svm))
print("SVM Testing F1 Score:", f1_score(y_test, y_test_pred_svm))