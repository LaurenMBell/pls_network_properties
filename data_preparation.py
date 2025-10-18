"""
THIS SCRIPT PERFORMS QC, FILTERING FOR CTRLs, AND NORMALIZING ALL THREE VARIATIONS
"""
import pandas as pd
import qnorm
import numpy as np


#ORIGINAL DATA
pls_og_data = pd.read_excel("data/CARLONI_ORIGINAL_DATA.xlsx", header=[0,1])

#======================================================================================================cd ..

#1) normalization -> median, quantile
#       (and turning nans into zeros)
def median_normalization():
    median = pd.read_csv("data/PLS_normalized.csv", header=[0,1])
    to_drop = [median.columns[1], "pls_76_13012022"]
    median.drop(columns=to_drop, inplace=True, errors='ignore')
    
    return median

def log_1(x):
    return np.log2(x + 1)

def quantile_normalization(data):
    numeric_cols = data.select_dtypes(include=[np.number]).columns
    non_numeric_cols = data.select_dtypes(exclude=[np.number])

    log_transformed = data[numeric_cols].apply(log_1)
    qnormed_values = qnorm.quantile_normalize(log_transformed, axis=1)
    qnormed = pd.DataFrame(qnormed_values, index=data.index, columns=numeric_cols)

    # Merge non-numeric columns back in (preserve metabolite names)
    qnormed_full = pd.concat([non_numeric_cols, qnormed], axis=1)

    # Reinsert NaNs where original zeros were
    qnormed_nans = qnormed_full.copy()
    qnormed_nans[log_transformed == 0] = np.nan
    qnormed_nans.to_csv("quantile/pls_converted_nans.csv", index=False, na_rep="NaN")
    return qnormed_nans

    
def median_without_zeros(median):
    median_wo = median.copy()
    median_wo[median == 0] = np.nan
    return median_wo

# NORMALIZE IN BOTH WAYS
quantile = quantile_normalization(pls_og_data)

median = median_normalization() #if median nomalization was actually implemented, pass it here
median_w = median
median_w.to_csv("median_w/median_w_preqc.csv", index=False)
median_wo = median_without_zeros(median)
median_wo.to_csv("median_wo/median_wo_preqc.csv", index=False)

#=====================================================================================================
#step 2: qc on each normalization variant
to_filter = pd.read_csv("data/Filtered Metabolomics - Plasma.csv",header=None)

def filter_metabolomics(data, qc, outfile):
    # First row of original data = feature_id
    #first_row = data.iloc[0, :]
    #print("Original # of columns:", len(first_row))
    
    # QC data - keep only rows marked "Keep"
    filtered_keep = qc[qc[3] == "Keep"]
    keep_ids = set(filtered_keep[0].astype(str))

    cols_to_keep = [col for col in data.columns if str(col[0]) in keep_ids]
    filtered = data.loc[:, cols_to_keep]
    filtered.to_csv(outfile, index=False)
    return filtered

filtered_quantile = filter_metabolomics(quantile, to_filter, "quantile/quantile_qc.csv")
filtered_median_w = filter_metabolomics(median_w, to_filter, "median_w/median_w_qc.csv")
filtered_median_wo = filter_metabolomics(median_wo, to_filter, "median_wo/median_wo_qc.csv")

#=========================================================================================================================================================

#STEP 3: FILTERNING OUT CONTROLS

def filter_by_group(data, outfile, label):
    label_row = data.columns.get_level_values(1)
    matching_cols = [col for col in data.columns if col[1] == label]

    ctrl_df = data.loc[:, matching_cols]
    ctrl_df.to_csv(outfile, index=False)
    print(f"[{label}] Saved {outfile} with {len(matching_cols)} data columns")
    return ctrl_df




filter_by_group(filtered_quantile, "quantile/quantile_ctrls.csv", "LPS_CTRL")
filter_by_group(filtered_median_w, "median_w/median_w_ctrls.csv", "VECPAC_CTRL")
filter_by_group(filtered_median_wo, "median_wo/median_wo_ctrls.csv", "DSS_T0")

#at this point this *should* be the clean data we want


