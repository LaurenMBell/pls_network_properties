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
    #you could put manual median normalization instead of just reading another file
    median = pd.read_csv("data/PLS_normalized.csv") 
    to_drop = [median.columns[1], "pls_76_13012022"]
    median.drop(columns=to_drop, inplace=True, errors='ignore')
    return median

def log_1(x):
    return np.log2(x + 1)

def quantile_normalization(data): 
    to_drop = [data.columns[1], "pls_76_13012022"]
    data.drop(columns=to_drop, inplace=True, errors='ignore')

    names = data.index if data.index.name else data.iloc[:, 0]
    numeric = data._get_numeric_data()
    log_transformed = numeric.map(log_1)

    named = log_transformed.copy()
    named.insert(0, "Names", names)
    named.to_csv(f"quantile/pls_log_transformed.csv", index=False)

    qnormed = qnorm.quantile_normalize(log_transformed, axis=1)
    qnormed.insert(0, "", names)
    qnormed.to_csv(f"quantile/pls_qnormed.csv", index=False)

    qnormed_nans = qnormed.copy()
    qnormed_nans[log_transformed==0] = np.nan
    qnormed_nans.to_csv(f"quantile/pls_converted_nans.csv", index=False, na_rep="NaN")
    return qnormed_nans
    
def median_without_zeros(median):
    median_wo = median

    median_wo[median==0] == np.nan

    return median_wo

# NORMALIZE IN BOTH WAYS
quantile = quantile_normalization(pls_og_data)

median = median_normalization() #if median nomalization was actually implemented, pass it here
median_w = median
median_wo = median_without_zeros(median)

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
    data_filtered = data.loc[:, [col for col in data.columns if col[0] in keep_ids]]

    data_filtered.to_csv(outfile, index=False)
    return data_filtered

filtered_quantile = filter_metabolomics(quantile, to_filter, "quantile/quantile_qc.csv")
filtered_median_w = filter_metabolomics(median_w, to_filter, "median_w/median_w_qc.csv")
filtered_median_wo = filter_metabolomics(median_wo, to_filter, "median_wo/median_wo_qc.csv")

#=========================================================================================================================================================

#STEP 3: FILTERNING OUT CONTROLS

def filter_by_group(data, outfile, label):
    ctrl_df = data.loc[:, [col for col in data.columns if col[1] == label]]
    ctrl_df.to_csv(outfile, index=False)
    print(f"[{label}] Saved {outfile} with {ctrl_df.shape[1]} columns")


filter_by_group(filtered_quantile, "quantile/quantile_ctrls.csv", "LPS_CTRL")
filter_by_group(filtered_median_w, "median_w/median_w_ctrls.csv", "VECPAC_CTRL")
filter_by_group(filtered_median_wo, "median_wo/median_wo_ctrls.csv", "DSS_T0")

#at this point this *should* be the clean data we want


