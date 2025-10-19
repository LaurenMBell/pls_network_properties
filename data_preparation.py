"""
THIS SCRIPT PERFORMS QC, FILTERING FOR CTRLs, AND NORMALIZING ALL THREE VARIATIONS
"""
import pandas as pd
import qnorm
import numpy as np


#ORIGINAL DATA
pls_og_data = pd.read_excel("data/CARLONI_ORIGINAL_DATA.xlsx", header=[0,1])

def fix_names(s):
    s = str(s).strip()
    s=s.replace(" ", "_")
    return s

pls_og_data.columns = pd.MultiIndex.from_tuples(
    [(fix_names(a), fix_names(b)) for a, b in pls_og_data.columns]
)

#======================================================================================================cd ..

#1) normalization -> median, quantile
#       (and turning nans into zeros)

def median_without_zeros(median, og_data):
    median_wo = median.copy()
    median_wo[og_data == 0] = np.nan
    return median_wo

def median_normalization():
    median = pd.read_csv("data/PLS_normalized.csv", header=[0,1])
    to_drop = [median.columns[1], "pls_76_13012022"]
    median.drop(columns=to_drop, inplace=True, errors='ignore')
    median.to_csv("median_w/median_w_preqc.csv", index=False)

    median_wo = median_without_zeros(median, pls_og_data)
    median_wo.to_csv("median_wo/median_wo_preqc.csv", index=False)
    
def log_1(x):
    return np.log2(x + 1)

def quantile_normalization():
    og_data = pls_og_data
    
    names = og_data.iloc[:,0]

    data_df = og_data.iloc[0:,1:]
    log_transformed = data_df.apply(log_1)

    named = log_transformed.copy()
    named.insert(0, "Names", names)
    named.to_csv(f"quantile/quantile_log_transformed.csv", index=False)

    qnormed = qnorm.quantile_normalize(log_transformed, axis=1)
    qnormed.insert(0, "", names)
    qnormed.to_csv(f"quantile/quantile_qnormed.csv", index=False)

    qnormed_nans = qnormed.copy()
    qnormed_nans[log_transformed==0] = np.nan
    qnormed_nans.to_csv(f"quantile/quantile_preqc.csv", index=False, na_rep="NaN")



def fix_quantile_format():
    df = pd.read_csv("quantile/pls_quantile_preqc.csv", header=None)
    df.iat[0, 0] = "Label"
    df.to_csv("quantile/pls_quantile_preqc.csv",index=False)

# NORMALIZE IN BOTH WAYS
#quantile_normalization()
#fix_quantile_format() #FOR NOW, I HAVE GIVEN UP AND FIXED THE FILES IN EXCEL
#median_normalization() 

#=====================================================================================================
#step 2: qc on each normalization variant
to_filter = pd.read_csv("data/Filtered Metabolomics - Plasma.csv",header=None)

def filter_metabolomics(data, filtered_qc, out_file): #, label=""):
    data = pd.read_csv(data, header=None)
    first_row = data.iloc[0, :]
    
    filtered_keep = filtered_qc[filtered_qc[3] == "Keep"]
    
    class_to_replace = dict(zip(filtered_keep[0], filtered_keep[5]))

    columns_to_keep = []
    for i, val in enumerate(first_row.values):
        if val in filtered_keep[0].values:
            columns_to_keep.append(i)

    # Keeps row 0 (samplie IDs)
    columns_to_keep = sorted(set([0] + columns_to_keep))
    # Copy Data
    data_filtered = data.iloc[:, columns_to_keep].copy()

    # Replace second row with new class if available
    replacements = 0
    for sid in data_filtered.columns:
        # ID in row 0
        feature_name = data_filtered.loc[0, sid] 
        # current value in row 1
        old_val = data_filtered.loc[1, sid] 
        new_val = class_to_replace.get(feature_name, old_val)
        if new_val != old_val:
            replacements += 1
        # update row 1
        data_filtered.loc[1, sid] = new_val 

    # Save result
    data_filtered.to_csv(out_file, index=False, header=False)
    return data_filtered

filter_metabolomics("quantile/pls_quantile_preqc.csv", to_filter, "quantile/quantile_qc.csv")
filter_metabolomics("median_w/median_w_preqc.csv", to_filter, "median_w/median_w_qc.csv")
filter_metabolomics("median_wo/median_wo_preqc.csv", to_filter, "median_wo/median_wo_qc.csv")

#=========================================================================================================================================================

#STEP 3: FILTERNING OUT CONTROLS

def filter_by_group(input_file, outdir, outfile, label):
    df = pd.read_csv(input_file, header=None)
    second_row = df.iloc[1]
    
    # Find columns that match the control group
    columns_to_keep = [i for i, val in enumerate(second_row) if val == label]
    
    # Always keep first column (metabolite/sample IDs)
    columns_to_keep = sorted(set([0] + columns_to_keep))
    
    # Save filtered dataset
    filtered_df = df.iloc[:, columns_to_keep]

    filtered_df.to_csv(f"{outdir}/{label}_{outfile}", index=False, header=False)
    print(f"[{label}] Saved {outfile} with {len(columns_to_keep)} columns")



#lauren fix this tomorow when youre less tired
filter_by_group("quantile/quantile_qc.csv", "quantile", "quantile.csv", "VECPAC_CTRL")
filter_by_group("quantile/quantile_qc.csv", "quantile", "quantile.csv", "LPS_CTRL")
filter_by_group("quantile/quantile_qc.csv", "quantile", "quantile.csv", "DSS_T0")

filter_by_group("median_w/median_w_qc.csv", "median_w", "median_w.csv", "VECPAC_CTRL")
filter_by_group("median_w/median_w_qc.csv", "median_w", "median_w.csv", "LPS_CTRL")
filter_by_group("median_w/median_w_qc.csv", "median_w", "median_w.csv", "DSS_T0")

filter_by_group("median_wo/median_wo_qc.csv", "median_wo", "median_wo.csv", "VECPAC_CTRL")
filter_by_group("median_wo/median_wo_qc.csv", "median_wo", "median_wo.csv", "LPS_CTRL")
filter_by_group("median_wo/median_wo_qc.csv", "median_wo", "median_wo.csv", "DSS_T0")

#at this point this *should* be the clean data we want


