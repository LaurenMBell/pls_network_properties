import os
import pandas as pd
import itertools
import csv
from scipy import stats
import chime
import numpy as np

def metabolite(file_name, output_dir):
    #Read data
    filtered_data = pd.read_csv(file_name, header=None, na_values=["", "NaN"], keep_default_na=True)
    #print("test")
    
    # Extract metabolite names. Skip the first 2 rows (headers + Sample IDs)
    names = filtered_data.iloc[2:, 0].tolist()
        
    # Generate all 2-element combinations
    pairs = list(itertools.combinations(names, 2))
    results = []

    #temp = 0

  
    for met1, met2 in pairs:
        # Find the row where the first column = met1
        l1 = filtered_data[filtered_data.iloc[:, 0] == met1].iloc[:, 1:]
        #list1 = l1.iloc[:, 1:]
        #print(l1)
        
        # Find the row where the first column = met2
        l2 = filtered_data[filtered_data.iloc[:, 0] == met2].iloc[:, 1:]
        #print(l2)
        #list2 = l2.iloc[:, 1:]
                 
        # Converts to float
        row1 = l1.astype(float).values.flatten().tolist()
        row2 = l2.astype(float).values.flatten().tolist()
            
        #print(f" {met1} before: \n{row1}\n")
        #print(f"{met2} before:\n {row2}\n")
            
        #remove any comparisons with NaN values
        df = pd.DataFrame({'x': row1, 'y': row2}).dropna()
        row1, row2 = df['x'], df['y']
        
            #mask  = ~np.isnan(row1) & ~np.isnan(row2)
            #r1_masked = row1[mask]
            #r2_masked = row2[mask]
            
            
        if len(row1) <3:
            corr, pval = np.nan, np.nan
            results.append((met1, met2, corr, pval))
            #print("CORR and PVAL: np.nan\n")
        
        else:
            #print(f"{met1}: \n{row1}\n")
            #print(f"{met2}:\n {row2}\n")
                
            # Calculate Spearman correlation
            corr, pval = stats.spearmanr(row1, row2)
        
            #print(f"CORR AFTER: {corr}\n")
            #print(f"PVAL AFTER: {pval}\n")
                    
            # Append metabolite names + correlation + p-value
            results.append((met1, met2, corr, pval))
        #print(f"done with {met1} and {met2}!")
        
    # Create output file name inside "Metabolite Pairs"
    parts = os.path.basename(file_name).split("_")
    prefix = "_".join(parts[:2])
    
    output_file = os.path.join(output_dir, f"COMBINATIONS_{prefix}.csv")
    
    with open(output_file, "w", newline='', encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Metabolite 1", "Metabolite 2", "Spearman Coefficient", "p-value"])
        writer.writerows(results)

    print(f"Saved {output_file} with {len(results)} metabolite pairs.")

#metabolite("quantile/VECPAC_CTRL_quantile_ctrls.csv", "quantile/met_pairs")
#metabolite("quantile/LPS_CTRL_quantile_ctrls.csv", "quantile/met_pairs")
#metabolite("quantile/DSS_T0_quantile_ctrls.csv", "quantile/met_pairs")

#metabolite("median_w/VECPAC_CTRL_median_w.csv", "median_w/met_pairs")
#metabolite("median_w/LPS_CTRL_median_w.csv", "median_w/met_pairs")
metabolite("median_w/DSS_T0_median_w.csv", "median_w/met_pairs")


#metabolite("median_wo/VECPAC_CTRL_median_wo.csv", "median_wo/met_pairs")
#metabolite("median_wo/LPS_CTRL_median_wo.csv", "median_wo/met_pairs")
#metabolite("median_wo/DSS_T0_median_wo.csv", "median_wo/met_pairs")

