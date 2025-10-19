import pandas as pd
import numpy as np

# FUNCS TO CALC NETWORK CHARACTERISTICS FOR EACH PARAM ####################

def filter_thresholds(ind, fdr, df):
    #takes total fdr table, returns filtered one

    filtered_df = df[df["FDR"] <= fdr]
    filtered_df = filtered_df[(filtered_df["VECPAC p-values"] <= ind) 
                              & (filtered_df["DSS p-values"] <= ind) 
                              & (filtered_df["LPS p-values"] <= ind)]

    return filtered_df #new dataframe

def nodes(df):
    # takes filtered df, counts uniqe metabolites, then finds + and - mets
    
    #to be changed later
    pos_nodes = np.nan
    neg_nodes = np.nan

    unique_mets = set() 

    unique_mets = set(df["Metabolite 1"]).union(set(df["Metabolite 2"]))
    total_nodes = len(unique_mets)
    
    # question: how to calculate positive and negative nodes, fix this lauren
    return total_nodes, pos_nodes, neg_nodes 

def edges(df):
    # takes filtered df, counts number of corr coefs
    #finds + and - correlations and the ratio btwn the two

    total_edges = len(df[df["Meta-analysis_validity"] == True])

    to_check = ["VECPAC r", "LPS r", "DSS r"]
    """
    if there are three values: 
        the sign of the edge will be the same as any of the CCs
    if there are two values: 
        make a new df and drop the nan value

        the sign of the edge will be the same as either of those CCs
    """ 

    if len(df[to_check] == 3):
        pos_edges = len(df[df["DSS r"] > 0])
        neg_edges = len(df[df["DSS r"] < 0])
    else:
        pos_edges = ((np.sign(df[to_check]).sum(axis=1)) == 2)
        pos_edges = ((np.sign(df[to_check]).sum(axis=1)) == -2)


    #pos_edges = len(df[df["edge_dir"] == 1])
    #neg_edges = len(df[df["edge_dir"] == -1])

    if neg_edges > 0:
        ratio = pos_edges/neg_edges 
    else:
        ratio = np.nan

    return total_edges, pos_edges, neg_edges, ratio

def mean_corr(df):
    # Only take correlation columns
    r_cols = ["VECPAC r", "DSS r", "LPS r"]

    edge_means = df[r_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    overall_mean = edge_means.mean()

    return overall_mean


# constructs the final table for each file !!
def table_making(file, name):
    final_table = pd.DataFrame()

    full_graph = file.dropna()
    #number of edges in the full graph = the number of consistent corr coeffs
    full_graph_edges = len(full_graph)

    for i in ind_p:
        for j in fdr_p:
            return_table = pd.DataFrame()
            
            for col in ["VECPAC r", "VECPAC p-values", "DSS r", "DSS p-values", "LPS r", "LPS p-values", "FDR"]:
                file[col] = pd.to_numeric(file[col], errors="coerce")
            file.columns = file.columns.str.strip()

        
            filtered = filter_thresholds(i, j, full_graph)


            total_nodes, pos_nodes, neg_nodes = nodes(filtered)
            total_edges, pos_edges, neg_edges, ratio = edges(filtered)
            full_edges = (total_nodes**2 - total_nodes)/2
            density = total_edges/full_edges
            mean_corr_coeff = mean_corr(filtered)

            return_table = pd.DataFrame([{
                "Individual_Pval": i,
                "FDR": j,
                "Total_Nodes": total_nodes,
                #"Positive_Nodes": pos_nodes,
                #"Negative_Nodes": neg_nodes,
                "Total_Edges": total_edges,
                "Positive_Edges": pos_edges,
                "Negative_Edges": neg_edges,
                "Pos_to_Neg": ratio,
                "Density": density,
                "Mean_Corr_Coeff": mean_corr_coeff
            }])

            final_table = pd.concat([final_table, return_table])

    final_table.to_csv(f"{name}_network_properties.csv", index=False)

    return final_table


# TABLE MAKING #########################################################################

#df_quant = pd.read_csv("")
df_median = pd.read_csv("median_w/fdr_table.csv")
#df_median_without_zeros = pd.read_csv("/Users/laurenbell/Desktop/Metabolite_Brain_Gut_MorgunLab/network_properties/median_norm_drop_zero_data/updated_metaanalysis_fdr_table.csv")


#thresholds to iterate over
ind_p = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1]
fdr_p = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1]

#table_making(df_quant, "quantile")
table_making(df_median, "median_w/median_w")
#table_making(df_median_without_zeros, "median_without_zeros")
