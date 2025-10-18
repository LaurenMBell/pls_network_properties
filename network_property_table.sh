echo "Preparing data..."
python data_preparation.py
echo "done!"


echo "Calculating metabolite correlations..."
python metabolite_correlations.py
echo "done!"

echo "making the FDR and meta-analysis table..."
python fdr_table.py
echo "done!"

echo "making the network properties table..."
python network_properties_table.py
echo "done with everything!"