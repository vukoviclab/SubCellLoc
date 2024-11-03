# SubCellLoc

To retrieve datasets from Pharos, ChEMBL (Release 32), and DrugCentral, run the `pharos_ex.py` script to fetch the data.

## Notebook Descriptions

1. **sublocation.ipynb**  
   This notebook reads CSV files and converts UniProt IDs to sublocations.

2. **prepare_labels**  
   This notebook introduces labels as presented in **S1** of the paper. It reads a CSV file and sets up 33 main sublocation columns, with around 400 sublocations as values under each column. This setup simplifies mapping the locations for data fetched from ChEMBL, Pharos, and DrugCentral.

3. **making_final_dataset_and_save_with_labels.ipynb**  
   This notebook creates directories to save CSV files with location labels mapped to each sublocation.

4. **deepchem.ipynb**  
   This notebook is used to train models with the DeepChem library.

5. **training_with_diff_models**  
   This notebook trains models using various scikit-learn algorithms.

