UKB CLI tutorial: https://community.ukbiobank.ac.uk/hc/en-gb/articles/26962205402525-Introduction-to-Command-Line-Interface-CLI


All files have been run with swiss army knife in UKB-RAP

Swiss army knife requires all files used in the script, as well as the script file itself to be added as inputs.




00_Preprocessing

Concatenate Proteomics.R
Used for: concatenating 6 files with proteomics data, when table exporter failed to export everything at once
Input parameters: None
Input files: prot1.csv, prot2.csv, prot3.csv, prot4.csv, prot5.csv, prot6.csv

Create Cohorts.R
Used for: generating training and test sets, ran locally
Input parameters: None
Input files: ids_95k_female.txt (female image IDs),
ids_95k_male.txt (male image IDs),
double_proteomics.csv (participants IDs with proteomics both instance 0 and 2),
proteomics_b0-6_eids.csv (participants IDs with proteomics instance 0, from batch 0-6),
metabolomics_ids.csv (participants IDs with metabolomics instance 0)

change col_names.R
Used for: changing column names in age_data
Input parameters: None
Input files: age_data.csv

adjust_age.R
Used for: adding difference between age_0 and age_2 to the age data
Input parameters: None
Input files: age_data.csv


01_Statistics

calc_adj_corr.R
Used for: calculating correlations with both pearson and spearman, and both with and without height adjustment
Input parameters: input_omic_path, input_image_path, output_name, cohort_eids_path, gender
Input files: age_data.csv, height_data.csv, omics and image data files, eid_file

generate_heatmap.R
Used for: generating heatmap of correlations, run locally
Input parameters: None
Input files: One correlation file

met_vs_img.R/prot_vs_img.R
Used for: performing the normalization, residualization and PLS
Input parameters:
    input_metabolomics_path/input_proteomics_path,
    input_image_path,
    train_eids_path,
    test_eids_path,
    img_type_string (not used anymore),
    interested_feature (name of image feature),
    height_included (boolean)
Input files:
    Omics file,
    image file,
    train eid file,
    test eid file,
    age_data,
    height_data,
    encoding file,
Outputs:
    xx_train_vs_test: Plot of heldout vs traing R2 scores
    xx_train: Plot of train score
    xx_heldout: Plot of heldout score
    xx_results_2: Results from CV folds
    xx_train_added_value: Plot of added value per PLS component
    xx_heldout_added_value: Plot of added value per PLS component
    xx_R2_results: Results from built_in tuning function (the ones used to create added_value plots)
    xx_normalizers/xx_residualizers: Model values for the normalization and residualization (for application on test set)
    xx_results: Final model on residualized data
    xx_sample_spread: Sample spread across first two components
    xx_Variance: Variance in the omics data across first two components
    xx_Loadings: Top 10 omics loadings for first component

02_Evaluation

test_met.R/test_prot.R
Used for: Evaluating model on test set
Input parameters:
    input_metabolomics_path/input_proteomics_path,
    input_image_path,
    train_eids_path,
    test_eids_path,
    img_type_string (not used anymore),
    interested_feature (name of image feature),
    height_included (boolean),
    ncomp (number of components to use)
Input files:
    Omics file,
    image file,
    train eid file,
    test eid file,
    age_data,
    height_data,
    encoding file,
    xx_results,
    xx_image_normalizer.rds,
    xx_image_residualizer.rds,
    xx_metabolite_normalizers.rds/xx_protein_normalizers.rds,
    xx_metabolite_normalizers.rds/xx_protein_normalizers.rds,
    (NB xx= image feature used)
Outputs:
    xx_train_vs_test: Plot of heldout vs traing R2 scores
    xx_train: Plot of train score
    xx_heldout: Plot of heldout score
    xx_results_2: Results from CV folds
    xx_train_added_value: Plot of added value per PLS component
    xx_heldout_added_value: Plot of added value per PLS component
    xx_R2_results: Results from built_in tuning function (the ones used to create added_value plots)
    xx_normalizers/xx_residualizers: Model values for the normalization and residualization (for application on test set)
    xx_results: Final model on residualized data
    xx_sample_spread: Sample spread across first two components
    xx_Variance: Variance in the omics data across first two components
    xx_Loadings: Top 10 omics loadings for first component
