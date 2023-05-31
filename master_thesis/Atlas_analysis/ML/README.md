# Description of files and notebooks

## Python script

* **utils_ML.py** <br>
In this folder multiple utility functions that are used during data analysis are implemented. These include the following:
    - filtering functions and classes
    - Custom data splitter class
    - Imputer classes
    - Feature selector class
    - scoring functions

* **grid_search50.py** <br>
Script to find the best hyperparameters of the models

* **norm_impute_pbkf.py** and **norm_impute.py**<br>
Scripts evaluating the most optimal normalisation-imputation combination with model performance. Either leave-one-project-out cross validation is performed (*norm_impute_pbkf*) or stratified k-fold split (*norm_impute*)

*  **opt_feature_selection.py** <br>
Script used to evaluate the most optimal feature selector. This is based on classification performance when only using the 100 selected features defined by each selector.

* **oversampler_comparison.py** <br>
Script evaluating classification performance when first oversampling on the fPEMatrix and then only retaining the selected features. In the *Oversampling* notebook, the script for the reverse (feature selection --> oversampling) is written. The reason why this difference is important to evaluate is that oversamplers can leverage the complete dataset in the former method in contrast to the latter, where only the selected features guide oversampling.

* **project_split_eval**<br>
Final model classification performance estimation with the leave-one-project-out cross validation for the optimised model and pipeline

* **vae_filtered_hyperparameter_opt.py**
Script used to perform parameter sweep for the VAE with feature selected dataset

* **vae_hyperparameter_opt.py**
Script used to perform parameter swee for the VAE on the fPEMatrix

* **vae_utils.py** <br>
Contains classes implementing the Variational Autoencoder network and functions to train, use and oversample from it.

## Notebooks

* **corr_prots.ipynb**<br>
Contains two parts: 
    - Evaluation of the predictive value of high pairwise Pearson correlations for biological associations as annotated by STRING
    - Interpretation of feature importances as described in *feature_importance_analysis.ipynb* in context of the complete dataset with correlation clustering

* **evaluate_unseen_data.ipynb** <br>
Here, a part of the left out data (due to having less than 10 samples in a class) are reused here to evaluate the final model. Additionally, some predictions are interpreted with SHAP.

* **feature_importance_analysis.ipynb** <br>
Notebook describing the SHAP-based feature analysis of the trained model

* **find_tissue_profile.ipynb** <br>
Notebook written by Tine Claeys, exploring the tissue specificity of the most important features for cell line prediction


## Files

* **all_feature_importances.csv** <br>
csv-file with all feature importances of the optimised Logistic Regression model. Feature importance estimation was performed with SHAP.

* **Most_important_feature.csv** <br>
Subset of *all_feature_importances.csv*.

* **GO_\*_uniprot_map.json** <br>
Mapped GO-terms with each uniprot identifier

* **string_interaction_scores_all.csv** <br>
STRING-score file for all pairwise proteins. Downloaded from "https://string-db.org/cgi/download.pl (version 11.5)"

* **ML** <br>
Contains notebooks and files related to machine learning modelling of the data. Also contains the folders with machine learning evaluation metrics and SHAP-values after model optimalisation

## Folders

* **shap_values**<br>
Contains the files with thecomputed SHAP-values

* **results** <br>
Folder with result files from the experiments

* **bin** <br>
Abandoned notebooks and scripts