import utils_ML as uml

import pandas as pd
import numpy as np

from lightgbm import LGBMClassifier as lgbm
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler, normalize
from sklearn.utils.class_weight import compute_class_weight
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.utils.fixes import loguniform

from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import SMOTE
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions

import vae_utils as vae_utils
device = 'cuda' if torch.cuda.is_available() else "cpu"

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from lightgbm import LGBMClassifier as lgbm


import utils_ML as uml

import pandas as pd
import numpy as np

from lightgbm import LGBMClassifier as lgbm
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler, normalize
from sklearn.utils.class_weight import compute_class_weight
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.utils.fixes import loguniform

from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import SMOTE
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions

import vae_utils as vae_utils
device = 'cuda' if torch.cuda.is_available() else "cpu"

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from lightgbm import LGBMClassifier as lgbm


def main():
    # Machine learning classification loop

    # read data
    skf = StratifiedKFold(n_splits=10, shuffle=True)

    print("Reading data...")
    data_quantile = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/quantile_norm_NSAF_50.csv", index_col = "assay_id")
    protein_columns = data_quantile.columns
    meta = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/unified_metadata.csv")
    meta = meta[meta.assay_id.isin(data_quantile.index)]
    groups = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/group_cells_annotation.csv", sep =";", index_col="Unnamed: 0")
    meta["Group"] = meta.cell_line.apply(lambda x: groups[groups.cell_line == x]["group"].values[0])
    meta = meta.set_index("assay_id")

    data_quantile.sort_index(inplace=True)
    meta.sort_index(inplace=True)

    with open("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/selected_features.txt", "r") as f:
        features = f.readlines()
        features = [x.strip() for x in features]
    data_quantile = data_quantile.loc[:, features]

    data_quantile = data_quantile.reset_index(drop=True).rename(columns={data_quantile.columns[x]:x for x in range(len(data_quantile.columns))})


    target_encoder = LabelEncoder()
    targets = target_encoder.fit_transform(meta.Group)
    unique_labels = pd.Series(targets).unique()
    class_weights = compute_class_weight(class_weight='balanced', classes=unique_labels, y=targets)

    weights = {unique_labels[i]: class_weights[i] for i in range(len(unique_labels))}
    print(weights)

    vae = vae_utils.VariationalAutoencoder(6, 161, 50)
    vae.load_state_dict(torch.load("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/VAE_model_filtered"))

    lr_clf = LogisticRegression(max_iter=10000)
    svm_clf = SVC()
    rf_clf = RandomForestClassifier()
    lgbm_clf = lgbm()
    models = [lr_clf, svm_clf, rf_clf, lgbm_clf]
    vae = vae_utils.VariationalAutoencoder(6, 161, 50)
    vae.load_state_dict(torch.load("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/VAE_model_filtered"))
    
    # Models to be used:
    print("Initializing models")
    svc_clf = SVC(C=15,kernel="rbf")
    rf_clf = RandomForestClassifier(n_estimators=150, criterion="entropy", max_depth=10)
    lr_clf = LogisticRegression(max_iter = 10000, C=50, solver="liblinear", penalty="l2")
    models = [svc_clf, rf_clf, lr_clf]

    skf = StratifiedKFold(n_splits=5, shuffle=True)
    pbkf = uml.ProjectBasedSplit(splits=None, metadata=meta, LOOCV=True, on="Group")
    fold=0

    print("Starting outer loop...")
    print("Data shape:", data_quantile.shape)

    for train, test in pbkf.split(dataset=data_quantile, metadata=meta, n_projects=35):
        fold += 1
        if fold < 5:
            continue
        print(fold)
        # Split data
        X_train_quantile = data_quantile.iloc[train,:].reset_index(drop=True)
        X_test_quantile = data_quantile.iloc[test,:].reset_index(drop=True)

        Y_train = targets[train]
        Y_test = targets[test]

        # Preprocess
        print("Imputing...")
        #imputer = uml.MNAR_MCAR_Imputer(max_iter=15)
        imputer = uml.LowestValueImputerGaussian()
        imputer.fit(X_train_quantile, Y_train)
        imputed_train = imputer.transform(X_train_quantile, Y_train)
        imputed_test = imputer.transform(X_test_quantile)

        quant_scaler = MinMaxScaler()
        scaled_train = quant_scaler.fit_transform(imputed_train)
        scaled_test = quant_scaler.transform(imputed_test)
        scaled_train = scaled_train.astype("float32")
        scaled_test = scaled_test.astype("float32")

        # Oversampling
        print("Oversampling...")
        smote = SMOTE()
        try:
            smote_quant, smote_y = smote.fit_resample(scaled_train, Y_train)
        except:
            smote = SMOTE(k_neighbors=3)
            print(pd.Series(Y_train).value_counts().min())
            smote_quant, smote_y = smote.fit_resample(scaled_train, Y_train)
            
        
        vae_quant, vae_y = vae_utils.resampleVAE(vae, scaled_train, Y_train, 10)

        print(fold, ': Oversampled')

        for model in models:

            model.fit(smote_quant, smote_y)
            y_pred_smote = model.predict(scaled_test)

            model.fit(vae_quant, vae_y)
            y_pred_vae = model.predict(scaled_test)

            # Compute scores and save
            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=y_pred_smote, 
                                                                        Y_test=Y_test,
                                                                        labels=unique_labels)

            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], 
                                            "micro_f1": [micro_f1], 
                                            "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm], 
                                            "oversampler": ["SMOTE"]})
                
            uml.save_results(results_df, "project_split_eval") 
            print("SMOTE: {:.2f} for {} samples".format(micro_f1, len(y_pred_vae)))

            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=y_pred_vae, 
                                                                        Y_test=Y_test,
                                                                        labels=unique_labels)

            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], 
                                            "micro_f1": [micro_f1], 
                                            "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm], 
                                            "oversampler": ["VAE"]})
                
            uml.save_results(results_df, "project_split_eval") 
            print("VAE: {:.2f} for {} samples".format(micro_f1, len(y_pred_vae)))

            print("Saved results")

if __name__ == '__main__':
    main()