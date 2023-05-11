from lightgbm import LGBMClassifier as lgbm
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

# preprocessors
from sklearn.preprocessing import LabelEncoder, MinMaxScaler

# Samplers
from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import SMOTE

# metrics and splitters
from sklearn.model_selection import StratifiedKFold

# utils
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.class_weight import compute_class_weight
import pandas as pd
import sys
import matplotlib.cm as cm


import vae_utils as vae_utils

sys.path.append("../")
import utils_ML as uml

import torch
import torch.utils
import torch.distributions


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

#with open("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/selected_features.txt", "r") as f:
#    features = f.readlines()
#    features = [x.strip() for x in features]
#data_quantile = data_quantile.loc[:, features]

data_quantile = data_quantile.reset_index(drop=True).rename(columns={data_quantile.columns[x]:x for x in range(len(data_quantile.columns))})

target_encoder = LabelEncoder()
targets = target_encoder.fit_transform(meta.Group)
unique_labels = pd.Series(targets).unique()
class_weights = compute_class_weight(class_weight='balanced', classes=unique_labels, y=targets)

weights = {unique_labels[i]: class_weights[i] for i in range(len(unique_labels))}
print(weights)

vae = vae_utils.VariationalAutoencoder(10,2615,500)
vae.load_state_dict(torch.load("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/VAE_model"))

lr_clf = LogisticRegression(max_iter=10000)
svm_clf = SVC()
rf_clf = RandomForestClassifier()
lgbm_clf = lgbm()
models = [lr_clf, svm_clf, rf_clf, lgbm_clf]

param = {"class_weight": weights}
no_param = {"class_weight": None}

fs = uml.FeatureSelector("anova MI LR SVC RF".split(), 300, .5)

fold = 0
for train, test in skf.split(X=data_quantile, y=targets):
    fold+=1
    print(fold,"/10")
    X_train_quantile = data_quantile.iloc[train,:].reset_index(drop=True)
    X_test_quantile = data_quantile.iloc[test,:].reset_index(drop=True)

    Y_train = targets[train]
    Y_test = targets[test]

    # Imputation
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
    smote_quant, smote_y = smote.fit_resample(scaled_train, Y_train)

    smotetomek = SMOTETomek()
    SMOTETomek_quant, SMOTETomek_y = smotetomek.fit_resample(scaled_train, Y_train)

    smoteenn = SMOTEENN()
    SMOTEENN_quant, SMOTEENN_y =smoteenn.fit_resample(scaled_train, Y_train)

    vae_quant, vae_y = vae_utils.resampleVAE(vae, scaled_train, Y_train, 1)

    annotation = ["weights", "SMOTE", "SMOTETomek", "SMOTEENN", "VAE"]

    X_ = [scaled_train, smote_quant, SMOTETomek_quant, SMOTEENN_quant, vae_quant]
    y = [Y_train, smote_y, SMOTETomek_y, SMOTEENN_y, vae_y]
    
    X_train = []
    X_test = []

    for i in range(len(X_)):
        X = pd.DataFrame(X_[i])
        subset = fs.fit_transform(quant_scaler.inverse_transform(X), y[i])
        subset_test = fs.transform(pd.DataFrame(scaled_test))
        X = X.loc[:, subset.columns]

        X_train.append(X)
        X_test.append(subset_test)
        print(f"Features selected {i+1}/{len(X_)}")


    for model in models:
        print("Training models...")
        for i, balancer in enumerate(annotation):

            if balancer == "weights":
                model.set_params(**param)
            else:
                model.set_params(**no_param)
            
            model.fit(X_train[i], y[i])
            y_pred = model.predict(X_test[i])

            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=y_pred, Y_test=Y_test, labels=unique_labels)
            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], "micro_f1": [micro_f1],
                                        "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm],
                                        "balancer": [balancer], "n_features": [X_train[i].shape[1]]})
            uml.save_results(results_df, "oversampling_evaluation_fs_after")