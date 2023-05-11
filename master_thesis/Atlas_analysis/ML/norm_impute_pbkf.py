# imports

# Models
from sklearn.base import BaseEstimator, TransformerMixin

import xgboost as xgb
from lightgbm import LGBMClassifier as lgbm
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier

# pipelines
from sklearn.pipeline import Pipeline
from sklearn.pipeline import make_pipeline

# preprocessors
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler, quantile_transform
from sklearn.feature_selection import chi2, f_classif, mutual_info_classif, SelectKBest, SelectFromModel, RFE
from sklearn.impute import SimpleImputer, KNNImputer


# metrics and splitters
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score

from sklearn.model_selection import StratifiedKFold, GridSearchCV, RandomizedSearchCV
from sklearn.model_selection import train_test_split

# utils
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.class_weight import compute_class_weight
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("../")

import utils_ML as uml


# Load datasets
data_combat = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/combat_NSAF_50.csv", index_col = "assay_id")
data_quantile = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/quantile_norm_NSAF_50.csv", index_col = "assay_id")
data_median_norm = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/median_scaling_50.csv", index_col = "Unnamed: 0")
data_nsaf = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/NSAF_50.csv", index_col = "assay_id")
data_nsaf = np.log2(data_nsaf)

meta = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/unified_metadata.csv")
meta = meta[meta.assay_id.isin(data_combat.index)]

groups = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/group_cells_annotation.csv", sep =";", index_col="Unnamed: 0")
meta["Group"] = meta.cell_line.apply(lambda x: groups[groups.cell_line == x]["group"].values[0])
meta = meta.set_index("assay_id")

data_combat.sort_index(inplace=True)
data_quantile.sort_index(inplace=True)
data_median_norm.sort_index(inplace=True)
data_nsaf.sort_index(inplace=True)
meta.sort_index(inplace=True)

missing_value_mask = data_nsaf.isna()
data_combat = data_combat.where(~missing_value_mask, other=np.nan)

data_combat = data_combat.reset_index(drop=True).rename(columns={data_combat.columns[x]:x for x in range(len(data_combat.columns))})
data_quantile = data_quantile.reset_index(drop=True).rename(columns={data_quantile.columns[x]:x for x in range(len(data_quantile.columns))})
data_median_norm = data_median_norm.reset_index(drop=True).rename(columns={data_median_norm.columns[x]:x for x in range(len(data_median_norm.columns))})
data_nsaf = data_nsaf.reset_index(drop=True).rename(columns={data_nsaf.columns[x]:x for x in range(len(data_nsaf.columns))})


target_encoder = LabelEncoder()
targets = target_encoder.fit_transform(meta.Group)
unique_labels = pd.Series(targets).unique()
class_weights = compute_class_weight(class_weight='balanced', classes=unique_labels, y=targets)

weights = {unique_labels[i]: class_weights[i] for i in range(len(unique_labels))}
print(weights)


lr_clf = LogisticRegression(max_iter=10000)
svm_clf = SVC()
rf_clf = RandomForestClassifier()
xgb_clf = xgb.XGBClassifier()
lgbm_clf = lgbm()
nb_clf = GaussianNB()
knn_clf = KNeighborsClassifier()

models = [lr_clf, svm_clf, rf_clf, lgbm_clf]

# Parameter grids
lr_grid = {"penalty" : ['l2', 'l1'],
            "dual": [False],
            "max_iter": [10000],
            "class_weight": [weights],
            "C": np.linspace(0.005, 15, 10),
            'solver': ['newton-cg', 'sag', 'lbfgs', "liblinear"]}

svc_grid = {'decision_function_shape': ["ovr"],
            "kernel": ['linear', 'poly', 'rbf'],
            "C": np.linspace(0.0005, 5, 10),
            "class_weight": [weights]}

rf_grid = {'n_estimators': np.linspace(10, 200, 4, dtype = int),
            "criterion": ["entropy", "gini"], 
            "max_depth": [10,20,40, None],
            "class_weight": [weights]}

xgb_grid = {"verbosity": [0],
            'eta': np.linspace(0.005,0.5,5),
            'gamma': np.linspace(0.005,10,5),
            'max_depth': [3,5,7,10]}

gnb_grid = {'var_smoothing': np.logspace(0,-15,15)}


grids = {"lr": lr_grid, "svc": svc_grid, "rf": rf_grid, "gnb": gnb_grid}


imputation_methods = [uml.MNAR_MCAR_Imputer(max_iter=15, MCAR_estimator='knn'), uml.MNAR_MCAR_Imputer(max_iter=15, MCAR_estimator='pca'), uml.LowestValueImputerGaussian(), KNNImputer(n_neighbors=10), uml.MNAR_MCAR_Imputer(missing_percentage=1, MCAR_estimator='pca', max_iter=15)]
imputation_ids = ["KNN+LOD", "PCA+LOD", "Shifted Gaussian", "KNNImpute", "PCA"]
dataset_ids = ["NSAF", "ComBat", "Quantile", "Median Norm"]

selector = uml.FeatureSelector(selectors=["anova", "MI", "LR", "SVC"], num_features=300, threshold=.75)

# Step one: Estimate performance of vanilla models on each type of normalization, each imputed with a certain imputation strategy, with feature selection
# Performance measured both with SKF and Projectbased splits
# Measured: Which model performs optimal overall? Which imputation method performs optimal overall? Which dataset performs optimal overall?

pbkf = uml.ProjectBasedSplit(splits=10, metadata=meta, on="Group", LOOCV=True)

fold = 0
for train, test in pbkf.split(dataset=data_nsaf, metadata=meta, n_projects=35):
    fold+=1
    # For each normalization method
    # Get train datasets
    X_train_nsaf = data_nsaf.iloc[train, :].reset_index(drop=True)
    X_train_combat = data_combat.iloc[train, :].reset_index(drop=True)
    X_train_quantile = data_quantile.iloc[train,:].reset_index(drop=True)
    X_train_med = data_median_norm.iloc[train, :].reset_index(drop=True)
    X_train = [X_train_nsaf, X_train_combat, X_train_quantile, X_train_med]

    # Get test datasets
    X_test_nsaf = data_nsaf.iloc[test, :].reset_index(drop=True)
    X_test_combat = data_combat.iloc[test, :].reset_index(drop=True)
    X_test_quantile = data_quantile.iloc[test,:].reset_index(drop=True)
    X_test_med = data_median_norm.iloc[test, :].reset_index(drop=True)
    X_test = [X_test_nsaf, X_test_combat, X_test_quantile, X_test_med]

    # Get label stratifications
    Y_train = targets[train]
    Y_test = targets[test]

    # For each imputation method
    for impute_i, imputation_method in enumerate(imputation_ids):
        # For each normalization method
        for i in range(4):
            
            # Impute
            if imputation_method == "KNNImpute":
                scaler = MinMaxScaler()
                scaled_train = scaler.fit_transform(X_train[i])
                scaled_test = scaler.fit_transform(X_test[i])

                imputer = KNNImputer(n_neighbors=10)
                imputer.fit(scaled_train)
                scaled_train = imputer.transform(scaled_train)
                scaled_test = imputer.transform(scaled_test)

            else:
                imputer = imputation_methods[impute_i]
                imputer.fit(X_train[i], Y_train)
                imputed_train = imputer.transform(X_train[i], Y_train)
                imputed_test = imputer.transform(X_test[i])

                # Standardize
                scaler = MinMaxScaler()
                scaled_train = scaler.fit_transform(imputed_train)
                scaled_test = scaler.transform(imputed_test)

            selector.fit(scaled_train, Y_train)
            X_train_selection = selector.transform(scaled_train)
            X_test_selection = selector.transform(scaled_test)

            # For each model
            print("Training models")
            for clf in models:
                #if type(clf).__name__ == 'XGBClassifier':
                    
                #    xgb_train = xgb.DMatrix(X_train_selection, label=Y_train, weight = [weights[x] for x in Y_train])
                #    xgb_test = xgb.DMatrix(X_test_selection)
                    
                #    clf.fit(xgb_train, Y_train)
                #    Y_pred = clf.predict(xgb_test)
                
                #else:
                clf.set_params(class_weight = weights)
                clf.fit(X_train_selection, Y_train)
                Y_pred = clf.predict(X_test_selection)

                train_pred_unique_labels = np.unique(np.append(Y_test, Y_pred))

                micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=Y_pred, Y_test=Y_test, labels=train_pred_unique_labels)
                results_df = pd.DataFrame({"model": [type(clf).__name__], "fold": [fold], "micro_f1": [micro_f1],
                                        "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm],
                                        "Imputation": [imputation_method], "Normalization": [dataset_ids[i]], "PXD": [pbkf.taken_PXD[-1]]})
                
                uml.save_results(results_df, "norm_imputer_evaluation_pbkf2")

                print("Fold, imputation: ",fold, imputation_method)