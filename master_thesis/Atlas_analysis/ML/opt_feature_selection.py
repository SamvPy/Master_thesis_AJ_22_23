# imports
from itertools import combinations
# Models
from sklearn.base import BaseEstimator, TransformerMixin

from xgboost import XGBClassifier
from lightgbm import LGBMClassifier as lgbm
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression, BayesianRidge, Ridge
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans

# pipelines
from sklearn.pipeline import Pipeline
from sklearn.pipeline import make_pipeline

# preprocessors
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler


from sklearn.model_selection import StratifiedKFold

# utils
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.class_weight import compute_class_weight
import pandas as pd
import numpy as np

import sys

from scipy.stats import pearsonr

sys.path.append("../")

import utils_ML as uml

from scipy.stats import pearsonr

# Load datasets

data_quantile = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/quantile_norm_NSAF_50.csv", index_col = "assay_id")
meta = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/unified_metadata.csv")
meta = meta[meta.assay_id.isin(data_quantile.index)]

groups = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/group_cells_annotation.csv", sep =";", index_col="Unnamed: 0")
meta["Group"] = meta.cell_line.apply(lambda x: groups[groups.cell_line == x]["group"].values[0])
meta = meta.set_index("assay_id")

data_quantile.sort_index(inplace=True)
meta.sort_index(inplace=True)

data_quantile = data_quantile.reset_index(drop=True).rename(columns={data_quantile.columns[x]:x for x in range(len(data_quantile.columns))})
target_encoder = LabelEncoder()
targets = target_encoder.fit_transform(meta.Group)
unique_labels = pd.Series(targets).unique()
class_weights = compute_class_weight(class_weight='balanced', classes=unique_labels, y=targets)
weights = {unique_labels[i]: class_weights[i] for i in range(len(unique_labels))}

labels = np.array(meta.Group)



skf = StratifiedKFold(n_splits=10, shuffle=True)

lr_clf = LogisticRegression(max_iter=10000, class_weight=weights)
svm_clf = SVC(class_weight=weights)
rf_clf = RandomForestClassifier(class_weight=weights)
lgbm_clf = lgbm(class_weight=weights)

models = [lr_clf, svm_clf, rf_clf, lgbm_clf]


selectors = "anova MI LR RF SVC".split()
model_names = "lr svm rf lgbm".split()
model_scores = {f:{x:[] for x in selectors} for f in model_names}

fold = 0
for train, test in skf.split(X=data_quantile, y=targets):
    print(fold)

    X_train = data_quantile.loc[train, :].reset_index(drop=True)
    X_test = data_quantile.loc[test,:].reset_index(drop=True)

    Y_train = targets[train]
    Y_test = targets[test]

    #imputer = uml.MNAR_MCAR_Imputer(max_iter=15)
    imputer = uml.LowestValueImputerGaussian()
    imputer.fit(X_train, Y_train)
    imputed_train = imputer.transform(X_train, Y_train)
    imputed_test = imputer.transform(X_test)

    for selector in selectors:
        fs = uml.FeatureSelector(selectors=[selector], num_features=100)
        fs.fit(imputed_train, Y_train)
        X_train_fs = fs.transform(imputed_train)
        X_test_fs = fs.transform(imputed_test)

        scaler = MinMaxScaler()
        X_train_scaled = scaler.fit_transform(X_train_fs)
        X_test_scaled = scaler.fit_transform(X_test_fs)
    
        print("Fitting models...")
        for model_i,model in enumerate(models):
                
                model.fit(X_train_scaled, Y_train)

                Y_pred = model.predict(X_test_scaled)
            
                micro_f1, macro_f1, weighted_f1, _ = uml.scoring_functions(Y_pred=Y_pred, Y_test=Y_test, labels=unique_labels)
                results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], "micro_f1": [micro_f1],
                                        "macro_f1": [macro_f1], "weighted_f1": [weighted_f1],
                                        "selector": [selector]})
                model_scores[model_names[model_i]][selector].append(macro_f1)
                uml.save_results(results_df, "fs_comparison_LVIG_100")
    fold+=1