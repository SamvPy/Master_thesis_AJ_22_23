import utils_ML as uml

import pandas as pd
import numpy as np

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler, normalize
from sklearn.utils.class_weight import compute_class_weight
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.utils.fixes import loguniform

from imblearn.combine import SMOTEENN, SMOTETomek

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions

device = 'cuda' if torch.cuda.is_available() else "cpu"

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import StratifiedKFold


def main():
    # Machine learning classification loop

    # read data
    print("Reading data...")
    data = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/PEMatrix/norm_NSAF_data2.csv", index_col = "assay_id")
    meta = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/unified_metadata.csv")
    meta = meta[meta.assay_id.isin(data.index)]

    groups = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/group_cells_annotation.csv", sep =";", index_col="Unnamed: 0")
    meta["Group"] = meta.cell_line.apply(lambda x: groups[groups.cell_line == x]["group"].values[0])
    meta = meta.set_index("assay_id")

    data.sort_index(inplace=True)
    meta.sort_index(inplace=True)

    target_encoder = LabelEncoder()
    targets = target_encoder.fit_transform(meta.Group)
    unique_labels = pd.Series(targets).unique()

    class_weights = compute_class_weight(class_weight='balanced', classes=unique_labels, y=targets)

    weights = {unique_labels[i]: class_weights[i] for i in range(len(unique_labels))}

    # Models to be used:
    print("Initializing models and grids...")
    svc_clf = SVC()
    rf_clf = RandomForestClassifier()
    lr_l1_clf = LogisticRegression(max_iter = 10000)
    lr_l2_clf = LogisticRegression(max_iter = 10000)
    gnb_clf = GaussianNB()

    # Parameter grids
    lr_l1_grid = {"penalty" : ['l2'],
                "dual": [False],
                "max_iter": [10000],
                "class_weight": [weights, None],
                "C": loguniform(.001, 1e2),
                'solver': ['newton-cg', 'sag', 'lbfgs', "liblinear"]}

    lr_l2_grid = {"penalty" : ['l1'],
                "dual": [False],
                "max_iter": [10000],
                "class_weight": [weights, None],
                "C": loguniform(.001, 1e2),
                'solver': ["liblinear"]}

    svc_grid = {'decision_function_shape': ["ovr"],
                "kernel": ['linear', 'poly', 'rbf'],
                "C": loguniform(.001, 1e2),
                "class_weight": [weights,None]}

    rf_grid = {'n_estimators': np.linspace(10, 200, 100, dtype = int),
                "criterion": ["entropy", "gini"], 
                "max_depth": [5,10,15,20,40, None],
                "class_weight": [weights,None]}

    gnb_grid = {'var_smoothing': np.logspace(0,-15,100)}

    model_grids = {"lr1": (lr_l1_clf ,lr_l1_grid), "lr2":(lr_l2_clf, lr_l2_grid), 
                "svc": (svc_clf, svc_grid), "rf": (rf_clf, rf_grid), "gnb": (gnb_clf, gnb_grid)}

    """
    preprocessor1 = Pipeline(steps=[
        ('filtering', uml.FilterByClass(keep=True)),
        ('imputation', uml.LowestValueImputer()),
        ('scaler', MinMaxScaler()),
    ])

    preprocessor2 = Pipeline(steps=[
        ('LowestValueImputer', uml.LowestValueImputer()),
        ('scaler', MinMaxScaler()),
        ('feature_selection', uml.FeatureSelector(selectors= ['LR', 'MI', 'SVC'], num_features=3000, threshold=.5))
    ])
    """

    preprocessor3 = Pipeline(steps=[
        ('filtering', uml.FilterByOccurence(0.5)),
        ('imputation', uml.LowestValueImputer()),
        ('scaler', MinMaxScaler()),
    ])

    skf = StratifiedKFold(n_splits=5, shuffle=True)

    fold=0

    print("Starting outer loop... (5 folds)")
    for train, test in skf.split(X=data, y=targets):
        fold += 1
        print(fold)
        # Split data
        X_train = data.iloc[train,:]
        Y_train = targets[train]
        X_test = data.iloc[test,:]
        Y_test = targets[test]

        # Preprocess
        """
        preprocessor1.fit(X_train, Y_train)
        X_train_preprocessed1 = preprocessor1.transform(X_train)
        X_test_preprocessed1 = preprocessor1.transform(X_test)
        print(fold, ': Preprocessed 1/3')

        preprocessor2.fit(X_train, Y_train)
        X_train_preprocessed2 = preprocessor2.transform(X_train)
        X_test_preprocessed2 = preprocessor2.transform(X_test)
        print(fold, ': Preprocessed 2/3')
        """

        preprocessor3.fit(X_train, Y_train)
        X_train_preprocessed3 = preprocessor3.transform(X_train)
        X_test_preprocessed3 = preprocessor3.transform(X_test)
        print(fold, ': Preprocessed 3/3')

        # Oversample on training set
        # X_train_oversampled1, Y_train_oversampled1 = SMOTETomek().fit_resample(X_train_preprocessed1, Y_train)
        # X_train_oversampled2, Y_train_oversampled2 = SMOTETomek().fit_resample(X_train_preprocessed2, Y_train)
        X_train_oversampled3, Y_train_oversampled3 = SMOTETomek().fit_resample(X_train_preprocessed3, Y_train)
        print(fold, ': Oversampled')

        for key, (model, grid) in model_grids.items():
            
            """
            print('Starting gridsearch for', key)
            # Hyperparameter tuning on the training set with nested cv
            clf1 = RandomizedSearchCV(model, grid, n_iter=10, scoring = 'f1_macro', verbose = 3)
            search1 = clf1.fit(X_train_oversampled1, Y_train_oversampled1)
            print(fold, ': grid_search 1/3; score: ', search1.best_score_)

            clf2 = RandomizedSearchCV(model, grid, n_iter=10, scoring = 'f1_macro', verbose = 3)
            search2 = clf2.fit(X_train_oversampled2, Y_train_oversampled2)
            print(fold, ': grid_search 2/3; score: ', search2.best_score_)
            """
            
            clf3 = RandomizedSearchCV(model, grid, n_iter=10, scoring = 'f1_macro', verbose = 3)
            search3 = clf3.fit(X_train_oversampled3, Y_train_oversampled3)
            print(fold, ': grid_search 3/3; score: ', search3.best_score_)
            
            # Predict with optimized hyperparameters
            # Y_pred1 = search1.best_estimator_.predict(X_test_preprocessed1)
            # Y_pred2 = search2.best_estimator_.predict(X_test_preprocessed2)
            Y_pred3 = search3.best_estimator_.predict(X_test_preprocessed3)


            # Compute scores and save
            '''
            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=Y_pred1, 
                                                                        Y_test=Y_test,
                                                                        labels=unique_labels)

            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], 
                                            "micro_f1": [micro_f1], "cv_macro_f1": [search1.best_score_],
                                            "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm], 
                                            "preprocessor": ["class_based"], "best_params": [search1.best_params_],
                                            "oversampler": ["SMOTETomek"]})
                
            uml.save_results(results_df, "grid_search") 

            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=Y_pred2, 
                                                                        Y_test=Y_test,
                                                                        labels=unique_labels)       
            
            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], 
                                            "micro_f1": [micro_f1], "cv_macro_f1": [search2.best_score_],
                                            "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm], 
                                            "preprocessor": ["cpu"], "best_params": [search2.best_params_],
                                            "oversampler": ["SMOTETomek"]})
                
            uml.save_results(results_df, "grid_search")     
            '''
            
            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=Y_pred3, 
                                                                        Y_test=Y_test,
                                                                        labels=unique_labels)       
            
            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], 
                                            "micro_f1": [micro_f1], 'cv_macro_f1': [search3.best_score_],
                                            "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm], 
                                            "preprocessor": ["global"], 'best_params': [search3.best_params_],
                                            "oversampler": ["SMOTETomek"]})
                
            uml.save_results(results_df, "grid_search")
            print("Saved results")

if __name__ == '__main__':
    main()