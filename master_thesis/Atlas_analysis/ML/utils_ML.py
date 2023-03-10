# Imports
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import chi2, f_classif, mutual_info_classif, SelectKBest, SelectFromModel, RFE
from sklearn.utils.validation import check_is_fitted
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score

import random

import numpy as np
import pandas as pd

# based on: https://towardsdatascience.com/coding-a-custom-imputer-in-scikit-learn-31bd68e541de

# Filtering reoccuring proteins functions
def _identify_global_reoccured_proteins(subset: pd.DataFrame, percentage_reoccurence):
    '''Returns protein_ids that reoccur in >= percentage_reoccurence of the subset
    
    Also returns protein_ids that are filtered out as second term'''

    # Get percentage of samples in the dataset a protein occurs in 
    reoccuring_proteins = (subset != 0).sum() / len(subset)

    # Get the proteins that do not meet the required percentage of occurence but are identified in at least 1 sample
    deleted_proteins = reoccuring_proteins[(reoccuring_proteins < percentage_reoccurence) & (reoccuring_proteins != 0)].index.tolist()

    # Get reoccuring proteins
    result = reoccuring_proteins[(reoccuring_proteins >= percentage_reoccurence) & (reoccuring_proteins != 0)].index.tolist()
    
    return result, deleted_proteins

# Splitter
# Index must be 0-n_samples bcz prob index is taken with iloc instead of loc
class ProjectBasedSplit():
    def __init__(self, splits: int, metadata: pd.DataFrame, on = str):
        """Called when training model and splitting_procedure is set as 'project'
        
        metadata: the metadata table that is used to generate splits
        
        on: the column name that represent the class column name"""

        self.splits = splits
        self.metadata = metadata.reset_index(drop=True)
        self.label = on
        self.label_indices = list(range(metadata[self.label].nunique()))
        self.dropped_pxds = []

    def split(self, dataset, metadata, groups = None):
        
        dataset = dataset.reset_index(drop=True)
        metadata = self.metadata.loc[self.metadata.index.isin(dataset.index),:]

        index_splits = []
        for split in range(self.splits):
            train_index, test_index, dropped_pxds = self.train_test_project_split(dataset, metadata=metadata)
            index_splits.append((train_index, test_index))
            self.dropped_pxds.append(dropped_pxds)

            yield train_index, test_index
        
    def train_test_project_split(self, dataset, metadata: pd.DataFrame, groups = None):

        indices = list(range(15))
        random.shuffle(indices)

        choosen_PXD = []

        for group, PXD in self.metadata.groupby(self.label).PXD_accession.unique().iloc[indices].iteritems():
            
            if True in [pxd in choosen_PXD for pxd in PXD]:
                continue
            if len(PXD) > 1:
                choosen_PXD.append(random.choice(PXD))
            if self.metadata[~self.metadata.PXD_accession.isin(choosen_PXD)].groupby(self.label).PXD_accession.nunique().shape[0] != 15:
                choosen_PXD = choosen_PXD[:-1]
            if len(choosen_PXD) == 5:
                break

        test_index = self.metadata[self.metadata.PXD_accession.isin(choosen_PXD)].index
        train_index = dataset.loc[~dataset.index.isin(test_index), :].index.to_numpy()

        return train_index, test_index, choosen_PXD   

    def get_n_splits(self, x, y, groups = None):
        return self.splits


# Imputer

class LowestValueImputer(BaseEstimator, TransformerMixin):
    '''
    Imputes the lowest value in that column
    '''

    def __init__(self):
        self.fitted = False

    def fit(self, X, y = None):
        self.impute_map = X.min()

        # If protein columns has all nans, impute the lowest value in the dataframe
        self.impute_map = self.impute_map.fillna(self.impute_map.min())
        return self

    def transform(self, X, y=None):
        check_is_fitted(self, 'impute_map')
        
        X = X.copy()

        for index, row in self.impute_map.iteritems():
            ind = X.index[X.loc[:,index].isna()]
            X.loc[ind, index] = X.loc[ind, index].fillna(self.impute_map[index])
        
        return X

#--------------------------------------------------------------------------------------
# Filters 

class FilterByOccurence(BaseEstimator, TransformerMixin):
    '''
    Filters the proteins that occur in 50% of the training samples
    '''
    def __init__(self, percentage=.5):
        self.percentage = percentage

    def fit(self, X: pd.DataFrame, y = None):
        proteins, _ = _identify_global_reoccured_proteins(X.fillna(0), self.percentage)
        
        self.filtered_proteins = proteins
        self.n_features = len(proteins)

        return self

    def transform(self, X: pd.DataFrame, y=None):
        check_is_fitted(self, 'filtered_proteins')

        return X.loc[:, self.filtered_proteins]

class FilterByClass(BaseEstimator, TransformerMixin):
    '''
    Filters the proteins by class that occur in 50% of the training samples

    keep: bool (keep the proteins of a sample in a class that that do not occur in 50% of the samples in that class yet do in other classes)
    '''

    def __init__(self, keep=False, percentage=0.5):
        self.keep = keep
        self.filter_per_class = {}
        self.percentage = percentage

    def fit(self, X: pd.DataFrame, y = False):
        assert y.any()

        reoccuring_proteins = []

        # Make subsets
        for cls in pd.Series(y).unique():
            ind = np.where(y == cls)[0]
            subset = X.reset_index(drop=True).loc[ind,:]

            proteins, deleted = _identify_global_reoccured_proteins(subset=subset.fillna(0), percentage_reoccurence=self.percentage)
            reoccuring_proteins += proteins
            
            self.filter_per_class[cls] = proteins

        self.filtered_proteins = set(reoccuring_proteins)

        return self

    def transform(self, X: pd.DataFrame, y=False):
        check_is_fitted(self, 'filtered_proteins')

        if self.keep:
            return X.loc[:, self.filtered_proteins]

        # Data leakage is necessary to execute this function. The effect of leakage is evaluated by shuffling the labels
        assert y.any()
        # Make the abundance values np.nan that do not occur > percentage in samples of a class
        X = X.copy()
        X = X.reset_index(drop=True)

        for cls in pd.Series(y).unique():
            ind_row = np.where(y == cls)[0]
            X.loc[ind_row, ~X.columns.isin(self.filter_per_class[cls])] = np.nan

        return X.loc[:, self.filtered_proteins]


class FeatureSelector(BaseEstimator, TransformerMixin):
    def __init__(self, selectors = ['anova', 'MI', 'LR'], num_features = 2000, threshold = 0.5):
        """
        Kind: Which feature selector to use? one/multiple of the following:
            - anova
            - chi2
            - MI
            - LR
            - SVC
            - all
        """

        self.selectors = selectors
        self.num_features = num_features
        self.threshold = threshold
        

    def fit(self, X, y):
    
        self.selector_models = []

        for selector in self.selectors:
            if selector == "anova":
                self.selector_models.append(SelectKBest(f_classif, k = self.num_features))
            if selector == "chi2":
                self.selector_models.append(SelectKBest(chi2, k = self.num_features))
            if selector == "MI":
                self.selector_models.append(SelectKBest(mutual_info_classif, k = self.num_features))
            if selector == "LR":
                self.selector_models.append(SelectFromModel(LogisticRegression(penalty="l1", solver = "liblinear", multi_class="ovr"), max_features=500))
            if selector == "SVC":
                self.selector_models.append(RFE(SVC(class_weight="balanced", kernel = "linear"), n_features_to_select=self.num_features, step=500))

        self.supports = []
        
        for selector in self.selector_models:
            selector.fit(X,y)
            self.supports.append(selector.get_support())
        
        votes = np.sum(self.supports, axis = 0) / len(self.selector_models)
        final_mask = votes >= self.threshold
        self.final_mask = final_mask

        print(f"Selected {sum(final_mask)} features.")

        return self


    def transform(self, X, y = None):
        return pd.DataFrame(X).loc[:,self.final_mask]


    
def scoring_functions(Y_pred, Y_test, labels):
    """Returns pxd_f1, micro_f1, macro_f1, cm"""
        
    micro_f1 = f1_score(Y_test, Y_pred, labels=labels, average = "micro")
    macro_f1 = f1_score(Y_test, Y_pred, labels=labels, average = "macro")
    weighted_f1 = f1_score(Y_test, Y_pred, labels=labels, average = "weighted")
    cm = confusion_matrix(Y_test, Y_pred, labels=labels)
    
    return micro_f1, macro_f1, weighted_f1, cm

def save_results(df: pd.DataFrame, name_file: str):
    """Writes results away to a file"""
    try:
        file = pd.read_csv(f"results/{name_file}.csv", sep = ";")
        df = pd.concat([file, df], ignore_index = True)
        df.to_csv(f"results/{name_file}.csv", sep = ";", index=False)
    except:
        df.to_csv(f"results/{name_file}.csv", sep = ";", index=False)

