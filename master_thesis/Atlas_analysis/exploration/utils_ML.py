# Imports
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import chi2, f_classif, mutual_info_classif, SelectKBest, SelectFromModel, RFE
from sklearn.utils.validation import check_is_fitted
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score, mean_squared_error
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_auc_score
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.decomposition import PCA

import random

import numpy as np
import pandas as pd

# based on: https://towardsdatascience.com/coding-a-custom-imputer-in-scikit-learn-31bd68e541de

# Filtering reoccuring proteins functions
def _identify_global_reoccured_proteins(subset, percentage_reoccurence):
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
    def __init__(self, splits, metadata, LOOCV=False, on = str):
        """Called when training model and splitting_procedure is set as 'project'
        
        metadata: the metadata table that is used to generate splits
        
        on: the column name that represent the class column name"""

        self.LOOCV = LOOCV
        self.splits = splits
        self.metadata = metadata.reset_index(drop=True)
        self.label = on
        self.label_indices = list(range(metadata[self.label].nunique()))
        self.dropped_pxds = []

    def split(self, dataset, metadata, n_projects=5, groups = None):
        self.taken_PXD = []
        self.n_projects = n_projects
        dataset = dataset.reset_index(drop=True)
        metadata = self.metadata.loc[self.metadata.index.isin(dataset.index),:]

        dataset = dataset.sort_index()
        metadata = metadata.sort_index()

        index_splits = []

        if self.LOOCV:
            for split in range(n_projects):
                train_index, test_index, dropped_pxd = self.LOOCV_splitter(dataset, metadata=metadata)
                index_splits.append((train_index,test_index))
                self.dropped_pxds.append(dropped_pxd)

                yield train_index,test_index

        else:
            for split in range(self.splits):
                train_index, test_index, dropped_pxds = self.train_test_project_split(dataset, metadata=metadata)
                index_splits.append((train_index, test_index))
                self.dropped_pxds.append(dropped_pxds)

                yield train_index, test_index

    def LOOCV_splitter(self, dataset, metadata, groups=None):
        
        for PXD in self.metadata.PXD_accession.unique():

            if PXD in self.taken_PXD:
                continue

            # All classes must be in the training set
            if self.metadata[~self.metadata.PXD_accession.isin([PXD])].groupby(self.label).PXD_accession.nunique().shape[0] != 15:
                continue

            # Condition needed to oversample (at least 3 neighbours of the class)
            if (self.metadata[~self.metadata.PXD_accession.isin([PXD])].Group.value_counts() < 4).sum() > 0:
                continue
            
            self.taken_PXD.append(PXD)
            break
        
        test_index = self.metadata[self.metadata.PXD_accession.isin([PXD])].index
        train_index = dataset.loc[~dataset.index.isin(test_index), :].index.to_numpy()

        return train_index, test_index, [PXD]
        
    def train_test_project_split(self, dataset, metadata, groups = None):

        indices = list(range(15))
        random.shuffle(indices)

        choosen_PXD = []

        for group, PXD in self.metadata.groupby(self.label).PXD_accession.unique().iloc[indices].iteritems():
            
            if True in [pxd in choosen_PXD for pxd in PXD]:
                continue
            if len(PXD) > 1:
                choosen_PXD.append(random.choice(PXD))

            # All classes must be in the training set
            if self.metadata[~self.metadata.PXD_accession.isin(choosen_PXD)].groupby(self.label).PXD_accession.nunique().shape[0] != 15:
                choosen_PXD = choosen_PXD[:-1]
            
            # Condition needed to oversample (at least 3 neighbours of the class)
            if (self.metadata[~self.metadata.PXD_accession.isin(choosen_PXD)].Group.value_counts() < 4).sum() > 0:
                choosen_PXD = choosen_PXD[:-1]
                
            if len(choosen_PXD) == self.n_projects:
                break

        test_index = self.metadata[self.metadata.PXD_accession.isin(choosen_PXD)].index
        train_index = dataset.loc[~dataset.index.isin(test_index), :].index.to_numpy()

        return train_index, test_index, choosen_PXD   

    def get_n_splits(self, x, y, groups = None):
        return self.splits


# ----------------------------------------------------------------------------
# Imputers

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

# Variant of zero-imputation: left shifted gaussian, similar to the one implemented in Perseus (https://www.nature.com/articles/nmeth.3901)
class LowestValueImputerGaussian(BaseEstimator, TransformerMixin):
    '''
    Imputes by sampling from left shifted Gaussian distribution. Imputes 0 if protein was not identified in 50% of samples belonging to the class
    '''

    def __init__(self, mean_shift=2, std_scale=0.3, lowest_val = False):
        self.fitted = False
        self.mean_shift = mean_shift
        self.std_scale = std_scale
        self.lowest_val = lowest_val
        self.imputed_values = {}
    
    def fit(self, X, y=None):
        
        X = pd.DataFrame(X)
        self.mean_std_df = pd.DataFrame({"mean": X.mean(), "std": X.std()})

        #imputed_data = X.apply(lambda x: self.impute(x), axis=1)

        return self

    def transform(self, X, y=None):
        check_is_fitted(self, 'mean_std_df')

        X_imputed = pd.DataFrame(X, columns=self.mean_std_df.index).copy()

        # Create imputation matrix
        imputation_matrix = []
        for i in self.mean_std_df.index:
            m, s = self.mean_std_df.loc[i]
            imputation_matrix.append(np.random.normal(m-self.mean_shift*s, s*self.std_scale, size=len(X)))
        imputation_matrix = pd.DataFrame(data=imputation_matrix, index=self.mean_std_df.index).T

        if self.lowest_val:
            imputation_matrix.where(imputation_matrix>self.lowest_val, other = self.lowest_val, inplace=True)

        for i in X_imputed.columns:
            X_imputed[i].fillna(imputation_matrix[i], inplace=True)
        
        return X_imputed

# Initialize random values
def PCA_imputation(dataset, explained_variance=.95, max_iter=50, keep_PC_constant = True, get_knn_corr=False):
    '''Accepts non minmax normalized dataset'''

    dataset = dataset.reset_index(drop=True).rename(columns={dataset.columns[x]:x for x in range(len(dataset.columns))})
    mm_scaler = MinMaxScaler()
    scaled_dataset = pd.DataFrame(mm_scaler.fit_transform(dataset))
    non_MV_indices = dataset.melt()["value"].dropna().index
    MV_indices = dataset.melt().index[dataset.melt().isna()["value"]]
    MSEs = []
    MSEs_unrescaled = []
    corr_PCA_KNN = []

    if get_knn_corr:
        knn_imputer = KNNImputer(n_neighbors=10, weights="distance")
        X_knn_imputed = knn_imputer.fit_transform(scaled_dataset)
        X_knn_imputed = pd.DataFrame(mm_scaler.inverse_transform(X_knn_imputed))
        mv_knn = X_knn_imputed.melt().loc[MV_indices, 'value']

    for i in range(max_iter):

        # Initialize missing values with mean imputation
        if i == 0:
            scaler = MinMaxScaler()
            dataset_scaled = scaler.fit_transform(dataset)
            mean_imputation = SimpleImputer(strategy='mean')
            X_imputed = mean_imputation.fit_transform(dataset_scaled)
            print("First iteration:", i)

        # Else use the calculated PC imputation and iterate
        else:
            # Scale data
            scaler = MinMaxScaler()
            X_imputed = scaler.fit_transform(X_imputed)

        
        print("iteration:", i)
        # Select component explaining 95% of the data 
        if not keep_PC_constant or (i ==0):  
            pca_model = PCA()
            embedded = pca_model.fit_transform(X_imputed)
            components_to_keep = sum(np.cumsum(pca_model.explained_variance_ratio_) < explained_variance)
            print('Components:', components_to_keep)

        # Refit and reconstruct
        pca_model = PCA(n_components=components_to_keep)
        embedded = pca_model.fit_transform(X_imputed)
        reconstructed = pca_model.inverse_transform(embedded)

        # Standardise back
        rescaled_reconstructed = scaler.inverse_transform(reconstructed)
        rescaled_reconstructed = pd.DataFrame(rescaled_reconstructed)

        # Fill in new dataframe with computed missing values retaining non-missing values
        X_imputed = dataset.copy()
        for i in X_imputed.columns:
            X_imputed[i].fillna(rescaled_reconstructed[i], inplace=True)
        
        # Compute reconstruction error with MSE of the observed values only. This is done by setting the missing values to same value in both dataframes (original and reconstructed)
        MSE_unrescaled = mean_squared_error(scaled_dataset.melt()["value"].dropna(), pd.DataFrame(reconstructed).melt()["value"][non_MV_indices])
        MSEs_unrescaled.append(MSE_unrescaled)

        MSE = mean_squared_error(dataset.melt()["value"].dropna(), rescaled_reconstructed.melt()["value"][non_MV_indices])
        MSEs.append(MSE)

        if get_knn_corr:
            mv_pca = rescaled_reconstructed.melt().loc[MV_indices, 'value']
            corr_PCA_KNN.append(pearsonr(mv_knn, mv_pca)[0])

    if get_knn_corr:
        return X_imputed, MSEs, MSEs_unrescaled, corr_PCA_KNN

    return X_imputed, MSEs, MSEs_unrescaled


class MNAR_MCAR_Imputer(BaseEstimator, TransformerMixin):
    '''
    Imputes values either with left-censored gaussian or PCA/KNN imputation.

    This delineation is only made during the fit phase for training data when labels are seen.
    For test data, PCA/KNN imputation is done dependant on which one was chosen during initialization ('pca', 'knn')

    MNAR is done when > a% of a protein in a group is missing.
    For example: A0AVT1 is only found in (a=.8) 20% of samples in group Glioblastoma --> perform MNAR method for this protein for samples in Glioblastoma

    '''

    def __init__(self, missing_percentage = .8, MCAR_estimator = "pca", mean_shift=2, 
                 std_scale=.3, lowest_val=False, explained_variance=.95, max_iter=50, keep_PC_constant=True,
                 n_neighbors = 10, weights='distance'):
        
        assert MCAR_estimator in ["pca", "knn"]
        self.missing_percentage = missing_percentage
        self.MCAR_estimator = MCAR_estimator
        self.fitted = False
        self.mean_shift = mean_shift
        self.std_scale = std_scale
        self.lowest_val = lowest_val
        self.explained_variance = explained_variance
        self.max_iter = max_iter
        self.keep_PC_constant = keep_PC_constant
        self.n_neighbors = n_neighbors
        self.weights = weights

    def fit(self, X, y):
        assert len(X) == len(y)
        X = pd.DataFrame(X).reset_index(drop=True).rename(columns={X.columns[x]:x for x in range(len(X.columns))})
        
        self.mean_std_df = pd.DataFrame({"mean": X.mean(), "std": X.std()})
        
        self.mv_LOD_group = {}
        for group in np.unique(y):

            mv_perc_group = X[y==group].isna().sum() / (y==group).sum()
            self.mv_LOD_group[group] = X.columns[mv_perc_group < self.missing_percentage]
            
        return self

    def transform(self, X, y=np.array([None])):
        check_is_fitted(self, 'mean_std_df')

        # If train set is imputed and LOD must be used
        if not (y == None).all():
            
            # 1. MNAR-imputation

            X_imputed = pd.DataFrame(X).reset_index(drop=True).rename(columns={X.columns[x]:x for x in range(len(X.columns))})
            melted_df = X_imputed.melt()
            all_mv_indices = melted_df.index[melted_df.isna()["value"]]

            # Create imputation matrix. This can be accessed to compare which values would be imputed if PCA/KNN would be used for the complete dataset
            imputation_matrix = []
            for i in self.mean_std_df.index:
                m, s = self.mean_std_df.loc[i]
                imputation_matrix.append(np.random.normal(m-self.mean_shift*s, s*self.std_scale, size=len(X)))
            self.imputation_matrix = pd.DataFrame(data=imputation_matrix, index=self.mean_std_df.index).T

            # Make sure no values are imputed below 0 (when minmax normalized) (this case, no minmax is assumed)
            #self.imputation_matrix.where(self.imputation_matrix>0, other = 0, inplace=True)

            # Set imputation columns that do not fullfill MNAR-imputation condition to np.nan. This prevents imputation for certain cells.
            for group in np.unique(y):
                not_impute_columns = self.mv_LOD_group[group]
                self.imputation_matrix.loc[y==group, not_impute_columns] = np.nan
            
            for i in X_imputed.columns:
                X_imputed[i].fillna(self.imputation_matrix[i], inplace=True)

            # 2. PCA/KNN imputation
            # Store how each MV is interpreted (either MNAR or MCAR) for the training data
            melted_df = X_imputed.melt()
            self.pca_knn_mv_indices = melted_df.index[melted_df.isna()["value"]]
            self.lod_mv_indices = [x for x in all_mv_indices if x not in self.pca_knn_mv_indices] 
            
            self.X_pca_imputed, self.MSE_rescaled, self.MSEs_scaled = PCA_imputation(X_imputed, explained_variance=self.explained_variance, max_iter=self.max_iter, keep_PC_constant=self.keep_PC_constant)

            knn_imputer = KNNImputer(n_neighbors=self.n_neighbors, weights=self.weights)
            self.X_knn_imputed = knn_imputer.fit_transform(X_imputed)

            mv_knn = pd.DataFrame(self.X_knn_imputed).melt().loc[self.pca_knn_mv_indices, 'value']
            mv_pca = pd.DataFrame(self.X_pca_imputed).melt().loc[self.pca_knn_mv_indices, 'value']
                
            self.corr_knn_PCA = pearsonr(mv_knn, mv_pca)
                
            if self.MCAR_estimator == "pca":
                return self.X_pca_imputed
            elif self.MCAR_estimator == "knn":
                return self.X_knn_imputed

        # else test set must be fitted. Retrain PCA/KNN together with imputed training data
        else:
            check_is_fitted(self, 'X_pca_imputed', msg="First transform a training set by providing an array of labels to y during transform")

            X_imputed = pd.DataFrame(X, columns=self.mean_std_df.index).reset_index(drop=True)
            train_and_test = pd.concat([self.X_pca_imputed, X_imputed]).reset_index(drop=True)
            test_indices = list(range(len(self.X_pca_imputed), len(train_and_test)))

            melted_train_and_test = train_and_test.melt()
            test_mv_indices = melted_train_and_test.index[melted_train_and_test.isna()["value"]]

            self.X_pca_imputed_test, self.MSE_rescaled_test, self.MSEs_scaled_test = PCA_imputation(train_and_test, explained_variance=self.explained_variance, max_iter=self.max_iter, keep_PC_constant=self.keep_PC_constant)
            self.X_pca_imputed_test = pd.DataFrame(self.X_pca_imputed_test).loc[test_indices,:]

            knn_imputer = KNNImputer(n_neighbors=self.n_neighbors, weights=self.weights)
            self.X_knn_imputed_test = knn_imputer.fit_transform(train_and_test)
            self.X_knn_imputed_test = pd.DataFrame(self.X_knn_imputed_test).loc[test_indices, :]
                
            if self.MCAR_estimator == "pca":
                return self.X_pca_imputed_test
            elif self.MCAR_estimator == "knn":
                return self.X_knn_imputed_test
            
#--------------------------------------------------------------------------------------
# Filters 

class FilterByOccurence(BaseEstimator, TransformerMixin):
    '''
    Filters the proteins that occur in 50% of the training samples
    '''
    def __init__(self, percentage=.5):
        self.percentage = percentage

    def fit(self, X, y = None):
        proteins, _ = _identify_global_reoccured_proteins(X.fillna(0), self.percentage)
        
        self.filtered_proteins = proteins
        self.n_features = len(proteins)

        return self

    def transform(self, X, y=None):
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

    def fit(self, X, y = False):
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

    def transform(self, X, y=False):
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
                self.selector_models.append(SelectFromModel(LogisticRegression(penalty="l1", solver = "liblinear", multi_class="ovr"), max_features=self.num_features))
            if selector == "SVC":
                self.selector_models.append(RFE(SVC(class_weight="balanced", kernel = "linear"), n_features_to_select=self.num_features, step=50))

        self.supports = []
        
        for i, selector in enumerate(self.selector_models):
            if self.selectors[i] in ["LR", "SVC"]:
                X_scaled = MinMaxScaler().fit_transform(X)
                selector.fit(X_scaled, y)
                self.supports.append(selector.get_support())
            selector.fit(X,y)
            self.supports.append(selector.get_support())
        
        votes = np.sum(self.supports, axis = 0) / len(self.selector_models)
        final_mask = votes >= self.threshold
        self.final_mask = final_mask

        # print(f"Selected {sum(final_mask)} features.")

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

def save_results(df, name_file):
    """Writes results away to a file"""
    try:
        file = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/ML/results/{}.csv".format(name_file), sep = ";")
        df = pd.concat([file, df], ignore_index = True)
        df.to_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/ML/results/{}.csv".format(name_file), sep = ";", index=False)
    except:
        df.to_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/ML/results/{}.csv".format(name_file), sep = ";", index=False)

def main():
    return

if __name__ == '__main__':
    main()