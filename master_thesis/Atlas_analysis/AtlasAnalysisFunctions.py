import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, GridSearchCV, ParameterGrid
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from xgboost import XGBClassifier
from sklearn.metrics import f1_score
from sklearn.linear_model import LogisticRegression
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import random

# Ignore following warning: 
# /home/compomics/miniconda3/envs/ionbot/lib/python3.7/site-packages/sklearn/metrics/_classification.py:1580: UndefinedMetricWarning: 
# F-score is ill-defined and being set to 0.0 in labels with no true nor predicted samples. Use `zero_division` parameter to control this behavior.
# _warn_prf(average, "true nor predicted", "F-score is", len(true_sum))

import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("always")

from typing import Tuple, Iterable, TypeVar
import os

# Filtering reoccuring proteins functions
def _identify_global_reoccured_proteins(subset: pd.DataFrame, percentage_reoccurence) -> Tuple[list, list]:
    '''Returns protein_ids that reoccur in >= percentage_reoccurence of the subset
    
    Also returns protein_ids that are filtered out as second term'''

    # Get percentage of samples in the dataset a protein occurs in 
    reoccuring_proteins = (subset != 0).sum() / len(subset)

    # Get the proteins that do not meet the required percentage of occurence but are identified in at least 1 sample
    deleted_proteins = reoccuring_proteins[(reoccuring_proteins < percentage_reoccurence) & (reoccuring_proteins != 0)].index.tolist()

    # Get reoccuring proteins
    result = reoccuring_proteins[(reoccuring_proteins >= percentage_reoccurence) & (reoccuring_proteins != 0)].index.tolist()
    
    return result, deleted_proteins

def _identify_local_reoccured_proteins(subset: pd.DataFrame, percentage_reoccurence) -> Tuple[list, list]:
    '''Modified protein_reoccurence filtering
    
    Idea: Many proteins are lost due to some samples identifying +1000 proteins more. These samples will
    certainly contain proteins that are not in the other samples simply due to their much higher identifications
    
    Therefor this new filtering process is proposed by filtering proteins only occuring in the samples identifying more proteins
    that occur in that top x percent of samples. This will drop less proteins than the global reoccured protein filtering.'''

    # Sort the samples according to how much proteins were identified
    subset = subset.reset_index(drop = True)
    sorted_count_subset = subset.sum(axis = 1).sort_values(ascending = False)
    #Get the indices in that order
    indices = sorted_count_subset.index.tolist()

    filter_indices = []
    reoccuring_proteins = []
 
    for index in indices:

        filter_indices = filter_indices + [index]
        
        # Only look at the samples most identified
        filter_subset = subset.iloc[filter_indices,:]
        
        # Get the proteins that are identified in the other samples that have less identifications
        less_quantified = subset.iloc[~subset.index.isin(filter_indices),:]
        less_quantified_proteins = (less_quantified != 0).sum(axis = 0)
        less_quantified_proteins = less_quantified_proteins[less_quantified_proteins != 0].index
        #print(len(less_quantified_proteins))
        
        # Get proteins ONLY present in the samples with most quantifications that occur in x percent or more of those samples
        proteins_to_keep, proteins_to_ignore = _identify_global_reoccured_proteins(filter_subset.loc[:,~subset.columns.isin(less_quantified_proteins)], percentage_reoccurence=percentage_reoccurence)
        reoccuring_proteins += proteins_to_keep

        #print(f"{filter_indices}\nreoccur: {len(proteins_to_keep)}\ndelete: {len(proteins_to_ignore)}")
    deleted_proteins = [x for x in subset.columns if x not in proteins_to_keep]

    return set(reoccuring_proteins), deleted_proteins

def reoccurence_filtering(dataset: pd.DataFrame, labels: pd.Series, percentage_reoccurence: float, how = str, drop_missing_protein_columns = True, info = False) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    '''
    Make sure the dataset indices match the label indices!

    Return a dataset with proteins filtered that do not occur in percentage_reoccurence of the samples corresponding to a class

    If info = True returns:
    
    - filtered_dataframe
    - dataframe containing protein count filtering info per sample
    - dictionary of deleted proteins per class
    '''
    # Define filtering procedure
    if how == "global":
        filter_procedure = _identify_global_reoccured_proteins
    if how == "local":
        filter_procedure = _identify_local_reoccured_proteins

    # Init variables
    counts_before_filtering = (dataset != 0 ).sum(axis = 1)
    dataset_copy = dataset.copy()
    amount_assays = len(dataset_copy)
    deleted_proteins_per_class = {}

    for label in labels.unique():
        # Get subset dataframe of the label
        label_indeces = labels[labels == label].index
        subset = dataset_copy.loc[dataset_copy.index.isin(label_indeces)]

        # Get the proteins that reoccur in that class and the proteins that were filtered out
        reoccuring_proteins, deleted_proteins = filter_procedure(subset=subset, percentage_reoccurence=percentage_reoccurence)
        
        # Set the abundances to 0 for the deleted proteins
        dataset_copy.loc[dataset_copy.index.isin(label_indeces), ~dataset_copy.columns.isin(reoccuring_proteins)] = 0
        deleted_proteins_per_class[label] = deleted_proteins

    # Protein identifications per sample after filtering and how many proteins were lost per sample
    counts_after_filtering = (dataset_copy != 0).sum(axis = 1)
    dropped = counts_before_filtering-counts_after_filtering
    merged_info = pd.concat([counts_before_filtering, counts_after_filtering, dropped], axis = 1)
    merged_info.rename(columns = {0: "before", 1: "after", 2: "proteins_lost"}, inplace = True)
    
    # Delete the protein columns for which the protein do not occur in any sample in the dataset
    if drop_missing_protein_columns:
        sieved_proteins = (dataset_copy == 0).sum()[(dataset_copy == 0).sum()==amount_assays].index
        dataset_copy = dataset_copy.loc[:,~dataset_copy.columns.isin(sieved_proteins)]
    
    if info:
        return dataset_copy, merged_info, deleted_proteins_per_class
    return dataset_copy, "", ""

def _conjunction(proteins: pd.Series, lowerbound, upperbound):
    upperbound = round(len(proteins) * upperbound)
    lowerbound = round(len(proteins) * lowerbound)
    filtered_proteins = proteins[0:upperbound].index.tolist()
    filtered_proteins += proteins[lowerbound:len(proteins)].index.tolist()
    return filtered_proteins

def _union(proteins: pd.Series, lowerbound, upperbound):
    upperbound = round(len(proteins) * upperbound)
    lowerbound = round(len(proteins)* lowerbound)
    return proteins[lowerbound:upperbound].index.tolist()

def quantification_filtering(dataset: pd.DataFrame, upperbound = 1, lowerbound = 1):
    '''Take only the proteins quantified in a certain abundance range.
    
    The ranges are calculated as percentages for each sample seperately and will be the conjunction of the two values. 
    
    Must be values in [0,1] range.
    - upperbound: Most x% quantified proteins
    - Lowerbound: Least x% quantified proteins
    
    Example: 
    
    - UNION (overlap): Upperbound (.8) and lowerbound (.3) --> [70% : 80%] 
    - CONJUNCTION (non-overlap) Upperbound (.6) and lowerbound (.3) --> [0% : 60%] and [70% : 100%]'''
    
    assert upperbound <= 1 and upperbound >= 0 and lowerbound <= 1 and lowerbound >= 0

    upperbound = upperbound
    lowerbound = 1-lowerbound

    if upperbound > lowerbound: # .8 and .3 -> Union
        lower_index = lowerbound # .7
        upper_index = upperbound # .8
        slicing = _union
    else:                       # .6 and .3 -> Conjunction
        lower_index = upperbound # 0 : .6
        upper_index = lowerbound # .7 : 1
        slicing = _conjunction

    dataset_copy = dataset.copy()

    for index, row in dataset.iterrows():
        proteins = row[row != 0].sort_values(ascending=False)
        proteins = slicing(proteins, lower_index, upper_index)
        dataset_copy.loc[index, ~dataset.columns.isin(proteins)] = 0

    return dataset_copy

def project_F1_score(Y_true: pd.Series, Y_pred: np.ndarray, project_ids: pd.Series, labels):
    '''Costum scoring function
        
    - Y_pred: Predictions of the test set
    - Y_test: The true values in the test set
    - project_ids: Series of project ids mapped by index to the test samples
    - labels: The class labels'''
    
    score = 0
    # Enumerate over the projects
    for count, project_id in enumerate(project_ids.unique()):
        project_indices = project_ids[project_ids == project_id].index
        subset_Y_pred = Y_pred[project_indices].tolist()
        subset_Y_test = Y_true.loc[project_indices].tolist()
                
        # Calculate the f1_score for each project separately (value between 0 and 1)    
        score += f1_score(subset_Y_test, subset_Y_pred, labels = labels, average = "micro")
        
    # Return average of all the f1 scores
    return score/(count+1) 



class GridSearchProjectF1():
    def __init__(self, clf: SVC, grid: dict, cv: StratifiedKFold):
        '''Will perform gridsearch with the adapted F1 score metric which normalizes the f1 scores for each project'''

        self.clf = clf
        self.grid = grid
        self.cv = cv

    def fit(self, X_train: pd.DataFrame, Y_train: pd.Series, project_labels: pd.Series):
        # Iterate over every combination of parameters
        self.labels = Y_train.unique()

        best_param = {}
        best_score = {}
        fold = 0
        for train, test in self.cv.split(X_train, Y_train):
            fold += 1

            # Split data in validation and test set.
            X_val = X_train.iloc[train,:].reset_index(drop = True)
            Y_val = Y_train.iloc[train].reset_index(drop = True)
            X_test = X_train.iloc[test,:].reset_index(drop = True)
            Y_test = Y_train.iloc[test].reset_index(drop = True)
            project_ids = project_labels.iloc[test].reset_index(drop = True)

            # Iterate over the parameters
            fold_scores = {}
            for params in ParameterGrid(self.grid):
                clf = self.clf.set_params(**params)
                clf.fit(X_val, Y_val)
                Y_pred = clf.predict(X_test)
                
                # Calculate costum f1-score for each parameter and save them in dictionary {parameter: score} 
                fold_scores[clf.get_params] = project_F1_score(Y_test, Y_pred, project_ids, self.labels)
            
            # Get the best parameters based on the costum f1-score and save both the parameters and the scores
            best_param[fold] = max(fold_scores, key = fold_scores.get)
            best_score[fold] = max(fold_scores.values())
        
        # save the parameters with best f1-score across the folds in the attribute 'best_param'
        print(f"Gridsearch scores: {best_score.values()}")
        best_fold = max(best_score, key = best_score.get)
        self.best_params_ = best_param[best_fold]
        self.clf = self.clf.set_params(**params).fit(X_train, Y_train)

    def predict(self, X_test):
        """Predict the label for an unknown test set. The model with the best parameters is used to predict the labels."""
        return self.clf.predict(X_test)




# Index must be 0-n_samples bcz prob index is taken with iloc instead of loc
class ProjectBasedSplit():
    def __init__(self, splits: int, metadata: pd.DataFrame, on = str):
        """Called when training model and splitting_procedure is set as 'project'
        
        metadata: the metadata table that is used to generate splits
        
        on: the column name that represent the class column name"""

        self.splits = splits
        self.metadata = metadata
        self.label = on
        self.dropped_pxds = []

    def split(self, dataset, metadata, groups = None):
        """Generate indices to split dataset into training and test set

        metadata should be a dataframe containing following columns:
            - PXD_accession
            - assay_id (unique)
            - label (given when initiating this splitting class)
        
        Splits are done with complete projects
        
        No unique splits are ensured.
        
        returns list of sets containing (train_index, test_index)"""
        
        metadata = self.metadata.loc[self.metadata.index.isin(dataset.index),:]

        index_splits = []
        for split in range(self.splits):
            train_index, test_index, dropped_pxds = self.train_test_project_split(dataset, metadata=metadata)
            index_splits.append((train_index, test_index))
            self.dropped_pxds.append(dropped_pxds)

            yield train_index, test_index
        
    def train_test_project_split(self, dataset, metadata: pd.DataFrame, groups = None):

        
        pxds = metadata.PXD_accession.unique()
        n_samples = len(dataset)
        random.shuffle(pxds)
        dropped_pxds = []
        indeces = np.array([],dtype = int)

        for pxd in pxds:
            
            pxd_indeces = np.concatenate((metadata.loc[metadata.PXD_accession == pxd,:].index.to_numpy(), indeces))
            counter = metadata.loc[~metadata.index.isin(pxd_indeces),:].groupby(self.label).nunique()
            assays = len(pxd_indeces)

            condition1 = (assays/n_samples) < 0.1 # Ensure no more than 10% of data in testset
            condition2 = not False in (counter.PXD_accession > 1).unique() # Ensure at least 2 projects per class in training set
            condition3 = not False in (counter.assay_id > 10).unique() # Ensure > 10 samples per class in training set for cv

            if condition1 and condition2 and condition3:
                dropped_pxds.append(pxd)
                indeces = pxd_indeces

        classes_in_test = metadata[metadata.PXD_accession.isin(dropped_pxds)][self.label].unique()
        #if the resulting test train split contains at most 4 projects and represents at most 3 classes, try again
        if len(dropped_pxds)<3 and len(classes_in_test) < 3:
            return self.train_test_project_split(dataset = dataset, metadata=metadata)

        train_index = dataset.loc[~dataset.index.isin(indeces), :].index.to_numpy()
        test_index = indeces 

        return train_index, test_index, dropped_pxds   

    def get_n_splits(self, x, y, groups = None):
        return self.splits





class ModelModule():
    def __init__(self, dataset: pd.DataFrame, metadata: pd.DataFrame):
        '''This class encapsulates every transformation that should be applied on the dataset and metadata including modeltraining.
        
        It is also possible to extract the training data after certain filtering procedures etcetera
        
        This class will ensure data and metadata integrity
        
        Insert the metadata table including project ids.'''

        if len(metadata) != len(dataset):
            raise Exception("Error: Ensure equal metadata and dataset length")

        # Make sure the indices match for metadata and dataset
        self.metadata = metadata.reset_index(drop = True)
        self.dataset = dataset.reset_index(drop = True)

        self.filtered_dataset = None
        self.labels = None
        self.models = {}
        self.grids = {}
        self.filtering_info = []

    def set_labels(self, label):
        '''Label requires a column name in the metadata. This will be perceived as the classes during modelling'''

        if label not in self.metadata.columns:
            raise Exception(f"Label '{label}' not in metadata. Metadata has the following columns:\n{self.metadata.columns}")

        assert label in self.metadata.columns

        self.labels = self.metadata.loc[:, label]

    def filter_dataset(self, method: str, dataset = "base", percentage_reoccurence = None, drop_missing_protein_columns = False, upperbound = None, lowerbound = None, info = True):
        """# Collection of filtering options
        
        filtered dataset is saved in attribute 'filtered_dataset'

        dataset must be either 'base' or 'filtered' to select either the base loaded dataset or the already filtered dataset to filter.
        
        method = the manner of filtering

        options: ('reoccur global', 'reoccur local', 'abundance')

        ## 'reoccur global' or 'reoccur local': 
        Filters the dataset on percentage reoccuring proteins belonging to a class.
        saves a dataset with proteins filtered that do not occur in percentage_reoccurence of the samples corresponding to a class
        
        Required parameters:
        
        - percentage_reoccurence -> how much samples should a protein be in
        - drop_missing_protein_columns -> mutate the dataframe by dropping protein colums with only 0 values
        - info = true -> returns:
            - filtered_dataframe (attribute: filtered_dataset)
            - dataframe containing protein count filtering info per sample (attribute: deletions_per_sample)
            - dictionary of deleted proteins per class (attribute: deletions_per_class)

        global --> perform filtering regardless of protein identifications of the individual samples

        local --> filtering while considering different amounts of protein identifications in the individual samples

        ## 'abundance': Filters dataset based on top or lowest quantified proteins.
        Required parameters:
        - upperbound: [0,1] percentage of top quantified proteins
        - lowerbound: [0,1] percentage of lowest quantified proteins"""
        
        # Select the dataset
        if dataset == "filtered":
            assert self.filtered_dataset != None
            dataset = self.filtered_dataset
        elif dataset == "base":
            dataset = self.dataset
        else: 
            raise Exception("'dataset' must be either 'filtered' or 'base'.")

        # Do the filtering as queried.
        if method == "abundance":
            assert upperbound != None and lowerbound != None
            self.filtered_dataset = quantification_filtering(dataset, upperbound, lowerbound)

        elif method == "reoccur global" or method == "reoccur local":
            self.filtered_dataset, deletions_per_sample, deletions_per_class = reoccurence_filtering(dataset, self.labels, 
                                            percentage_reoccurence, method.split()[1], drop_missing_protein_columns, info)
            if info == True:
                self.deletions_per_sample = deletions_per_sample
                self.deletions_per_class = deletions_per_class
        else:
            raise Exception('Method of filtering was ill defined.')

        print(f"Filtered dataset on {method} stored in 'filtered_dataset'")        
        self.filtering_info += [method]

    def initialize_models(self, model):
        """Pass a dictionary of models where the key is the name of the algorithm and the value is the model function        
        
        Remembers previously initialized models

        Overwrites the previous initialized models"""
        for model_name, model_function in model.items():
            self.models[model_name] = model_function
        print(f"Models initialized.\nModels: {self.models.keys()}")

    def initialize_grid(self, grid):
        """Pass a dictionary of grids where the key is the name of the algorithm and the value is the dictionary being the grid
        
        Remembers previously initialized grids

        Overwrites the previous initialized grids"""
        
        for key, values in grid.items():
            self.grids[key] = values
        print("Grids initialized.")

        for key in self.models.keys():
            if key not in self.grids.keys():
                print(f"Warning: The model for '{key}' grid has not been initialized yet.")

    def scoring_functions(self, Y_pred, Y_test, labels, project_ids = None):
        """Returns pxd_f1, micro_f1, macro_f1, cm"""
        
        pxd_f1 = project_F1_score(Y_test, Y_pred, project_ids, labels)
        micro_f1 = f1_score(Y_test, Y_pred, labels=labels, average = "micro")
        macro_f1 = f1_score(Y_test, Y_pred, labels=labels, average = "macro")
        cm = confusion_matrix(Y_test, Y_pred, labels = labels)
        return pxd_f1, micro_f1, macro_f1, cm

    def return_feature_importance(self): # Use separate class to interpret/evaluate metrics and models 
        pass

    def nested_grid_search(self, filtered: bool, model_names: list, splitting_procedure: set, gridsearch_scoring: str, name_results_file: str,
                    outerloop_scoring: Tuple[str] = None, outerloop_cv = 5, innerloop_cv = 5):
        '''Nested cross validation grid search on dataset

        - model_names: The models that should be used with the grids that were initiated earlier.
        ---
        - splitting procedure (str,str): the way the test and train set should be split in outer and innerloop respectively
            - provide splitting class itself from sklearn, if second split == "grid", normal gridsearch is done
            - project: on project level -> ensures 3 conditions
                - the test set is < 10 % of the dataset
                - the train set must contain at least 2 projects
                - each class must contain > 10 samples
        ---
        - gridsearch_scoring: scoring function that should be used during gridsearch to optimize hyperparameters
            - project_norm: costum scoring function that first calculates the percentage of well predicted samples in 1 project 
            and then takes an average of that score across all tested projects. 
            This prevents one project with many samples to dominate the end score as these
            will most likely be very similar and thus will all predict either right or wrong.
            NOTE! Provide smaller grid as this gridsearch is not as optimized as the scikit version which is incompatible with the scoring function
            - any scoring function found on https://scikit-learn.org/stable/modules/model_evaluation.html that 
            is compatible with GridSearchCV
        '''
        # Make sure the models and grids queried are initialized
        for model in model_names:
            if model not in self.models.keys() or model not in self.grids.keys():
                raise Exception(f"Initialize {model} first.")

        # Define the dataset that will be used for gridsearch
        data = self.dataset
        labels = self.labels
        if filtered:
            try:
                data = self.filtered_dataset
            except:
                raise Exception("No filter was applied yet. Call 'filter_dataset' to generate a filtered dataset.")

        # Initiate the splitters for the nested cross validation
        if "project" == splitting_procedure[0]:
            # Make sure everything is added to start project based split
            outer_splitter = ProjectBasedSplit(outerloop_cv, self.metadata, "tissue_type")
        else:
            outer_splitter = splitting_procedure[0]
            outer_splitter.n_splits = outerloop_cv


        # Iterate over every model queried
        for model_name in model_names:

            clf = self.models[model_name]
            grid = self.grids[model_name]

            outerloop_count = 1
            print(f"Start evaluation: {model_name}")

            # Start outer loop
            for train, test in outer_splitter.split(data, labels):
                print(f"CV: {outerloop_count}/10")

                # Split data in train + validation and test set.
                X_train = data.iloc[train,:].reset_index(drop = True)
                Y_train = labels.iloc[train].reset_index(drop = True)
                train_project_labels = self.metadata.iloc[train,:].PXD_accession.reset_index(drop = True)
                
                X_test = data.iloc[test,:].reset_index(drop = True)
                Y_test = labels.iloc[test].reset_index(drop = True)
                test_project_labels = self.metadata.iloc[test,:].PXD_accession.reset_index(drop = True)


                if "project" == splitting_procedure[1]:
                    inner_splitter = ProjectBasedSplit(innerloop_cv, metadata = self.metadata.iloc[train,:].reset_index(), on = "tissue_type")
                else:
                    inner_splitter = splitting_procedure[1]
                    inner_splitter.n_splits = innerloop_cv

                # Start Gridsearch to finetune hyperparameters
                # project_norm scoring uses extra info on project level, which would be incompatible with gridsearch
                if gridsearch_scoring != 'project_norm':
                    gridsearch = GridSearchCV(clf, param_grid=grid, scoring = gridsearch_scoring, cv = inner_splitter)
                    gridsearch.fit(X_train, Y_train)

                else:
                    gridsearch = GridSearchProjectF1(clf, grid = grid, cv = inner_splitter)
                    gridsearch.fit(X_train, Y_train, train_project_labels)

                Y_pred = gridsearch.predict(X_test)

                # Scores are calculated in different class. The Y_test and Y_pred is instead written to a file as are the labels
                pxd_f1, micro_f1, macro_f1, cm = self.scoring_functions(Y_pred, Y_test, labels.unique(), test_project_labels)

                # Save metrics and model to file
                results_df = pd.DataFrame({"model": model_name,"fold" : outerloop_count, "pxd_f1": pxd_f1, "micro_f1": micro_f1,
                "macro_f1": macro_f1, "cm": [cm], "parameters": [gridsearch.best_params_]})
                
                self._save_results(results_df, f"gs_{name_results_file}")
                outerloop_count += 1
        
        print(f"Grid search done. See the results in results/gs_{name_results_file}.csv")

    def compare_models(self, dataset_type: str, models: list, splitting_procedure: StratifiedKFold, name_results_file: str, filter_trainset_only = False, percentage_reoccurence: int = None):
        """Models: Algorithms to be used provided in a list, with all parameters provided
        
        if models = "init", use the initialized models instead.

        dataset type: "filtered" or "base"

        splitting_procedure: accepts a splitting class like StratifiedKFold"""
        
        if models == "init":
            models = self.models.values()
        
        if dataset_type == "filtered":
            dataset = self.filtered_dataset
        elif dataset_type == "base":
            dataset = self.dataset
        else:
            raise Exception("Dataset type was not well specified.")

        # Loop over the models
        for clf in models:
            fold = 0
            print(type(clf).__name__)
            for train, test in splitting_procedure.split(dataset, self.labels):
                fold += 1
                print(fold)
                X_train = dataset.iloc[train,:].reset_index(drop = True)
                Y_train = self.labels.iloc[train].reset_index(drop = True)
                train_project_labels = self.metadata.iloc[train,:].PXD_accession.reset_index(drop = True)
                
                X_test = dataset.iloc[test,:].reset_index(drop = True)
                Y_test = self.labels.iloc[test].reset_index(drop = True)
                test_project_labels = self.metadata.iloc[test,:].PXD_accession.reset_index(drop = True)

                if filter_trainset_only:
                    if percentage_reoccurence == None:
                        raise Exception("Define percentage of reoccurence for filtering")
                    
                    X_train, _, _ = reoccurence_filtering(X_train, Y_train, percentage_reoccurence=percentage_reoccurence, how = "global",
                                drop_missing_protein_columns=True)
                    X_test = X_test.loc[:, X_test.columns.isin(X_train.columns)]

                clf.fit(X_train, Y_train)

                Y_pred = clf.predict(X_test)

                pxd_f1, micro_f1, macro_f1, cm = self.scoring_functions(Y_pred, Y_test, self.labels.unique(), test_project_labels)
                
                results_df = pd.DataFrame({"model": type(clf).__name__, "fold": fold, "pxd_f1": pxd_f1, "micro_f1": micro_f1,
                                            "macro_f1": macro_f1, "cm": [cm]})

                self._save_results(results_df, f"comparison_{name_results_file}")

            print(f"Comparing results done. See the results in results/comparison_{name_results_file}.csv")



    def _save_results(self, df: pd.DataFrame, name_file: str):
        """Writes results away to a file"""
        try:
            file = pd.read_csv(f"results/{name_file}.csv", sep = ";")
            df = pd.concat([file, df], ignore_index = True)
            df.to_csv(f"results/{name_file}.csv", sep = ";", index=False)
        except:
            df.to_csv(f"results/{name_file}.csv", sep = ";", index=False)

