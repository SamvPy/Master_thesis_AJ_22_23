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
from sklearn.model_selection import StratifiedKFold


# model architecture and train function

class VariationalEncoder(nn.Module):
    def __init__(self, latent_dims, input_shape):
        super(VariationalEncoder, self).__init__()

        self.linear1 = nn.Linear(input_shape, 500)
        self.linear2 = nn.Linear(500, latent_dims)
        self.linear3 = nn.Linear(500, latent_dims)

        self.N = torch.distributions.Normal(0,1)
        self.kl = 0

        self.dropout = nn.Dropout(.15)

    def forward(self, x):
        x = torch.flatten(x, start_dim = 1)
        x = self.linear1(x)
        x = F.relu(x)
        x = self.dropout(x)

        mu = self.linear2(x)
        sigma = torch.exp(self.linear3(x))

        z = mu + sigma * self.N.sample(mu.shape)
        self.kl = (sigma**2 + mu**2 - torch.log(sigma) - 1/2).sum()

        return z

class Decoder(nn.Module):
    def __init__(self, latent_dims, input_shape):
        super(Decoder, self).__init__()
        
        self.linear1 = nn.Linear(latent_dims, 500)
        self.linear2 = nn.Linear(500, input_shape)
        self.dropout = nn.Dropout(.15)
        
    def forward(self, z):
        z = self.linear1(z)
        z = F.relu(z)
        z = self.dropout(z)

        z = self.linear2(z)
        z = torch.sigmoid(z)

        return z
    
class VariationalAutoencoder(nn.Module):
    def __init__(self, latent_dims, input_shape):
        super(VariationalAutoencoder, self).__init__()

        self.encoder = VariationalEncoder(latent_dims, input_shape)
        self.decoder = Decoder(latent_dims, input_shape)

    def forward(self, x):
        z = self.encoder(x)
        reconstruction = self.decoder(z)

        return reconstruction

def trainDVAE2(autoencoder, data, validation_data=False, epochs = 20):
    optimization = torch.optim.Adam(autoencoder.parameters())

    training_losses = []
    validation_losses = []

    running_loss = 0
    last_loss = 0

    for epoch in range(epochs):

        for label, (x_1, x_2) in data:

            x_1 = x_1.to(device)
            optimization.zero_grad()

            x1_hat = autoencoder(x_1)
            # x2_hat = autoencoder(x_2)

            reconstruction_error = ((x_1-x1_hat)**2).sum() # Reconstructed datapoint must be similar
            intra_class_difference_error = ((x_2-x1_hat)**2).sum() # Incentivize class identity reconstruction by using other samples from the same class
            kullback_leibler = autoencoder.encoder.kl # Error term to ensure a smooth latent space

            loss = .6* reconstruction_error + .4* intra_class_difference_error + kullback_leibler
            loss.backward()
            optimization.step()

            running_loss += loss.item()

        last_loss = running_loss / len(data)
        #print('  epoch {} loss: {}'.format(epoch, last_loss))
        training_losses.append(last_loss)
        
        running_loss = 0
        valid_loss = 0

        print(f"\t Epoch: ({epoch}/{epochs}); Loss: {last_loss}")

        if validation_data == False:
            continue

        for label, (x_1, x_2) in validation_data:

            x1_hat = autoencoder(x_1)
            x2_hat = autoencoder(x_2)
            
            reconstruction_error = ((x_1-x1_hat)**2).sum() 
            intra_class_difference_error = ((x1_hat-x2_hat)**2).sum()
            kullback_leibler = autoencoder.encoder.kl

            loss = reconstruction_error + intra_class_difference_error + 2*kullback_leibler

            valid_loss += loss.item()
        
        valid_loss = valid_loss/len(validation_data)
        print('     Validation loss: {}\n'.format(valid_loss))
        validation_losses.append(valid_loss)

    return autoencoder, training_losses, validation_losses


# Utility functions for VAE
def decodeBatch(vae, encodings):
    '''Return reconstructed samples given some encodings'''
    
    encodings = torch.utils.data.DataLoader(encodings, batch_size=1, shuffle=False)

    flag = True
    for encoding in encodings:
        if flag:
            reconstructed = vae.decoder(encoding).detach().numpy()[0]
            flag = False
            continue
        reconstructed = np.vstack([reconstructed, vae.decoder(encoding).detach().numpy()[0]])

    return reconstructed


def useVAE(vae, data):
    """
    Returns encodings and reconstructed data

    Accepts a numpy array as data and an autoencoder
    """
    input_data = torch.utils.data.DataLoader(data, batch_size=1, shuffle=False)

    flag = True
    for datapoint in input_data:
        
        encoding = vae.encoder(datapoint)
        if flag:
            encodings = encoding.detach().numpy()[0]
            reconstructed = vae.decoder(encoding).detach().numpy()[0]
            flag = False

            continue

        encodings = np.vstack([encodings, encoding.detach().numpy()[0]])
        reconstructed = np.vstack([reconstructed, vae.decoder(encoding).detach().numpy()[0]])

    return encodings, reconstructed

class PairedDataset(torch.utils.data.Dataset):
    def __init__(self, dataset, targets):
        self.dataset = dataset
        self.targets = targets

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, index):
        image1 = self.dataset[index]
        label1 = self.targets[index]
        
        while True:
            index2 = torch.randint(len(self.dataset), size=(1,)).item()
            image2, label2 = self.dataset[index2], self.targets[index2]
            if label2 == label1 and index != index2:
                break

        return label1, (image1, image2)


def main():
    # Machine learning classification loop
    # read data

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
    
    svc_grid = {'decision_function_shape': ["ovr"],
                "kernel": ['linear', 'poly', 'rbf'],
                "C": loguniform(.001, 1e2),
                "class_weight": [weights,None]}
    rf_grid = {'n_estimators': np.linspace(10, 200, 100, dtype = int),
                "criterion": ["entropy", "gini"], 
                "max_depth": [5,10,15,20,40, None],
                "class_weight": [weights,None]}
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
    
    model_grids = [(LogisticRegression(max_iter=10000), lr_l1_grid), 
                   (LogisticRegression(max_iter=10000), lr_l2_grid),
                   (SVC(), svc_grid), 
                   (RandomForestClassifier(), rf_grid),
                   ]

    preprocessor = Pipeline(steps=[
        ('filtering', uml.FilterByOccurence(.5)),
        ('imputation', uml.LowestValueImputer()),
        ('scaler', MinMaxScaler()),
    ])

    skf = StratifiedKFold(n_splits=5, shuffle=True)

    fold=0
    for train, test in skf.split(X=data, y=targets):
        fold += 1
        print(fold)
        # Split data
        X_train = data.iloc[train,:]
        Y_train = targets[train]
        X_test = data.iloc[test,:]
        Y_test = targets[test]

        preprocessor.fit(X_train, Y_train)
        X_train_preprocessed = preprocessor.transform(X_train)
        X_test_preprocessed = preprocessor.transform(X_test)

        X_train_preprocessed = X_train_preprocessed.astype('float32')
        X_test_preprocessed = X_test_preprocessed.astype('float32')

        input_shape = X_train_preprocessed.shape[1]

        # Train autoencoder on train data
        train_paired_dataset = PairedDataset(X_train_preprocessed, Y_train)
        train_loader = torch.utils.data.DataLoader(train_paired_dataset, batch_size=10, shuffle=True)

        print('training VAE...')
        vae_test = VariationalAutoencoder(32, input_shape)
        vae_test, _, _ = trainDVAE2(vae_test, train_loader, epochs=50)

        print('Oversample datapoints from latent space and reconstruct...')
        # Oversample on training set latent space
        train_encoded, train_reconstruction = useVAE(vae_test, X_train_preprocessed)
        train_oversample_encoded, y_train_oversample = SMOTETomek().fit_resample(train_encoded, Y_train)

        # Reconstruct datapoints from the oversampled encodings
        train_reconstructed_smote = decodeBatch(vae_test, train_oversample_encoded)

        for (model, grid) in model_grids:

            print("Starting gridsearch...")
            clf = RandomizedSearchCV(model, grid, n_iter = 10, scoring='f1_macro', verbose = 2)
            search = clf.fit(train_reconstructed_smote, y_train_oversample)
            
            # Get prediction for original and reconstructed test samples
            test_encoding, test_reconstructed = useVAE(vae_test, X_test_preprocessed)
            Y_pred_original = search.predict(X_test_preprocessed)
            Y_pred_reconstruction = search.predict(test_reconstructed)

            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=Y_pred_original, 
                                                                        Y_test=Y_test,
                                                                        labels=unique_labels)

            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], "micro_f1": [micro_f1], 'cv_macro_f1': [search.best_score_],
                                            "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm], "prediction": ["original"],
                                            "oversampler": ["VAE/SMOTETomek"], 'best_params': [search.best_params_]})
                
            uml.save_results(results_df, "VAE_evaluation_grid_search") 

            micro_f1, macro_f1, weighted_f1, cm = uml.scoring_functions(Y_pred=Y_pred_reconstruction, 
                                                                        Y_test=Y_test,
                                                                        labels=unique_labels)       
            
            results_df = pd.DataFrame({"model": [type(model).__name__], "fold": [fold], "micro_f1": [micro_f1], 'cv_macro_f1': [search.best_score_],
                                            "macro_f1": [macro_f1], "weighted_f1": [weighted_f1] ,"cm": [cm], "prediction": ["reconstruction"],
                                            "oversampler": ["VAE/SMOTETomek"], 'best_params': [search.best_params_]})
                
            uml.save_results(results_df, "VAE_evaluation_grid_search") 
            print('Saved results')    

if __name__ == '__main__':
    main()