import utils_ML as uml

import pandas as pd
import numpy as np

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler, normalize
from sklearn.utils.class_weight import compute_class_weight
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split, sRandomizedSearchCV
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



splitter = uml.ProjectBasedSplit(10, metadata=meta, on="Group")

fold=0

vaes = []
X_train_pxds = []
X_test_pxds = []

Y_train_pxds = []
Y_test_pxds = []

train_losses = []
test_losses = []

for train, test in splitter.split(data, None):
    
    fold += 1
    print(fold)
    # Split data
    X_train_pxd = data.iloc[train,:]
    Y_train_pxd = targets[train]
    X_test_pxd = data.iloc[test,:]
    Y_test_pxd = targets[test]

    preprocessor.fit(X_train_pxd, Y_train_pxd)
    X_train_preprocessed_pxd = preprocessor.transform(X_train_pxd)
    X_test_preprocessed_pxd = preprocessor.transform(X_test_pxd)

    X_train_preprocessed_pxd = X_train_preprocessed_pxd.astype('float32')
    X_test_preprocessed_pxd = X_test_preprocessed_pxd.astype('float32')

    train_paired_dataset = PairedDataset(X_train, y_train)
    test_paired_dataset = PairedDataset(X_test, y_test)

    train_loader = torch.utils.data.DataLoader(train_paired_dataset, batch_size=10, shuffle=True)
    test_loader = torch.utils.data.DataLoader(test_paired_dataset, batch_size=5, shuffle=True)

    input_shape = X_train_preprocessed_pxd.shape[1]
    vae = VariationalAutoencoder(latent_dims=50, num_features=input_shape).to(device)

    pxd_denoising_vae, train_loss_pxd, test_loss_pxd = trainDVAE2(vae, train_loader, validation_data=test_loader, epochs = 100)

    vaes.append(pxd_denoising_vae)
    train_losses.append(train_loss_pxd)
    test_losses.append(test_loss_pxd)

    X_train_pxds.append(X_train_pxd)
    X_test_pxds.append(X_test_pxd)

    Y_train_pxds.append(Y_train_pxd)
    Y_test_pxds.append(Y_test_pxd)