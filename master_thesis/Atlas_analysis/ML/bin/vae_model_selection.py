import utils_ML as uml

import pandas as pd
import numpy as np

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler, normalize
from sklearn.utils.class_weight import compute_class_weight
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split

from imblearn.combine import SMOTEENN, SMOTETomek

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions

class VariationalEncoder(nn.Module):
    def __init__(self, latent_dims, num_features, hidden_features):
        super(VariationalEncoder, self).__init__()

        self.linear1 = nn.Linear(num_features, hidden_features)
        self.linear2 = nn.Linear(hidden_features, latent_dims)
        self.linear3 = nn.Linear(hidden_features, latent_dims)

        self.N = torch.distributions.Normal(0,1)
        self.kl = 0

        # Add regularization with dropouts
        #self.dropout = nn.Dropout(0.25) # approximately 250 neurons will be deactivated

    def forward(self, x):
        x = torch.flatten(x, start_dim = 1)
        x = self.linear1(x)
        x = F.relu(x)
        #x = self.dropout(x)

        mu = self.linear2(x)
        sigma = torch.exp(self.linear3(x))

        z = mu + sigma * self.N.sample(mu.shape)
        self.kl = (sigma**2 + mu**2 - torch.log(sigma) - 1/2).sum()

        return z

class Decoder(nn.Module):
    def __init__(self, latent_dims, num_features, hidden_features):
        super(Decoder, self).__init__()
        
        self.linear1 = nn.Linear(latent_dims, hidden_features)
        self.linear2 = nn.Linear(hidden_features, num_features)
        #self.dropout = nn.Dropout(0.25)

    def forward(self, z):
        z = self.linear1(z)
        z = F.relu(z)

        #z = self.dropout(z)

        z = self.linear2(z)
        z = torch.sigmoid(z)

        return z
    
class VariationalAutoencoder(nn.Module):
    def __init__(self, latent_dims, num_features, hidden_features):
        super(VariationalAutoencoder, self).__init__()

        self.encoder = VariationalEncoder(latent_dims, num_features, hidden_features)
        self.decoder = Decoder(latent_dims, num_features, hidden_features)

    def forward(self, x):
        z = self.encoder(x)
        reconstruction = self.decoder(z)

        return reconstruction
    
def trainVAE(autoencoder: VariationalAutoencoder, data, epochs = 20):
    optimization = torch.optim.Adam(autoencoder.parameters())

    running_loss = 0
    last_loss = 0

    for epoch in range(epochs):
        print(epoch)

        for i, x in enumerate(data):
            x = x.to(device)
            optimization.zero_grad()
            x_hat = autoencoder(x)

            loss = ((x-x_hat)**2).sum() + autoencoder.encoder.kl
            loss.backward()
            optimization.step()
        
            running_loss += loss.item()

            if i == 129:
                last_loss = running_loss/129
                print('batch {} loss: {}'.format(i + 1, last_loss))

                running_loss = 0
        print(loss)
    return autoencoder

def trainDVAE2(autoencoder: VariationalAutoencoder, data, validation_data=False, epochs = 20, lr=.001):
    optimization = torch.optim.Adam(autoencoder.parameters(), lr = lr)

    training_losses = []
    validation_losses = []

    running_loss = 0
    last_loss = 0

    for epoch in range(epochs):

        for label, (x_1, x_2) in data:

            x_1 = x_1.to(device)
            optimization.zero_grad()

            x1_hat = autoencoder(x_1)
            x2_hat = autoencoder(x_2)

            reconstruction_error = ((x_1-x1_hat)**2).sum() # Reconstructed datapoint must be similar
            intra_class_difference_error = ((x1_hat-x2_hat)**2).sum() # Incentivize class identity reconstruction by using other samples from the same class
            kullback_leibler = autoencoder.encoder.kl # Error term to ensure a smooth latent space

            loss = reconstruction_error + kullback_leibler
            loss.backward()
            optimization.step()

            running_loss += loss.item()

        last_loss = running_loss / len(data)
        print('  epoch {} loss: {}'.format(epoch, last_loss))
        training_losses.append(last_loss)
        
        running_loss = 0
        valid_loss = 0

        if validation_data == False:
            continue

        for label, (x_1, x_2) in validation_data:

            x1_hat = autoencoder(x_1)
            x2_hat = autoencoder(x_2)
            
            reconstruction_error = ((x_1-x1_hat)**2).sum() 
            intra_class_difference_error = ((x1_hat-x2_hat)**2).sum()
            kullback_leibler = autoencoder.encoder.kl

            loss = reconstruction_error + kullback_leibler

            valid_loss += loss.item()
        
        valid_loss = valid_loss/len(validation_data)
        print('     Validation loss: {}\n'.format(valid_loss))
        validation_losses.append(valid_loss)

    return autoencoder, training_losses, validation_losses

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


device = 'cuda' if torch.cuda.is_available() else "cpu"

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


imputer = uml.MNAR_MCAR_Imputer(max_iter=15)
imputer.fit(data_quantile, targets)
imputed_combat = imputer.transform(data_quantile, targets)

imputed_combat_ = imputed_combat.astype('float32').to_numpy()
X_train, X_test, y_train, y_test = train_test_split(imputed_combat_, targets, test_size=0.20, random_state=42, stratify=targets)

scaler = MinMaxScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)


train_paired_dataset = PairedDataset(X_train, y_train)
test_paired_dataset = PairedDataset(X_test, y_test)

train_loader = torch.utils.data.DataLoader(train_paired_dataset, batch_size=5, shuffle=True)
test_loader = torch.utils.data.DataLoader(test_paired_dataset, batch_size=5, shuffle=True)



latent_dims_k = [10, 30, 50, 75]
hidden_dims_k = [100,300,500,700]
learning_rates = [.0001, .001, .01]

VAEs = {}
results = {}
counts = 0
for latent_dim in latent_dims_k:
    for hidden_dim in hidden_dims_k:
        vae = VariationalAutoencoder(latent_dims = latent_dim, num_features=2615, hidden_features=hidden_dim)
        for lr in learning_rates:
            counts+= 1
            print(counts, ":", latent_dim, hidden_dim, lr)

            fitted_vae, train_loss, validation_loss = trainDVAE2(vae, train_loader, validation_data=test_loader, epochs = 50, lr = lr)
            results[counts] = {"latent_dim": latent_dim, "hidden_dim": hidden_dim, "lr": lr, "train_loss": train_loss, "validation_loss": validation_loss}

            VAEs[counts] = fitted_vae

pd.DataFrame(results).to_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/ML/results/VAE_architecture_search.csv")