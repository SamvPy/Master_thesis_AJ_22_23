# preprocessors
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.model_selection import train_test_split

# utils

import pandas as pd
import numpy as np
import sys

import vae_utils as vae_utils

sys.path.append("../")

import utils_ML as uml

import torch
import torch.utils
import torch.distributions


data_quantile = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/quantile_norm_NSAF_50.csv", index_col = "assay_id")
meta = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/unified_metadata.csv")
meta = meta[meta.assay_id.isin(data_quantile.index)]

groups = pd.read_csv("/home/compomics/Sam/git/python/master_thesis/Metadata/group_cells_annotation.csv", sep =";", index_col="Unnamed: 0")
meta["Group"] = meta.cell_line.apply(lambda x: groups[groups.cell_line == x]["group"].values[0])
meta = meta.set_index("assay_id")

data_quantile.sort_index(inplace=True)
meta.sort_index(inplace=True)

with open("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/preprocessing/selected_features.txt", "r") as f:
    features = f.readlines()
    features = [x.strip() for x in features]

data_quantile = data_quantile.loc[:,features]

data_quantile = data_quantile.reset_index(drop=True).rename(columns={data_quantile.columns[x]:x for x in range(len(data_quantile.columns))})
target_encoder = LabelEncoder()
targets = target_encoder.fit_transform(meta.Group)
unique_labels = pd.Series(targets).unique()

labels = np.array(meta.Group)
#mmi_quant = uml.MNAR_MCAR_Imputer(max_iter=15)
mmi_quant = uml.LowestValueImputerGaussian()
mmi_quant.fit(data_quantile, labels)
imputed_df_quant = mmi_quant.transform(data_quantile, labels)

hyperparameter_sweep = {'epoch': [],
                        'learning_rate': [],
                        'batch_size': [],
                        'latent_space': [],
                        'kappa': [],
                        'train_loss': [],
                        'test_loss': [],
                        'train_loss_kl': [],
                        'train_loss_mse': [],
                        'test_loss_kl': [],
                        'test_loss_mse': []}

def update_hyperparameter_dict(hyperparameter_dict, epochs, lr, bs, ls, kappa, train_loss, val_loss):
    
    hyperparameter_dict['epoch'] += list(range(1,epochs+1))
    hyperparameter_dict['learning_rate'] += [lr]*epochs
    hyperparameter_dict['batch_size'] += [bs]*epochs
    hyperparameter_dict['latent_space'] += [ls]*epochs
    hyperparameter_dict['kappa'] += [kappa]*epochs


    avg_train_loss, train_kl, train_mse = [value[0] for value in train_loss.values()], [value[2] for value in train_loss.values()], [value[1] for value in train_loss.values()]
    avg_val_loss, val_kl, val_mse = [value[0] for value in val_loss.values()], [value[2] for value in val_loss.values()], [value[1] for value in val_loss.values()]

    hyperparameter_dict['train_loss'] += avg_train_loss
    hyperparameter_dict['test_loss'] += avg_val_loss
    hyperparameter_dict['train_loss_kl'] += train_kl
    hyperparameter_dict['train_loss_mse'] += train_mse
    hyperparameter_dict['test_loss_kl'] += val_kl
    hyperparameter_dict['test_loss_mse'] += val_mse

    return hyperparameter_dict

imputed_df_quant = imputed_df_quant.astype("float32")
X_train, X_test, y_train, y_test = train_test_split(imputed_df_quant, targets, test_size=0.20, random_state=42, stratify=targets)

mm_scaler = MinMaxScaler()
X_train = mm_scaler.fit_transform(X_train)
X_test = mm_scaler.transform(X_test)

total_count = 5*2*7*2
count=0

print(f"Starting paramater sweep... ({total_count})")

for latent_dim in [2,5,8,10,15]:
    for batch_size in [5, 10]:

        train_loader = torch.utils.data.DataLoader(X_train, batch_size=batch_size, shuffle=True)
        test_loader = torch.utils.data.DataLoader(X_test, batch_size=batch_size, shuffle=True)

        for lr in [.0001, .001]:
            for kappa in [.001, .01, .05, .1, .3, .5, .75]:
                vae = vae_utils.VariationalAutoencoder(latent_dims=latent_dim, num_features=161, hidden_features=50)
                trained_vae, train_loss, val_loss = vae_utils.trainVAE(vae, train_loader, test_loader, epochs = 100, lr = lr, kappa=kappa, verbose=True)

                hyperparameter_sweep = update_hyperparameter_dict(hyperparameter_sweep, 100,
                                                                lr, batch_size, latent_dim, kappa, train_loss, val_loss)

                count+= 1
                print(f"{count}/{total_count}")

pd.DataFrame(hyperparameter_sweep).to_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/ML/results/VAE_hyperparameter_opt_filtered_LWVGI.csv")