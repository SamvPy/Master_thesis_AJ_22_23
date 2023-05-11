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

data_quantile = data_quantile.reset_index(drop=True).rename(columns={data_quantile.columns[x]:x for x in range(len(data_quantile.columns))})
target_encoder = LabelEncoder()
targets = target_encoder.fit_transform(meta.Group)
unique_labels = pd.Series(targets).unique()

labels = np.array(meta.Group)
#mmi_quant = uml.MNAR_MCAR_Imputer(max_iter=15)
mmi_quant = uml.LowestValueImputerGaussian()
mmi_quant.fit(data_quantile, labels)
imputed_df_quant = mmi_quant.transform(data_quantile, labels)
mm_scaler = MinMaxScaler()
minmax_df_quant = mm_scaler.fit_transform(imputed_df_quant)

hyperparameter_sweep = {'epoch': [],
                        'learning_rate': [],
                        'batch_size': [],
                        'latent_space': [],
                        'train_loss': [],
                        'test_loss': [],
                        'train_loss_kl': [],
                        'train_loss_mse': [],
                        'test_loss_kl': [],
                        'test_loss_mse': []}

def update_hyperparameter_dict(hyperparameter_dict, epochs, lr, bs, ls, train_loss, val_loss):
    
    hyperparameter_dict['epoch'] += list(range(1,epochs+1))
    hyperparameter_dict['learning_rate'] += [lr]*epochs
    hyperparameter_dict['batch_size'] += [bs]*epochs
    hyperparameter_dict['latent_space'] += [ls]*epochs


    avg_train_loss, train_kl, train_mse = [value[0] for value in train_loss.values()], [value[2] for value in train_loss.values()], [value[1] for value in train_loss.values()]
    avg_val_loss, val_kl, val_mse = [value[0] for value in val_loss.values()], [value[2] for value in val_loss.values()], [value[1] for value in val_loss.values()]

    hyperparameter_dict['train_loss'] += avg_train_loss
    hyperparameter_dict['test_loss'] += avg_val_loss
    hyperparameter_dict['train_loss_kl'] += train_kl
    hyperparameter_dict['train_loss_mse'] += train_mse
    hyperparameter_dict['test_loss_kl'] += val_kl
    hyperparameter_dict['test_loss_mse'] += val_mse

    return hyperparameter_dict


minmax_df_quant = minmax_df_quant.astype('float32')
X_train, X_test, y_train, y_test = train_test_split(minmax_df_quant, targets, test_size=0.20, random_state=42, stratify=targets)

total_count = 6*4*5
count=0
for latent_dim in [5, 10, 15, 20, 30, 50]:
    for batch_size in [5, 10, 15, 20]:

        train_loader = torch.utils.data.DataLoader(X_train, batch_size=batch_size, shuffle=True)
        test_loader = torch.utils.data.DataLoader(X_test, batch_size=batch_size, shuffle=True)

        

        for lr in [.0005, .001, .0025, .005, .01]:
            vae = vae_utils.VariationalAutoencoder(latent_dims=latent_dim, num_features=2615, hidden_features=500)
            trained_vae, train_loss, val_loss = vae_utils.trainVAE(vae, train_loader, test_loader, epochs = 100, lr = lr, verbose=True)

            hyperparameter_sweep = update_hyperparameter_dict(hyperparameter_sweep, 100,
                                                              lr, batch_size, latent_dim, train_loss, val_loss)

            count+= 1
            print(f"{count}/{total_count}")

pd.DataFrame(hyperparameter_sweep).to_csv("/home/compomics/Sam/git/python/master_thesis/Atlas_analysis/ML/results/VAE_hyperparameter_opt_LWVGI.csv")