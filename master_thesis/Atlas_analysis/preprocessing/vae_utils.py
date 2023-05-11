import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils
import torch.distributions
import numpy as np
device = 'cuda' if torch.cuda.is_available() else "cpu"

class VariationalEncoder(nn.Module):
    def __init__(self, latent_dims, num_features, hidden_features):
        super(VariationalEncoder, self).__init__()

        self.linear1 = nn.Linear(num_features, hidden_features)
        self.linear2_mu = nn.Linear(hidden_features, latent_dims)
        self.linear2_log_var = nn.Linear(hidden_features, latent_dims)

        self.N = torch.distributions.Normal(0,1)
        self.kl = 0

        # Add regularization with dropouts
        #self.dropout = nn.Dropout(0.25) # approximately 250 neurons will be deactivated

    def forward(self, x):
        x = torch.flatten(x, start_dim = 1)
        x = self.linear1(x)
        x = F.relu(x)
        #x = self.dropout(x)

        mu = self.linear2_mu(x)
        log_var = self.linear2_log_var(x)

        return mu, log_var
        

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
        self.N = torch.distributions.Normal(0,1)

    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5*log_var)
        eps = self.N.sample(mu.shape)
        z = mu + std*eps
        return z

    def forward(self, x):
        mu, log_var = self.encoder(x)
        z = self.reparameterize(mu, log_var)
        reconstruction = self.decoder(z)

        return (reconstruction, mu, log_var)

    def get_encoding_reconstruction(self, x):
        mu, log_var = self.encoder(x)
        z = self.reparameterize(mu, log_var)
        reconstruction = self.decoder(z)

        return (z, reconstruction)
    
def loss_function_annealing(x, reconstruction, mu, log_var, kappa, annealing, epoch):
    
    mse = F.mse_loss(reconstruction, x, reduction="sum") / len(x)

    kl = (-.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())) / len(x)

    if annealing:
        (np.array([np.log10(x)/2*1 if np.log10(x)/2*1<1 else 1 for x  in range(1,500)]) ==1).sum()
        if epoch == 0:
            weight = 0
        else:
            weight = np.log10(epoch)/2 * kappa
            if weight > kappa:
                weight = kappa

    else:
        weight = kappa
    
    kl = kl*weight
    total_loss = mse + kl

    return total_loss, mse, kl

def loss_function(x, reconstruction, mu, log_var):
    
    mse = F.mse_loss(reconstruction, x, reduction="sum") / len(x)

    kl = (-.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())) / len(x)
    

    total_loss = mse + kl

    return total_loss, mse, kl

def trainVAE(autoencoder: VariationalAutoencoder, data, validation_data=False, epochs = 20, lr=.001, kappa=1, annealing=False, verbose=False):
    optimization = torch.optim.Adam(autoencoder.parameters())

    # Stores losses per epoch [total, mse, kl]
    train_losses_dict = {}
    valid_losses_dict = {}

    for epoch in range(epochs):

        running_loss = 0
        mse_loss = 0
        kl_loss = 0
        valid_loss = 0
        valid_mse_loss = 0
        valid_kl_loss = 0

        for i, x in enumerate(data):

            x = x.to(device)
            optimization.zero_grad()
            x_hat, mu, log_var = autoencoder(x)

            if annealing:
                total_loss, mse, kl = loss_function_annealing(x, x_hat, mu, log_var, kappa, annealing, epoch)
            else:
                total_loss, mse, kl = loss_function(x, x_hat, mu, log_var)
            
            loss = mse+kappa*kl
            loss.backward()
            optimization.step()
        
            # Print statistics
            running_loss += total_loss.item()
            mse_loss += mse.item()
            kl_loss += kl.item()

        
        avg_loss = running_loss/(i+1)
        avg_mse_loss = mse_loss/(i+1)
        avg_kl_loss = kl_loss/(i+1)
  
        train_losses_dict[epoch] = [avg_loss, avg_mse_loss, avg_kl_loss]

        if verbose:
            print('Epoch {} loss: {}'.format(epoch+1, avg_loss))
            print("\tMSE: {}\tKL: {}".format(avg_mse_loss, avg_kl_loss))

        if validation_data == False:
            continue
    
        for i, x in enumerate(data):

            x_hat, mu, log_var = autoencoder(x)
            
            if annealing:
                total_loss, mse, kl = loss_function_annealing(x, x_hat, mu, log_var, kappa, annealing, epoch)
            else:
                total_loss, mse, kl = loss_function(x, x_hat, mu, log_var)


            valid_loss += total_loss.item()
            valid_mse_loss += mse.item()
            valid_kl_loss += kl.item()
        
        
        valid_loss = valid_loss/(i+1)
        valid_mse_loss = valid_mse_loss/(i+1)
        valid_kl_loss = valid_kl_loss/(i+1)
        
        if verbose:
            print("(VAL)\tMSE: {}\tKL: {}\t({})\n".format(valid_mse_loss, valid_kl_loss, valid_loss))
        
        valid_losses_dict[epoch] = [valid_loss, valid_mse_loss, valid_kl_loss]

    return autoencoder, valid_losses_dict, train_losses_dict

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

            loss = (reconstruction_error + intra_class_difference_error)/2 + kullback_leibler
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

            loss = (reconstruction_error + intra_class_difference_error)/2 + kullback_leibler

            valid_loss += loss.item()
        
        valid_loss = valid_loss/len(validation_data)
        print('     Validation loss: {}\n'.format(valid_loss))
        validation_losses.append(valid_loss)

    return autoencoder, training_losses, validation_losses

def useVAE(vae, data):
    """
    Returns encodings and reconstructed data

    Accepts a numpy array as data and an autoencoder
    """
    print(f"Reconstructing {len(data)} datapoints...")
    input_data = torch.utils.data.DataLoader(data, batch_size=1, shuffle=False)

    flag = True
    reconstruction_cycle=0
    for datapoint in input_data:
        reconstruction_cycle += 1
        if reconstruction_cycle %500 == 0:
            print(f"Reconstructed {reconstruction_cycle}/{len(data)}")

        encoding, reconstruction = vae.get_encoding_reconstruction(datapoint)
        if flag:
            encodings = encoding.detach().numpy()[0]
            reconstructed = reconstruction.detach().numpy()[0]
            flag = False
            continue
        
        encodings = np.vstack([encodings, encoding.detach().numpy()[0]])
        reconstructed = np.vstack([reconstructed, reconstruction.detach().numpy()[0]])

    return encodings, reconstructed

def resampleVAE(vae: VariationalAutoencoder, data: np.array, labels: np.array, n:int):
    
    # Count class amounts:
    class_count = {}
    for label in np.unique(labels):
        class_count[label] = np.sum(labels==label)
    
    # Get the sample count goal for all classes
    class_count_goal = max(class_count.values())*n
    to_resample = {"label": [], "data": []}
    
    flag = True
    
    # Resample by passing random samples through VAE
    cycle=0
    while flag:
        cycle+=1
        if cycle % 10 == 0:
            print(f"Minority: {min(class_count.values())}; Goal: {class_count_goal}")
            
        flag=False
        indices = list(range(0,len(labels)))
        np.random.shuffle(indices)

        for i in indices:
            label = labels[i]
            datapoint = data[i]
            if class_count[label] < class_count_goal:
                flag=True
                to_resample["label"].append(label)
                to_resample["data"].append(datapoint)
                class_count[label] += 1

    _, reconstructions = useVAE(vae, np.array(to_resample["data"]))
    print("Resampled from VAE")
    return reconstructions, to_resample["label"]