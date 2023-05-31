# Description of notebooks and files

## Python file

* **vae_utils.py**<br>
Contains classes implementing the Variational Autoencoder network and functions to train, use and oversample from it.

## Notebooks

* **normalization_methods.ipynb**<br>
This notebook describes several aspects:
    - Filtering sparse samples
    - Explore the characteristics of the missingness
    - Evaluate whether this affects the NSAF
    - Try 3 additional normalisation methods
    - Check which proteins are least variable
    - Evaluate whether variability and group-specific identification are related

* **imputation.ipynb**<br>
This notebook explores multiple methods of imputation

* **feature_selection.ipynb** <br>
This notebook evaluates five feature selectors and performs a correlation clustering analysis on the feature-selected proteins

* **Oversampling** <br>
This notebook evaluates four oversampling methods

* **VAE_optimization.ipynb**<br>
In this notebook, the optimal parameters of the VAE are chosen based on a small hyperparameter search. Additionally, the latent distributions are plotted as well as pairplots of the encodings

## Files

* **\*_NSAF_50.csv**<br>
Normalised datasets

* **enrichment.all.tsv**<br>
STRING-enrichment file from the interaction network of the least variable proteins, as described at the end of *normalization_methods.ipynb*

* **selected_features.txt**<br>
Features selected in *feature_selection.ipynb*

* **VAE_model(_filtered)**<br>
The trained VAE models from *VAE_optimization.ipynb*

* **THPA_\*.csv** <br>
Downloaded from the Human Protein Atlas (https://www.proteinatlas.org/about/download)