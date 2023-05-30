# Description of files and notebooks

## Python script
* **master_functions.py**<br>
A bundle of functions used to populate the mySQL database as well parsing the ionbot result files.

* **utils_ML.py** <br>
In this folder only the FilterByOccurence function is used: Filtering proteins that occur in >x of samples

## Notebooks

* **NSAF_atlas.ipynb**<br>
Proteins from each RAW-file are quantified with the Normalised Spectral Abundance Factor (NSAF). Also includes pre-processing steps such as filtering proteotypic peptides and contaminant removal.
Quantification with pooled and unpooled strategy are described here.

* **PoolingEvaluation.ipynb** <br>
Contains an evaluation of the combining fractions (pooling) strategy.

## Files
* **pooling.log** <br>
Contains part of the fraction identifier annotation

* **proteome_nsaf_3.h5**<br>
The NSAF-atlas without pooling

* **proteome_nsaf_pooled_3.h5**<br>
The NSAF-atlas with the pooling strategy

* **CRAP.tsv**<br>
Contains proteins that are likely contaminants as described in "https://www.thegpm.org/crap/"