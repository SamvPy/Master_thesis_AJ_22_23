# Description of files and notebooks

## Notebooks
* **check_PXD_reprocessed.ipynb** <br>
Extracts the reprocessed filenames of the projects deemed useable from the compomics internal server and creates annotation files for them which were manually filled out.

* **unify_metadata.ipynb**<br>
Notebook used to aggregate all the disparate annotation files that were created throughout the thesis at multiple timepoints. For the sake of completeness, all the original files are kept and this notebook is shown.


## Files
* **unified_metadata.csv** <br>
Main annotation file, which groups all annotation files described below. For the sake of completeness, all the original annotation files are kept as well as the notebook *unify_metadata.ipynb* which describes how the files are bundled together.

---

* **annotation_excel4.xlsx** <br>
Large annotation file containing the following annotation labels:
    - useability
    - cell line name
    - disease type
    - tissue type
    - treatment
    - sub cell name
    - fraction identifier

* **annotation_excel_use2.csv**,  **file_anntoation_update.csv**, **annotation_PXD018450.xlsx**<br>
Similar as *annotation_excel4.xlsx*  but with other raw files.

* **cellosuaurs_webscraping_filtered_update.csv** <br>
Manual evaluation of useability of the PRIDE projects extracted with the webscraper in the *Scraper* folder.

* **filter_found_pxd.csv** <br>
More PRIDE projects that were found through searching in PRIDE with keywords or by searching ProteomeCentral at "https://proteomecentral.proteomexchange.org/cgi/GetDataset"

* **project_metadata.csv** <br>
Project level annotation of the following aspects
    - Digestion method
    - Instrument name
    - PubMed identifier

* **Group_cells_annotation.csv** <br>
The used *cell_line* annotation from previous annotation files is expanded to allow grouping. This grouping is later used as a *class* in the machine learning workflow. The file contains the following rows:
    - cell_line annotation
    - count: number of files of this cell line
    - description
    - group: the class name