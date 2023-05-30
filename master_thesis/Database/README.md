# Description of files and notebooks

## Python script
* **master_functions.py**<br>
A bundle of functions used to populate the mySQL database as well parsing the ionbot result files.

## Notebooks
* **database_parser.ipynb** <br>
Notebook containing the scripts to populate the mySQL database.

* **ionbot_version_comparison.ipynb**<br>
Notebook exploring the difference between ionbot versions. This is thought to be mainly caused by differences in the format of the result files used between ionbot versions. Indeed, some ionbot result file parsers are probably not equally adequate for all versions and should be investigated further.

* **summary_stats_data.ipynb** <br>
Notebook exploring differences in PSMs, peptide identifications and protein identifications for each project

## Files
* **protein_table_uniprot.csv** <br>
File used to fill in the protein table of the mySQL database.

* **CRAP.tsv**<br>
Contains proteins that are likely contaminants as described in "https://www.thegpm.org/crap/"

* **ionbot_assays.log**, **ionbot_assays1.log**<br>
log-file providing information why certain ionbot result files are not parsed and number of proteins extracted when succesfully parsed.

* **ionbot_version_comparison3.json**<br>
Containing comparative statistics when parsing different ionbot version result files coming from the same raw-file

---

NOTE: Most files below are bundled in the *unified_metadata.csv* file.

* **assays_not_found.csv** <br>
Contains the assays that were not present on the compomics server at the time when ionbot results were put in the database.

* **new_assay_file_paths.csv** <br>
the assays in *assays_not_found.csv* were on 13/02/2023 once again searched on the compomics server and the filepaths are added

* **parsed_manual_meta2.csv** <br>
Annotation file with the file paths on the compomics server added.

* **parsedAvailableFiles.csv** **parsedFailedFiles.csv** **parsedNewFiles.csv**<br>
Files with raw-file annotation overwritten with raw-file level statistics (PSM, peptide, protein counts).

* **parser_failed_ae4.csv** <br>
Contains the filenames of the assays that could not be parsed with the ionbot_parser function in *master_functions.py*

* **parser_tested_new_files.csv** <br>
Contains information on which files can be parsed with the ionbot_parser function. These are overwritten with information on PSM, peptide and protein numbers present in the raw-files as performed in *summary_stats_data.ipynb*