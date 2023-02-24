import pandas as pd
from collections import defaultdict
import mysql.connector
import glob
import os
import csv
import numpy as np
import re
import logging
class dbf():
    '''A selection of functions used to maintain the expression_atlas_cells database.
    '''

    def __init__(self):
        self.conn = mysql.connector.connect(user='root', password='password', host='127.0.0.1', port='3306',
                               database='expression_atlas_cells')
        self.mycursor = self.conn.cursor(buffered=True)

        # check the connection
        if self.conn.is_connected():
            pass
            print("connection succesfull")
        else:
            print("no connection")
        

    def check_new_project(self, new_projects: list):
        '''Accept a list of PXD-accession codes and checks whether they are already registered in the database. Returns "[]" if no duplicates have been found.'''
        for x in new_projects:
            if x[0:3] != "PXD":
                print("new_projects must be list of PXDxxxxxx")
                return False
                
        query = "SELECT PXD_accession FROM project"
        old_projects = pd.read_sql(query, self.conn).PXD_accession.tolist()
        
        duplicates = []
        for project in new_projects:
            if project in old_projects:
                print(project, "is already in the database")
                duplicates.append(project)
        
        print("Projects checked.")
        return duplicates

    def input_cell_info(self, excel_annotation_file: pd.DataFrame):
        '''Inputs disease and tissue_type info for given cell lines if they are already represented in the sql database. Returns annotation_file'''
        query = "SELECT cell_line, disease, tissue_type, treatment, sub_cell FROM cell"
        db_df = pd.read_sql(query, self.conn)

        for index, row in db_df.iterrows():
            excel_annotation_file.loc[excel_annotation_file.cell_line == row["cell_line"], ["disease", "tissue_type"]] = (row["disease"], row["tissue_type"])
        
        no_entry = excel_annotation_file.loc[excel_annotation_file.disease == '/', "cell_line"].unique().tolist()
        print("Following cell lines still need annotation: ", no_entry)
        return excel_annotation_file

    def find_pxd_path(self, pxds: list, return_path = False):
        '''Given a list of pxds, this function searches the compomics directories for the pxd directory and returns the pxds it finds
        
        If return path is set to True, returns (x,y,z)
        
        - x = list of pxd that were found
        - y = dictionary of path: amount of files
        - z = ionbot version paths'''
        
        path_found = []
        version_path = []
        file_output = {}
        pxd_found = []
        for pxd in pxds:
            for path in glob.glob("/home/compomics/mounts/*/*/PRIDE_DATA/" + str(pxd)):
                if pxd not in pxd_found:
                    pxd_found.append(pxd)
                path_found.append(path)
        
        flag = False
        for path in path_found:
            count = 0
            for version in glob.glob(path + '/IONBOT_v*/*'):
                if count == 0:
                    version_path.append("/".join(version.split('/')[0:-1]))
                count += 1
            
            file_output[path] = count

        print(f"Found {len(pxd_found)} out of {len(pxds)}.")
        if not return_path:
            return pxd_found
        return pxd_found, file_output, version_path

    def find_file_path(self, pxd, raw):
        '''Used in pd.apply as function to return the path of a given RAW-file. If more paths are present, returns the string "multiple_paths"'''

        files = []
        
        for file in glob.glob("/home/compomics/mounts/*/*/PRIDE_DATA/" + str(pxd) + '/IONBOT_v*/*/' + str(raw) + '.mgf.gzip.ionbot.csv'):
            files.append(file)
        
        if files == []:
            for file in glob.glob("/home/compomics/mounts/*/*/PRIDE_DATA/" + str(pxd) + '/IONBOT_v*/' + str(raw) + '.mgf.ionbot.csv'):
                files.append(file)

        if len(files) == 1:
            return files[0]
        
        if len(files) > 1:
            versions = [file.split("/")[8] for file in files]
            versions = [[int(version.split('.')[1]), int(version.split('.')[2])] for version in versions]
            best = [0,0]
            for version in versions:
                if version[0] > best[0]:
                    best[0], best[1] = version[0], version[1]

                if version[0] == best[0]:
                    if version[1] > best[1]:
                        best[0], best[1] = version[0], version[1]
            return files[versions.index(best)]
            
        return np.nan

    def build_project_table(self, meta_df, list_of_pxds):
        '''meta_df must be a dataframe with columns: accession, digestion, instrumentNames, PMID.
        
        list_of_pxds are the pxds that will be added to the table.'''

        count = 0
        meta_df = meta_df[meta_df['accession'].isin(list_of_pxds)]
        
        check = self.check_new_project(meta_df.accession.unique().tolist())
        if check != []:
            print(f"Duplicates detected: {check}. \nNo entries have been added")
            return

        meta_df = meta_df[['accession', 'digestion', 'instrumentNames', 'PMID']]
        meta_df = meta_df.astype(str)
        meta_tuples = list(meta_df.to_records(index=False)) #a list of tuples is easily iteratible and easy to store in the database
        for i in meta_tuples:
            count += 1
            project = "INSERT INTO project(PXD_accession, experiment_type, instrument, pmid) VALUES (%s, %s, %s, %s)"
            i = list(i)
            self.mycursor.execute(project, i)
            self.conn.commit()
            
        print(f"{count} projects added in table 'project'.")
            
    def build_cell_table(self, cell_df: pd.DataFrame):
        '''Inserting non-duplicate entries in the cell database. 
        
        cell_df must have following columns:
        - cell_line
        - disease
        - tissue_type
        - treatment
        - sub_cell'''

        cell_df = cell_df["cell_line disease tissue_type treatment sub_cell".split()]
        cell_df = cell_df.drop_duplicates()
        print(f"{len(cell_df)} cell entries in file.")
        count = 0
        cell_tuples = list(cell_df.to_records(index=False)) #a list of tuples is easily iteratible and easy to store in the database
        for i in cell_tuples:
            count += 1
            cell = "INSERT IGNORE INTO cell(cell_line, disease, tissue_type, treatment, sub_cell) VALUES (%s, %s, %s, %s, %s)"
            i = list(i)
            self.mycursor.execute(cell, i)
            self.conn.commit()
        print(f"{count} entries added in table 'cell'.")

    def build_assay_cell_table(self, assay_df):
        '''assay_df must have columns of PXD, RAW, Useable, cell_line, disease, tissue_type, treatment, sub_cell'''
        count = 0
        assay_df = assay_df["PXD RAW Useable cell_line disease tissue_type treatment sub_cell".split()]
        assay_tuples = list(assay_df.to_records(index = False))
        for i in assay_tuples:
            (accession, filename, useable, cell_line, disease, tissue_type, treatment, sub_cell) = i
            
            #filename = filename.split(".")[0]
            
            #select project_id
            self.mycursor.execute("SELECT project_id FROM project where PXD_accession = %s", (accession,))
            projectID_tup = self.mycursor.fetchone()
            (projectID,) = projectID_tup
            #insert into assay table
            assay = "insert into assay(project_id, filename) VALUES(%s, %s)"
            projectID_filename = (projectID, filename)
            self.mycursor.execute(assay, projectID_filename)
            self.conn.commit()
            #store this automatically generated assay ID for the cell_to_assay table
            assayID = self.mycursor.lastrowid
            #select cellID
            self.mycursor.execute("SELECT cell_id FROM cell WHERE cell_line = %s AND treatment = %s AND disease = %s AND sub_cell = %s", (cell_line, treatment, disease, sub_cell))
            cellID_tup = self.mycursor.fetchone()
            (cellID,) = cellID_tup
            #insert cellID and assayID in cell_to_assay
            cell_to_assay = "INSERT INTO cell_to_assay(assay_id, cell_id) VALUES(%s, %s)"
            assayID_cellID = (assayID, cellID)
            self.mycursor.execute(cell_to_assay, assayID_cellID)
            self.conn.commit()
            count += 1
        print(count)

    def id_regex(self, string):
        uniprot_regex = re.compile('\(\([a-zA-Z0-9]+\)\)')
        alt_uniprot_regex = re.compile('\(\([a-zA-Z0-9\|_-]+\)\)')
        protID = uniprot_regex.search(string)
        
        if protID == None: 
            protID = alt_uniprot_regex.findall(string)
            
            if len(protID) != 2:
                logging.warning(f"{string} does not match regex.")
                return np.nan

            protID = protID[1]
            
            if protID[2:5] == "sp|":
                return protID[5:-2]
            else:
                return protID[2:-2]
        
        return protID.group()[2:-2]

    def ionbot_parse(self, file, version):
        df = pd.read_csv(file, sep=',')
        
        if df.empty:
            logging.debug(f"{file.split('/')[-1]} is empty.")

        # best_psm is equal to 1
        df = df.loc[df['best_psm'] == 1]
        #  q-value-best <= 0.01
        df = df.loc[df['q_value'] <= 0.01]
        # DB column needs to contain 'T' (otherwise decoy hit) +  extra check: only retain swissprot entries (start with sp)
        df = df.loc[df['DB'] == 'T']

        if df.empty:
            logging.debug(f"{file.split('/')[-1]} lowest q-value: {pd.read_csv(file, sep = ',').q_value.min()}")
            return False

        #Only supports versions 
        versions = "IONBOT_v0.6.2 IONBOT_v0.6.3 IONBOT_v0.7.0 IONBOT_v0.8.0".split()
        if version not in versions:
            logging.debug(f"{file} not supported by versions {versions}")
            return False

        if version in ["IONBOT_v0.6.2", "IONBOT_v0.6.3"]:

            df_validated = df[df['proteins'].astype(str).str.startswith('sp')]
            # remove peptides that are not uniquely identified and are linked to multiple proteins = containing || in proteins
            x = '||'
            # regex is False otherwise it also detects a single | which is in every protein present
            df_validated = df_validated[~df_validated['proteins'].str.contains(x, regex=False)]
            #Reformat proteins to uniprotID

            if df_validated.empty:
                logging.debug(f"{file.split('/')[-1]}: no proteins after excluding duplicates.")
                return False

            df_validated["proteins"] = df_validated.apply(lambda x: x["proteins"].split('|')[1], axis = 1)
        
        elif version in ["IONBOT_v0.7.0", "IONBOT_v0.8.0"]:
            x = "||"
            df_validated = df[~df['proteins'].str.contains(x, regex = False)]

            if df_validated.empty:
                logging.debug(f"{file.split('/')[-1]}: no proteins after excluding duplicates.")
                return False

            df_validated["proteins"] = df_validated.apply(lambda x: self.id_regex(x["proteins"]), axis = 1)
            df_validated = df_validated[df_validated["proteins"].notna()]

        # check not all entries were removed
        if df_validated.empty:
            return False
        # modifications can be linked to unimod id: peptide_modifications: unimod ID vs peptide
        # calculte the spectral counts from each peptide: dict: count
        peptides = df_validated['matched_peptide'].tolist()
        spectral_counts = defaultdict(int)
        for pep in peptides:
            spectral_counts[pep] += 1
        spectral_counts = dict(sorted(spectral_counts.items(), key=lambda item: item[1], reverse=True))
        logging.info(f'{file.split("/")[-1]}\n#proteins: {len(df_validated.proteins.unique())}')
        #logging.info('parsed check')
        return df_validated, spectral_counts

    def ionbot_store(self, file, filename):
        #check if the assay isn't already in the assay table
        version = filename.split('/')[8]
        filename = filename.split('/')[-1].split('.')[0]

        self.mycursor.execute("SELECT assay_id FROM assay WHERE filename = %s", (filename,))
        assayIDtup = self.mycursor.fetchone()
        if assayIDtup is None:
            logging.info(f'{filename} is not in assays.')
            return
        (assayID,) = assayIDtup

        parser = self.ionbot_parse(file, version=version)
        if parser is False:
            logging.warning(f"parser failed for {filename}.")
            return
        df_validated, spectral_counts = parser

        # use the pandeylines in assay format
        # pandeylines resulted in a pd dataframe with all the proteins and sequences of validated peptides
        # loop over all rows/peptides present in the file (pandey_validated dataframe)
        df_validated_store = df_validated[['proteins', 'matched_peptide', 'modifications']]
        df_validated_tuples = [tuple(x) for x in df_validated_store.to_numpy()]
        for t in df_validated_tuples:
            protID = (t[0])
            pepseq = tuple((t[1],))
            mod = list((t[2],))

            # peptide storage - peptide ID
            sequence = "INSERT INTO peptide(peptide_sequence) VALUES (%s) " \
                        "ON DUPLICATE KEY UPDATE peptide_sequence=peptide_sequence"
            self.mycursor.execute(sequence, pepseq)
            self.conn.commit()

            # retrieve peptide_id, do not generate a new id each time!
            self.mycursor.execute("SELECT peptide_id FROM peptide WHERE peptide_sequence = %s", (pepseq))
            pepIDtup = self.mycursor.fetchone()
            (pepID,) = pepIDtup

            # link uniProtID = protein in assay to peptide
            proteinID = "INSERT INTO protein(uniprot_id) VALUES (%s) ON DUPLICATE KEY UPDATE uniprot_id=uniprot_id"
            uniprotID = (protID,)
            self.mycursor.execute(proteinID, uniprotID)
            self.conn.commit()

            # relation peptide to protein
            pepToProt = "INSERT INTO peptide_to_protein(uniprot_id, peptide_id) VALUES (%s,%s) " \
                        "ON DUPLICATE KEY UPDATE peptide_id=peptide_id, uniprot_id=uniprot_id"
            uniprotIDstr = ''.join(uniprotID)
            uniprotID_peptideID = (uniprotIDstr, pepID)
            self.mycursor.execute(pepToProt, uniprotID_peptideID)
            self.conn.commit()

            for i in mod:
                if pd.isnull(i):
                    break
                else:
                    # retrieve modification id, peptide id is present
                    location = (i.split('|')[0],)
                    id = (i[i.find("[")+1:i.find("]")],)

                    #retrieve modID
                    self.mycursor.execute("SELECT mod_id FROM modifications WHERE mod_id = %s", (id))
                    modIDtup = self.mycursor.fetchone()
                    if modIDtup is None:
                        break
                    (modID,) = modIDtup
                    # relation peptide_to_modification
                    peptoMod = "INSERT INTO peptide_modifications(peptide_id, location, mod_id, assay_id) VALUES (%s, %s, %s, %s)" \
                                "ON DUPLICATE KEY UPDATE peptide_id = peptide_id, mod_id = mod_id, assay_id=assay_id"
                    peptoModvalues = pepIDtup + location + modIDtup + assayIDtup
                    self.mycursor.execute(peptoMod, peptoModvalues)
                    self.conn.commit()

            # spectral count for peptide
            count = float('inf')
            for k, v in spectral_counts.items():
                if k == (''.join(pepseq)):
                    count = v
                    break
            peptideToAssay = "INSERT INTO peptide_to_assay(peptide_id, assay_id, quantification) VALUES (%s, %s, %s) " \
                                "ON DUPLICATE KEY UPDATE peptide_id=peptide_id, assay_id=assay_id"
            peptideID_assayID_count = (pepID, assayID, count)
            self.mycursor.execute(peptideToAssay, peptideID_assayID_count)
            self.conn.commit()
        logging.info(f'{filename} was stored')

    def find_ionbot_files(self, df: pd.DataFrame):
        """df must contain file_path column with the paths to the .ionbot.csv file"""
        
        logging.basicConfig(filename='ionbot_assays.log', level=logging.DEBUG)

        file_paths = df.file_path.tolist()
        number_of_files = 0
        for file in file_paths:
            number_of_files += 1

            if os.path.getsize(file) != 0:
                filename = str(file)
                self.ionbot_store(file, filename)
                
        print(number_of_files)

