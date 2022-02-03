class NSAF():
    import pandas as pd
    import mysql.connector
    import math
    def __init__(self, db, tissue, percentage, disease_status):
        """Start building the atlas. Provide the database connection, level of the atlas (tissue_name or cell_type), the percentage of filtering and the disease_status which is desired (Healthy, None, Cancer,...)"""
        self.db = db
        self.percentage = percentage
        self.conn = self.mysql.connector.connect(user='root', password='password', host='127.0.0.1', port='3306',database=self.db, auth_plugin='mysql_native_password')
        mycursor = self.conn.cursor()
        if tissue in ('tissue_name', 'cell_type'):
            self.tissue = tissue
        else:
            print('tissue type not defined')
            raise ValueError
        self.disease_status = disease_status

    def check_connection(self):
        if self.conn.is_connected():
            print("connection succesfull")
        else:
            print("no connection")
    
    def get_tissue_assay_quantification(self):
        """ Build the table containing assay - peptide - quantification - tissue information"""
        assayprojectsql = "SELECT * FROM assay where assay_id in (SELECT assay_id FROM tissue_to_assay);"
        assayprojectData = self.pd.read_sql_query(assayprojectsql, self.conn)
        
        assaysql = "SELECT assay_id, peptide_id, quantification FROM peptide_to_assay"
        assayData = self.pd.read_sql_query(assaysql, self.conn)
        
        assaytissuesql = "SELECT assay_id, tissue_id FROM tissue_to_assay"
        assaytissueData = self.pd.read_sql_query(assaytissuesql, self.conn)
        
        tissuesql = "SELECT tissue_id, tissue_name,cell_type, disease_status FROM tissue"
        tissueData = self.pd.read_sql_query(tissuesql, self.conn)
        
        tissue_assay = self.pd.merge(assaytissueData, tissueData, on='tissue_id', how='left')
        tissue_assay = self.pd.merge(assayData, tissue_assay, on='assay_id', how='left')
        return tissue_assay
    
    def get_sequence_data(self):
        """ select all the proteins from the database and their length """
        seqsql = "SELECT uniprot_id, length FROM protein WHERE length IS NOT NULL"
        seqData = self.pd.read_sql_query(seqsql, self.conn)
        seqData['length'] = self.pd.to_numeric(seqData['length'], errors='coerce')
        return seqData

    def get_proteins(self):
        "Get table containing the proteins and their quantification, cleaned for proteotypic data"
        pepsql = "SELECT peptide_to_protein.peptide_id, peptide_to_protein.uniprot_id FROM peptide_to_protein"
        pepData = self.pd.read_sql_query(pepsql, self.conn)    
   
        proteotypicData = pepData.groupby("peptide_id").filter(lambda x: len(x) == 1)
        proteins = proteotypicData.groupby("uniprot_id").filter(lambda x: len(x) > 2)
        non_human_proteins = ['TRYP_PIG', 'TRY2_BOVIN','TRY1_BOVIN','SSPA_STAAU','SRPP_HEVBR','REF_HEVBR', 'ADH1_YEAST', 'ALBU_BOVIN', 'CAS1_BOVIN', 'CAS2_BOVIN', 'CASK_BOVIN', 'CASB_BOVIN', 'OVAL_CHICK', 'ALDOA_RABIT', 'BGAL_ECOLI', 'CAH2_BOVIN', 'CTRA_BOVIN', 'CTRB_BOVIN', 'CYC_HORSE', 'DHE3_BOVIN', 'GAG_SCVLA', 'GFP_AEQVI', 'K1C15_SHEEP', 'K1M1_SHEEP', 'K2M2_SHEEP', 'K2M3_SHEEP', 'KRA3A_SHEEP', 'KRA3_SHEEP', 'KRA61_SHEEP', 'LALBA_BOVIN', 'LYSC_CHICK', 'LYSC_LYSEN', 'MYG_HORSE', 'K1M2_SHEEP', 'K2M1_SHEEP']
        proteins = proteins[~proteins['uniprot_id'].isin(non_human_proteins)]
        return proteins

    def merge_protein_and_tissue_assay(self):
        "merge the information from the proteins with the assay and tissue information"
        tissue_assay = self.get_tissue_assay_quantification()
        proteins = self.get_proteins()
        protData = self.pd.merge(tissue_assay, proteins, on = 'peptide_id').sort_values(['assay_id','uniprot_id'])
        if self.tissue == 'tissue_name':
            del protData['cell_type']
        if self.tissue == 'cell_type':
            del protData['tissue_name']
        del protData['peptide_id']
        del protData['tissue_id']
        del protData['disease_status']
        return protData

    def filter_protein_data(self):
        """filter the proteins per tissue based on the percentage given in the initialisation"""
        print('started protein filtering')
        protData = self.merge_protein_and_tissue_assay()
        assays = protData[self.tissue].unique()
        reduction = []
        DataFrameDict = {elem : self.pd.DataFrame for elem in assays}
        for key in DataFrameDict.keys():
            DataFrameDict[key] = protData[:][protData[self.tissue] == key]
            perc = self.math.floor(self.percentage * len(self.pd.unique(DataFrameDict[key]['assay_id'])))
            before= DataFrameDict[key]['uniprot_id'].nunique()
            DataFrameDict[key] = DataFrameDict[key].groupby('uniprot_id').filter(lambda x : len(x)>perc)
            after= DataFrameDict[key]['uniprot_id'].nunique()
            reduction.append(before-after)
            print('{} is done'.format(key))
        filteredData = self.pd.DataFrame()
        for key in DataFrameDict.keys():
            filteredData = filteredData.append(DataFrameDict[key])
        del filteredData[self.tissue]
        filteredData = filteredData
        reduction = sum(reduction)
        print('finished protein filtering')
        print('first there were {} proteins and now there are {}'.format(before, after))
        return filteredData

    def calculate_NSAF(self):
        """calculate the NSAF from the filtered protein data"""
        seqData = self.get_sequence_data()
        filteredData = self.filter_protein_data()
        assays = filteredData['assay_id'].unique()
        print('started NSAF calculations')
        DataFrameDict3 = {elem : self.pd.DataFrame for elem in assays}
        counter = 0
        for key in DataFrameDict3.keys():
            DataFrameDict3[key] = filteredData[:][filteredData['assay_id'] == key]
            sumSaf = 0
            counter += 1
            df = DataFrameDict3[key]
            df = df.drop(columns=['assay_id'])
            # calculate sum of spectral counts for each protein
            df_seq = self.pd.merge(df.groupby('uniprot_id').sum().reset_index(), seqData, on = 'uniprot_id')
            df_seq.insert(loc = 2, column = 'SAF', value = 0)
            df_seq.insert(loc = 3, column = 'NSAF', value = 0)
            # calculate SAF score for each protein by dividing sum of spectral counts by protein length
            df_seq['SAF'] =  df_seq['quantification']/df_seq['length']
            # calculate sum of SAF scores in assay
            sumSaf = df_seq['SAF'].sum()
            # Calculate NSAF score by normalizing each SAF score
            df_seq['NSAF'] = df_seq['SAF']/ sumSaf
            del df_seq['length']
            del df_seq['quantification']
            del df_seq['SAF']
            df_seq.insert(loc = 0, column = 'assay_id', value = key)
            DataFrameDict3[key] = df_seq
            print(counter)
        proteinData = self.pd.DataFrame()
        for key in DataFrameDict3.keys():
            proteinData = proteinData.append(DataFrameDict3[key])
        self.proteinData = proteinData
        print('finished NSAF calculations')
        return self.proteinData

class Atlas():
    import pandas as pd
    import mysql.connector
    def __init__(self, db, tissue, disease_status, df_nsaf):
            self.db = db
            self.conn = self.mysql.connector.connect(user='root', password='password', host='127.0.0.1', port='3306',database=self.db, auth_plugin='mysql_native_password')
            mycursor = self.conn.cursor()
            if tissue in ('tissue_name', 'cell_type'):
                self.tissue = tissue
            else:
                print('tissue type not defined')
                raise ValueError
            self.disease_status = disease_status
            self.nsaf = df_nsaf
    
    def get_tissue(self):
        """ get all the tissue information with the assay ids and select based on the disease status"""
        tissuesql = """SELECT tissue_to_assay.assay_id, tissue.cell_type, tissue.tissue_name, tissue.disease_status, tissue.organ_id, tissue.fluid FROM tissue_to_assay JOIN tissue ON tissue_to_assay.tissue_id = tissue.tissue_id"""
        tissueData = self.pd.read_sql_query(tissuesql,self.conn)
        if self.disease_status != None:
            tissueData = tissueData[tissueData['disease_status']==self.disease_status]
        return tissueData

    def create_atlas(self):
        """ create the atlas using the NSAF proteome and the assay-tissue values"""
        tissueData = self.get_tissue()
        atlas = self.pd.merge(self.nsaf, tissueData, on='assay_id')
        return atlas

class Predictor():
    import pandas as pd
    import mysql.connector
    def __init__(self, db, tissue, disease_status, df_nsaf):
            self.db = db
            self.conn = self.mysql.connector.connect(user='root', password='password', host='127.0.0.1', port='3306',database=self.db, auth_plugin='mysql_native_password')
            mycursor = self.conn.cursor()
            if tissue in ('tissue_name', 'cell_type'):
                self.tissue = tissue
            else:
                print('tissue type not defined')
                raise ValueError
            self.disease_status = disease_status
            self.nsaf = df_nsaf
    
    def get_tissue(self):
        """get the tissue and assay information from the database"""
        tissuesql = """SELECT tissue_to_assay.assay_id, tissue.cell_type, tissue.tissue_name, tissue.disease_status, tissue.fluid FROM tissue_to_assay JOIN tissue ON tissue_to_assay.tissue_id = tissue.tissue_id"""
        tissueData = self.pd.read_sql_query(tissuesql,self.conn)
        if self.disease_status != None:
            tissueData = tissueData[tissueData['disease_status']==self.disease_status]
        if self.tissue == 'tissue':
            tissueData = tissueData.drop(['cell_type', 'disease_status', 'fluid'], axis=1)
        elif self.tissue == 'cell_type':
            tissueData = tissueData.drop(['tissue_name', 'disease_status', 'fluid'], axis=1)
        return tissueData
    
    def get_assay_atlas(self):
        """merge the NSAF values with the tissue data on assay to get the final predictor atlas"""
        tissueData = self.get_tissue()
        assay_atlas = self.pd.pivot_table(self.nsaf, values = 'NSAF', index = 'assay_id', columns = 'uniprot_id').fillna(0).reset_index()
        atlas = self.pd.merge(assay_atlas, tissueData, on='assay_id')
        atlas = atlas.drop(columns=['assay_id'])
        return atlas

