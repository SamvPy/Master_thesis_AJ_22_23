from pathlib import Path
import os
os.chdir(Path.home() / 'master_thesis_folder')

tissues = ['Adipose tissue', 'Adrenal gland', 'B-cells', 'Bone',
'Brain', 'Cartilage', 'Cervix', 'Colon', 'Dental plaque', 'Esophagus',
'Eye', 'Heart', 'Kidney', 'Liver', 'Lung', 'Monocytes', 'NK-cells',
'Nasal polyps', 'Ovary', 'PBMC', 'Palatine tonsils', 'Pancreas', 'Parotid gland',
'Retina', 'Skeletal muscle', 'Skin', 'Small intestine', 'T-cells', 'Testis', 'Tooth',
'Umbilical cord', 'Ureter', 'Urinary bladder']

text_string = ""
for count, i in enumerate(tissues):
    text_string += f'{count+1}. {i}:\n'

pxd_list_file = open('pxd_list.txt', 'w')
pxd_list_file.write(text_string)
pxd_list_file.close()