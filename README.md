# Master thesis 2022-2023
Ghent University<br>
Sam van Puyenbroeck <br>
---
This repository functions as lab notebook for the work done during the master thesis with the following topic title: <br>
**"USING REPROCESSED PUBLIC PROTEOMIC DATA TO PREDICT PROTEIN PATTERNS IN HUMAN CELL LINES"** <br> <br>
In the *Reference* folder, all material created by Tine Claeys is located. Some parts on building the database, parsing ionbot result files and quantification is based on this material. <br>
All analyses performed during the master's thesis is located in the *master_thesis* folder. <br>
<br>

## Abstract
Cancer cell lines are widely used in cancer research as a model system to study the aberrant 
pathways that give rise to cancer and to test the efficacy of cancer treatments. In this project, 43 
PRIDE-projects are reprocessed and combined to build a model capable of accurately classifying 
cell line groups which through feature importance analysis can provide insights into which 
proteins are most discriminative for a group of cell lines. To build such a model, for each pre-processing step consisting of: (i) normalisation; (ii) imputation; (iii) feature selection; and (iv) 
oversampling; several methods were implemented and evaluated. Additionally, the systematic 
difference in protein identifications due to sample preparation was explored and whether or not 
correlations in our dataset can be predictive of functional associations. Our findings suggest that 
by performing in-gel digestion a more hydrophobic part of the proteome is identified than 
compared to in-solution digestion. Additionally, we found pairwise protein Pearson correlations 
between protein abundances to have a predictive value for functional associations. Based on our 
evaluations of the pre-processing methods, the optimal pre-processing pipeline included (i) 
quantile normalisation; (ii) limit-of-detection imputation; (iii) an ensemble of feature selection 
methods; and (iv) SMOTE. Subsequently, a Logistic Ridge Regression was trained and reached 
93.7% classification performance. A preliminary biological exploration of the model showed slight 
concordance with the annotations made by the Human Protein Atlas. By exploring the model 
further and leveraging the correlations within the dataset, more insights can be gained into what 
makes cell line groups unique.

## Layman summary with societal impact
Cancer continuous to be a prominent global cause of mortality, underscoring the importance of 
early detection and effective treatment. However, due to the intrinsic heterogeneity of cancer 
across patients, the same treatment does not work for everyone and a highly personalised 
treatment strategy is necessary. 
In the realm of cancer research, surrogate model systems known as cell lines are frequently used 
to study the effectiveness of various compounds in treating cancer. However, due to the fact that 
cell lines are intrinsically different from the cancer subtypes of patients, the translation of 
compound effectiveness on cell lines to the clinic is highly inefficient.
In our study, we have tried to capture the differences between cell lines. This approach can be 
beneficial in two ways. Firstly, it can help to understand the biological diversity between cell lines
on a systematic level. As cell lines model cancer subtypes, this could help to identify new 
biomarkers able to distinguish between cancer subtypes and thus make patient treatment more 
precise and lower the devastating effects of unresponsive therapy. Additionally, our approach
could match a patient’s cancer subtype to the most closely related cancer cell line. By capitalizing 
on this resource, researchers could conduct more precise assessments of the efficacy of potential 
treatments and interventions on cell lines that closely mirror the patient's specific cancer subtype.
This can enhance the predictive value of preclinical studies and lower the high costs associated 
with developing unsuccessful drug candidates. 

