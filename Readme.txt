
1. Overview
This repository provides scripts for metabolomics data processing, statistical analysis, and figure generation associated with the study, A Brain-Region-Resolved Untargeted Metabolomics Dataset of Prolonged Sleep Disruption in Mice.

2. Software requirements
MATLAB R2024b
Required toolboxes
Statistics and Machine Learning Toolbox
Bioinformatics Toolbox

3.  Directory structure
Code/
(1)  step1_Mergemetabolites.m
Merge positive- and negative-ion metabolic dataset.

(2)  step2_datafiltering.m
Perform data filtering, missing value imputation, normalization and log2 transformation, PCA analysis and RSD analysis.

(3) step3_SuperClass_piechart
Plot the superclass composition of identified metabolites. 

(4) step4_mergeexcel.m
Merge metabolites with VIP score 

(5) step5_vocanolplot.m
Identify significantly altered metabolites and generate volcano plots.

(6) step6_KEGGbubbleplot.m
Generate KEGG pathway enrichment bubble plots.

(7) step7_Clustering.m
Generate hierarchical clustering heatmaps.

4. Input files

NEG_MS2score_cutoff.xlsx
Negative ion metabolite dataset.

POS_MS2score_cutoff.xlsx
Positive ion metabolite dataset.

Merged_Unique_Metabolites.xlsx
Metabolite dataset with unique MS2 name, MS2 score,superclass and relative abundance.

oplsda_vip.csv
Results of  the OPLS-DA analysis,including metabolite lists and VIP scores.

Metabolites_with_statistics.xlsx
Metabolite dataset with metabolites name, mean relative abundance of metabolites in control and CSD groups,statistical analysis and Superclass.

pathway_results.csv
KEGG pathway enrichment results.

Sigdata folder
Significant altered metabolites in each brain region.

5. Output files
(1)  Processed metabolite matrices

• Merged_Unique_Metabolites.xlsx
• filtered_metabolites.xlsx
• filtered_metabolites_imputed.xlsx
• filtered_metabolites_normalized.xlsx
• filtered_metabolites_log2_transformation.xlsx

(2)  Figures

• PCA plots
• QC RSD distribution
• Volcano plots
• KEGG enrichment bubble plots
• Hierarchical clustering heatmaps

(3)  Intermediate data

• total_metabolites.mat
• up_regulated.mat
• down_regulated.mat
• Significant_Metabolites.xlsx


(4)  Shared metabolite tables

• Shared_Significant_Metabolites_4Regions.xlsx
• Shared_Significant_Metabolites_in_3or4Regions.xlsx
• Shared_Metabolites_Zscore_Data.xlsx

6. How to run
Run the script in the following order:
1. step1_Mergemetabolites.m

2. step2_datafiltering.m

3. step3_Superclass_piechart.m

4. step4_mergeexcel.m

5. step5_vocanolplot.m

6. step6_KEGGbubbleplot.m

7. step7_Clustering.m

7. Notes
• All scripts were developed and tested using MATLAB R2024b.
• Input files should be placed in the same directory as the scripts or the file paths should be modified accordingly.

