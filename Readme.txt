
1. Overview
This repository provides scripts for metabolomics data processing, statistical analysis, and figure generation associated with the study, A Brain-Region-Resolved Untargeted Metabolomics Dataset of Prolonged Sleep Disruption in Mice.

2. Software requirements
MATLAB R2024b
Required toolboxes
Statistics and Machine Learning Toolbox
Bioinformatics Toolbox

3.  Directory structure
Code/
(1) Step1_QC_RSD_distribution.m
1.1 Plot the distribution of QC-RSD values for each feature.
1.2 Plot the cross-correlation of QC samples under both POS and NEG modes.
1.3 Perform data filtering and imputation under both POS and NEG modes.

(2) Step2_Data_filtering.m
1.1 Integrate metabolites from both POS and NEG modes based on QC-RSD and MS2 score.
1.2 Plot the distribution of MS2 score values for all filtered metabolites.

(3) Step3_Supreclass_piechart.m
Plot the superclass composition of identified metabolites. 

(4) Step4_Vocanolplot.m
Identify significantly altered metabolites and generate volcano plots.

(5) Step5_KEGGbubbleplot.m
Generate KEGG pathway enrichment bubble plots.

(6) Step6_Upsetplot.m
Generate UpSet plots to visualize the intersection of significantly altered metabolites across brain regions.

(7) Step7_Heatmap_plot
Generate Z-score heatmaps to visualize commonly altered and region-specific metabolites.

4. How to run
Run the scripts sequentially from Step 1 through Step 7.

5. Notes
5.1 All scripts were developed and tested using MATLAB R2024b.
5.2 Input files should be placed in the same directory as the scripts or the file paths should be modified according
