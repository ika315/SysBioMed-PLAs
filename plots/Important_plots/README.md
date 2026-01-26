# PLA Identification Pipeline – Results Overview

This folder contains the final benchmarking and identification results for PLAs. 

## Folder Structure

### Reference_and_Overview/
1_Reference_UMAP.png: Displays the ground truth cell type distribution based on manual gating. This serves as the spatial reference for all subsequent classification attempts.

### Global_Benchmarking/
Contains global comparison plots across all tested signatures and methods:

Roadmap_Extension_Benefit.png: Quantifies the F1-score change  by using "Extended Gene Lists" (ExtTRUE) compared to baseline signatures.

PR_Comparison_Incl_Youden.png: A global Precision-Recall overview comparing supervised and unsupervised strategies.

Ultimate_Approach_Comparison.png: Summarizes the F1-score evolution, showing the impact of extension and thresholding modes.

Method_Runtime_Efficiency_Log.png: Benchmarks the computational time per enrichment method.

### Detailed_Performance_Heatmaps/

Heatmaps showing cell-type specific performance metrics (Balanced Accuracy, FP count, Mean_Z). 
These plots are faceted by the Threshold Mode to show how specific lineages (e.g., Monocytes) are affected by different gating strategies.

### Individual_Method_Deep_Dives/

Detailed analysis of selected representative runs using the AUCell + ExtTRUE configuration.
Note: For these deep dives, specific gene signatures were chosen that best represent the strengths of each thresholding variant:

AUCell_ExtTRUE_youden/ (Signature: REACTOME):

	- Approach: Supervised (utilizes ground truth)
	- Purpose: Serves as the "Gold Standard" reference to show the maximum theoretical performance of the signature

AUCell_ExtTRUE_gmm_dist_platelet/ (Signature: GOBP_REG):

	- Approach: Unsupervised 1D-Gating
	- Purpose: Demonstrates automated identification using only the platelet signal distribution

AUCell_ExtTRUE_gmm_dist_dual/ (Signature: MANNE_DN):

	- Approach: Unsupervised 2D-Gating (Platelet + Immune Activation)
	- Purpose: It combines platelet signals with leukocyte activation markers to minimize false positives in highly active cell types

## Thresholding Logic Explained

The pipeline uses Gaussian Mixture Models (GMM) via the mclust package to identify PLAs without relying on manual labels:

1. gmm_dist_platelet (1D-Clustering): This variant performs statistical clustering solely on the Platelet Z-Score. It identifies a "high" and a "low" population. While automated, it can sometimes misclassify highly active leukocytes as PLAs due to non-specific background signals.
2. gmm_dist_dual (2D-Clustering): It first identifies cells with high platelet scores and then performs a second clustering step on the Leukocyte Activation Score. A cell is only classified as a "Positive PLA" if it belongs to the "High" cluster in both dimensions. This significantly increases precision by filtering out technical noise and non-aggregated platelets.


