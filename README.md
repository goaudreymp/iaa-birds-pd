# iaa-birds-pd
code and data for chapter 2 of thesis

# Updated supertree for passerines from the Indo-Australian Archipelago highlights superendemism in New Guinea.
This folder includes data and code to run the analysis for the manuscript in progress/thesis Chapter 2 

The associated code files includes:
- Phylo folder: all scripts and data required for the phylogenetic inference part of the analyses
- R_PD folder: all R code and input required for the phylogenetic diversity part of the analyses

Below is a description for each data file included in this folder that is required to run the above code:

## Phylogenetic Inference
- calib folder - subtrees that are used with calibration
- list folder - list of each taxa in each subtree
- mcmc folder - scripts and ctl files used for the baseml and mcmc dating
- `concat_090.fasta` - full concatenated alignments of all taxa
- `final_grafting.py` - final python code for grafting subtrees to backbone tree
- `finalbirdsp.tree` - final bird tree result
 
## Phylogenetic Diversity & Endemism
- `iaa-tree-pd.Rproj` - R project file for easy management
- `mrcanode_tree.R` - used to identify MRCA nodes in the tree for visualisation
- `cleanscript.R` - includes all R code required to load trees and run the phylogenetic diversity analyses
