# svmGNT workflow script

R code for pre-processing of Illumina methylation array data, consensus clustering, t-SNE plots, and SVM glioneuronal classification model training. Used in "DNA methylation-based classification of glioneuronal tumours synergises with histology and radiology to refine accurate molecular stratification." by Stone et al.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Hardware requirements](#hardware-requirements)
- [Software Requirements](#software-requirements)

# Overview
The workflow script `svmGNT.R` contains the step-by-step commands used to read in and pre-process Illumina 450K/EPIC methylation array data, in addition to subsequent consensus clustering analysis, t-SNE plotting, and construction of a simple support vector machine classification model. 

`850Kremove_probes.txt` is provided to document probes removed following pre-processing. These are derived from Pidsley et al (2016) and correspond to probes within 50bp of an SNP, probes with at least 1 cross reactive target, and probes with a MAF >5%.

# System Requirements
## Hardware Requirements
The script provided should be compatible with any standard system that supports *R 4.1.1* with enough RAM to support the constituent operations. This script was implemented originally on a system with 8Gb RAM.

## Software Requirements
The script provided depends on a number of R packages. Functionality has been tested on *R 4.1.1* operating under *macOS Catalina 10.15.7*
### R Dependencies

```
 [1] IlluminaHumanMethylation450kanno.ilmn12.hg19
 [2] IlluminaHumanMethylation450kmanifest      
 [3] MLmetrics                                  
 [4] doParallel                                
 [5] caret                                    
 [6] lattice                                  
 [7] ggplot2                                    
 [8] Rtsne                                       
 [9] ConsensusClusterPlus                      
[10] minfi                                    
[11] bumphunter                                
[12] locfit                                   
[13] iterators                                
[14] foreach                                   
[15] Biostrings                                
[16] XVector                                  
[17] SummarizedExperiment                      
[18] Biobase                                   
[19] MatrixGenerics                             
[20] matrixStats                              
[21] GenomicRanges                           
[22] GenomeInfoDb                              
[23] IRanges                                   
[24] S4Vectors                                 
[25] BiocGenerics
```
