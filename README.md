# Focused BH

A novel FDR-controlling method for hierarchically structured hypotheses that inputs any pre-specified filter and outputs a non-redundant rejection set with respect to that filter. 

Accompanying paper:
> *FDR control following filtering for hierarchically structured hypotheses* <br />
> E. Katsevich, C. Sabatti, M. Bogomolov <br />
> arXiv 1809.01792: [https://arxiv.org/abs/1809.01792](https://arxiv.org/abs/1809.01792)

## Overview

Focused BH is a multiple testing procedure designed for a broad range of applications where structured hypotheses arise. It is mainly motivated by problems with hierarchical structure, such as (1) phenome-wide association studies, with tree structured diseases based on the [International Classification of Diseases](https://icd.who.int/browse10/2016/en), with the outer nodes filter and (2) Gene Ontology enrichment analysis, with directed acyclic graph structured biological processes based on the [Gene Ontology](http://geneontology.org/), with the [REVIGO](http://revigo.irb.hr/) filter. Focused BH also extends beyond hierarchically structured applications, such as to spatially structured applications like genome-wide association studies or neuroimaging analysis.

This repository provides software implementing Focused BH, as well as code to reproduce all numerical simulations and data analysis in the paper. 

## Repository structure

* [data/](https://github.com/ekatsevi/Focused-BH/tree/master/data): contains raw and processed data from the PheWAS analysis of the [UK Biobank](https://www.ukbiobank.ac.uk/) data[^1] and the GO enrichment analysis of the [breast cancer outcome data](https://www.ncbi.nlm.nih.gov/pubmed/11823860). 

[^1] Note that some of the UK Biobank data is not publicly accessible. 



## Files 

* FocusedBH.R:      implementation of Focused BH for the outer nodes filter
* Test_FocusedBH.R: demonstrates the use of the Focused BH code
* aux.R:            auxiliary functions (necessary only for Test_FocusedBH.R)
* coding6.tsv:      tree-structured self-reported phenotype coding from UK Biobank

## Author

* [Eugene Katsevich](http://www.andrew.cmu.edu/user/ekatsevi/)