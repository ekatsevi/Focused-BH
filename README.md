# Focused BH

A novel FDR-controlling method for structured hypotheses that inputs any pre-specified filter and outputs a non-redundant rejection set with respect to that filter. 

Accompanying paper:
> *FDR control following filtering for hierarchically structured hypotheses* <br />
> E. Katsevich, C. Sabatti, M. Bogomolov <br />
> arXiv 1809.01792: [https://arxiv.org/abs/1809.01792](https://arxiv.org/abs/1809.01792)

## Overview

Focused BH is a multiple testing procedure designed for a broad range of applications where structured hypotheses arise. It is mainly motivated by problems with hierarchical structure, such as (1) phenome-wide association studies, with tree structured diseases based on the [International Classification of Diseases](https://icd.who.int/browse10/2016/en), with the outer nodes filter and (2) Gene Ontology enrichment analysis, with directed acyclic graph structured biological processes based on the [Gene Ontology](http://geneontology.org/), with the [REVIGO](http://revigo.irb.hr/) filter. Focused BH also extends beyond hierarchically structured applications, such as to spatially structured applications like genome-wide association studies or neuroimaging analysis.

This repository provides software implementing Focused BH, as well as code to reproduce all numerical simulations and data analysis in the paper. 

## Repository structure

* [data/](https://github.com/ekatsevi/Focused-BH/tree/master/data): contains raw and processed data from the PheWAS analysis of the [UK Biobank](https://www.ukbiobank.ac.uk/) data<sup>1</sup> and the GO enrichment analysis of the [breast cancer outcome data](https://www.ncbi.nlm.nih.gov/pubmed/11823860). 
* [precomp/](https://github.com/ekatsevi/Focused-BH/tree/master/precomp): contains precomputation results for our numerical simulations.
* [results/](https://github.com/ekatsevi/Focused-BH/tree/master/results): contains final results for our numerical simulations.
* [src/](https://github.com/ekatsevi/Focused-BH/tree/master/src): contains all source code for methodology, numerical simulations, data analysis, and plotting.

<sup>1</sup> Note that our analysis of the UK Biobank data (application number 27837) is based on data fields [41202](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41202) (Diagnoses - main ICD10) and [22182](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=22182) (HLA imputation values), which are not publicly accessible. If you do not have access to these fields, we have provided summary statistics in [data/processed/biobank](https://github.com/ekatsevi/Focused-BH/tree/master/data/processed/biobank) so that the multiple testing portion of our analysis can be reproduced. If you do have access to these fields, place them into a comma-separated file called "ukb25261.csv" in the directory [data/raw/biobank](https://github.com/ekatsevi/Focused-BH/tree/master/data/raw/biobank) to reproduce our full data analysis.

## Author

* [Eugene Katsevich](http://www.andrew.cmu.edu/user/ekatsevi/)