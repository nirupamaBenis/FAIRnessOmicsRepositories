# FAIRnessOmicsRepositories

This repository contains supplementary information for the manuscript titled **"FAIRness assessment of metadata of omics datasets in online repositories"**. A Python script and a pickle file are in this repository.  
  
**Python script: apiSearch.py**  
This script searches on the follwoing 7 omics repositories using their APIs. Excuting this script will allow you to replicate the search we have done. The search terms are hard coded into the script and will have to be modified for each repository for searches with different keywords.  
* ArrayExpress
* DbGaP
* ENA
* GEO
* Metabolomics Workbench
* OmicsDI
* PRIDE  

**Pickle file: globalsaveAllAPIsAllTerms.pkl**  
We used the package `dill` to obtain a data dump of the Python session after the repositories were queried. This data dump contains all the raw results from the repositories and the processed information that was ultimately used in the manuscript. We considered this snapshot important to include since omics repositories are updated often and the same search done at a different time can yield different results. 
