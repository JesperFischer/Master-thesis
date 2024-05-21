# Uncertainty In Cognitive Science
Public repository of the Master thesis "Uncertainty in Cognitive Science" containing data, source code visualizations and the manuscript used in this master thesis.

## Table of Contents
1. [Abstract](#abstract)
2. [Introduction of repository](#introduction)
3. [Directory Structure](#directory-structure)
4. [Access](#access)
5. [Reproducibility](#reproducibility)
6. [Usage](#usage)

## Abstract


## Introduction of repository


This repository contains all scripts necessary to execute and generate the results, plot, analyses, and the manuscript for the Master thesis. To ensure complete reproducibility, the project is integrated with the Open Science Framework [OSF](https://osf.io/qrsc3/) to facilitate either the rerunning of analyses or the retrieval of results which are located on the [OSF](https://osf.io/qrsc3/). 



Additionally, the computational environment required for these analyses can be precisely replicated using the [renv](https://rstudio.github.io/renv/articles/renv.html) package, which is specifically designed for managing R environments. This setup guarantees that other researchers can accurately reproduce and build upon the work presented in this study.


## Directory Structure

The repository is structured in the following way:

```         
Master-thesis/
│
├── Analyses/                   # Directory containing Directories for all analyses in the manuscript.
│   ├── ICC analysis                              # Directory for all parameter recovery analyses as well as uncertainty minimization.
│   ├── Legrand reanalysis                        # Directory for all models and analyses of the reanalysis of Legrand (2021).
│   └── Power analysis                            # Directory for the scripts and analysis of the power analysis.
│
├── data ICC/                   # Directory containing the simulated data for the parameter recovery analysis.
│   └── ... 
│
├── Figures/                    # Figures generated from code to the from the final manuscript.
│   └── ... 
│
├── manuscript/                 # Folder containing everything needed to recreate the final manuscript in pdf,docx and html.
│   ├── 1_Introduction.Rmd                        # Rmarkdown for the introduction.
│   ├── 2_Measurement uncertainty.Rmd             # Rmarkdown for the section on measurement uncertainty.
│   ├── 3_Modeling definitions.Rmd                # Rmarkdown for the section on modeling definitions.
│   ├── 4_uncertainty minimazation.Rmd            # Rmarkdown for the section on uncertainty minimazation.
│   ├── 5_Real data.Rmd                           # Rmarkdown for the section on real data.
│   ├── 6_Power analysis.Rmd                      # Rmarkdown for the section on power analysis.
│   ├── 7_discussion.Rmd                          # Rmarkdown for the section on discussion.
│   ├── Knitting files                            # Directory containing files for knitting the manuscript. Including templates and citations.
│   │    └── ... 
│   ├── Manuscript.Rmd                            # Rmarkdown for the final manuscript. Containing all of the above sections (used for knitting).
│   └── Supplementary materia                     # Directory for all supplementary materials. Including analyses and materials.
│       └── ... 
│
├── osf/                        # Directory containing a script for getting the workspace for the analysis.
│   └── ... 
│
├── plot scripts/               # Directory containing scripts for generating all figures in the Thesis.
│   └── ... 
│
├── python/                     # Directory for the single python script that excutes the PSI-adaptive design optimization procedure.
│   └── ... 
│
├── README.md                   # Overview of the project.
│
├── Stanmodels/                 # Directory containing Stan scripts for particular Stan-models used in the Thesis.
│   └── ... 
│
├── Supplementary figures/      # Directory with all supplementary figures.
│   └── ... 
│
├── Supplementary tables/       # Directory with all supplementary tables.
│   └── ... 
│
└── tables/                     # Directory with all main tables.
    └── ... 


```

## Access

To get access to the already run workspace for the a OSF-token is needed, as the workspaces are stored in the following [OSF-project](https://osf.io/7pu6a/) due to space limitations of github. 


It is recommended to make an osf folder that contains an osf.txt file which on the first line contains the OSF-token, however make sure that this token is not shared or pushed to github. In the current gitignore an osf folder will be ignored.
To get access to the repository users are recommended to clone the respository with the following command in the terminal





## Reproducibility

To enhance reproducibility this project is setup with R-package "renv"", ensuring the same packages and versions of these are loaded. This means that users should install the "renv" package and after cloning the repository, run the following code in the console (Not terminal).

```r
renv::restore()
```

This will download and install all the needed packages to run the analysis from scratch with the same version of packages used to generate the manuscript. 


## Usage

You can use this repository in various ways. Choose the one that suits your needs: Note you have to complete the [Access](#access) step to get get the data from OSF. 

1. **Method 1 - Re-create manuscript without rerunning analysis**:
   -  Go to the terminal and write ```bash bash run.sh```, which will rewrite the Manuscript in all 3 document types.
   -  Open the Manuscript folder and knit the [Manuscript.Rmd](./Manuscripts/Manuscript.Rmd) file to the desired format (.docx, .pdf or .html). 

2. **Method 2 - Running the main analysis line by line**:
   - Open the Markdown folder and go through the [Analysis.Rmd](./Markdowns/Analysis.Rmd) markdown.

3. **Method 3 - Rerunning the plots line by line**:
   - Open the Markdown folder and go through the [Plots.Rmd](./Markdowns/Plots.Rmd) markdown.
   - To view the plotting functions users are encouraged to check out the [plot.R](./scripts/plots.R) script inside the scripts folder.

