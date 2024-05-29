# Uncertainty In Cognitive Science
Public repository of the Master thesis "Uncertainty in Cognitive Science" containing data, source code visualizations and the manuscript used in this master thesis.

## Table of Contents
1. [Abstract](#abstract)
2. [Introduction of repository](#introduction)
3. [Directory Structure](#directory-structure)
4. [Access](#access)
5. [Usage](#usage)

## Abstract

Understanding human cognition and behavior is the primiary aim of Cognitive Science. To acheive this quantitative methods are frequently used. When these quantitative methods are employed assumptions and simplifications must be made, which are embedded in the models used, which try to make sense of the quantitative data. This thesis explores some of these assumptions herein, propagation of uncertainty and validation of the models themselves in a simulated setting.  propagation of uncertainty stems from the fact that when quantitative data is collected uncertainty is embedded in these measurements. A failure to account for these uncertainties can have a substantial impact on the inferences made based on these measures. With a focus in how uncertainties in statistical and cognitive models propagate, the thesis will further investigate how the validation process of cognitive models can be improved, by embracing and quantifying the inevitable uncertainties associated with our cognitive models. The framework proposed revolves around simulating agents with known properties which are then fitted to the cognitive model to asses the model's ability to detect these simulated properties while accounting for uncertainties embedded in these estimations. The thesis will explore these considerations through the three-parameter psychometric function, a widely used cognitive model.

It will be shown that a correlational approach to determine internal model validity is at best quite insensible compared to a more sophisticated approach based on the intra class coefficient. The thesis will then demonstrate how uncertainties in parameter estimates and in these two metrics can be minimized through more sophisticated methods. This will be done without a need for increasing the number of trials or subjects, which is the standard approach. In this regard the thesis highlights two important methods for minimizing uncertainty. Firstly, optimizing the design of the experiment such that each trial will contain the most information possible. Secondly, incorporating already collected data, such as reaction times, into the cognitive model as a means of decreasing the uncertainties in the measures of interest. The thesis goes on to explore and re-analyses published data using the psychometric function. Here it is demonstrated that incorporating structural assumptions of how the data was collected, as well as incorporating reaction times, does not only decrease uncertainty in the reliability, but also well describes the data. Lastly the thesis highlights and demonstrates novel opportunities for conducting power analyses using. Here it is demonstrated, based on the re-analysis of the published data, that by using simulations its possible to build predictive-models that accurately estimate the number of trials, subjects and effect-size needed for the psychometric function to find group differences in a particular parameter estimates. This highlights an avenue for researchers building cognitive models to inform others, about their models' strengths and weaknesses in estimating parameters of interest. Lastly with this novel way of generalizing power analyses it is shown that the number of trials in a cognitive science experiment is highly relevant in estimating the psychometric models ability to pick up on group differences in parameters, which is completely neglected by commonly used power analysis soft-wares.

## Introduction of repository


This repository contains all scripts necessary to execute and generate the results, plot, analyses, and the manuscript for the Master thesis. To ensure complete reproducibility, the project is integrated with the Open Science Framework [OSF](https://osf.io/uebmj/) to facilitate either the rerunning of analyses or the retrieval of results, which are located on the [OSF](https://osf.io/uebmj/). 

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
├── Shiny app/                  # Directory containing the scripts for the shiny app, that visualizes the joint modeling interactively.
│   └── ... 
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

To get access to the already run workspace a OSF-token is needed, as the workspaces are stored in the following on OSF due to space limitations of github. 

It is recommended to make an osf folder that contains an osf.txt file which on the first line contains the OSF-token, however make sure that this token is not shared or pushed to github. In the current gitignore an osf folder will be ignored.

## Usage

You can use this repository in various ways. Choose the one that suits your needs: 

1. **Method 1 - Re-create manuscript without rerunning analysis**:
   -  Open the Manuscript folder and knit the [Manuscript.Rmd](./manuscripts/Manuscript.Rmd) file to the desired format (.docx, .pdf or .html). 

2. **Method 2 - Running the Shiny app to better understand the joint modeling introducted in the Uncertainty minimization section**:
   - Open the Shiny app folder and run (knit) the [Joint_model.Rmd](./Shiny app/Joint_model.Rmd) markdown.

3. **Method 3 - Rerunning code and scripts for making the plots of the manuscript**:
   - Navigate to plot scripts directory and go through the different markdowns for re-producing the plots used in the manuscript.
   
4. **Method 4 - Investigating the analysis piple, parallelization and stan models**:
   - Navigate to Analyses directory and choose the section you want to investigate (power analysis, ICC Analysis (uncertainty minmization) and the reanalysis of Legrand.
