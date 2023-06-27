<p align="center">
<img src="https://user-images.githubusercontent.com/47888790/184152615-5b94905a-9ddc-4d8f-8f27-3b10a2fda858.png" width="300" />
</p>

# Ant_Tracking 
### ⚠️ Warning: Ant_Tracking (Work in Progress)
This is a re-organized version of [PhD-exp1-data-analysis](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis). Please note that this repository is still being updated. Older versions or missing files can be found in the original repository.

<strong>Lasius niger Colony Maturation Experiment</strong> :ant: 🦠

Testing Social Immunity role in disease hindering dynamics in societies, using ant colonies interaction networks and behaviour | experimental epidemiology



<strong>R and C++ scripts for:</strong>
- Tracking Data Post-processing
- Geometrical capsules assignment
- Interactions detection
- Machine Learning Behavioural Classification
- Personal immunity investment and pathogen load quantification
- RTqPCR data pre-processing pipeline with MNAR data simulations
- extra related material

<em>Code written by Adriano Wanderlingh, Nathalie Stroeymeyt and Tom Richardson with contributions by Enrico Gavagnign</em>

Ant Epidemiology Lab


#

Ant_Tracking is a GitHub repository containing code and data for tracking and analyzing ant colonies. The repository is organized into two main folders: `Scripts` and `Data`. The `Scripts` folder contains R and C++ code for data processing, analysis, and visualization. The `Data` folder contains raw and processed data files.

## Table of Contents

- [Folder Structure](#folder-structure)
- [Scripts](#scripts)
  - [PhD-Ant_Colonies_Tracking_Analysis](#phd-ant_colonies_tracking_analysis)
    - [Automated_Behavioural_Inference](#automated_behavioural_inference)
    - [molecular_bio_analysis](#molecular_bio_analysis)
  - [code_Social_Network_Plasticity_Exp_2018_AW](#code_social_network_plasticity_exp_2018_aw)
    - [1_data_post_processing](#1_data_post_processing)
    - [2_statistics_and_plotting](#2_statistics_and_plotting)
- [Data](#data)
  - [PhD-Ant_Colonies_Tracking_Analysis](#phd-ant_colonies_tracking_analysis-data)
    - [Automated_Behavioural_Inference](#automated_behavioural_inference-data)
    - [molecular_bio_analysis](#molecular_bio_analysis-data)

## Folder Structure

```
.
├── Scripts
│   ├── PhD-Ant_Colonies_Tracking_Analysis
│   │   ├── Automated_Behavioural_Inference
│   │   └── molecular_bio_analysis
│   └── code_Social_Network_Plasticity_Exp_2018_AW
│       ├── 1_data_post_processing
│       │   └── source
│       └── 2_statistics_and_plotting
│           └── source
└── Data
    └── PhD-Ant_Colonies_Tracking_Analysis
        ├── Automated_Behavioural_Inference
        └── molecular_bio_analysis
```

## Scripts

### PhD-Ant_Colonies_Tracking_Analysis

#### Automated_Behavioural_Inference (UPDATE!)

This folder contains scripts for automated behavioural inference, including interaction detection, trajectory extraction, and classifier evaluation. Key scripts include:

- `BEH_self_defined_functions.R`: Contains several functions for angle calculations, classification performance, output text file generation, and filtering short interactions.
- `interaction_detection.R`: Detects and filters interactions between ants based on time, distance, angle, and ant IDs.
- `trajectory_extraction.R`: Extracts ant trajectories from a dataset and stores them in successive segments.
- `STEP2_evaluating_classifiers.R`: Evaluates the best parameter values for training and test using Tukey post-hoc comparisons.
- `add_angles.cpp`: Calculates various angle measures for a given trajectory.
- `BEH_libraries.R`: Loads multiple R packages/libraries for various purposes.
- `BEH_Auto_Man_agreement_matrix_fort081_FUNCTIONS.R`: Calculates the agreement between automatic and manual behavioural annotations.
- `ClassifierScores_check_plots.R`: Generates visualizations comparing data from two different sources and scatterplots of candidate grooming events.
- `Classifier_preformance_and_comparison_plots.R`: Fits the selected model, plots the model's error rate and variable importance, and creates a precision versus sensitivity plot for all methods.

#### Data_preparation_after_postprocessing

see [here](https://github.com/AdrianoWanderlingh/Ant_Tracking/tree/main/Scripts) and the [pre-processing_Adriano_June2022](https://uob.sharepoint.com/:w:/r/teams/grp-AntsEpidemiologyLab/_layouts/15/Doc.aspx?sourcedoc=%7B2562631B-A6E5-4289-907F-89502F6C27E6%7D&file=pre-processing_Adriano_June2022.docx&action=default&mobileredirect=true&DefaultItemOpen=1) guide.

#### molecular_bio_analysis

This folder contains scripts for molecular biology analysis, including qPCR and RT-qPCR analysis.

### code_Social_Network_Plasticity_Exp_2018_AW

This folder contains scripts for data post-processing and statistical analysis of ant behavioural interactions, adapted from [Stroeymeyt et al., 2018](https://www.science.org/doi/epdf/10.1126/science.aat4793).

#### 1_data_post_processing

- `A_Main_experiment.R`: Main script for data processing and analysis of ant behavioural interactions.
- `source`: Contains various R and C++ scripts for data processing, including trajectory analysis, interaction detection, transmission simulation, network analysis, and heatmap generation.

#### 2_statistics_and_plotting

- `Statistics_and_plots.R`: the original file
- `Statistics_and_plots_new.R`: a modified version to work with Adriano's data structure
- `source`: Contains various R scripts for statistical analysis and plotting, including libraries, functions, analysis parameters, and plotting parameters.

## Data

### PhD-Ant_Colonies_Tracking_Analysis

#### Automated_Behavioural_Inference

- `Mean_ant_length_per_TrackingSystem.txt`: [Information needed]

#### molecular_bio_analysis

- `Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv`: [Information needed]
- `Adriano_qPCR_pathogen_load_MASTER_REPORT.csv`: [Information needed]

#
<small><i>This file was generated by extracting information on the repo with [README-AI](https://github.com/eli64s/README-AI) and using `find . d > dirs.txt` to get the folder structure, then feeding both to GPT-4, which was prompted to write a README.md file.</i></small>
