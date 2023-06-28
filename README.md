
<div align="center">
<h1 align="center">
<img src="https://user-images.githubusercontent.com/47888790/184152615-5b94905a-9ddc-4d8f-8f27-3b10a2fda858.png" width="300" />
<br>
Ant Tracking
</h1>

<p align="center">
<img src="https://img.shields.io/badge/R-276DC3.svg?style&logo=R&logoColor=white" alt="R" />
<img src="https://img.shields.io/badge/C++-00599C?style=flat-square&logo=C%2B%2B&logoColor=white" alt="C++" />
<!--- <img src="https://img.shields.io/badge/FORmicidae%20Tracker-8A2BE2" alt="formicidae-tracker" /> --->
</p>

[![formicidae-tracker](https://img.shields.io/badge/FORmicidae%20Tracker-8A2BE2)](https://github.com/formicidae-tracker/myrmidon)


![GitHub top language](https://img.shields.io/github/languages/top/AdrianoWanderlingh/Ant_Tracking?style&color=5D6D7E)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/AdrianoWanderlingh/Ant_Tracking?style&color=5D6D7E)
![GitHub license](https://img.shields.io/github/license/AdrianoWanderlingh/Ant_Tracking?style&color=5D6D7E)
</div>







### ⚠️ Warning: Work in Progress
>This is a re-organized version of [PhD-exp1-data-analysis](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis). Please note that this repository is still being updated. Older versions or missing files can be found in the original repository.


## 📍 Overview
Code used for Adriano Wanderlingh's Lasius niger Colony Maturation Experiment.  :ant: 🦠

Testing Social Immunity role in disease hindering dynamics in societies, using ant colonies interaction networks and behaviour | experimental epidemiology.



<strong>R and C++ scripts for:</strong>
- Tracking Data Post-processing
- Geometrical capsules assignment
- Interactions detection
- Machine Learning Behavioural Classification
- Personal immunity investment and pathogen load quantification
- RTqPCR data pre-processing pipeline with MNAR data simulations
- extra related material


The repository is organized into two main folders: `Scripts` and `Data`. The `Scripts` folder contains R and C++ code for data processing, analysis, and visualization. The `Data` folder contains raw and processed data files.

## Table of Contents

- [Folder Structure](#folder-structure)
- [Scripts](#scripts)
  - [PhD-Ant_Colonies_Tracking_Analysis](#phd-ant_colonies_tracking_analysis)
    - [Automated_Behavioural_Inference](#automated_behavioural_inference)
    - [molecular_bio_analysis](#molecular_bio_analysis)
  - [code_Social_Network_Plasticity_Exp_2018_AW](#code_social_network_plasticity_exp_2018_aw)
- [Data](#data)
  - [PhD-Ant_Colonies_Tracking_Analysis](#phd-ant_colonies_tracking_analysis-data)
    - [Automated_Behavioural_Inference](#automated_behavioural_inference-data)
    - [molecular_bio_analysis](#molecular_bio_analysis-data)

## 📂 Folder Structure

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


#### Data_preparation_after_postprocessing

For more information, read the folder [README](https://github.com/AdrianoWanderlingh/Ant_Tracking/tree/main/Scripts) and the [pre-processing_Adriano_June2022](https://uob.sharepoint.com/:w:/r/teams/grp-AntsEpidemiologyLab/_layouts/15/Doc.aspx?sourcedoc=%7B2562631B-A6E5-4289-907F-89502F6C27E6%7D&file=pre-processing_Adriano_June2022.docx&action=default&mobileredirect=true&DefaultItemOpen=1) guide.


#### Automated_Behavioural_Inference

For more information, read the folder [README](https://github.com/AdrianoWanderlingh/Ant_Tracking/tree/main/Scripts/PhD-Ant_Colonies_Tracking_Analysis/Automated_Behavioural_Inference)


#### molecular_bio_analysis

This folder contains scripts for molecular biology analysis, including qPCR and RT-qPCR analysis.


### code_Social_Network_Plasticity_Exp_2018_AW

This folder contains scripts for data post-processing and statistical analysis of ant behavioural interactions, adapted from [Stroeymeyt et al., 2018](https://www.science.org/doi/epdf/10.1126/science.aat4793).


## Data

### PhD-Ant_Colonies_Tracking_Analysis

#### Automated_Behavioural_Inference

- `Mean_ant_length_per_TrackingSystem.txt`: [Information needed]

#### molecular_bio_analysis

- `Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv`: [Information needed]
- `Adriano_qPCR_pathogen_load_MASTER_REPORT.csv`: [Information needed]
- 


## 👏 Acknowledgments

> - Code written by Adriano Wanderlingh, Nathalie Stroeymeyt and Tom Richardson with contributions by Enrico Gavagnign. Ant Epidemiology Stroeymeyt Lab, University of Bristol, UK



#
<small><i>This file was generated by extracting information on the repo with [README-AI](https://github.com/eli64s/README-AI) and using `find . d > dirs.txt` to get the folder structure, then feeding both to GPT-4, which was prompted to write a README.md file.</i></small>
