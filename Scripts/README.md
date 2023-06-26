# Tracking Analysis Scripts Order 📊

This file will help you to run the analyses on the tracking data from start to finish. The pre-processing steps have been written for tracking data produced by Bristol's new tracking system, the post-processing has been adapted from [Stroeymeyt et al., 2018](https://www.science.org/doi/epdf/10.1126/science.aat4793). The final statistical analysis steps will change according to your experiment structure and, currently, they are designed to fit Adriano's data (2 size treatments x 2 exposure treatments).

  **Disclaimer**: Most of the base stuff you will be looking for is in [here](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/tree/main/scriptsR/EXP1_base_analysis), for example: create metadata files, define ant tasks, etc. This pipeline is currently quite messy, but at least it exists and is mostly correct. All the needed scripts are present but some paths may need to be changed. I never had the time to polish it but please, feel free to do it if you have the guts!

  **Disclaimer 2**: Please read the available guides and the comments in the scripts to understand how things work. Always check where things are stored before running scripts blindly. If there is no guide, write one! 📜

:octocat: To download the content of a git subfolder instead of a full Repository, either use [gitdir](https://github.com/sdushantha/gitdir) (linux terminal) or this [easy web tool](https://download-directory.github.io/).

## A. Pre-processing

- Follow the [pre-processing_Adriano_June2022](https://uob.sharepoint.com/:w:/r/teams/grp-AntsEpidemiologyLab/_layouts/15/Doc.aspx?sourcedoc=%7B2562631B-A6E5-4289-907F-89502F6C27E6%7D&file=pre-processing_Adriano_June2022.docx&action=default&mobileredirect=true&DefaultItemOpen=1) guide. Please update this file as it is in a terrible state. 
- Example files are in `Seagate\ Portable\ Disk` or in `DISK4` (Ask Adriano)

1. Create myrmidon files and create the ants

   **Source**: [Define_Ant_Identifications_standalone.R](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/tree/main/scriptsR)

2. Auto-orient, assign capsules

   **Source**: [Data_preparation_after_postprocessing](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/tree/main/scriptsR/EXP1_base_analysis/Data_preparation_after_postprocessing)

    **Files needed**: Capsule files

 The capsules defined by Adriano are in either in `Seagate\ Portable\ Disk` or `DISK4` `/ADRIANO/EXPERIMENT_DATA/CURRENT_CAPS_DEF_BASE_FILES_MAN_ORIENTED`. The .myrmidon files need to be moved at the same folder level as their relative tracking data, in this case either `/ADRIANO/EXPERIMENT_DATA` `/REP3` or `/REP9`. The currently present capsules are for the general definition of interactions for all analytical purposes (CapsuleDef2018) or specifically designed to test the Ant Behaviours Classifier (CapsuleDef 3 to 12).


## B. Extract the Interactions and Space Use

 To ensure that the outputs of this step are stored correctly and in the same format as Stroeymeyt et al. (2018) for the following analyses, it is necessary to use the same folder structure. I created  the [folder structure](https://github.com/AdrianoWanderlingh/Ant_Tracking/blob/main/Scripts/code_Social_Network_Plasticity_Exp_2018_AW/dirs.txt) with `find . -type d > dirs.txt` and you can recreate it in your destination of choice with `xargs mkdir -p < dirs.txt` (in the linux terminal).

  **Source**: [EXP1_base_analysis.R](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/tree/main/scriptsR/EXP1_base_analysis)

Dependencies, in the same folder:
- SpaceUse_v082.R
- NetworkProperties_v082.R _(note: Net properties now calculated with Stroeymeyt et al. 2018 pipeline, name obsolete but it contains the interaction detection call)_

These scripts use the `interaction_detection` produced for the Ant Behaviours Classifier.

  **Source**: [Automated_Behavioural_Inference](https://github.com/AdrianoWanderlingh/Ant_Tracking/tree/main/Scripts/PhD-Ant_Colonies_Tracking_Analysis/Automated_Behavioural_Inference)

_Note: I'm currently moving to a new, clean, GitHub repo so part of the code is found in `/Ant_Tracking`_


   **Files needed**
   | File | Source |
|------|--------|
| [Metadata summary file](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/blob/main/EXP_summary_data/Metadata_Exp1_2021_2023-02-27.txt) | [metadata extraction script](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/blob/main/scriptsR/EXP1_base_analysis/Extract_Metadata_v082.R) |
| [Mean ant lenght per TS](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/blob/main/EXP_summary_data/Mean_ant_length_per_TrackingSystem.txt) | [Mean ant length script](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/blob/main/scriptsR/EXP1_base_analysis/Data_preparation_after_postprocessing/MOD_OF_DataPrep2_FOR_Mean_ant_length_per_TrackingSystem_v082.R) |
| [Return times after treatment](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/blob/main/EXP_summary_data/Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU_All_colonies.csv) | Compile to have, at least, the columns REP_treat and ReturnExposed_time |

## C. Adapt Outputs to Plug it into Stroeymeyt et al., 2018 Modified Pipeline 🔌

- Example files are in `Seagate\ Portable\ Disk` or in `DISK4` (Ask Adriano) or, for the original files, in `New Passport` (Ask Nathalie)

  **Source**: [EXP1_base_analysis](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/blob/main/scriptsR/EXP1_base_analysis/)

- Tables_to_Match_Stroeymeyt2018_Pipeline.R
- Convert_grooming_to_Match_Stroeymeyt2018_Pipeline.R _(note: needed only if you used the Ant Behaviours Classifier)_

## D. Run Post-processing and Analysis Pipeline

[Code_Social_Network_Plasticity_Exp_2018_AW](https://github.com/AdrianoWanderlingh/Ant_Tracking/tree/main/Scripts/code_Social_Network_Plasticity_Exp_2018_AW)


### Processing Time ⏲️

Sample processing time for Adriano's Experiment (44 colonies, 48h of tracking data per colony):

| Task                                      | Processing Time                                         |
|:------------------------------------------|:--------------------------------------------------------|
| Capsules copy                             | Few hours                                               |
| individual_behavioural_data.txt (from Exp1_base_analysis.R : SpaceUse) | 20 hours                 |
| Produce interactions list (from Exp1_base_analysis.R : Interactions) | 13 hours                  |
| 11_randomise_interactions (if 6/7 instances are run in parallel) | ~ 8 hours                 |
| 12_simulate_transmission (4 seed files)   | ~ 50 hours                                              |
| 13_network_analysis                       | Estimated: 20 hours (15 mins if performed only on observed) |
| 14_summarise_interactions.R               | Estimated: 40 hours (15 mins if performed only on observed) |




### Log of differences of the code/files from Science2018 Data Analysis Pipeline

#### General

- The pipeline now includes some adjustments and extra measures to handle grooming data
- Following the naming convention of Science2018, the reported "Tag" is Bristol's tracking system's "antID" (which are used to identify an individual even after retag). I did not report my TagIDs.
- Changed colony naming convention, from “001” to “01B” to describe my R1BS, retaining the colony size information in the filename. This is the only way of preventing that size colonies Big and Small subject to the same treatment will overwrite each other.
- In the METADATA FILES, treatment is not just pathogen or control but "pathogen.big", "pathogen.small","control.big","control.small"

#### _Files produced by EXP1_base_analysis.R_

##### Interactions

- Extra Interaction files columns not present in Science files: (added at the end of the file output)
  - "REP_treat": unique replicate*colonySize*treatment identifier
  - "period": "pre", "post"
  - "ant1.zones": nest zone (nest=1,foraging=2) for AntID1 (now: Tag1) at the time of the interaction between a pair of ants, will always equal to ant2.zones
  - "ant2.zones": see above
  - "duration": interaction duration

##### Space Use

In our case, **individual_behavioural_data.txt** is generated by the EXP1_base_analysis.R script (SpaceUse), which is much richer than the same file produced by Nathalie in 6_time_investment.R as it already includes info that in her pipeline is generated by 8_process_trajectory_files.R. **network_position_vs_time_outside.dat** is obtained from Individual_behavioural_data.txt keeping only period = pre. 

#### _Files produced by Tables_to_Match_Stroeymeyt2018_Pipeline.R_

##### Task_groups.txt

- Added extra column "treatment" to identify colony, as the colony code alone is not enough. Added REP_treat

##### Treated_worker_list.txt

- Added "treatment" and "REP_treat"

##### Info.txt

- Added "REP_treat"
- Set as NA: Ynestmin, Ynestmax, nb_larvae, nb_pupae, colony_age

##### QPCR_file.txt

- Added "REP_treat"
- Set as NA: age

##### Individual_behavioural_data.txt

- Tab-separated instead of space separated




**Enjoy the ride!** 🚀
