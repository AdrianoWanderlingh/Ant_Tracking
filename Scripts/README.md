# Tracking Analysis Scripts Order 📊

**Disclaimer**: Most of the stuff you will be looking for is in `PhD-Ant_Colonies_Tracking_Analysis/scriptsR/EXP1_base_analysis/` for example: create metadata files, define ant tasks, etc. This pipeline is currently quite messy, but at least it exists and is mostly correct. All the needed scripts are present but some paths may need to be changed. I never had the time to polish it but please, feel free to do it if you have the guts!

**Disclaimer 2**: Please read the available guides and the comments in the scripts to understand how things work. If there is no guide, write one! 📜

## A. Pre-processing

- Follow the `pre-processing_Adriano_June2022` guide (find it in the AEL sharepoint)
- Example files are in `Seagate\ Portable\ Disk` or in `DISK4` (Ask Adriano)

1. Create myrmidon files and create the ants

   **Source**: [Define_Ant_Identifications_standalone.R](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/tree/main/scriptsR)

2. Auto-orient, assign capsules

   **Source**: [Data_preparation_after_postprocessing](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/tree/main/scriptsR/EXP1_base_analysis/Data_preparation_after_postprocessing)

## B. Extract the Interactions and Space Use

**Source**: [EXP1_base_analysis.R](https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/tree/main/scriptsR/EXP1_base_analysis)

Dependencies, in the same folder:
- SpaceUse_v082.R
- NetworkProperties_v082.R _(note: Net properties now calculated with Stroeymeyt et al. 2018 pipeline, name obsolete but it contains the interaction detection call)_

These scripts use the `interaction_detection` produced for the Ant Behaviours Classifier.

**Source**: [Automated_Behavioural_Inference](https://github.com/AdrianoWanderlingh/Ant_Tracking/tree/main/Scripts/PhD-Ant_Colonies_Tracking_Analysis/Automated_Behavioural_Inference)

_Note: I'm currently moving to a new, clean, GitHub repo so part of the code is found in `/Ant_Tracking`_

## C. Adapt Outputs to Plug it into Nathalie Stroeymeyt Modified Post-processing and Analysis Pipeline 🔌

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
- The reported "Tag1","tag2" are the new tracking system's "antID1","antID2" (which are used to identify an individual even after retag). I did not report my TagIDs.
- Changed colony naming convention, from “001” to “01B” to describe my R1BS, retaining the colony size information in the filename. This is the only way of preventing that size colonies Big and Small subject to the same treatment will overwrite each other.
- In the METADATA FILES, treatment is not just pathogen or control but "pathogen.big", "pathogen.small","control.big","control.small"

#### Interactions

- Extra Interaction files columns not present in Science files: (added at the end of the file output)
  - "REP_treat": unique replicate*colonySize*treatment identifier
  - "period": "pre", "post"
  - "ant1.zones": nest zone (nest=1,foraging=2) at the time of the interaction, should be always equal to ant2.zones
  - "ant2.zones": see above
  - "duration": interaction duration

#### Task_groups.txt

- Added extra column "treatment" to identify colony, as the colony code alone is not enough. Added REP_treat

#### Treated_worker_list.txt

- Added "treatment" and "REP_treat"

#### Info.txt

- Added "REP_treat"
- Set as NA: Ynestmin, Ynestmax, nb_larvae, nb_pupae, colony_age

#### QPCR_file.txt

- Added "REP_treat"
- Set as NA: age

#### Individual_behavioural_data.txt

- Tab separated instead of space separated




**Enjoy the ride!** 🚀
