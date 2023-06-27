# Ant Behaviours Classifier
### Automated_Behavioural_Inference

The Automated_Behavioural_Inference project is a comprehensive suite of scripts and functions for analyzing and inferring behaviours from tracking data of ant colonies. It provides functionalities for extracting movement variables, detecting and filtering interactions, predicting grooming behaviour, and evaluating the performance of machine learning classifiers. The project aims to automate the analysis process and provide insights into ant behaviours, improving efficiency and accuracy in studying ant colonies.

<p align="center">
<img src="https://github.com/AdrianoWanderlingh/Ant_Tracking/assets/47888790/bbb71348-f135-4b52-a287-846a728d4f8b" />
</p>




_Code written by Adriano Wanderlingh, Nathalie Stroeymeyt and Tom Richardson with contributions by Enrico Gavagnign_

## Main Scripts

### STEP1_BEH_MAIN_behaviours_analysis_fort081.R

This script performs behaviour analysis on ant interactions. It includes functionalities for loading libraries and functions, setting parameters, loading and modifying annotation data, splitting the data into training and test sets, extracting movement variables, performing manual and automatic agreement analysis, and fitting classifiers using LDA/SVM/RF. The script also includes loops to vary parameters such as capsule shapes, maximum interaction gaps, distance thresholds, frame thresholds, and trim lengths, and saves the output in separate directories for each parameter combination. The script aims to evaluate the quality of interaction parameters and classify interactions based on precision and sensitivity.

### STEP2_evaluating_classifiers.R

This script reads quality score files, filters the data based on specified parameters, evaluates the best parameter values using statistical modelling and Tukey posthoc comparisons, and selects the best-performing methods based on specified measures. The script is designed to assist in selecting the best methods for machine learning based on quality scores.

### STEP2b_Final_Classifier_Validation.R

This script performs various tasks related to data analysis and prediction of interactions in a behavioural dataset. It includes functions for cleaning and modifying the dataset, as well as for extracting movement variables and predicting interactions using a chosen classifier. The script also splits the dataset into training and test sets and evaluates the quality of the predictions using different measures. Overall, the script aims to provide insights into the accuracy of the interaction predictions and the performance of the chosen classifier.

### STEP3_detect_grooming_on_all_data.R

This script is for analyzing and inferring ant behaviours from myrmidon tracking data. It includes functions for extracting movement variables, detecting interactions, predicting grooming behaviour using a classification model, and saving the inferred grooming events. The script allows for the customization of parameters such as the behaviour of interest, frame rate, and measurement criteria. It also includes the sourcing of various function scripts and libraries necessary for the analysis.

## Functions and Libraries

### BEH_self_defined_functions.R

This script contains several functions that serve different purposes, ranging from mathematical calculations to data visualization and file manipulation.

### interaction_detection.R

This script contains several functions related to detecting and filtering interactions between objects in a given experiment.

### BEH_PCA_fort081_FUNCTIONS.R

This script contains functions for performing classification analysis on a dataset, predicting class labels, and calculating quality metrics.

### trajectory_extraction.R

This script consists of two functions: "merge_trajectory_segments" and "extract_trajectories", which are used for merging trajectory segments and extracting trajectory segments from a given start to end time.

### add_angles.cpp

This C++ function takes a DataFrame object as input and performs calculations on the data to generate several derived variables.

### BEH_libraries.R

This script loads multiple libraries in R that are used for various purposes.

### BEH_Auto_Man_agreement_matrix_fort081_FUNCTIONS.R

This script consists of two functions: "make_binary_interaction_matrix" and "auto_manual_agreement", which are used for generating a binary interaction matrix and calculating the agreement between automatic and manual interactions.

### merge_interactions.cpp

This C++ implementation merges and summarizes interactions between pairs of entities (e.g., ants) based on a given set of collisions.

###
