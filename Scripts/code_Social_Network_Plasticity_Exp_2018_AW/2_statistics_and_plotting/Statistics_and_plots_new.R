rm(list=ls())

USER <- "2A13_Office" # Nath_office 

if (USER == "2A13_Office") {
  usr <- "cf19810"
  source_path <- paste0("/home/",usr,"/Ant_Tracking/Scripts/code_Social_Network_Plasticity_Exp_2018_AW/2_statistics_and_plotting/source")
} else {
  usr <- "bzniks"
  warning("the up to date code is in https://github.com/AdrianoWanderlingh/PhD-Ant_Colonies_Tracking_Analysis/. \n\nPlease use git clone <URL> to download")
  #source_path  <- "~/Dropbox/Papers/4_Lausanne/5_OI/2_Science/SHARED_DATA_SCIENCE_REPOSITORIES/code/2_statistics_and_plotting_Adriano/source"
}


#####Overall parameters and functions ##########
####define folders
figurefolder <- "~/figures"
disk_path    <- paste0("/media/",usr,"/DISK4/Lasius-Bristol_pathogen_experiment")

treatment_order <- c("control.small","control.big","pathogen.small","pathogen.big")
#treatment_labs_order <- c("control small","control large","pathogen small","pathogen large")
exposure_order  <- c("control","pathogen")
size_order      <- c("small","big")
period_order    <- c("pre","post")
task_group_order <- c("queen","nurse","forager","untreated","treated")

####source programs
source(paste(source_path,"/libraries.R",sep=""))
source(paste(source_path,"/plotting_parameters.R",sep=""))
source(paste(source_path,"/functions_new.R",sep=""))
source(paste(source_path,"/analysis_parameters.R",sep=""))

# NOTES
# plot_untransformed is always TRUE as the the variable fed to the plotting is transformed beforehand (see section: transform variable)

###################################################################
###I - Example collective analysis - without rescaling ############
###################################################################

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="colony_data.txt"
variable_list <- c("modularity","clustering","task_assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
names(variable_list) <- c("modularity","clustering","task assortativity","efficiency","mean degree","degree maximum","density","diameter")
transf_variable_list <- c("log"      ,"none"      ,"none"         ,"log"      ,"log"       ,"log"          ,"none"   ,"log")   ######"none", "sqrt" "log","power2"

coll_no_rescal_net <- collective_analysis_no_rescal(data_path)

# to print a plot
#coll_no_rescal$barplot_delta_collective_list$modularity

# Reshape the data
#stats_outcomes_reshaped <- reshape(stats_outcomes[,which(names(stats_outcomes)!="df")], idvar = "predictor", timevar = "variable", direction = "wide")
#colnames(reshaped_data)[-1] <- gsub("pval.", "", colnames(reshaped_data)[-1])


#######################################################################################################
###II - Example collective analysis - with rescaling by pre-exposure mean for each size##############
#######################################################################################################

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="colony_data.txt"
variable_list <- c("modularity","clustering","task_assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
names(variable_list) <- c("modularity","clustering","task assortativity","efficiency","mean degree","degree maximum","density","diameter")
transf_variable_list <- c("log"       ,"sqrt"      ,"none"          ,"log"       ,"log"        ,"log"           ,"log"   ,"power0.01")   ######"none", "sqrt" "log","power2"
# NOTE: task_assortativity is hard to normalise because of the bimodal distribution (small vs big cols)

coll_rescal_net <- collective_analysis_rescal(data_path)

# Reshape the data
#stats_outcomes_reshaped <- reshape(stats_outcomes[,which(names(stats_outcomes)!="df")], idvar = "predictor", timevar = "variable", direction = "wide")

###################################################################################################################################
###III - Example individual analysis - for ONE type of individuals (e.g. treated workers only, or queens only)#####################
###################################################################################################################################

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("degree","aggregated_distance_to_queen") 
names(variable_list) <- c("degree","aggregated distance to queen")
transf_variable_list <- c("none"  ,"log")  ######"none", "sqrt" "log","power2"


ind_treated_net <- individual_ONE_analysis_rescal(data_path,which_individuals="treated") # "treated","queen","nurse","forager"


### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside","proportion_time_active", "average_bout_speed_pixpersec" ,"total_distance_travelled_pix") #inter_caste_contact_duration?
names(variable_list) <- c("prop. time outside","prop. time active", "average bout speed pixpersec" ,"total distance travelled pix")
transf_variable_list <- c("power0.01"        ,"none"                  ,"none"                          ,"none"                        )  ######"none", "sqrt" "log","power2"

ind_treated_beh <- individual_ONE_analysis_rescal(data_path,which_individuals="treated") # "treated","queen","nurse","forager"



#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("inter_caste_contact_duration", "duration_grooming_received_min_zone1","duration_grooming_received_min_zone2 ") #ZONE GROOMING HAVE TO BE COMBINED SOMEHOW (MAYBE EASIST IS PROP_REC_OUTISE = ZONE2/(ZONE1+ZONE2) )
names(variable_list) <- c("prop. time outside","prop. time active", "average bout speed pixpersec" ,"total distance travelled pix","inter caste contact_duration")
transf_variable_list <- c("power0.01"        ,"none"                  ,"none"                          ,"none"                        ,"none"      )   ######"none", "sqrt" "log","power2"


#ind_treated_grooming <- individual_ONE_analysis_rescal(data_path,which_individuals="treated") # "treated","queen","nurse","forager"



###################################################################################################################################
###IV - Example individual analysis - for TWO types of individuals (e.g. nurses and foragers; or "untreated" and "treated") #######
###################################################################################################################################

#### ALL INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside","inter_caste_contact_duration","duration_of_contact_with_treated_min")
names(variable_list) <- c("prop_time_outside","inter_caste_contact_duration","duration_of_contact_with_treated_min")
transf_variable_list <- c("power0.01"        , "none"                       , "none"            )   ######"none", "sqrt" "log","power2"


ind_TWO_analysis <- individual_TWO_analysis_rescal(data_path,which_individuals= c("nurse","forager")) ## "treated","queen","nurse","forager"

# data: main_experiment/network_properties/individual_data
# Mean distance to treated



#### GROOMING INTERACTIONS ####

# data: main_experiment_grooming/individual_behaviour/individual_behavioural_data 
# duration_grooming_given_to_treated_min 
