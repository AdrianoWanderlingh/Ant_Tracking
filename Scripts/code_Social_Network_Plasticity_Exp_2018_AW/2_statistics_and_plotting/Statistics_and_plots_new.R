#### Statistic and plots #####

#### code adapted from Stroeymeyt et al. 2018

### written by Nathalie Stroeymeyt and Adriano Wanderlingh

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

RUN_UNSCALED_NETS <- F

# NOTES
# plot_untransformed is always TRUE as the the variable fed to the plotting is transformed beforehand (see section: transform variable)

###################################################################
### Calculate missing variables ###################################
###################################################################

### Grooming zone 
root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"

setwd(data_path)
file_list <- list.files(pattern=pattern)
print(file_list)
data <- NULL
for (file in file_list){
  data <- read.table(file,header=T,stringsAsFactors = F)
  #calculate new var and save in the data only if it does not exist
  if(is.null(data$prop_duration_grooming_received_outside_min)){
    data$prop_duration_grooming_received_outside_min <- with(data,duration_grooming_received_min_zone2/(duration_grooming_received_min_zone1+duration_grooming_received_min_zone2) )
    write.table(data, file,col.names=T,row.names=F,quote=F,append=F)
  }
}


###################################################################
###I - Example collective analysis - without rescaling ############
###################################################################

if(RUN_UNSCALED_NETS){

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
}

#######################################################################################################
###II - collective analysis - with rescaling by pre-exposure mean for each size #######################
#######################################################################################################

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="colony_data.txt"
variable_list <- c("modularity","clustering","task_assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
names(variable_list) <- c("modularity","clustering","task assortativity","efficiency","mean degree","degree maximum","density","diameter")
transf_variable_list <- c("log"       ,"sqrt"      ,"none"              ,"log"       ,"power0.01"        ,"log"           ,"log"   ,"power0.01")   ######"none", "sqrt" "log","power2"
# TRANSFORMATION NOTE: task_assortativity is hard to normalise (no transformation has the best result)

coll_rescal_net <- collective_analysis_rescal(data_path)

# Reshape the data
#stats_outcomes_reshaped <- reshape(stats_outcomes[,which(names(stats_outcomes)!="df")], idvar = "predictor", timevar = "variable", direction = "wide")

###################################################################################################################################
###III - individual analysis - for ONE type of individuals (e.g. treated workers only, or queens only) ############################
###################################################################################################################################

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("degree","aggregated_distance_to_queen") 
names(variable_list) <- c("degree","aggregated distance to queen")
transf_variable_list <- c("none"  ,"log")  ######"none", "sqrt" "log","power2"


ind_treated_net <- individual_ONE_analysis(data_path,which_individuals="treated") # "treated","queen","nurse","forager"


### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside","proportion_time_active", "average_bout_speed_pixpersec" ,"total_distance_travelled_pix") #inter_caste_contact_duration?
names(variable_list) <- c("prop. time outside","prop. time active", "average bout speed pixpersec" ,"total distance travelled pix")
transf_variable_list <- c("power0.01"        ,"none"                  ,"log"                          ,"log"                        )  ######"none", "sqrt" "log","power2"

ind_treated_beh <- individual_ONE_analysis(data_path,which_individuals="treated") # "treated","queen","nurse","forager"



#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_received_min", "inter_caste_contact_duration") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2"
names(variable_list) <- c("duration grooming received min", "inter caste grooming duration") # , "prop duration grooming received outside min","duration grooming received outside min"
transf_variable_list <- c("none"                          ,"log"                         )   ######"none", "sqrt" "log","power2"


ind_treated_grooming <- individual_ONE_analysis(data_path,which_individuals="treated") # "treated","queen","nurse","forager"



###################################################################################################################################
###IV - individual analysis - for TWO types of individuals (e.g. nurses and foragers; or "untreated" and "treated") ###############
###################################################################################################################################

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("mean_aggregated_distance_to_treated")
names(variable_list) <- c("mean aggregated distance to treated")
transf_variable_list <- c("power0.01")  ######"none", "sqrt" "log","power2"

warning("Error in get(which_levels) : object 'level_names' not found")
ind_TWO_net <- individual_TWO_analysis(data_path,which_individuals= c("nurse","forager")) ## "treated","queen","nurse","forager"


### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside","inter_caste_contact_duration","duration_of_contact_with_treated_min")
names(variable_list) <- c("prop time outside","inter caste contact duration","duration of contact with treated min")
transf_variable_list <- c("power0.01"        , "power0.01"                       , "power0.01"            )   ######"none", "sqrt" "log","power2"

warning("Error in get(which_levels) : object 'level_names' not found WHEN inter_caste_contact_duration and duration_of_contact_with_treated_min are transformed")
ind_TWO_beh <- individual_TWO_analysis(data_path,which_individuals= c("nurse","forager")) ## "treated","queen","nurse","forager"
# to fix: double letters assigned when cols are split by task

# data: main_experiment/network_properties/individual_data
# Mean distance to treated



#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_given_to_treated_min")
names(variable_list) <- c("duration grooming given to treated min")
transf_variable_list <- c("none"        )   ######"none", "sqrt" "log","power2"

warning("Error in get(which_levels) : object 'level_names' not found")
ind_TWO_beh <- individual_TWO_analysis(data_path,which_individuals= c("nurse","forager")) ## "treated","queen","nurse","forager"











warning("plot hist total and split by size!!!!!")






#TEMP (The green line highlights the value at which there is an inversion in the sign of the density difference. This value was used as athreshold to distinguish high from low simulated loads in allsubsequent analyses)
high_threshold <- 0.0411
  
# logging not enough for qPCR data, as it has kurtois = 27.28894
# sqrt                                                  27.46589
# power0.01                                             70.34206
# power0.001                                            26.71343
# power0.0001                                           26.64026

## issue in functions_new.R, line 2321

######Open pdf Figure 2 #####
#pdf(file=paste(figurefolder,"/Figure2.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=page_height*0.25,pointsize=pointsize_less_than_2row2col)
######Plot qPCR data #####
full_statuses_names_ori <- full_statuses_names
full_statuses_names[full_statuses_names%in%c("Nurses","Untreated\nnurses")] <- "Nurses\n"

statuses_colours_ori <- statuses_colours
statuses_colours[names(statuses_colours)%in%c("queen","forager","nurse")] <- "black"

par(pars)
par_mar_ori <- par()$mar
par(mar=par_mar_ori+c(1,0,1,1))
widz <- c(2,1.5)
layout(matrix(c(1,2),nrow=1),widths=widz)
translated_high_threshold <- plot_qpcr(experiments=c("main_experiment"))
to_keep <- c(to_keep,"translated_high_threshold")
####Add letters ####
par(xpd=NA)
x_text1 <- grconvertX(1/80, from='ndc');x_text2 <- grconvertX(widz[1]/sum(widz)+1/80, from='ndc')
y_text <- grconvertY(0.97, from='ndc')
text(x_text1,y_text,labels=panel_casse("a"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
text(x_text2,y_text,labels=panel_casse("b"),font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
par(xpd=F)
full_statuses_names <- full_statuses_names_ori
statuses_colours <- statuses_colours_ori
par(mar=par_mar_ori)
####Close figure 2#####
#dev.off()
######## clean before next step###
clean();
Sys.sleep(2)


