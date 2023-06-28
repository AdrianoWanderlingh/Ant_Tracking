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
variable_list <- c("modularity","clustering","task_assortativity","efficiency","degree_mean","density") #,"degree_maximum","diameter"
names(variable_list) <- c("modularity","clustering","task assortativity","efficiency","mean degree","density") # ,"degree maximum","diameter"
transf_variable_list <- c("log"       ,"sqrt"      ,"Box_Cox"              ,"log"       ,"power0.01"        ,"log" ) # ,"log"   ,"power0.01"  ######"none", "sqrt" "log","power2"
# TRANSFORMATION NOTE: task_assortativity is hard to normalise (no transformation has the best result) - used Box_Cox

coll_rescal_net <- collective_analysis_rescal(data_path)

# Reshape the data
#stats_outcomes_reshaped <- reshape(stats_outcomes[,which(names(stats_outcomes)!="df")], idvar = "predictor", timevar = "variable", direction = "wide")

###################################################################################################################################
###III - individual analysis - for ONE type of individuals (e.g. treated workers only, or queens only) ############################
###################################################################################################################################

################    FOR TREATED INDIVIDUALS   ################

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("degree")#,"aggregated_distance_to_queen") 
names(variable_list) <- c("degree")#,"aggregated distance to queen")
transf_variable_list <- c("none"  )#,"log")  ######"none", "sqrt" "log","power2"


ind_treated_net <- individual_ONE_analysis(data_path,which_individuals="treated") # "treated","queen","nurse","forager"


### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside") #,"proportion_time_active", "average_bout_speed_pixpersec" ,"total_distance_travelled_pix", "inter_caste_contact_duration") #inter_caste_contact_duration?
names(variable_list) <- c("prop. time outside") # ,"prop. time active", "average bout speed pixpersec" ,"total distance travelled pix", "inter caste contact duration")
transf_variable_list <- c("power0.01"        )#,"none"                  ,"log"                          ,"log"                      ,"sqrt"                   )  ######"none", "sqrt" "log","power2"

ind_treated_beh <- individual_ONE_analysis(data_path,which_individuals="treated") # "treated","queen","nurse","forager"

line_plot(data_path,which_individuals="treated")


#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_received_min", "inter_caste_contact_duration","prop_duration_grooming_received_outside_min") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2"
names(variable_list) <- c("duration grooming received min", "inter caste grooming duration","prop duration grooming received outside_min") # , "prop duration grooming received outside min","duration grooming received outside min"
transf_variable_list <- c("log"                          ,"log"                         , "Box_Cox")   ######"none", "sqrt" "log","power2"


ind_treated_grooming <- individual_ONE_analysis(data_path,which_individuals="treated") # "treated","queen","nurse","forager"

line_plot(data_path,which_individuals="treated")

################    FOR UNTREATED INDIVIDUALS    ################    


#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("degree")#, "mean_aggregated_distance_to_treated")
names(variable_list) <- c("degree")#,"mean aggregated distance to treated")
transf_variable_list <- c("none")#,"power0.01")  ######"none", "sqrt" "log","power2"

ind_untreated_net_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse") ## "treated","queen","nurse","forager"
ind_untreated_net_forag <- individual_ONE_analysis(data_path,which_individuals="forager") ## "treated","queen","nurse","forager"



### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside","inter_caste_contact_duration","duration_of_contact_with_treated_min")
names(variable_list) <- c("prop time outside","inter caste contact duration","duration of contact with treated min")
transf_variable_list <- c("power0.01"        , "Box_Cox"                    , "Box_Cox"            )   ######"none", "sqrt" "log","power2"

ind_untreated_beh_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse") ## "treated","queen","nurse","forager"
ind_untreated_beh_forag <- individual_ONE_analysis(data_path,which_individuals="forager") ## "treated","queen","nurse","forager"


#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_given_to_treated_min")
names(variable_list) <- c("duration grooming given to treated min")
transf_variable_list <- c("log"        )   ######"none", "sqrt" "log","power2"

ind_untreated_grooming_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse") ## "treated","queen","nurse","forager"
ind_untreated_grooming_forag <- individual_ONE_analysis(data_path,which_individuals="forager") ## "treated","queen","nurse","forager"


###################################################################################################################################
###IV - individual analysis - for TWO types of individuals (e.g. nurses and foragers; or "untreated" and "treated") ###############
###################################################################################################################################

# #### ALL INTERACTIONS ####
# 
# ### network properties
# root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
# data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep="")
# pattern="individual_data"
# variable_list <-        c("degree")#, "mean_aggregated_distance_to_treated")
# names(variable_list) <- c("degree")#,"mean aggregated distance to treated")
# transf_variable_list <- c("none")#,"power0.01")  ######"none", "sqrt" "log","power2"
# 
# ind_TWO_net <- individual_TWO_analysis(data_path,which_individuals= c("nurse","forager")) ## "treated","queen","nurse","forager"
# 
# 
# ### behavioural data
# root_path <- paste(disk_path,"/main_experiment",sep="")
# data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
# pattern="individual_behavioural_data"
# variable_list <-        c("prop_time_outside","inter_caste_contact_duration","duration_of_contact_with_treated_min")
# names(variable_list) <- c("prop time outside","inter caste contact duration","duration of contact with treated min")
# transf_variable_list <- c("power0.01"        , "power0.01"                       , "power0.01"            )   ######"none", "sqrt" "log","power2"
# 
# ind_TWO_beh <- individual_TWO_analysis(data_path,which_individuals= c("nurse","forager")) ## "treated","queen","nurse","forager"
# # to fix: double letters assigned when cols are split by task
# 
# 
# #### GROOMING INTERACTIONS ####
# root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
# data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
# pattern="individual_behavioural_data"
# variable_list <-        c("duration_grooming_given_to_treated_min")
# names(variable_list) <- c("duration grooming given to treated min")
# transf_variable_list <- c("none"        )   ######"none", "sqrt" "log","power2"
# 
# ind_TWO_beh <- individual_TWO_analysis(data_path,which_individuals= c("nurse","forager")) ## "treated","queen","nurse","forager"
# 




###################################################################################################################################
### PLOT GRIDS ####################################################################################################################
###################################################################################################################################

warning("\n-LEGEND SHOULD HAVE NEW treatment LABEL (control small, not control.small), TO BE UPDATED AS IT WAS DONE FOR THE X-ASIS
        \n-WRAP LEGEND TO BE 2 COLS https://stackoverflow.com/questions/39552682/base-r-horizontal-legend-with-multiple-rows" )

warning("ISSUE: SOME POST-HOCS OF IND_ONE_ANALYSIS LOOK WRONG (SEE DURATION_GROOMING_RECEIVED)")

### ind_net_properties ### 
## degree

plot_list <- list(ind_treated_net$barplot_delta_period_list$degree,
                  ind_untreated_net_nurse$barplot_delta_period_list$degree,
                  ind_untreated_net_forag$barplot_delta_period_list$degree)

# Set the same y-axis limits for all plots
YLIM_extra <- 1

plot_comps <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)

allplots <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps$y_limits) + ggtitle("treated nurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]] + ylim(plot_comps$y_limits) + ggtitle("untreated nurses")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 plot_list[[3]] + ylim(plot_comps$y_limits) + ggtitle("untreated foragers")+ fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 align="h")

ind_net_properties <- cowplot::plot_grid(
  cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],
                     ncol=3, rel_widths = c(0.28,0.24,0.24))
  , plot_comps$leg, ncol=1, rel_heights = c(0.9, 0.1))


ind_net_properties

  
### ind_beh_measures ### 3 panels

## prop_time_outside
plot_list <- list(ind_treated_beh$barplot_delta_period_list$prop_time_outside,
                  ind_untreated_beh_nurse$barplot_delta_period_list$prop_time_outside,
                  ind_untreated_beh_forag$barplot_delta_period_list$prop_time_outside)

YLIM_extra <- 0.01


plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)

allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + ggtitle("treated nurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]] + ylim(plot_comps1$y_limits) + ggtitle("untreated nurses")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 plot_list[[3]] + ylim(plot_comps1$y_limits) + ggtitle("untreated foragers")+ fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 align="h")

prop_time_outside <- cowplot::plot_grid(
  cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]],
                     ncol=3, rel_widths = c(0.28,0.24,0.24))
  , plot_comps1$leg, ncol=1, rel_heights = c(0.9, 0.1))



### ind_beh_measures### 2 panels


## duration_of_contact_with_treated_min

plot_list <- list(ind_untreated_beh_nurse$barplot_delta_period_list$duration_of_contact_with_treated_min,
                  ind_untreated_beh_forag$barplot_delta_period_list$duration_of_contact_with_treated_min)

YLIM_extra <- 0.05


plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)

allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + ggtitle("untreated nurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                  plot_list[[2]] + ylim(plot_comps1$y_limits) + ggtitle("untreated foragers")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")

duration_of_contact_with_treated <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                   ncol=2, rel_widths = c(0.28,0.24))
## inter_caste_contact_duration

plot_list <- list(ind_untreated_beh_nurse$barplot_delta_period_list$inter_caste_contact_duration,
                  ind_untreated_beh_forag$barplot_delta_period_list$inter_caste_contact_duration)

YLIM_extra <- 0.5


plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)

allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + ggtitle("untreated nurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                  plot_list[[2]] + ylim(plot_comps1$y_limits) + ggtitle("untreated foragers")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")


inter_caste_contact_duration <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                       ncol=2, rel_widths = c(0.28,0.24))


## duration_grooming_given_to_treated_min

plot_list <- list(ind_untreated_grooming_nurse$barplot_delta_period_list$duration_grooming_given_to_treated_min,
                  ind_untreated_grooming_forag$barplot_delta_period_list$duration_grooming_given_to_treated_min)

#YLIM_extra <- 0


plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)

allplots1 <- cowplot::align_plots(plot_list[[1]] + ggtitle("untreated nurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") , # ylim(plot_comps1$y_limits)
                                  plot_list[[2]] + ggtitle("untreated foragers")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs, # ylim(plot_comps1$y_limits)
                                  align="h")


duration_grooming_given_to_treated_min <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                   ncol=2, rel_widths = c(0.28,0.24))


### COMBINE ind_beh_measures 2 panels together

cowplot::plot_grid(
  duration_of_contact_with_treated,
  inter_caste_contact_duration,
  duration_grooming_given_to_treated_min,
  plot_comps1$leg, ncol=1, rel_heights = c(0.30,0.30,0.30, 0.1))


# ### ind_grooming_received ### 3 panels
# 
# plot_list <- list(ind_treated_grooming$barplot_delta_period_list$duration_grooming_received_min,
#                   ind_treated_grooming$barplot_delta_period_list$inter_caste_contact_duration,
#                   ind_treated_grooming$barplot_delta_period_list$prop_duration_grooming_received_outside_min)
# 
# YLIM_extra <- 1
# 
# plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
# 
# allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + ggtitle("treated nurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
#                                   plot_list[[2]] + ylim(plot_comps1$y_limits) + ggtitle("untreated nurses")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
#                                   plot_list[[3]] + ylim(plot_comps1$y_limits) + ggtitle("untreated foragers")+ fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
#                                   align="h")
# 
# treated_grooming <- cowplot::plot_grid(
#   cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]],
#                      ncol=3, rel_widths = c(0.28,0.24,0.24))
#   , plot_comps1$leg, ncol=1, rel_heights = c(0.9, 0.1))
# 
# 












### collective_net_properties ### 
collective_net_properties <- cowplot::plot_grid(
  coll_rescal_net$barplot_delta_period_list$modularity,
  coll_rescal_net$barplot_delta_period_list$clustering,
  coll_rescal_net$barplot_delta_period_list$task_assortativity,
  coll_rescal_net$barplot_delta_period_list$density,
  coll_rescal_net$barplot_delta_period_list$efficiency,
  coll_rescal_net$barplot_delta_period_list$degree_mean,
  labels=c("", "","",""), nrow = 2)


collective_net_properties



###################################################################################################################################
### Simulated Load plots ##########################################################################################################
###################################################################################################################################


### Experimentally measured and simulated M.brunneum transmission in pathogen-exposed colonies ###

#TEMP (The green line highlights the value at which there is an inversion in the sign of the density difference. This value was used as athreshold to distinguish high from low simulated loads in allsubsequent analyses)
high_threshold <- 0.0258
  
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


### Pathogen-induced changes in simulated diseasetransmission ###

######Open pdf Figure 3 #####
#pdf(file=paste(figurefolder,"/Figure3.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=page_height*0.45,pointsize=pointsize_more_than_2row2col)
##### Set-up layout and plot parameters #####
par(pars)
ncoli <- 4
heits <- c(5/9,0.05,4/9)
widz <- c(0.075,0.45,0,0.45)
layout(matrix(c(rep(1,ncoli/2),rep(4,ncoli/2),
                rep(5,ncoli),
                rep(5,ncoli/4),rep(2,ncoli/4),rep(5,ncoli/4),rep(3,ncoli/4)
), 3, ncoli, byrow = TRUE),heights=heits,widths = widz)
######First, plot distribution and threshold identification#######
plot_distribution(experiments="main_experiment",desired_treatments=c("pathogen"))
####Second, individual simulation results - comparison #######
root_path <- paste(disk_path,"/main_experiment",sep="")
queen <- T; treated <- F; nurses <- T; foragers <- T; 
unit_ori <- unit; unit <- 24
time_window <- 24

variable_list <- c("probability_high_level","probability_low_level")
names(variable_list) <- c("prob. receiving high load","prob. receiving low load")
transf_variable_list <- c("power3","sqrt")
predictor_list <- c("task_group","task_group")
analysis <- list(variable_list=variable_list,
                 transf_variable_list=transf_variable_list,
                 predictor_list=predictor_list,
                 violin_plot_param = list(c(1.5,0,-0.02,0.35,0.11),c(1.5,-0.02,-0.02,0.35,0.11)))
plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_simulation_results_observed",plot_untransformed=T,aligned=T)
unit <- unit_ori
# ####Third, add survival curve  #######
# par(pars)
# survival_analysis(experiment="survival_experiment",which_to_plot="second_only")
####Fourth, add letters  #########
par(xpd=NA)
##LETTERS
x_text1 <- grconvertX(0+1/80, from='ndc'); x_text2 <- grconvertX((sum(widz[1:2]))/(sum(widz))+1/80, from='ndc')
y_text1 <- grconvertY((1-1/80), from='ndc')
y_text2 <- grconvertY((1-heits[1]-1/80), from='ndc')
text(x_text1,y_text1,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text1,y_text2,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
text(x_text2,y_text1,labels=panel_casse("c"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
par(xpd=F)
####Fifth, Close figure 3 #######
#dev.off()