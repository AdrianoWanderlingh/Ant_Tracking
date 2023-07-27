#### Statistic and plots #####

#### code adapted from Stroeymeyt et al. 2018

### written by Nathalie Stroeymeyt and Adriano Wanderlingh

rm(list=ls())

USER <- "2A13_Office" # Nath_office 
DISK <-  "DISK4" #"Seagate Portable Drive"

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
disk_path    <- paste0("/media/",usr,"/",DISK,"/Lasius-Bristol_pathogen_experiment")
figurefolder <- paste0(disk_path,"/figures/") #"~/figures"
  
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

fixed_aspect_theme <- theme(aspect.ratio = 2) #move to plotting_parmas.R

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

coll_no_rescal_net <- collective_analysis_no_rescal(data_path,showPlot=F)

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
transf_variable_list <- c("log"       ,"sqrt"      ,"Box_Cox"              ,"log"       ,"log"        ,"log" ) # ,"log"   ,"power0.01"  ######"none", "sqrt" "log","power2"
# TRANSFORMATION NOTE: task_assortativity is hard to normalise (no transformation has the best result) - used Box_Cox

coll_rescal_net <- collective_analysis_rescal(data_path,showPlot=F)

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


ind_treated_net <- individual_ONE_analysis(data_path,which_individuals="treated",showPlot=F) # "treated","queen","nurse","forager"


### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside") #,"proportion_time_active", "average_bout_speed_pixpersec" ,"total_distance_travelled_pix", "inter_caste_contact_duration") #inter_caste_contact_duration?
names(variable_list) <- c("prop. time outside") # ,"prop. time active", "average bout speed pixpersec" ,"total distance travelled pix", "inter caste contact duration")
transf_variable_list <- c("power0.01"        )#,"none"                  ,"log"                          ,"log"                      ,"sqrt"                   )  ######"none", "sqrt" "log","power2"

ind_treated_beh <- individual_ONE_analysis(data_path,which_individuals="treated",showPlot=F) # "treated","queen","nurse","forager"

ind_treated_beh_lineplot <- line_plot(data_path,which_individuals="treated",showPlot=F)


#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_received_min", "inter_caste_contact_duration","prop_duration_grooming_received_outside_min") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2"
names(variable_list) <- c("duration grooming received (min)", "inter caste grooming duration","prop. duration grooming received outside (min)") # , "prop duration grooming received outside min","duration grooming received outside min"
transf_variable_list <- c("log"                          ,"log"                         , "Box_Cox")   ######"none", "sqrt" "log","power2"


ind_treated_grooming <- individual_ONE_analysis(data_path,which_individuals="treated",showPlot=F) # "treated","queen","nurse","forager"

ind_treated_grooming_lineplot <- line_plot(data_path,which_individuals="treated",showPlot=F)

################    FOR UNTREATED INDIVIDUALS    ################    


#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("degree")#, "mean_aggregated_distance_to_treated")
names(variable_list) <- c("degree")#,"mean aggregated distance to treated")
transf_variable_list <- c("none")#,"power0.01")  ######"none", "sqrt" "log","power2"

ind_untreated_net_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_net_forag <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"


### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside","inter_caste_contact_duration","duration_of_contact_with_treated_min")
names(variable_list) <- c("prop. time outside","inter caste contact duration","duration of contact with treated (min)")
transf_variable_list <- c("Box_Cox"        , "Box_Cox"                    , "Box_Cox"            )   ######"none", "sqrt" "log","power2"

ind_untreated_beh_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_beh_forag <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"


#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_given_to_treated_min")
names(variable_list) <- c("duration grooming given to treated (min)")
transf_variable_list <- c("log"        )   ######"none", "sqrt" "log","power2"

ind_untreated_grooming_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_grooming_forag <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"

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
# names(variable_list) <- c("prop. time outside","inter caste contact duration","duration of contact with treated min")
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
### COMPARE  grooming and time outside ############################################################################################
###################################################################################################################################

root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_received_min") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2"
names(variable_list) <- c("duration grooming received (min)") # , "prop duration grooming received outside min","duration grooming received outside min"
#transf_variable_list <- c("log"                           )   ######"none", "sqrt" "log","power2"

dur_groom_rec_data <- read_data(data_path,which_individuals="treated")


root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside") #,"proportion_time_active", "average_bout_speed_pixpersec" ,"total_distance_travelled_pix", "inter_caste_contact_duration") #inter_caste_contact_duration?
names(variable_list) <- c("prop. time outside") # ,"prop. time active", "average bout speed pixpersec" ,"total distance travelled pix", "inter caste contact duration")
#transf_variable_list <- c("power0.01"        )#,"none"                  ,"log"                          ,"log"                      ,"sqrt"                   )  ######"none", "sqrt" "log","power2"

prop_time_out_data <- read_data(data_path,which_individuals="treated")

#merge
common_col_names <- intersect(names(dur_groom_rec_data), names(prop_time_out_data))
#common_col_names <- common_col_names_NetSpace[!common_col_names_NetSpace %in% c("age")]
CompareBehavs <- dplyr::left_join(dur_groom_rec_data, prop_time_out_data, by = common_col_names[])

# remove extra cols
CompareBehavs <- CompareBehavs %>%
  dplyr::select(colony, tag, antID, time_hours, colony_size, treatment, age, period, time_of_day, 
                prop_time_outside, duration_grooming_received_min
  )

# melt the dataframe to long format
CompareBehavs <- reshape(CompareBehavs,
                         varying = c("prop_time_outside", "duration_grooming_received_min"),
                         v.names = "measure",
                         times = c("prop_time_outside", "duration_grooming_received_min"),
                         timevar = "variables",
                         direction = "long")
CompareBehavs$variables <- as.factor(CompareBehavs$variables)
CompareBehavs$variables <- gsub("_", " ", CompareBehavs$variables)

###### MODEL # only post exposure

#scaling for model (for plot, performed only after that the vars means are calculated)
CompareBehavs<- CompareBehavs %>%
  group_by(variables)  %>%
  dplyr::mutate(measure_scaled = scale(measure)) %>%
  ungroup()

mod1 <- lmer(log_transf(measure_scaled) ~ variables*time_hours + (1|colony) + (1|antID), data = CompareBehavs[which(CompareBehavs$period=="post"),])
#output_lmer(mod1)
anov  <- anova(mod1)
p_interaction_vars_time <- anov["variables:time_hours","Pr(>F)"]


###### PLOT

# 1. Calculate the mean of the variable
mean_data <- aggregate(measure ~ period + time_hours + variables + colony,
                       FUN = mean, na.rm = T, na.action = na.pass, CompareBehavs)

# 2. Calculate the grand mean and standard error dropping the colony AND TREATMENT factors
grand_mean_data <- mean_data %>%
  group_by(period, time_hours, variables) %>% #treatment
  summarise(grand_mean = mean(measure),
            standard_error = sd(measure) / sqrt(n()))

# Add NA values at time_hours == -3
unique_treatments <- unique(grand_mean_data$variables)
unique_periods <- unique(grand_mean_data$period)
na_rows <- expand.grid(period = unique_periods,
                       time_hours = -3,
                       variables = unique_treatments,
                       grand_mean = NA,
                       standard_error = NA)

grand_mean_data <- rbind(grand_mean_data, na_rows) %>%
  arrange( period, time_hours, variables) #treatment,


# Perform scaling by group # scale() function standardizes a vector to have mean 0 and standard deviation 1
grand_mean_data_scaled <- grand_mean_data %>%
  group_by(variables)  %>%
  dplyr::mutate(scaled_grand_mean = scale(grand_mean),scaled_standard_error = scale(standard_error)) %>%
  ungroup()


# Wrap the legend labels
grand_mean_data_scaled$variables <- str_wrap(grand_mean_data_scaled$variables, width = 20)  # Adjust the width as needed

#plot fit
GroomingVsTimeOutside <- ggplot(grand_mean_data_scaled, aes(x = time_hours, y = scaled_grand_mean,  fill = variables, color = variables, group = variables)) +
  geom_smooth(data = subset(grand_mean_data_scaled, period == "pre"), method = "lm", se = T, linetype = "dashed") +
  geom_smooth(data = subset(grand_mean_data_scaled, period == "post"), method = "lm", se = T, linetype = "solid") +
  scale_x_continuous(limits = c(min(grand_mean_data_scaled$time_hours), max(grand_mean_data_scaled$time_hours)), expand = c(0, 0)) +
  labs(#title = "Prop Time Outside by Time Hours and Treatment",
    x = "Time Hours since treatment exposure",
    #y = "scaled variables" #names(variable_list[i])
    y = "\nScaled\nvariables"
  ) +
  #geom_text(aes(x = 10, label = from_p_to_ptext(p_interaction_vars_time))) 
  annotate("text", x = 10, y = 3, label = from_p_to_ptext(p_interaction_vars_time)) +
  STYLE +
  theme(legend.position = c(.05, .95), # Position legend inside plot area
        legend.justification = c(0, 1), # Justify legend at top left
        legend.box.just = "left",
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent'),
        legend.title = element_blank()) + 
  
  guides(fill = guide_legend(nrow = 2,ncol = 1)) 


###################################################################################################################################
### PLOT GRIDS ####################################################################################################################
###################################################################################################################################

warning("ISSUE: SOME POST-HOCS OF IND_ONE_ANALYSIS LOOK WRONG (SEE DURATION_GROOMING_RECEIVED)")

###################################################################################################################################
### ind_net_properties ### 
## degree
plot_list <- list(ind_treated_net$barplot_delta_period_list$degree,
                  ind_untreated_net_nurse$barplot_delta_period_list$degree,
                  ind_untreated_net_forag$barplot_delta_period_list$degree)

# Set the same y-axis limits for all plots
YLIM_extra <- 1
plot_comps <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps$y_limits) + ggtitle("treated\nnurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]] + ylim(plot_comps$y_limits) + ggtitle("untreated\nnurses")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 plot_list[[3]] + ylim(plot_comps$y_limits) + ggtitle("untreated\nforagers")+ fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 align="h")
ind_net_degree <- cowplot::plot_grid(
  cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],
                     ncol=3, rel_widths = c(0.28,0.24,0.24))
  , plot_comps$leg, ncol=1, rel_heights = c(0.9, 0.1))


###################################################################################################################################  
### ind_beh_measures ### 3 panels

## prop_time_outside
plot_list <- list(ind_treated_beh$barplot_delta_period_list$prop_time_outside,
                  ind_untreated_beh_nurse$barplot_delta_period_list$prop_time_outside,
                  ind_untreated_beh_forag$barplot_delta_period_list$prop_time_outside)
YLIM_extra <- 0.01
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + ggtitle("treated\nnurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]] + ylim(plot_comps1$y_limits) + ggtitle("untreated\nnurses")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 plot_list[[3]] + ylim(plot_comps1$y_limits) + ggtitle("untreated\nforagers")+ fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 align="h")
prop_time_outside <- cowplot::plot_grid(
  cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]],
                     ncol=3, rel_widths = c(0.28,0.24,0.24))
  , plot_comps1$leg, ncol=1, rel_heights = c(0.9, 0.1))


### ind_beh_measures### 2 panels
# create text boxes for titles (to ensure equal size of all plotted objects)
titlePlots <- cowplot::align_plots(ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "untreated\nnurses", family = "Liberation Serif",  size = 4, hjust = 0.5), #fontface = "bold",
                                   ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "untreated\nforagers",family = "Liberation Serif", size = 4, hjust = 0.5), #fontface = "bold",
                                   align="h")
titlePlots_untreated <- cowplot::plot_grid(titlePlots[[1]], titlePlots[[2]], ncol=2, rel_widths = c(0.28,0.24))


## duration_of_contact_with_treated_min
plot_list <- list(ind_untreated_beh_nurse$barplot_delta_period_list$duration_of_contact_with_treated_min,
                  ind_untreated_beh_forag$barplot_delta_period_list$duration_of_contact_with_treated_min)
YLIM_extra <- 0.05
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits)    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits)   + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")
duration_of_contact_with_treated <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                   ncol=2, rel_widths = c(0.28,0.24))
## inter_caste_contact_duration
plot_list <- list(ind_untreated_beh_nurse$barplot_delta_period_list$inter_caste_contact_duration,
                  ind_untreated_beh_forag$barplot_delta_period_list$inter_caste_contact_duration)
YLIM_extra <- 0.5
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")
inter_caste_contact_duration <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                       ncol=2, rel_widths = c(0.28,0.24))

## duration_grooming_given_to_treated_min
plot_list <- list(ind_untreated_grooming_nurse$barplot_delta_period_list$duration_grooming_given_to_treated_min,
                  ind_untreated_grooming_forag$barplot_delta_period_list$duration_grooming_given_to_treated_min)
#YLIM_extra <- 0
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots1 <- cowplot::align_plots(plot_list[[1]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)), # ylim(plot_comps1$y_limits)
                                  plot_list[[2]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs, # ylim(plot_comps1$y_limits)
                                  align="h")
duration_grooming_given_to_treated_min <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                   ncol=2, rel_widths = c(0.28,0.24))



### COMBINE ind_beh_measures 2 panels together
ind_beh_measures <- cowplot::plot_grid(
  titlePlots_untreated,
  duration_of_contact_with_treated,
  inter_caste_contact_duration,
  duration_grooming_given_to_treated_min,
  plot_comps1$leg, ncol=1, rel_heights = c(0.08,0.30,0.30,0.30, 0.05))


###################################################################################################################################
### ind_grooming_received ### 3 panels
plot_list <- list(ind_treated_grooming$barplot_delta_period_list$duration_grooming_received_min,
                  ind_treated_grooming$barplot_delta_period_list$inter_caste_contact_duration,
                  ind_treated_grooming$barplot_delta_period_list$prop_duration_grooming_received_outside_min)
YLIM_extra <- 0
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots1 <- cowplot::align_plots(plot_list[[1]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[2]]$labels$y)),
                                  plot_list[[3]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[3]]$labels$y)),
                                  align="h")
treated_grooming <- cowplot::plot_grid(
  cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]],
                     ncol=3, rel_widths = c(0.24,0.25,0.265))
  , plot_comps1$leg, ncol=1, rel_heights = c(0.9, 0.1))

###################################################################################################################################
### comparing timeline of grooming and time_outside line_plots  ### 3 panels
plot_list <- list(ind_treated_beh_lineplot$prop_time_outside,
                  ind_treated_grooming_lineplot$duration_grooming_received_min,
                  ind_treated_grooming_lineplot$prop_duration_grooming_received_outside_min,
                  GroomingVsTimeOutside)
YLIM_extra <- 0
#plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
#leg <-cowplot::get_legend(plot_list[[1]] + theme(legend.direction="horizontal"))

#inherit legend from elsewhere
allplots1 <- cowplot::align_plots(plot_list[[1]] + theme(aspect.ratio = 0.5)  + remove_x_labs + guides(fill = "none") + theme(legend.position = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + theme(aspect.ratio = 0.5)  + remove_x_labs + guides(fill = "none") + theme(legend.position = "none") + labs(y = split_title(plot_list[[2]]$labels$y)),
                                  plot_list[[3]] + theme(aspect.ratio = 0.5)                  + guides(fill = "none") + theme(legend.position = "none") + labs(y = split_title(plot_list[[3]]$labels$y)),
                                  plot_list[[4]] + theme(aspect.ratio = 0.55)                  + guides(fill = "none")                                   ,
                                  align="v")
GroomVSTimeOut <- cowplot::plot_grid(
  cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]],plot_list[[4]],
                     ncol=2, rel_widths = c(0.24,0.24,0.24,0.20))
  , plot_comps1$leg, ncol=1, rel_heights = c(0.9, 0.15))


###################################################################################################################################
### collective_net_properties ### 
plot_list <- list(coll_rescal_net$barplot_delta_period_list$modularity,
                  coll_rescal_net$barplot_delta_period_list$clustering,
                  coll_rescal_net$barplot_delta_period_list$task_assortativity,
                  coll_rescal_net$barplot_delta_period_list$density,
                  coll_rescal_net$barplot_delta_period_list$efficiency,
                  coll_rescal_net$barplot_delta_period_list$degree_mean)
# Set the same y-axis limits for all plots
YLIM_extra <- 0.01
#have 2 scales: 1 for top row (measures expected to increase), 1 for bottom row (measures expected to decrease), 
plot_compsA <- multi_plot_comps(plot_list[1:3],ylim_extra=YLIM_extra)
plot_compsB <- multi_plot_comps(plot_list[4:6],ylim_extra=YLIM_extra)
allplots <- cowplot::align_plots(plot_list[[1]] + ylim(plot_compsA$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]] + ylim(plot_compsA$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[3]] + ylim(plot_compsA$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[4]] + ylim(plot_compsB$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[5]] + ylim(plot_compsB$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[6]] + ylim(plot_compsB$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 align="h")
collective_net_properties <- cowplot::plot_grid(
  cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],allplots[[4]], allplots[[5]], allplots[[6]],
                     ncol=3, rel_widths = c(0.1,0.1,0.1,0.1,0.1,0.1))
  , plot_compsB$leg, ncol=1, rel_heights = c(0.9, 0.1))
# width_pixels <- 600 
# height_pixels <- 600



###################################################################################################################################
# GRID PLOTS!! (add here the outputs of the models in them! stat_outcomes$formatted)

#size for a good font size, for a 3 horizontal panels plot
# width_pixels <- 450 (150 x panel)
# height_pixels <- 300 (2/3 of width)

SavePrint_plot(
  plot_obj = ind_net_degree,
  plot_name = "ind_net_degree",
  plot_size = c(430/ppi, 300/ppi),
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)

SavePrint_plot(
  plot_obj = prop_time_outside,
  plot_name = "prop_time_outside",
  plot_size = c(430/ppi, 300/ppi),
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)

# ISSUE WITH PVAL IN MODEL!!!!!!!!!!!!!!!!!!!!!!!!!
SavePrint_plot(
  plot_obj = ind_beh_measures, 
  plot_name = "ind_beh_measures",
  plot_size = c(330/ppi, 900/ppi), #extra length required to get y-axis labeled and unlabeled of same size
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)

# WHY NOW DIFF IS SIGNIFICANT? WHICH TRANSF IS MORE ADAPT? IT SHOULD NOT BE! LETTERS ARE WEIRD TOO
# remove inter_Caste_grooming_duration
SavePrint_plot(
  plot_obj = treated_grooming, 
  plot_name = "treated_grooming",
  plot_size = c(430/ppi, 280/ppi), #extra length required to fix the different digits on y-axis
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)


SavePrint_plot(
  plot_obj = GroomVSTimeOut, 
  plot_name = "GroomVSTimeOut",
  plot_size = c(520/ppi, 300/ppi), #extra length required to fix the different digits on y-axis
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)
 

SavePrint_plot(
  plot_obj = collective_net_properties, 
  plot_name = "collective_net_properties",
  plot_size = c(450/ppi,480/ppi), #extra length required to get y-axis labeled and unlabeled of same size
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)









warning("MAKE SURE THAT FOR THE INDIVIDUAL MEASURES OF UNTREATED WORKERS, THE DATASET RETAINED IS RIGHT")
warning("inter_caste_contact_duration REMOVE FROM EVERYWHERE, or add nurses if can be added to the story!")
ind_treated_grooming$barplot_delta_period_list$inter_caste_contact_duration #ADD PLOT



## posthoc letters looking weird:
# - prop. time outside (treated nurses)
# - duration grooming received (treated nurses) # why now significant?
# - duration grooming given to treated (unt. nurses)
# - simulated_load (nurses - is there an issue with the matrix?)
# - transmission_latency (foragers)

# stars not assigned (occasionally)
# - inter caste contact duration (forager)
# - duration grooming given to treated (forager) - ensure it is directed to treated only!!!



###################################################################################################################################
### Simulated Load plots ##########################################################################################################
###################################################################################################################################

root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="")
pattern="collective_simulation_results_observed"
variable_list <-  c("Prevalence", "Mean_load", "Load_skewness", "Queen_load", "logistic_r")
names(variable_list) <-  c("Prevalence", "Mean load", "Load skewness", "Queen load", "logistic r")
transf_variable_list <- c("none"       ,"none"        ,"none"           ,"log"      ,"log" )   ######"none", "sqrt" "log","power2"

coll_no_rescal_sim <- collective_analysis_no_rescal(data_path,showPlot=F)
warning(paste("-Prevalence", "-Mean load", "-Load skewness","Can't be modeled",sep= "\n"))

#https://stackoverflow.com/questions/68915173/how-do-i-fit-a-quasi-poisson-model-with-lme4-or-glmmtmb
#model <- glmer(Prevalence ~ period*treatment + (1|colony) ,data=data, family = poisson)



root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="")
pattern="individual_simulation_results_observed"
variable_list <-  c("simulated_load", "transmission_latency", "transmission_rank")
names(variable_list) <-  c("simulated load", "transmission latency", "transmission rank")
transf_variable_list <- c("none"           ,"log"                ,"none")   ######"none", "sqrt" "log","power2"

ind_untreated_sim_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_sim_forag <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"

#probability_of_transmission for nurses and foragers is almost always 1, can't really be modeled








###################################################################################################################################
### PLOT GRIDS ####################################################################################################################
###################################################################################################################################

###################################################################################################################################
### collective_sim_properties ### 
plot_list <- list(coll_no_rescal_sim$barplot_delta_period_list$Prevalence,
                  coll_no_rescal_sim$barplot_delta_period_list$Mean_load,
                  coll_no_rescal_sim$barplot_delta_period_list$Load_skewness,
                  coll_no_rescal_sim$barplot_delta_period_list$Queen_load,
                  coll_no_rescal_sim$barplot_delta_period_list$logistic_r)
# Set the same y-axis limits for all plots
YLIM_extra <- 0.001
# #have 2 scales: 1 for top row (measures expected to increase), 1 for bottom row (measures expected to decrease), 
# plot_compsA <- multi_plot_comps(plot_list[1:2],ylim_extra=YLIM_extra)
# plot_compsB <- multi_plot_comps(plot_list[3:4],ylim_extra=YLIM_extra)
allplots <- cowplot::align_plots(plot_list[[1]]  + ylim(-0.009,0.009) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]]  + ylim(-0.007,0.007) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[3]]  + ylim(-0.4,0.4) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[4]]  + ylim(-0.07,0.07)           + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[5]]                               + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 align="h")
collective_sim_properties <- cowplot::plot_grid(
  cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],allplots[[4]], allplots[[5]],
                     ncol=5, rel_widths = c(0.1,0.1,0.1,0.1,0.1))
  #, plot_compsB$leg
  ,ncol=1, rel_heights = c(0.9, 0.1))

# width_pixels <- 900 
# height_pixels <- 250

###################################################################################################################################
### individual_sim_properties ###


## simulated_load
plot_list <- list(ind_untreated_sim_nurse$barplot_delta_period_list$simulated_load,
                  ind_untreated_sim_forag$barplot_delta_period_list$simulated_load)
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=0.005)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) +  fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")
individual_sim_SL <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                       ncol=2, rel_widths = c(0.28,0.24))
## transmission_latency
plot_list <- list(ind_untreated_sim_nurse$barplot_delta_period_list$transmission_latency,
                  ind_untreated_sim_forag$barplot_delta_period_list$transmission_latency)
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=0.5)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")
individual_sim_TL <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                   ncol=2, rel_widths = c(0.28,0.24))

## transmission_rank
plot_list <- list(ind_untreated_sim_nurse$barplot_delta_period_list$transmission_rank,
                  ind_untreated_sim_forag$barplot_delta_period_list$transmission_rank)
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=0.8)
allplots1 <- cowplot::align_plots(plot_list[[1]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)), # ylim(plot_comps1$y_limits)
                                  plot_list[[2]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs, # ylim(plot_comps1$y_limits)
                                  align="h")
individual_sim_TR <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                             ncol=2, rel_widths = c(0.28,0.24))

### COMBINE ind_beh_measures 2 panels together
individual_sim_properties <- cowplot::plot_grid(
  titlePlots_untreated,
  individual_sim_SL,
  individual_sim_TL,
  individual_sim_TR,
  plot_comps1$leg, ncol=1, rel_heights = c(0.05,0.30,0.30,0.30, 0.05))


#still lightly offset
SavePrint_plot(
  plot_obj = individual_sim_properties, 
  plot_name = "individual_sim_properties",
  plot_size = c(350/ppi, 900/ppi), #extra length required to get y-axis labeled and unlabeled of same size
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)

SavePrint_plot(
  plot_obj = collective_sim_properties, 
  plot_name = "collective_sim_properties",
  plot_size = c(520/ppi, 120/ppi), #extra length required to get y-axis labeled and unlabeled of same size
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)




warning(paste("-Prevalence", "-Mean load", "-Load skewness","Can't be modeled",sep= "\n"))

###################################################################################################################################
### Experimentally measured and simulated M.brunneum transmission in pathogen-exposed colonies ###

#TEMP (The green line highlights the value at which there is an inversion in the sign of the density difference. This value was used as athreshold to distinguish high from low simulated loads in allsubsequent analyses)
# high_threshold <- 0.0258
 
# logging not enough for qPCR data, as it has kurtois = 27.28894
# sqrt                                                  27.46589
# power0.01                                             70.34206
# power0.001                                            26.71343
# power0.0001                                           26.64026
# Box_Cox fails
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
translated_high_threshold <- plot_qpcr(experiments=c("main_experiment")) #disabled the predict=high_threshold in plot_regression
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
#clean();
Sys.sleep(2)


### Pathogen-induced changes in simulated disease transmission ###

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
plot_distribution(experiments="main_experiment",desired_treatments=c("pathogen"),seeds="treated_workers")
plot_distribution(experiments="main_experiment",desired_treatments=c("pathogen"),seeds="random_workers") #desired_treatments may be changed to all
plot_distribution(experiments="main_experiment",desired_treatments=c("pathogen"),seeds="nurses") 
plot_distribution(experiments="main_experiment",desired_treatments=c("pathogen"),seeds="foragers")
#"foragers", "nurses", "random_workers", "treated_workers"

#dev.off()







# THIS PART IS NOT USED AS THE HIGH-LOW LOAD PROBS ARE NOT CALCULATEDA
# # ####Second, individual simulation results - comparison #######
# # root_path <- paste(disk_path,"/main_experiment",sep="")
# # queen <- T; treated <- F; nurses <- T; foragers <- T;
# # unit_ori <- unit; unit <- 24
# # time_window <- 24
# # 
# # variable_list <- c("probability_high_level","probability_low_level")
# # names(variable_list) <- c("prob. receiving high load","prob. receiving low load")
# # transf_variable_list <- c("power3","sqrt")
# # predictor_list <- c("task_group","task_group")
# # analysis <- list(variable_list=variable_list,
# #                  transf_variable_list=transf_variable_list,
# #                  predictor_list=predictor_list,
# #                  violin_plot_param = list(c(1.5,0,-0.02,0.35,0.11),c(1.5,-0.02,-0.02,0.35,0.11)))
# # plot_before_vs_after(root_path=root_path, data_path=paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep=""),analysis=analysis,status="all",collective=F,pool_plot=F,pattern="individual_simulation_results_observed",plot_untransformed=T,aligned=T)
# # unit <- unit_ori
# # # ####Third, add survival curve  #######
# # # par(pars)
# # # survival_analysis(experiment="survival_experiment",which_to_plot="second_only")
# # ####Fourth, add letters  #########
# # par(xpd=NA)
# # ##LETTERS
# # x_text1 <- grconvertX(0+1/80, from='ndc'); x_text2 <- grconvertX((sum(widz[1:2]))/(sum(widz))+1/80, from='ndc')
# # y_text1 <- grconvertY((1-1/80), from='ndc')
# # y_text2 <- grconvertY((1-heits[1]-1/80), from='ndc')
# # text(x_text1,y_text1,labels=panel_casse("a"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
# # text(x_text1,y_text2,labels=panel_casse("b"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
# # text(x_text2,y_text1,labels=panel_casse("c"),font=panel_font, cex=panel_cex,adj=c(0,1),xpd=NA)
# # par(xpd=F)
# # ####Fifth, Close figure 3 #######
# # #dev.off()


