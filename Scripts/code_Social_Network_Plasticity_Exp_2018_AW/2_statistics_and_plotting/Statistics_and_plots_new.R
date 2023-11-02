#### Statistic and plots #####

#### code adapted from Stroeymeyt et al. 2018

### written by Nathalie Stroeymeyt and Adriano Wanderlingh

rm(list=ls())

USER <- "2A13_Office" # Nath_office 
DISK <-  "Seagate Portable Drive" #"DISK4"


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
fixed_aspect_theme_PRE <- theme(aspect.ratio = 4) #move to plotting_parmas.R


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
  str(data)
}


###################################################################
###I - Example collective analysis - without rescaling ############
###################################################################

if(RUN_UNSCALED_NETS){

#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/network_properties_edge_weights_duration/pre_vs_post_treatment/all_workers",sep="")
pattern="colony_data.txt"
variable_list <- c("modularity_FacetNet","clustering","task_assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
names(variable_list) <- c("modularity_FacetNet","clustering","task assortativity","efficiency","mean degree","degree maximum","density","diameter")
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
data_path <- paste(root_path,"/network_properties_edge_weights_duration/pre_vs_post_treatment/all_workers",sep="")
pattern="colony_data.txt"
variable_list <- c("modularity_FacetNet","clustering","task_assortativity","efficiency","degree_mean","density") #,"degree_maximum","diameter"
names(variable_list) <- c("modularity (FacetNet)","clustering","task assortativity","efficiency","mean degree","density") # ,"degree maximum","diameter"
transf_variable_list <- c("sqrt"       ,"sqrt"      ,"Box_Cox"              ,"log"       ,"log"        ,"log" ) # ,"log"   ,"power0.01"  ######"none", "sqrt" "log","power2"
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
data_path=paste(root_path,"/processed_data/network_properties_edge_weights_duration/pre_vs_post_treatment/all_workers",sep="")
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
transf_variable_list <- c("power0.1"        )#,"none"                  ,"log"                          ,"log"                      ,"sqrt"                   )  ######"none", "sqrt" "log","power2"

ind_treated_beh <- individual_ONE_analysis(data_path,which_individuals="treated",showPlot=F) # "treated","queen","nurse","forager"

#timeline
ind_treated_beh_lineplot <- line_plot(data_path,which_individuals="treated",showPlot=F)


#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_received_min","number_contacts_received") #,"prop_duration_grooming_received_outside_min") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2", , "inter_caste_contact_duration"
names(variable_list) <- c("duration grooming received (min)","number grooming contacts received") #,"prop. duration grooming received outside (min)") # , "prop duration grooming received outside min","duration grooming received outside min"
transf_variable_list <- c("log"                             ,"log"                     ) #, "Box_Cox")   ######"none", "sqrt" "log","power2"


ind_treated_grooming <- individual_ONE_analysis(data_path,which_individuals="treated",showPlot=F) # "treated","queen","nurse","forager"
#timeline
ind_treated_grooming_lineplot <- line_plot(data_path,which_individuals="treated",showPlot=F)

################    FOR UNTREATED INDIVIDUALS    ################    


#### ALL INTERACTIONS ####

### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties_edge_weights_duration/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("degree")#, "mean_aggregated_distance_to_treated")
names(variable_list) <- c("degree")#,"mean aggregated distance to treated")
transf_variable_list <- c("none")#,"power0.01")  ######"none", "sqrt" "log","power2"

ind_untreated_net_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_net_forag <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"

# PRE period only (CONSTITUTIVE ORGANISATION)
ind_untreated_net_nurse_PRE <- individual_ONE_analysis(data_path,which_individuals="nurse",pre_only=T,showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_net_forag_PRE <- individual_ONE_analysis(data_path,which_individuals="forager",pre_only=T,showPlot=F) ## "treated","queen","nurse","forager"


### behavioural data
root_path <- paste(disk_path,"/main_experiment",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside","inter_caste_contact_duration","duration_of_contact_with_treated_min")
names(variable_list) <- c("prop. time outside","inter caste contact duration","duration of contact with treated (min)")
transf_variable_list <- c("power0.1"        , "Box_Cox"                    , "log"            )   ######"none", "sqrt" "log","power2"

ind_untreated_beh_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_beh_forag <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"
print("note: prop. time outside for foragers has a high kurtosis (8) but this transformation works well for treated and nurses and the posthoc matches the plot so I'll keep it")

# PRE period only (CONSTITUTIVE ORGANISATION)
ind_untreated_beh_nurse_PRE <- individual_ONE_analysis(data_path,which_individuals="nurse",pre_only=T,showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_beh_forag_PRE <- individual_ONE_analysis(data_path,which_individuals="forager",pre_only=T,showPlot=F) ## "treated","queen","nurse","forager"


#### GROOMING INTERACTIONS ####
root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_given_to_treated_min")
names(variable_list) <- c("duration grooming given to treated (min)")
transf_variable_list <- c("power0.1"        )   ######"none", "sqrt" "log","power2"

ind_untreated_grooming_nurse <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_grooming_forag <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"

# PRE period only (CONSTITUTIVE ORGANISATION)
ind_untreated_grooming_nurse_PRE <- individual_ONE_analysis(data_path,which_individuals="nurse",pre_only=T,showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_grooming_forag_PRE <- individual_ONE_analysis(data_path,which_individuals="forager",pre_only=T,showPlot=F) ## "treated","queen","nurse","forager"

#timeline
ind_untreated_grooming_lineplot_nurse <- line_plot(data_path,which_individuals="nurse",showPlot=T)
ind_untreated_grooming_lineplot_forag <- line_plot(data_path,which_individuals="forager",showPlot=T)
warning("it shows if treated are groomed by foragers or nurses, DON'T ADD IT TO PLOT GRID BUT ADD THE INFO UNDER Treated nurses grooming VS isolation ")

### Entropy measure
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"

Entropy_size <- calculate_entropy(data_path,which_individuals=all_workers,number_permutations=5,showPlot=F, pre_only = T) # 500 PERMS  # "treated","queen","nurse","forager"

Entropy_size$Dip_plot
Entropy_size$entropy_plot


################ check bimodality
# copied from 19_Facetnet
task_groups_A    <- read.table(paste0(disk_path,"/main_experiment/original_data/task_groups.txt"),header=T,stringsAsFactors = F)

size_order      <- c("small","big")
task_groups_A$size     <- unlist(lapply( task_groups_A$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))
task_groups_A$size      <- factor(task_groups_A$size    , levels=size_order   [which(size_order%in%task_groups_A$size )])

#task_groups_A[!is.na(task_groups_A$Forager_score),"Forager_score_log"] <- log_transf(task_groups_A[!is.na(task_groups_A$Forager_score),"Forager_score"] )

# Filter out NA values
task_groups_A <- task_groups_A %>%
  filter(!is.na(Forager_score_contact_count))

# A function to get density data frame along with confidence intervals
get_density_df <- function(data) {
  d <- density(data$Forager_score_contact_count, adjust = 1/1.2) # adjust is 1 for normal smoothing
  data.frame(x = d$x, y = d$y, ci = qnorm(0.975) * d$bw * sqrt(d$n) / sqrt(d$n - 1))
}

# Calculate density and confidence interval
density_data <- task_groups_A %>%
  group_by(size) %>%
  group_modify(~get_density_df(.)) %>%
  ungroup()


# Plot
SocialMaturity <- ggplot(density_data, aes(x = x, y = y, fill = size)) +
  geom_ribbon(aes(ymin = y - ci, ymax = y + ci), alpha = 0.3) +
  geom_line(aes(color=size),size=1, alpha=1) +
  #geom_line(aes(color=size), stat="density", size=2, alpha=1, adjust=1/1.2) +
  labs(
    x = "Social maturity",
    y = "Density"
  ) +
  scale_colour_manual(
    values = c("small" = pooled_small_cols_colour, "big" = pooled_large_cols_colour),
    labels = function(x) str_to_title(gsub("big", "large", gsub("\\.", " ", x)))
  ) +
  scale_fill_manual(
    values = c("small" = pooled_small_cols_colour, "big" = pooled_large_cols_colour),
    labels = function(x) str_to_title(gsub("big", "large", gsub("\\.", " ", x)))
  ) + # you can choose your own colors
  scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(ylim = c(0, max(density_data$y + density_data$ci)), xlim = c(0, 1)) +
  theme_classic() +
  theme(
    #legend.position = "top",   # Place the legend on top
    #legend.direction = "horizontal",  # Display legend items in a single line
    #legend.box = "horizontal",  # Arrange legend items horizontally
    text = element_text(family = "Liberation Serif") 

  )

# # Generate plots of Social Maturity
# ggplot(task_groups_A, aes(x = Forager_score_contact_count, colour = size)) + 
#   #geom_line(aes(color=size, group=colony), stat="density", size=1, alpha=0.2, adjust=1/1.2) +
#   geom_line(aes(color=size), stat="density", size=2, alpha=1, adjust=1/1.2) +
#   theme_minimal() +
#   xlab("Social maturity")





#bimodality comparison #########################################
# warning("RUN ON THE SOCIAL MAT. SCORES -FORAG SCORES-.\nMODIFY TO SPECIFY THE VARIABLE OUTSIDE, AS IN THE OTHER FUNCTIONS (NOW PICKS THE OBJ CALLED $variable")
# p_value <- bimodality_permutation_test(data,n_permutations = 1000)
# print(paste("Permutation test p-value:", p_value))

task_groups <- read.table(paste(root_path,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
task_groups$size     <- unlist(lapply( task_groups$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))

hist(task_groups[which(task_groups$size=="big"),"Forager_score"])
hist(task_groups[which(task_groups$size=="small"),"Forager_score"])

dip_results <- by(task_groups$Forager_score, task_groups$size, perform_dip_test)

bandwidth_permutation_test(task_groups, variable = "Forager_score")
warning("there is no proper way to statistically test which distribution is more bimodal")




###################################################################################################################################
###IV - individual analysis - for TWO types of individuals (e.g. nurses and foragers; or "untreated" and "treated") ###############
###################################################################################################################################

# #### ALL INTERACTIONS ####
# 
# ### network properties
# root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
# data_path=paste(root_path,"/processed_data/network_properties_edge_weights_duration/pre_vs_post_treatment/all_workers",sep="")
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


### check decline with time per HOUR ()see time 0, 3, etc
  # Filter data for specific times and variable
  filtered_data <- CompareBehavs %>%
    filter(time_hours %in% c(0, 3, 6) & variables == "duration grooming received min") %>%
    group_by(time_hours) %>%
    summarise(mean_measure = round(mean(measure, na.rm = TRUE)/3, 3), #divide by 3 to get the per hour
              sd_measure = round(sd(measure, na.rm = TRUE)/3, 3))
  
  # Extract values and format them into a string
  values <- paste0(filtered_data$time_hours, "h: ", filtered_data$mean_measure, " ± ", filtered_data$sd_measure)
  formatted_values <- paste(values, collapse = ", ")
  
  # Insert the formatted values into the sentence
  sentence <- paste0("It is of note that the duration of the grooming received by treated ants sharply decreased after 3h from the treatment (", formatted_values, ").")
  
  print(sentence)
  


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
mean_data <- aggregate(measure ~ period + time_hours + variables + colony + treatment,
                       FUN = mean, na.rm = T, na.action = na.pass, CompareBehavs)

# 2. Calculate the grand mean and standard error dropping the colony AND TREATMENT factors
grand_mean_data <- mean_data %>%
  group_by(period, time_hours, variables, treatment) %>% #treatment
  summarise(grand_mean = mean(measure),
            standard_error = sd(measure) / sqrt(n()))

# Add NA values at time_hours == -3
unique_treatments <- unique(grand_mean_data$treatment)
unique_vars       <- unique(grand_mean_data$variables)
unique_periods    <- unique(grand_mean_data$period)
na_rows <- expand.grid(period = unique_periods,
                       time_hours = -3,
                       variables = unique_vars,
                       treatment = unique_treatments,
                       grand_mean = NA,
                       standard_error = NA)

grand_mean_data <- rbind(grand_mean_data, na_rows) %>%
  arrange( period, time_hours, variables, treatment) #treatment,


# Perform scaling by group # scale() function standardizes a vector to have mean 0 and standard deviation 1
grand_mean_data_scaled <- grand_mean_data %>%
  group_by(variables)  %>%
  dplyr::mutate(scaled_grand_mean = scale(grand_mean),scaled_standard_error = scale(standard_error)) %>%
  ungroup()

# overall means
overall_grand_mean <- grand_mean_data_scaled %>%
  group_by(time_hours, period, variables) %>%
  summarise(scaled_grand_mean = mean(scaled_grand_mean))

#plot fit
GroomingVsTimeOutside <- ggplot(grand_mean_data_scaled, aes(x = time_hours, y = scaled_grand_mean, fill = variables, color= variables, group = variables)) +
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
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) + #inverting the order of scale fill and scale color breaks the plot
  theme(legend.position = c(.05, .95), # Position legend inside plot area
        legend.justification = c(0, 1), # Justify legend at top left
        legend.box.just = "left",
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent'),
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2,ncol = 1))

# Wrap the legend labels
# grand_mean_data_scaled$variables <- str_wrap(grand_mean_data_scaled$variables, width = 20)  # Adjust the width as needed
# #add colours for vars
# colors <- brewer.pal(length(unique(grand_mean_data_scaled$variables)), "Set2")

ggplot(grand_mean_data_scaled, aes(x = time_hours, y = scaled_grand_mean, color = treatment, group = treatment)) +
  # Add lines for each treatment
  geom_smooth(data = subset(grand_mean_data_scaled, period == "pre"), method = "lm", linetype = "dashed", aes(linetype = "pre")) +
  geom_smooth(data = subset(grand_mean_data_scaled, period == "post"),method = "lm", linetype = "solid", aes(linetype = "post")) +
  geom_ribbon(data = subset(grand_mean_data_scaled, period == "post"), stat = "smooth", method = "lm", fill = NA) +
  # Add overall mean line
  #geom_smooth(data = subset(overall_grand_mean, period == "pre"),  method = "lm", linetype = "dashed", aes(linetype = "pre", group = variables), color = "black") +
  #geom_smooth(data = subset(overall_grand_mean, period == "post"), method = "lm", linetype = "solid", aes(linetype = "post", group = variables),color = "black") + 
  STYLE +
  colFill_treatment +
  colScale_treatment +
  theme(legend.position = c(.05, .95), # Position legend inside plot area
        legend.justification = c(0, 1), # Justify legend at top left
        legend.box.just = "left",
        legend.direction = "horizontal",
        legend.background = element_rect(fill='transparent'),
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2,ncol = 1)) +
  facet_grid(.~variables) 


###################################################################################################################################
### ASSESS degree and time outside ############################################################################################
###################################################################################################################################


### network properties
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/network_properties_edge_weights_duration/pre_vs_post_treatment/all_workers",sep="")
pattern="individual_data"
variable_list <-        c("degree")#,"aggregated_distance_to_queen") 
names(variable_list) <- c("degree")#,"aggregated distance to queen")
#transf_variable_list <- c("none"  )#,"log")  ######"none", "sqrt" "log","power2"

degree_data <- read_data(data_path,which_individuals="treated")


root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("prop_time_outside") #,"proportion_time_active", "average_bout_speed_pixpersec" ,"total_distance_travelled_pix", "inter_caste_contact_duration") #inter_caste_contact_duration?
names(variable_list) <- c("prop. time outside") # ,"prop. time active", "average bout speed pixpersec" ,"total distance travelled pix", "inter caste contact duration")
#transf_variable_list <- c("power0.01"        )#,"none"                  ,"log"                          ,"log"                      ,"sqrt"                   )  ######"none", "sqrt" "log","power2"

prop_time_out_data <- read_data(data_path,which_individuals="treated")

root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_received_min") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2"
names(variable_list) <- c("duration grooming received (min)") # , "prop duration grooming received outside min","duration grooming received outside min"
#transf_variable_list <- c("log"                           )   ######"none", "sqrt" "log","power2"

dur_groom_rec_data <- read_data(data_path,which_individuals="treated")


#merge
common_col_names <- intersect(names(degree_data), names(prop_time_out_data))
common_col_names <-  common_col_names[common_col_names %in% c("colony", "tag","antID", "treatment", "time_hours","period")]
#common_col_names <- common_col_names_NetSpace[!common_col_names_NetSpace %in% c("age")]
CompareBehavs <- dplyr::left_join(prop_time_out_data,degree_data, by = common_col_names[])

#merge
common_col_names <- intersect(names(CompareBehavs), names(dur_groom_rec_data))
common_col_names <-  common_col_names[common_col_names %in% c("colony", "tag","antID", "treatment", "time_hours","period")]
#common_col_names <- common_col_names_NetSpace[!common_col_names_NetSpace %in% c("age")]
CompareBehavs <- dplyr::left_join(CompareBehavs,dur_groom_rec_data, by = common_col_names[])

# remove extra cols
CompareBehavs <- CompareBehavs %>%
  dplyr::select(colony, tag, antID, time_hours, treatment, period, 
                prop_time_outside.x, degree , duration_grooming_received_min
  )

CompareBehavs$prop_time_outside <- CompareBehavs$prop_time_outside.x; CompareBehavs$prop_time_outside.x <- NULL


# melt the dataframe to long format
CompareBehavs <- reshape(CompareBehavs,
                         varying = c("prop_time_outside", "degree","duration_grooming_received_min"),
                         v.names = "measure",
                         times = c("prop_time_outside", "degree","duration_grooming_received_min"),
                         timevar = "variables",
                         direction = "long")
CompareBehavs$variables <- as.factor(CompareBehavs$variables)
CompareBehavs$variables <- gsub("_", " ", CompareBehavs$variables)


###### MODEL # only post exposure

# #scaling for model (for plot, performed only after that the vars means are calculated)
# CompareBehavs<- CompareBehavs %>%
#   group_by(variables)  %>%
#   dplyr::mutate(measure_scaled = scale(measure)) %>%
#   ungroup()
# 
# mod1 <- lmer(log_transf(measure_scaled) ~ variables*time_hours + (1|colony) + (1|antID), data = CompareBehavs[which(CompareBehavs$period=="post"),])
# #output_lmer(mod1)
# anov  <- anova(mod1)
# p_interaction_vars_time <- anov["variables:time_hours","Pr(>F)"]
# 

###### PLOT

# 1. Calculate the mean of the variable
mean_data <- aggregate(measure ~ period + time_hours + variables + colony + treatment,
                       FUN = mean, na.rm = T, na.action = na.pass, CompareBehavs)

# 2. Calculate the grand mean and standard error dropping the colony AND TREATMENT factors
grand_mean_data <- mean_data %>%
  group_by(period, time_hours, variables, treatment) %>% #treatment
  summarise(grand_mean = mean(measure),
            standard_error = sd(measure) / sqrt(n()))

# Add NA values at time_hours == -3
unique_treatments <- unique(grand_mean_data$treatment)
unique_vars       <- unique(grand_mean_data$variables)
unique_periods    <- unique(grand_mean_data$period)
na_rows <- expand.grid(period = unique_periods,
                       time_hours = -3,
                       variables = unique_vars,
                       treatment = unique_treatments,
                       grand_mean = NA,
                       standard_error = NA)

grand_mean_data <- rbind(grand_mean_data, na_rows) %>%
  arrange( period, time_hours, variables, treatment) #treatment,


# # Perform scaling by group # scale() function standardizes a vector to have mean 0 and standard deviation 1
# grand_mean_data_scaled <- grand_mean_data %>%
#   group_by(variables)  %>%
#   dplyr::mutate(scaled_grand_mean = scale(grand_mean),scaled_standard_error = scale(standard_error)) %>%
#   ungroup()

# overall means
overall_grand_mean <- grand_mean_data %>%
  group_by(time_hours, period, variables) %>%
  summarise(grand_mean = mean(grand_mean))

# #plot fit
# DegreeVsTimeOutside <- ggplot(grand_mean_data, aes(x = time_hours, y = scaled_grand_mean, fill = variables, color= variables, group = variables)) +
#   geom_smooth(data = subset(grand_mean_data, period == "pre"), method = "lm", se = T, linetype = "dashed") +
#   geom_smooth(data = subset(grand_mean_data, period == "post"), method = "lm", se = T, linetype = "solid") +
#   scale_x_continuous(limits = c(min(grand_mean_data$time_hours), max(grand_mean_data$time_hours)), expand = c(0, 0)) +
#   labs(#title = "Prop Time Outside by Time Hours and Treatment",
#     x = "Time Hours since treatment exposure",
#     #y = "scaled variables" #names(variable_list[i])
#     y = "\nScaled\nvariables"
#   ) +
#   #geom_text(aes(x = 10, label = from_p_to_ptext(p_interaction_vars_time))) 
#   STYLE +
#   scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
#   scale_color_manual(values = c("#E69F00", "#56B4E9")) + #inverting the order of scale fill and scale color breaks the plot
#   theme(legend.position = c(.05, .95), # Position legend inside plot area
#         legend.justification = c(0, 1), # Justify legend at top left
#         legend.box.just = "left",
#         legend.direction = "horizontal",
#         legend.background = element_rect(fill='transparent'),
#         legend.title = element_blank()) +
#   guides(fill = guide_legend(nrow = 2,ncol = 1))



#library(tidyr)
grand_mean_data <- subset(grand_mean_data, time_hours != -3)


wide_data <- as.data.frame(grand_mean_data %>%
  tidyr::pivot_wider(id_cols = c(period, time_hours, treatment),
              names_from = variables, 
              values_from = c(grand_mean, standard_error)) %>%
  rename_with(~ gsub("grand_mean_", "", .x), starts_with("grand_mean")) %>%
  rename_with(~ gsub("standard_error_", "se_", .x), starts_with("standard_error")) )


DegreeLine <- ggplot(wide_data, aes(x = time_hours, y = degree, color = treatment, group = treatment)) +
  # Add lines for each treatment
  geom_smooth(data = subset(wide_data, period == "pre"), method = "lm", linetype = "dashed", aes(linetype = "pre")) +
  geom_smooth(data = subset(wide_data, period == "post"),method = "lm", linetype = "solid", aes(linetype = "post")) +
  geom_ribbon(data = subset(wide_data, period == "pre"), stat = "smooth", method = "lm", fill = NA) +
  geom_ribbon(data = subset(wide_data, period == "post"), stat = "smooth", method = "lm", fill = NA) +
  STYLE +
  colFill_treatment +
  colScale_treatment +
  labs(x = "hours since exposure to the treatment") +
  theme(legend.position = "bottom",      # Position legend at the bottom
        legend.justification = "center", # Center the legend
        legend.box.just = "center",      # Center the legend box
        legend.direction = "horizontal", # Make legend horizontal
        legend.background = element_rect(fill='transparent'), # Transparent legend background
        legend.title = element_blank(),  # Remove legend title
        legend.margin = margin(t = 10, b = 10)) # Add some margin to the top and bottom of the legend

PropOutLine <- ggplot(wide_data, aes(x = time_hours, y = `prop time outside`, color = treatment, group = treatment)) +
  # Add lines for each treatment
  geom_smooth(data = subset(wide_data, period == "pre"), method = "lm", linetype = "dashed", aes(linetype = "pre")) +
  geom_smooth(data = subset(wide_data, period == "post"),method = "lm", linetype = "solid", aes(linetype = "post")) +
  geom_ribbon(data = subset(wide_data, period == "pre"), stat = "smooth", method = "lm", fill = NA) +
  geom_ribbon(data = subset(wide_data, period == "post"), stat = "smooth", method = "lm", fill = NA) +
  STYLE +
  colFill_treatment +
  colScale_treatment +
  labs(x = "hours since exposure to the treatment",
       y= "prop. time outside") +
  theme(legend.position = "bottom",      # Position legend at the bottom
        legend.justification = "center", # Center the legend
        legend.box.just = "center",      # Center the legend box
        legend.direction = "horizontal", # Make legend horizontal
        legend.background = element_rect(fill='transparent'), # Transparent legend background
        legend.title = element_blank(),  # Remove legend title
        legend.margin = margin(t = 10, b = 10)) # Add some margin to the top and bottom of the legend

GroomLine <- ggplot(wide_data, aes(x = time_hours, y = `duration grooming received min`, color = treatment, group = treatment)) +
  # Add lines for each treatment
  geom_smooth(data = subset(wide_data, period == "pre"), method = "lm", linetype = "dashed", aes(linetype = "pre")) +
  geom_smooth(data = subset(wide_data, period == "post"),method = "lm", linetype = "solid", aes(linetype = "post")) +
  geom_ribbon(data = subset(wide_data, period == "pre"), stat = "smooth", method = "lm", fill = NA) +
  geom_ribbon(data = subset(wide_data, period == "post"), stat = "smooth", method = "lm", fill = NA) +
  STYLE +
  colFill_treatment +
  colScale_treatment +
  labs(x = "hours since exposure to the treatment",
       y = "duration grooming received (min)") +
  theme(legend.position = "bottom",      # Position legend at the bottom
        legend.justification = "center", # Center the legend
        legend.box.just = "center",      # Center the legend box
        legend.direction = "horizontal", # Make legend horizontal
        legend.background = element_rect(fill='transparent'), # Transparent legend background
        legend.title = element_blank(),  # Remove legend title
        legend.margin = margin(t = 10, b = 10)) # Add some margin to the top and bottom of the legend

###################################################################################################################################
### COMPARE  grooming given and duration of contacts with treated #################################################################
###################################################################################################################################

root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_given_to_treated_min") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2"
names(variable_list) <- c("duration grooming given to treated (min)") # , "prop duration grooming received outside min","duration grooming received outside min"
#transf_variable_list <- c("log"                           )   ######"none", "sqrt" "log","power2"

dur_groom_given_data <- read_data(data_path,which_individuals=c("nurse","forager"))
dur_groom_given_data$duration_of_contact_with_treated_min <- NULL #make sure there are no conflicts


root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_of_contact_with_treated_min") #,"proportion_time_active", "average_bout_speed_pixpersec" ,"total_distance_travelled_pix", "inter_caste_contact_duration") #inter_caste_contact_duration?
names(variable_list) <- c("duration of contact with treated (min)") # ,"prop. time active", "average bout speed pixpersec" ,"total distance travelled pix", "inter caste contact duration")
#transf_variable_list <- c("power0.01"        )#,"none"                  ,"log"                          ,"log"                      ,"sqrt"                   )  ######"none", "sqrt" "log","power2"

dur_contact_treated_data <- read_data(data_path,which_individuals=c("nurse","forager"))

#merge
common_col_names <- intersect(names(dur_groom_given_data), names(dur_contact_treated_data))
common_col_names <- common_col_names[!common_col_names %in% c("nb_frames_outside","nb_frames_inside","prop_time_outside", "proportion_time_active","average_bout_speed_pixpersec","total_distance_travelled_pix","inter_caste_contact_duration","inter_caste_contact_number"  )]
CompareBehavs2 <- dplyr::left_join(dur_groom_given_data, dur_contact_treated_data, by = common_col_names[])

# remove extra cols
CompareBehavs2 <- CompareBehavs2 %>%
  dplyr::select(colony, tag, antID, time_hours, colony_size, treatment, age, period, time_of_day, 
                duration_grooming_given_to_treated_min, duration_of_contact_with_treated_min, task_group
  )

#transform data as above for comparability
CompareBehavs2[!is.na(CompareBehavs2$duration_grooming_given_to_treated_min),"dur_groom_given_to_treat_min_log"]  <- log_transf(CompareBehavs2[!is.na(CompareBehavs2$duration_grooming_given_to_treated_min),"duration_grooming_given_to_treated_min"] )
CompareBehavs2[!is.na(CompareBehavs2$duration_of_contact_with_treated_min),"dur_contact_with_treated_min_log"]  <- log_transf(CompareBehavs2[!is.na(CompareBehavs2$duration_of_contact_with_treated_min),"duration_of_contact_with_treated_min"])

CompareBehavs2$duration_grooming_given_to_treated_min <- NULL
CompareBehavs2$duration_of_contact_with_treated_min <- NULL

#check normality
test_norm(CompareBehavs2$dur_groom_given_to_treat_min_log)
test_norm(CompareBehavs2$dur_contact_with_treated_min_log)

#check PEARSON correlation of vars, by group and period
# Calculate correlation, slope, and mean positions for each combination
calculated_values <- CompareBehavs2 %>%
  group_by(task_group, period) %>%
  summarise(correlation = cor(dur_contact_with_treated_min_log, dur_groom_given_to_treat_min_log),
            slope = coef(lm(dur_groom_given_to_treat_min_log ~ dur_contact_with_treated_min_log))[2], #for plotting
            mean_x = mean(dur_contact_with_treated_min_log, na.rm = TRUE),
            mean_y = mean(dur_groom_given_to_treat_min_log, na.rm = TRUE)) %>%
  ungroup()

  # melt the dataframe to long format
  CompareBehavs2_melt <- reshape(CompareBehavs2,
                           varying = c("dur_groom_given_to_treat_min_log", "dur_contact_with_treated_min_log"),
                           v.names = "measure",
                           times =  c("dur_groom_given_to_treat_min_log", "dur_contact_with_treated_min_log"),
                           timevar = "variables",
                           direction = "long")
  CompareBehavs2_melt$variables <- as.factor(CompareBehavs2_melt$variables)
  #CompareBehavs2_melt$variables <- gsub("_", " ", CompareBehavs2_melt$variables)
  
  ###### MODEL post vs pre
  mod1 <- lmer(measure ~ variables*period + (1|colony) + (1|antID) + (1|time_hours), data = CompareBehavs2_melt[which(CompareBehavs2_melt$task_group=="nurse"),])
  mod2 <- lmer(measure ~ variables*period + (1|colony) + (1|antID) + (1|time_hours), data = CompareBehavs2_melt[which(CompareBehavs2_melt$task_group=="forager"),], control = lmerControl(optimizer ="Nelder_Mead"))
  #output_lmer(mod1)

  p_interaction_vars_nuse <- anova(mod1)["variables:period","Pr(>F)"]
  p_interaction_vars_forag <- anova(mod2)["variables:period","Pr(>F)"]
  
  
  CompareBehavs2$period <- fct_rev(CompareBehavs2$period)
  
  # Create the plot
  #plotting by treatment is meaningless as they all overlap
  GroomingVsContacts <- ggplot(CompareBehavs2, aes(x=dur_contact_with_treated_min_log, 
                             y=dur_groom_given_to_treat_min_log, 
                             color=period)) + 
    # geom_point() +   # Scatter plot of points
    geom_smooth(method='lm', se=T) +   # Linear regression line, without standard error shading
    labs(title=" ", 
         x="Duration of contact (min) (log)", 
         y="Duration of grooming \ngiven (min) (log)") +
    theme_minimal() +
    annotate("text", x = -1, y = -1.5, label =  from_p_to_ptext(max(p_interaction_vars_nuse,p_interaction_vars_forag))) + #as they are the same, create stars once
    geom_text(data=calculated_values, 
              aes(x=mean_x+0.8, y=mean_y+0.1, label=sprintf("r = %.2f", correlation), angle=atan(slope)*180/pi),
              hjust=-0.1, vjust=-0.1, check_overlap = TRUE, size=3,show.legend = FALSE) +
    facet_grid(. ~ task_group) +
    coord_fixed(ratio = 1) +
    STYLE_generic 
  # STYLE_CONT
  

  

###################################################################################################################################
############### CH5: PRE-POST DIFFERENCES PLOTS ###################################################################################
###################################################################################################################################


### PLOT GRIDS ####################################################################################################################


#TO DO's  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# stars not assigned (occasionally)


###################################################################################################################################
### ind_net_properties ### 
## degree
plot_list <- list(ind_treated_net$barplot_delta_period_list$degree,
                  ind_untreated_net_nurse$barplot_delta_period_list$degree,
                  ind_untreated_net_forag$barplot_delta_period_list$degree)

# Set the same y-axis limits for all plots
YLIM_extra <- 1
plot_comps <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps$y_limits) + ggtitle("treated\nnurses")    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = "\n\nΔ degree"),
                                 plot_list[[2]] + ylim(plot_comps$y_limits) + ggtitle("untreated\nnurses")  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 plot_list[[3]] + ylim(plot_comps$y_limits) + ggtitle("untreated\nforagers")+ fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                 align="h")

ind_net_degree <-  cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],
                     ncol=3, rel_widths = c(0.3,0.24,0.24))
# ind_net_degree <- cowplot::plot_grid(
#   cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],
#                      ncol=3, rel_widths = c(0.28,0.24,0.24))
#   , plot_comps$leg, ncol=1, rel_heights = c(0.9, 0.1))


###################################################################################################################################  
### ind_beh_measures ### 3 panels

## prop_time_outside
plot_list <- list(ind_treated_beh$barplot_delta_period_list$prop_time_outside,
                  ind_untreated_beh_nurse$barplot_delta_period_list$prop_time_outside,
                  ind_untreated_beh_forag$barplot_delta_period_list$prop_time_outside)
YLIM_extra <- 0.01
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits)  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = "Δ prop. time\noutside (power0.1)") , # + ggtitle("treated\nnurses")   
                                 plot_list[[2]] + ylim(plot_comps1$y_limits)  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs, #+ ggtitle("untreated\nnurses") 
                                 plot_list[[3]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs, #+ ggtitle("untreated\nforagers")
                                 align="h")
prop_time_outside <- 
  cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]],
                     ncol=3, rel_widths = c(0.30,0.24,0.24))


### ind_beh_measures### 2 panels
# create text boxes for titles (to ensure equal size of all plotted objects)
titlePlots <- cowplot::align_plots(ggplot() + theme_void() + annotate("text", x = 0.2, y = 0.5, label = "untreated\nnurses", family = "Liberation Serif",  size = 4.7, hjust = 0.05, lineheight = 0.8), #fontface = "bold",
                                   ggplot() + theme_void() + annotate("text", x = 0.2, y = 0.5, label = "untreated\nforagers",family = "Liberation Serif", size = 4.7, hjust = 0.05, lineheight = 0.8), #fontface = "bold",
                                   align="v")
titlePlots_untreated <- cowplot::plot_grid(titlePlots[[1]], titlePlots[[2]], ncol=2, rel_widths = c(0.23,0.26))

spacingPlot <- ggplot() + theme_void() + annotate("text", x = 0.2, y = 0.5, label = "", family = "Liberation Serif",  size = 4, hjust = 0.05)
                                  


## duration_of_contact_with_treated_min
plot_list <- list(ind_untreated_beh_nurse$barplot_delta_period_list$duration_of_contact_with_treated_min,
                  ind_untreated_beh_forag$barplot_delta_period_list$duration_of_contact_with_treated_min)
YLIM_extra <- 0.02
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits)    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits)    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
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
YLIM_extra <- 0.03
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
plot_comps1$y_limits[1] <- 0
allplots1 <- cowplot::align_plots(plot_list[[1]]  + scale_y_continuous(limits = c(0, plot_comps1$y_limits[2]), expand = c(0, 0))  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)), # ylim(plot_comps1$y_limits)
                                  plot_list[[2]]  + scale_y_continuous(limits = c(0, plot_comps1$y_limits[2]), expand = c(0, 0))  + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs, # ylim(plot_comps1$y_limits)
                                  align="h")
duration_grooming_given_to_treated_min <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                   ncol=2, rel_widths = c(0.28,0.24))



### COMBINE ind_beh_measures 2 panels together
ind_beh_measures <- cowplot::plot_grid(
  titlePlots_untreated,
  duration_of_contact_with_treated,
  inter_caste_contact_duration,
  duration_grooming_given_to_treated_min,
  spacingPlot,
   ncol=1, rel_heights = c(0.04,0.30,0.30,0.30,0.01), labels = c("D","","E","F"))


###################################################################################################################################
### ind_grooming_received ### 2 panels
#titlePlots_treated <- ggplot() + theme_void() + annotate("text", x = 0.2, y = 0.5, label = "treated nurses", family = "Liberation Serif",  size = 4.7, lineheight = 0.8)
titlePlotsUnt <- cowplot::align_plots(ggplot() + theme_void() + annotate("text", x = 0.2, y = 0.5, label = "treated\nnurses", family = "Liberation Serif",  size = 4.7, lineheight = 0.8), #fontface = "bold",
                                   ggplot() + theme_void() + annotate("text", x = 0.2, y = 0.5, label = " ",family = "Liberation Serif", size = 4.7,  lineheight = 0.8), #fontface = "bold",
                                   align="v")
titlePlots_treated <- cowplot::plot_grid(titlePlotsUnt[[1]], titlePlotsUnt[[1]], titlePlotsUnt[[2]], ncol=3, rel_widths = c(0.2,0.2,0.2))



plot_list <- list(ind_treated_grooming$barplot_delta_period_list$duration_grooming_received_min,
                  ind_treated_grooming$barplot_delta_period_list$number_contacts_received)
YLIM_extra <- 0.18
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
plot_comps1$y_limits[1] <- 0
allplots1 <- cowplot::align_plots(plot_list[[1]] + scale_y_continuous(limits = c(0, plot_comps1$y_limits[2]), expand = c(0, 0)) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0))                    + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[2]]$labels$y)),
                                  align="h")
treated_grooming <- cowplot::plot_grid( titlePlots_treated,
  cowplot::plot_grid(allplots1[[1]], allplots1[[2]],spacingPlot,
                     ncol=3, rel_widths = c(0.24,0.24,0.24)) # extra widths for centering
  , ncol=1, rel_heights = c(0.1, 0.9))



###################################################################################################################################
#MASTER BEHAVIOURAL PLOT
plot_list <- list(ind_net_degree,
                  prop_time_outside,
                  treated_grooming)

# Set the same y-axis limits for all plots
YLIM_extra <- 1
plot_comps <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
allplots <- cowplot::align_plots(plot_list[[1]],
                                 plot_list[[2]],
                                 plot_list[[3]],
                                 align="v")

leftPanelBeh <-  cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],
                                      ncol=1, rel_heights  = c(0.26,0.21,0.24), 
                                    labels = c('A', 'B',"C"))
FullBehGrid <- cowplot::plot_grid(
cowplot::plot_grid(leftPanelBeh,ind_beh_measures, ncol=2, rel_widths = c(0.9,0.6)),
plot_comps1$leg,
ncol=1, rel_heights = c(0.9, 0.1))
# ind_net_degree <- cowplot::plot_grid(
#   cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],
#                      ncol=3, rel_widths = c(0.28,0.24,0.24))
#   , plot_comps$leg, ncol=1, rel_heights = c(0.9, 0.1))


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
### comparing timeline of degree and time_outside line_plots  ### 2 panels
plot_list <- list(DegreeLine, PropOutLine, GroomLine)
YLIM_extra <- 0
#plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
leg <-cowplot::get_legend(plot_list[[1]] + guides(fill = guide_legend(nrow = 2,ncol = 2)))

#inherit legend from elsewhere
allplots1 <- cowplot::align_plots(plot_list[[1]] + theme(aspect.ratio = 0.5) +  guides(fill = "none") + theme(legend.position = "none") ,
                                  plot_list[[2]] + theme(aspect.ratio = 0.5) + guides(fill = "none") + theme(legend.position = "none") ,
                                  plot_list[[3]] + theme(aspect.ratio = 0.5) +  guides(fill = "none") + theme(legend.position = "none") + labs(y = split_title(plot_list[[3]]$labels$y)),
                                  align="v")
DegreeVSPropOut <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]], leg,
                     ncol=2, rel_widths = c(0.24,0.24))

# DegreeVSPropOut <- cowplot::plot_grid(
#   cowplot::plot_grid(allplots1[[1]], allplots1[[2]], allplots1[[3]], 
#                      ncol=2, rel_widths = c(0.24,0.24,0.20))
#   , leg, ncol=1, rel_heights = c(0.9, 0.15))

###################################################################################################################################
### collective_net_properties ### 
plot_list <- list(coll_rescal_net$barplot_delta_period_list$modularity_FacetNet,
                  coll_rescal_net$barplot_delta_period_list$clustering,
                  coll_rescal_net$barplot_delta_period_list$task_assortativity,
                  coll_rescal_net$barplot_delta_period_list$density,
                  coll_rescal_net$barplot_delta_period_list$efficiency,
                  coll_rescal_net$barplot_delta_period_list$degree_mean)
# Set the same y-axis limits for all plots
YLIM_extra <- 0.01
#have 2 scales: 1 for top row (measures expected to increase), 1 for bottom row (measures expected to decrease), 
plot_compsA <- multi_plot_comps(plot_list[3],ylim_extra=YLIM_extra)
plot_compsA$y_limits[1] <- -plot_compsA$y_limits[2]  #minor adjustments to ensure alignment
plot_compsB <- multi_plot_comps(plot_list[4:6],ylim_extra=YLIM_extra,legsize=10)
plot_compsC <- multi_plot_comps(plot_list[1:2],ylim_extra=YLIM_extra)
plot_compsC$y_limits[1] <- -plot_compsC$y_limits[2]  #minor adjustments to ensure alignment

allplots <- cowplot::align_plots(plot_list[[1]] + ylim(plot_compsC$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]] + ylim(plot_compsC$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[3]] + ylim(plot_compsA$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = paste0("\n",plot_list[[3]]$labels$y)),
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

# SavePrint_plot(
#   plot_obj = ind_net_degree,
#   plot_name = "ind_net_degree",
#   plot_size = c(430/ppi, 300/ppi),
#   # font_size_factor = 4,
#   dataset_name = "Grid",
#   save_dir = figurefolder
# )
# 
# SavePrint_plot(
#   plot_obj = prop_time_outside,
#   plot_name = "prop_time_outside",
#   plot_size = c(430/ppi, 300/ppi),
#   # font_size_factor = 4,
#   dataset_name = "Grid",
#   save_dir = figurefolder
# )
# 
# # ISSUE WITH PVAL IN MODEL!!!!!!!!!!!!!!!!!!!!!!!!!
# SavePrint_plot(
#   plot_obj = ind_beh_measures,
#   plot_name = "ind_beh_measures",
#   plot_size = c(330/ppi, 900/ppi), #extra length required to get y-axis labeled and unlabeled of same size
#   # font_size_factor = 4,
#   dataset_name = "Grid",
#   save_dir = figurefolder
# )

# # remove inter_Caste_grooming_duration
# SavePrint_plot(
#   plot_obj = treated_grooming, 
#   plot_name = "treated_grooming",
#   plot_size = c(430/ppi, 230/ppi), #extra length required to fix the different digits on y-axis
#   # font_size_factor = 4,
#   dataset_name = "Grid",
#   save_dir = figurefolder
# )


# FullBehGrid
SavePrint_plot(
  plot_obj = FullBehGrid,
  plot_name = "FullBehGrid",
  plot_size = c(650/ppi, 660/ppi), #extra length required to fix the different digits on y-axis
  # font_size_factor = 4,
  dataset_name = "",
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
  plot_obj = DegreeVSPropOut, 
  plot_name = "DegreeVSPropOut",
  plot_size = c(520/ppi, 280/ppi), #extra length required to fix the different digits on y-axis
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)
 
#it looks terrible
SavePrint_plot(
  plot_obj = GroomingVsContacts, 
  plot_name = "GroomingVsContacts",
  plot_size = c(700/ppi, 300/ppi),
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


###################################################################################################################################
### Simulated Load plots ##########################################################################################################
###################################################################################################################################


###################################################################
###            collective analysis - without rescaling ############
###################################################################

#### experimentally exposed ####

root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="")
pattern="collective_simulation_results_observed"
variable_list <-  c("Prevalence", "Mean_load", "Load_skewness", "Queen_load", "logistic_r","Prop_high_level","Prop_low_level")
names(variable_list) <-  c("Prevalence", "Mean load", "Load skewness", "Queen load", "Transmission rate","Prop. high level","Prop. low level")
transf_variable_list <- c("none"       ,"none"        ,"none"           ,"log"      ,"log", "none"       ,"none"   )   ######"none", "sqrt" "log","power2"

coll_no_rescal_sim_EXP_seed <- collective_analysis_no_rescal(data_path,showPlot=F)
warning(paste("-Prevalence", "-Mean load", "-Load skewness","-Prop. high level","-Prop. low level",
              "don't meet normality regardless of transformation",sep= "\n"))
# but "-Load skewness" "-Prop. high level","-Prop. low level" have a high Shapiro-Wilk normality test statistic val (>0.9)

#https://stackoverflow.com/questions/68915173/how-do-i-fit-a-quasi-poisson-model-with-lme4-or-glmmtmb
#model <- glmer(Prevalence ~ period*treatment + (1|colony) ,data=data, family = poisson)

#### random_workers_seeds #### 
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/random_workers_seeds",sep="")

coll_no_rescal_sim_RAN_seed <- collective_analysis_no_rescal(data_path,showPlot=F)
warning(paste("-Prevalence", "-Mean load", "-Load skewness",
              "don't meet normality regardless of transformation",sep= "\n"))


#### nurses_seeds #### 
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/nurses_seeds",sep="")

coll_no_rescal_sim_NUR_seed <- collective_analysis_no_rescal(data_path,showPlot=F)
warning(paste("Prevalence", "Mean load","Prop. high level","Prop. low level",
              "don't meet normality regardless of transformation",sep= "\n"))


#### foragers_seeds #### 
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/foragers_seeds",sep="")

coll_no_rescal_sim_FOR_seed <- collective_analysis_no_rescal(data_path,showPlot=F)
warning(paste("-Prevalence", "-Load skewness","Queen_load",
              "don't meet normality regardless of transformation",sep= "\n"))


###################################################################
###            individual analysis                     ############
###################################################################

#### experimentally exposed ####

root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="")
pattern="individual_simulation_results_observed"
variable_list <-  c("simulated_load", "transmission_latency", "transmission_rank")
names(variable_list) <-  c("simulated load", "transmission latency", "transmission rank")
transf_variable_list <- c("none"           ,"log"                ,"none")   ######"none", "sqrt" "log","power2"

ind_untreated_sim_nurse_EXP_seed <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_sim_forag_EXP_seed <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"
#probability_of_transmission for nurses and foragers is almost always 1, can't really be modeled

#### random_workers_seeds #### 
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/random_workers_seeds",sep="")

ind_untreated_sim_nurse_RAN_seed <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_sim_forag_RAN_seed <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"

#### nurses_seeds #### 
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/nurses_seeds",sep="")

ind_untreated_sim_nurse_NUR_seed <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_sim_forag_NUR_seed <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"

#### foragers_seeds #### 
data_path <- paste(root_path,"/transmission_simulations/pre_vs_post_treatment/foragers_seeds",sep="")

ind_untreated_sim_nurse_FOR_seed <- individual_ONE_analysis(data_path,which_individuals="nurse",showPlot=F) ## "treated","queen","nurse","forager"
ind_untreated_sim_forag_FOR_seed <- individual_ONE_analysis(data_path,which_individuals="forager",showPlot=F) ## "treated","queen","nurse","forager"

##"transmission rank" is a measure of the average order in which an individual becomes contaminated in the simulations, with a lower rank indicating earlier contamination. It can be interpreted as a measure of how quickly or easily an individual tends to become infected in the simulated disease transmission scenario.

###################################################################################################################################
### PLOT GRIDS ####################################################################################################################
###################################################################################################################################


###################################################################################################################################
### individual_sim_properties ###

## make all ind simulation plots with loop by SEED
GROUP <- c("EXP","RAN","NUR","FOR")
names(GROUP) <- c("experimentally exposed","random workers","nurses","forager")

for ( SEED in GROUP) {
  
  SEED_group <- paste(SEED,"seed",sep="_")
  
  # create text boxes for seed
  titlePlotsSEED <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste("SEED",names(GROUP[which(GROUP %in% SEED)]),sep = " : "), family = "Liberation Serif",  size = 4, hjust = 0.5)
  #titlePlots_untreated <- cowplot::plot_grid(titlePlots[[1]], titlePlots[[2]], ncol=2, rel_widths = c(0.28,0.24))
  
  #data
  IND_UNTREATED_SIM_nurse <- get(paste0("ind_untreated_sim_nurse_",SEED_group))
  IND_UNTREATED_SIM_forag <- get(paste0("ind_untreated_sim_forag_",SEED_group))

## simulated_load
plot_list <- list(IND_UNTREATED_SIM_nurse$barplot_delta_period_list$simulated_load,
                  IND_UNTREATED_SIM_forag$barplot_delta_period_list$simulated_load)
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=0.005)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) +  fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")
individual_sim_SL <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                       ncol=2, rel_widths = c(0.28,0.24))
## transmission_latency
plot_list <- list(IND_UNTREATED_SIM_nurse$barplot_delta_period_list$transmission_latency,
                  IND_UNTREATED_SIM_forag$barplot_delta_period_list$transmission_latency)
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=0.5)
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs,
                                  align="h")
individual_sim_TL <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                   ncol=2, rel_widths = c(0.28,0.24))

## transmission_rank
plot_list <- list(IND_UNTREATED_SIM_nurse$barplot_delta_period_list$transmission_rank,
                  IND_UNTREATED_SIM_forag$barplot_delta_period_list$transmission_rank)
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=0.8)
allplots1 <- cowplot::align_plots(plot_list[[1]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)), # ylim(plot_comps1$y_limits)
                                  plot_list[[2]] + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") + remove_y_labs, # ylim(plot_comps1$y_limits)
                                  align="h")
individual_sim_TR <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                             ncol=2, rel_widths = c(0.28,0.24))


### COMBINE ind_beh_measures 2 panels together
individual_sim_properties <- cowplot::plot_grid(
  titlePlotsSEED,
  titlePlots_untreated,
  individual_sim_SL,
  individual_sim_TL,
  individual_sim_TR,
  plot_comps1$leg, ncol=1, rel_heights = c(0.05,0.05,0.30,0.30,0.30, 0.05))


#still lightly offset
SavePrint_plot(
  plot_obj = individual_sim_properties, 
  plot_name = paste("individual_sim_properties",SEED,sep="_"),
  plot_size = c(350/ppi, 900/ppi), #extra length required to get y-axis labeled and unlabeled of same size
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)


###################################################################################################################################
### collective_sim_properties ### 

#data
COLL_NO_RESCAL_SIM <- get(paste0("coll_no_rescal_sim_",SEED_group))

plot_list <- list(COLL_NO_RESCAL_SIM$barplot_delta_period_list$Prevalence,
                  COLL_NO_RESCAL_SIM$barplot_delta_period_list$Mean_load,
                  COLL_NO_RESCAL_SIM$barplot_delta_period_list$Load_skewness,
                  COLL_NO_RESCAL_SIM$barplot_delta_period_list$Queen_load,
                  COLL_NO_RESCAL_SIM$barplot_delta_period_list$logistic_r,
                  COLL_NO_RESCAL_SIM$barplot_delta_period_list$Prop_high_level,
                  COLL_NO_RESCAL_SIM$barplot_delta_period_list$Prop_low_level)
# Set the same y-axis limits for all plots
YLIM_extra <- 0.001
# #have 2 scales: 1 for top row (measures expected to increase), 1 for bottom row (measures expected to decrease), 
# plot_compsA <- multi_plot_comps(plot_list[1:2],ylim_extra=YLIM_extra)
# plot_compsB <- multi_plot_comps(plot_list[3:4],ylim_extra=YLIM_extra)
allplots <- cowplot::align_plots(plot_list[[1]]  + ylim(-0.009,0.009) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[2]]  + ylim(-0.007,0.007) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[3]]  + ylim(-0.4,0.4) + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[4]]  + ylim(-0.15,0.15)           + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[5]]                               + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[6]]                               + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 plot_list[[7]]                               + fixed_aspect_theme  + remove_x_labs + guides(fill = "none") ,
                                 align="h")
collective_sim_properties <- cowplot::plot_grid(
  cowplot::plot_grid(titlePlotsSEED, allplots[[1]], allplots[[2]], allplots[[3]],allplots[[4]], allplots[[5]],allplots[[6]], allplots[[7]],
                     ncol=1, rel_widths = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1))
  #, plot_compsB$leg
  ,ncol=1, rel_heights = c(0.9, 0.1))
# width_pixels <- 200 
# height_pixels <- 900

SavePrint_plot(
  plot_obj = collective_sim_properties, 
  plot_name = paste("collective_sim_properties",SEED,sep="_"),
  plot_size = c(200/ppi, 900/ppi), #extra length required to get y-axis labeled and unlabeled of same size
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)

}











###################################################################################################################################
### Experimentally measured and simulated M.brunneum transmission in pathogen-exposed colonies ###

#TEMP (The green line highlights the value at which there is an inversion in the sign of the density difference. This value was used as athreshold to distinguish high from low simulated loads in allsubsequent analyses)
high_threshold <- 0.0258
 
# logging not enough for qPCR data, as it has kurtois = 27.28894
# sqrt                                                  27.46589
# power0.01                                             70.34206
# power0.001                                            26.71343
# power0.0001                                           26.64026
# Box_Cox fails

######Open pdf Figure 2 #####
pdf(file=paste(figurefolder,"/Figure2.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col*0.5,height=double_col,pointsize=pointsize_less_than_2row2col) #height=page_height*0.25
######Plot qPCR data #####
full_statuses_names_ori <- full_statuses_names
full_statuses_names[full_statuses_names%in%c("Nurses","Untreated\nnurses")] <- "Nurses\n"

statuses_colours_ori <- statuses_colours
statuses_colours[names(statuses_colours)%in%c("queen","forager","nurse")] <- "black"

par(pars)
par_mar_ori <- par()$mar
par(mar=par_mar_ori+c(1,0,1,1))
widz <- c(2,1.5)
# layout(matrix(c(1,2),nrow=1),widths=widz)
translated_high_threshold <- plot_qpcr(experiments=c("main_experiment")) #thickness of violplots controlled by "wex" in ViolPlot function
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
dev.off()
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










#####################################################################################
##### CH4: 
#####################################################################################

## start with individual level measure as they don't scale.
## large vs small in pre (2 treatments) - sum exposure treatments
## subjects: nurses and foragers (individually)

## this has no post-hocs to be performed 

# Constitutive organisation compared to random (Treatment x Obs-Ran):
# - intra_caste_over_inter_caste_WW_contact_duration,  QNurse_over_QForager_contact_duration 
# 
# Proxy of DOL (Treatment):
# - task_assortativity, clustering, degree_mean, degree_maximum, density, diameter, efficiency, modularity (how to normalise given that it is only pre?) 
# 
# Individual behaviour 
# - Prop. Time outside (foragers)
# - Cross-worker variation in prop. Time outside 
# 
# Simulations outcomes (Treatment, Obs-Ran in Big, Obs-Ran in Small):
# - Big vs Small, Random vs observed in big, Random vs observed in small
# - Actual treated workers as seeds, Do results hold when using foragers as seeds?



###################################################################################################################################
############### CH4: PRE DIFFERENCES PLOTS ########################################################################################
###################################################################################################################################


### PLOT GRIDS ####################################################################################################################


### ind_net_properties ### 
## degree
plot_list <- list(ind_untreated_net_nurse_PRE$barplot_delta_period_list$degree,
                  ind_untreated_net_forag_PRE$barplot_delta_period_list$degree)

YLIM_extra <- 20
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
plot_comps1$y_limits <- c(0,plot_comps1$y_limits[2])
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits)   + fixed_aspect_theme_PRE  +  ggtitle("nurses") + guides(fill = "none"),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits)   + fixed_aspect_theme_PRE  +  ggtitle("foragers") + guides(fill = "none") + remove_y_labs, #ylim(plot_comps1$y_limits) +
                                  align="h")
ind_net_degree_PRE <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                       ncol=2, rel_widths = c(0.28,0.24))

### ind_beh_measures### 2 panels
# create text boxes for titles (to ensure equal size of all plotted objects)
titlePlots <- cowplot::align_plots(ggplot() + theme_void(), + annotate("text", x = 0.5, y = 0.5, label = "nurses", family = "Liberation Serif",  size = 4, hjust = 0.5), #fontface = "bold",
                                   ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "foragers",family = "Liberation Serif", size = 4, hjust = 0.5), #fontface = "bold",
                                   align="h")
titlePlots_untreated <- cowplot::plot_grid(titlePlots[[1]], titlePlots[[2]], ncol=2, rel_widths = c(0.28,0.24))


## prop_time_outside
plot_list <- list(ind_untreated_beh_nurse_PRE$barplot_delta_period_list$prop_time_outside,
                  ind_untreated_beh_forag_PRE$barplot_delta_period_list$prop_time_outside)

YLIM_extra <- 0.3
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
plot_comps1$y_limits[1] <- 0
#plot_comps1$y_limits <- c(0,plot_comps1$y_limits[2])
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits)   + fixed_aspect_theme_PRE  + guides(fill = "none"),
                                  plot_list[[2]] + ylim(plot_comps1$y_limits)   + fixed_aspect_theme_PRE  + guides(fill = "none") + remove_y_labs, #ylim(plot_comps1$y_limits) +
                                  align="h")
prop_time_outside_PRE <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                            ncol=2, rel_widths = c(0.28,0.24))


warning("- SIG STARS LOST BY NURSES AS TOO HIGH UP")


# ## duration_of_contact_with_treated_min
# plot_list <- list(ind_untreated_beh_nurse_PRE$barplot_delta_period_list$duration_of_contact_with_treated_min,
#                   ind_untreated_beh_forag_PRE$barplot_delta_period_list$duration_of_contact_with_treated_min)
# YLIM_extra <- 0.5
# plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
# #plot_comps1$y_limits <- c(0,plot_comps1$y_limits[2])
# allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits)    + fixed_aspect_theme_PRE  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),
#                                   plot_list[[2]] + ylim(plot_comps1$y_limits)    + fixed_aspect_theme_PRE  + remove_x_labs + guides(fill = "none") + remove_y_labs,
#                                   align="h")
# duration_of_contact_with_treated_PRE <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
#                                                        ncol=2, rel_widths = c(0.28,0.24))
## inter_caste_contact_duration
plot_list <- list(ind_untreated_beh_nurse_PRE$barplot_delta_period_list$inter_caste_contact_duration,
                  ind_untreated_beh_forag_PRE$barplot_delta_period_list$inter_caste_contact_duration)
YLIM_extra <- 0.8
plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
plot_comps1$y_limits <- c(0,plot_comps1$y_limits[2])
allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits)  + fixed_aspect_theme_PRE  + remove_x_labs + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)),# + ylim(plot_comps1$y_limits) 
                                  plot_list[[2]] + ylim(plot_comps1$y_limits)  + fixed_aspect_theme_PRE  + remove_x_labs + guides(fill = "none") + remove_y_labs,# + ylim(plot_comps1$y_limits) 
                                  align="h")
inter_caste_contact_duration_PRE <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
                                                   ncol=2, rel_widths = c(0.28,0.24))

# ## duration_grooming_given_to_treated_min
# plot_list <- list(ind_untreated_grooming_nurse_PRE$barplot_delta_period_list$duration_grooming_given_to_treated_min,
#                   ind_untreated_grooming_forag_PRE$barplot_delta_period_list$duration_grooming_given_to_treated_min)
# YLIM_extra <- 0.1
# plot_comps1 <- multi_plot_comps(plot_list,ylim_extra=YLIM_extra)
# plot_comps1$y_limits[1] <- 0
# allplots1 <- cowplot::align_plots(plot_list[[1]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme_PRE  + guides(fill = "none") + labs(y = split_title(plot_list[[1]]$labels$y)), # ylim(plot_comps1$y_limits)
#                                   plot_list[[2]] + ylim(plot_comps1$y_limits) + fixed_aspect_theme_PRE  + guides(fill = "none") + remove_y_labs, # ylim(plot_comps1$y_limits)
#                                   align="h")
# duration_grooming_given_to_treated_min_PRE <- cowplot::plot_grid(allplots1[[1]], allplots1[[2]],
#                                                              ncol=2, rel_widths = c(0.28,0.24))



### COMBINE ind_beh_measures 2 panels together
ind_beh_measures_PRE <- cowplot::plot_grid(
  titlePlots_untreated,
  prop_time_outside_PRE,
  #duration_of_contact_with_treated_PRE,
  inter_caste_contact_duration_PRE,
  #duration_grooming_given_to_treated_min_PRE,
  ncol=1, rel_heights = c(0.08,0.30,0.30))

warning("- FIX OFFSET OF VARS BY FIXING PLOT PROPORTIONS  rel_widths = c(0.28,0.24)
        \n- Occasionally the sig stars are not plotted!!!! ")






#done

SavePrint_plot(
  plot_obj = ind_net_degree_PRE,
  plot_name = "ind_net_degree_PRE",
  plot_size = c(230/ppi, 300/ppi),
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)

# SavePrint_plot(
#   plot_obj = prop_time_outside_PRE,
#   plot_name = "prop_time_outside_PRE",
#   plot_size = c(230/ppi, 300/ppi),
#   # font_size_factor = 4,
#   dataset_name = "Grid",
#   save_dir = figurefolder
# )

SavePrint_plot(
  plot_obj = ind_beh_measures_PRE, 
  plot_name = "ind_beh_measures_PRE",
  plot_size = c(250/ppi, 600/ppi), #extra length required to get y-axis labeled and unlabeled of same size
  # font_size_factor = 4,
  dataset_name = "Grid",
  save_dir = figurefolder
)



########################################################################

## RANDOM VS OBSERVED (split by size)

### DOL MEASURES
data_path <- paste(disk_path,"/main_experiment/processed_data/collective_behaviour/random_vs_observed",sep="")
pattern="interactions.dat"
variable_list        <- c("intra_caste_over_inter_caste_WW_contact_duration","QNurse_over_QForager_contact_duration")
names(variable_list) <- c("Within/Between-task contact","Q-N/Q-F contacts")
# DOL_plots <- plot_age_dol(data_path,experiments=c("main_experiment"))

starting_data <- NULL; warning("I've no idea what I did wrong but, for the love of god,without 'starting_data' in the global env the function crashes ")
dol_RandObs <- plot_observed_vs_random(data_path,experiments=c("main_experiment"))

# solved by calculating cohen's d differences
# # BIG VS SMALL (no random, no plot)
# data_path <- paste(disk_path,"/main_experiment/processed_data/collective_behaviour/random_vs_observed",sep="")
# pattern="interactions.dat"
# variable_list <- c("intra_caste_over_inter_caste_WW_contact_duration","QNurse_over_QForager_contact_duration")
# names(variable_list) <- c("Within/Between-task contact","Q-N/Q-F contacts")
# transf_variable_list <- c("log"      ,"log"   )   ######"none", "sqrt" "log","power2"
# 
# Constitutive_organisation_ratios <- collective_analysis_no_rescal(data_path,showPlot=F)
# Constitutive_organisation_ratios$stats_outcomes

### NETWORK MEASURES
# Science 2018: S2. Constitutive properties of the ant social networks compared to randomized networks
data_path <- paste(disk_path,"/main_experiment/processed_data/network_properties_edge_weights_duration/random_vs_observed",sep="")
pattern="network_properties"
variable_list        <- c("modularity_FacetNet","clustering","task_assortativity","efficiency","degree_mean","density")
names(variable_list) <- c("modularity","clustering","task assortativity","efficiency","mean degree","density")

starting_data <- NULL; warning("I've no idea what I did wrong but, for the love of god,without 'starting_data' in the global env the function crashes ")
net_RandObs <- plot_observed_vs_random(data_path,experiments=c("main_experiment"))

# PER TASK GROUP network measures 
#efficiency_FOR efficiency_NUR net_mean_dist_FOR net_mean_dist_NUR density_FOR density_NUR

data_path <- paste(disk_path,"/main_experiment/processed_data/network_properties_edge_weights_duration/random_vs_observed",sep="")
pattern="network_properties"
variable_list        <- c("efficiency_FOR","efficiency_NUR","net_mean_dist_FOR","net_mean_dist_NUR","density_FOR","density_NUR")
names(variable_list) <- c("efficiency_FOR","efficiency_NUR","net_mean_dist_FOR","net_mean_dist_NUR","density_FOR","density_NUR")

starting_data <- NULL; warning("I've no idea what I did wrong but, for the love of god,without 'starting_data' in the global env the function crashes ")

net_RandObs_TASK <- plot_observed_vs_random_TASK(data_path,experiments=c("main_experiment"))


######Extended data 3: simulation results from different seeds, obs vs. random, COLLECTIVE #####
# Science 2018: Fig. S3. Disease spread simulations over pre-treatment observed (Obs) and randomized (Rd) networks: colony-level results
#SET seeds
seeds <- c("random_workers","nurses","foragers")
names(seeds) <- c("Seeds = random workers","Seeds = nurses","Seeds = foragers" )

sim_RandObs_plots <- list()
sim_RandObs_formatted_outputs <- list() 
for (seed in seeds) {
data_path <- paste(disk_path,"/main_experiment/transmission_simulations/random_vs_observed/",seed,sep="")
pattern="collective_simulation_results_"
variable_list  <- c("logistic_r","Mean_load","Queen_load") #,"Load_skewness","Prop_high_level","Prop_low_level" 
names(variable_list) <- c("Transmission rate","Mean simulated load (W)","Simulated load (Q)") # ,"Simulated load skewness (W)","Prop. high level","Prop. low level"

sim_RandObs <- plot_observed_vs_random(data_path,experiments=c("main_experiment"),seedTitle=T)

sim_RandObs_plots[seed] <- sim_RandObs$saved_plot
sim_RandObs_formatted_outputs[seed] <- list(sim_RandObs$formatted_outcomes)
}



######Extended data 5: simulation results from different seeds, observed #####
# Science 2018: Fig. S5. Disease spread simulations run over pre-treatment observed networks: effect of disease origin.
seeds <- c("random_workers","nurses","foragers")
names(seeds) <- c("R","N","F")
#seeds <- seeds[c(2,3,4,5,6,1)]
#color_pal <- c(GetColorHex("grey50"),GetColorHex("grey70"),GetColorHex("grey20"),statuses_colours["forager"],statuses_colours["occasional_forager"],statuses_colours["nurse"])
color_pal <- c(statuses_colours["untreated"],statuses_colours["forager"],statuses_colours["nurse"])
#color_pal <- color_pal[c(2,3,4,5,6,1)]
names(color_pal) <- seeds

variables <- c("logistic_r","Prevalence","Mean_load","Load_skewness","Prop_high_level","Prop_low_level" ,"Queen_load")
names(variables) <- c("Transmission rate","Prevalence","Mean simulated load (W)","Simulated load skewness (W)","Prop. high level","Prop. low level","Simulated load (Q)")
transf <- c("none","none","none","none","none","none","none")
names(transf) <- variables

#pdf(file=paste(figurefolder,"/Extended_data_5.pdf",sep=""),family=text_font,font=text_font,bg="white",width=three_col,height=page_height*0.49,pointsize=pointsize_more_than_2row2col)
par(pars)
par(mar=c(2.6,2,2.8,1))
layout(matrix(c(1:(2*ceiling(length(variables)/2))), ceiling(length(variables)/2), 2, byrow = T))
plot_seeds(experiments="main_experiment",seeds=seeds,variables=variables,transf=transf,color_pal=color_pal)
#dev.off()


###########################################
###SAVE THESE PLOTSSSSS

############# DOL MEASURES
# plot_list <- list(Entropy_size$entropy_plot,
#                   dol_RandObs$saved_plot,
#                   SocialMaturity)
# 
# allplots <- cowplot::align_plots(plot_list[[1]]  ,
#                                  plot_list[[2]]    ,
#                                  plot_list[[3]] ,
#                                  align="h")
# 
#   cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]],
#                      ncol=3, rel_widths = c(0.2,0.2,0.6))
# 
# 

  # Convert base R plot to a grid object
  grid.newpage()
  dol_RandObs$saved_plot#plot(1:10) # Replace this with your base R plot: dol_RandObs$saved_plot
  grid.echo()
  base_plot <- grid.grab()
  
  library(cowplot)
  library(ggpubr)
  
  
  # Reduce the font size of the top row plots
  Entropy_size$entropy_plot <- Entropy_size$entropy_plot +
    theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6))
  
  SocialMaturity <- SocialMaturity +
    theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6))
  
  # Arrange ggplot objects horizontally
  top_row <- ggarrange(Entropy_size$entropy_plot, SocialMaturity, 
                       ncol = 2, widths = c(0.4, 0.6),labels = c("A", "B"))
  
  # Use plot_grid to arrange the top row and the base plot in two rows
  pre_DOL_plot <- cowplot::plot_grid(top_row, base_plot,  nrow = 2, rel_heights = c(0.3, 0.6),
                     labels = c("C"))
  
  
  #DON'T TOUCH IT, IT IS VERY DELICATE
  # THERE ARE SOME GRAPHICAL ABBERATIONS ON THE TOP ROW (RAN MISSING and top trimmed) BUT THAT CAN BE FIXED LATER
  SavePrint_plot(
    plot_obj = pre_DOL_plot,
    plot_name = "pre_DOL_plot",
    plot_size = c(480/ppi, 450/ppi),
    # font_size_factor = 4,
    dataset_name = "Grid",
    save_dir = figurefolder
  )
  
  # Entropy_size$entropy_plot <- Entropy_size$entropy_plot + theme(plot.margin = unit(c(0.5, 1.5, 1.5, 0.5), "cm"))
  # 
  # # Now use the 'grid' graphical object with cowplot::plot_grid along with other ggplots
  # allplots <- cowplot::align_plots(Entropy_size$entropy_plot, base_plot, SocialMaturity, align="h")
  # 
  # cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]], ncol=3, rel_widths = c(0.2,0.4,0.4))
  # 

  
  # # Assuming Entropy_size$entropy_plot and SocialMaturity are ggplot objects
  # # and base_plot is the grid object converted from a base R plot
  # 
  # # Customize the margin of the plot
  # Entropy_size$entropy_plot <- Entropy_size$entropy_plot + 
  #   theme(plot.margin = unit(c(0.5, 1.5, 1.5, 0.5), "cm"))
  # 
  # # Now use the 'grid' graphical object with cowplot::plot_grid along with other ggplots
  # allplots <- cowplot::align_plots(Entropy_size$entropy_plot, base_plot, SocialMaturity, align="h")
  # 
  # # Define the layout. 1:2 implies the first plot spans 2 columns, the second plot 1 column, and the third plot 3 columns.
  # layout_matrix <- rbind(c(1, 2, 2),
  #                        c(3, 3, 3))
  # 
  # cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]], 
  #                    layout_matrix = layout_matrix, 
  #                    rel_heights = c(0.5, 0.5))  # Modify the rel_heights to adjust the relative heights of the rows
  # 
  

# pdf(file=paste(figurefolder,"/Figure_DOL_RandObs_plots.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col*1.6,height=double_col*0.7,pointsize=pointsize_less_than_2row2col) #height=page_height*0.25
# dol_RandObs$saved_plot
# dev.off()




### NETWORK MEASURES
  pdf(file=paste(figurefolder,"/Figure_sim_NET_RandObs_plots.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col*1.6,height=double_col*0.7,pointsize=pointsize_less_than_2row2col) #height=page_height*0.25
  replayPlot(net_RandObs$saved_plot$density)
  dev.off()

############# SIMULATION
  # grid.newpage()
  # sim_RandObs_plots$foragers
  # grid.echo()
  # base_plot1 <- grid.grab()
  # 
  # grid.newpage()
  # sim_RandObs_plots$random_workers
  # grid.echo()
  # base_plot2 <- grid.grab()
  # 
  # grid.newpage()
  # sim_RandObs_plots$nurses
  # grid.echo()
  # base_plot3 <- grid.grab()
  # 
  # cowplot::plot_grid(base_plot1, base_plot2, base_plot3,  nrow = 2, rel_heights = c(0.3, 0.3))
  
pdf(file=paste(figurefolder,"/Figure_sim_RandObs_plots_ALL.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col*0.8,height=double_col*0.7,pointsize=pointsize_less_than_2row2col) #height=page_height*0.25
replayPlot(sim_RandObs_plots$foragers)
replayPlot(sim_RandObs_plots$random_workers)
replayPlot(sim_RandObs_plots$nurses)
dev.off()

pdf(file=paste(figurefolder,"/Figure_sim_RandObs_plots_TASK.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col*1.6,height=double_col*0.7,pointsize=pointsize_less_than_2row2col) #height=page_height*0.25
replayPlot(net_RandObs_TASK$saved_plot$net_mean_dist_small)
dev.off()



# ################ EXAMPLE CODE TO BUILD FUNCTION TO LATER INSERT IN THE RAN VS OBS SCRIPT
# 
# ### DOL MEASURES
# # data_path <- paste(disk_path,"/main_experiment/processed_data/collective_behaviour/random_vs_observed",sep="")
# # pattern="interactions.dat"
# # variable_list        <- c("intra_caste_over_inter_caste_WW_contact_duration","QNurse_over_QForager_contact_duration")
# # names(variable_list) <- c("Within/Between-task contact","Q-N/Q-F contacts")
# # experiments=c("main_experiment")
# 
# data_path <- paste(disk_path,"/main_experiment/processed_data/network_properties_edge_weights_duration/random_vs_observed",sep="")
# pattern="network_properties"
# variable_list        <- c("modularity","clustering","task_assortativity","efficiency","degree_mean","density")
# names(variable_list) <- c("modularity","clustering","task assortativity","efficiency","mean degree","density")
# 
# 
# setwd(data_path)
# file_list <- list.files(pattern=pattern)
# 
# data_input <- NULL
# for (file in file_list){
#   data_input <- rbind(data_input,read.table(file,header=T,stringsAsFactors = F))}
# 
# 
# if (is.null(data_input)) {
#   ###read-in data###
#   starting_data <- NULL
#   for (experiment in experiments){
#     print(experiment)
#     ### data files
#     setwd(data_path)
#     file_list <- list.files(pattern=pattern)
#     temp <- NULL
#     for (file in file_list){
#       dat <- read.table(file,header=T,stringsAsFactors=F)
#       dat <- dat[,which(names(dat)%in%c("randy","colony","treatment","period","time_hours","time_of_day",variable_list))]
#       temp <- rbind(temp,dat)
#       rm(list=c("dat"))
#     }
#     temp <- temp[,order(names(temp))]
#     temp <- data.frame(experiment=experiment,temp,stringsAsFactors = F)
#     if (!is.null(starting_data)){
#       if (!all(names(starting_data)%in%names(temp))){
#         temp[names(starting_data)[which(!names(starting_data)%in%names(temp))]] <- NA
#         temp <- temp[,names(starting_data)]
#       }
#     }
#     
#     starting_data <- rbind(starting_data,temp)
#     rm(list=c("temp"))
#     
#   }
#   
#   ####modify period values to be simple and match what scatterplot function expects
#   starting_data["colony"] <- as.character(interaction(starting_data$experiment,starting_data$colony))
#   
# }else{starting_data <- data_input}
# 
# 
# 
# starting_data$variable <-  starting_data[,which(names(starting_data)== variable_list[5])]

#result <- compare_cohens_d(values=starting_data$variable, orig_perm_factor= as.factor(starting_data$randy), group_factor=starting_data$size)

result <- compare_cohens_d(data=starting_data)



# # Calculate differences within each group
# #library(dplyr)
# differences <- starting_data %>%
#   group_by(size) %>%
#   summarize(mean_diff = mean(variable[randy == "observed"] - variable[randy == "random"]),
#             std_diff = sd(variable[randy == "observed"] - variable[randy == "random"]))
# 
#   
# 
# # Calculate pooled standard deviation of differences
# pooled_std_diff <- sqrt(weighted.mean(differences$std_diff^2, differences$n - 1))
# 
# # Calculate Cohen's d for each group
# differences <- differences %>%
#   mutate(cohens_d = mean_diff / pooled_std_diff)
# 
# # Print the results
# print(differences)













####################### GRID RESULTS


#2 plots
dol_RandObs
# 6 plots
net_RandObs

sim_RandObs_EXP_seed









# USED 24H FILES.
# picked a random plot. to go through all networks, change the selection system of "network_file" by removing the need for a renaming of the file with "observed" or "random"but using the name of the directory and loop through all original list files (pick only 1 random).

##################
#### plot networks PRE - POST #######
#good vertex options for plotting: 9,11,21
vertexi <- 21
pdf(file=paste(figurefolder,"/Pre-Post_Networks_example.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=double_col,pointsize=pointsize_less_than_2row2col)
root_path <- paste(disk_path,"/main_experiment",sep="")######linux laptop
layout(t(matrix(c(1:6),nrow=2)),heights = c(5,0.2,5)) 
plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_observed","PostTreatment_observed"), size="small",vertexi=vertexi)# ### clean before next step
par(mar = c(0, 0, 0, 0))  # Set margins to zero for empty plots
plot.new();plot.new() #empty spacer
par(mar = c(2, 2, 2, 2))
plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_observed","PostTreatment_observed"), size="big",vertexi=vertexi)# ### clean before next step
Sys.sleep(2)
dev.off()

#### plot networks PRE - RANDOM #######
#good vertex options for plotting: 9,11,21
vertexi <- 21
pdf(file=paste(figurefolder,"/Pre-Ran_Networks_example.pdf",sep=""),family=text_font,font=text_font,bg="white",width=double_col,height=double_col,pointsize=pointsize_less_than_2row2col)
root_path <- paste(disk_path,"/main_experiment",sep="")######linux laptop
layout(t(matrix(c(1:6),nrow=2)),heights = c(5,0.7,5)) 
plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_observed","PreTreatment_random"), size="small",vertexi=vertexi)# ### clean before next step
par(mar = c(0, 0, 0, 0))  # Set margins to zero for empty plots
plot.new();plot.new() #empty spacer
par(mar = c(2, 2, 2, 2))
plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_observed","PreTreatment_random"), size="big",vertexi=vertexi)# ### clean before next step
Sys.sleep(2)
dev.off()


# layout(t(matrix(c(1:10),nrow=2)))
# for (vertexi in c(8:15)) {
#   plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_observed","PostTreatment_observed"), size="big",vertexi=vertexi)# ### clean before next step
#   #plot_network(case="topology_comparison",which_to_draw=c("PreTreatment_observed","PostTreatment_observed"), size="small",vertexi=vertexi)# ### clean before next step
#   
# }



######### Results report


# List objects matching the specified patterns from the environment
patterns <- c("ind_untreated", "coll_rescal_net","coll_no_rescal_sim", "ind_treated","Constitutive_organisation_ratios")
matching_objects <- ls(pattern = paste(patterns, collapse = "|"))

# Define a custom sorting function
custom_sort <- function(x) {
  order_val <- c("PRE", "seed", "treated", "untreated")
  order_idx <- sapply(order_val, function(val) grepl(val, x))
  order_val[which(order_idx)][1]
}

# Sort matching_objects using the custom sorting function
sorted_objects <- matching_objects[order(sapply(matching_objects, custom_sort))]

vec <- NULL



transform_vector <- function(input_vector) {
  output_vector <- gsub("_", " ", input_vector)  # Remove underscores
  output_vector <- gsub("ind", "individual", output_vector)
  output_vector <- gsub("coll", "collective", output_vector)
  output_vector <- gsub("sim", "simulation", output_vector)
  output_vector <- gsub("net", "(network)", output_vector)
  output_vector <- gsub("beh", "(behaviour)", output_vector)
  output_vector <- gsub("forag", "forager", output_vector)
  output_vector <- gsub("PRE", "period:PRE", output_vector)
  output_vector <- gsub("EXP seed", "| seed:treated", output_vector)
  output_vector <- gsub("RAN seed", "| seed:random", output_vector)
  output_vector <- gsub("FOR seed", "| seed:forager", output_vector)
  output_vector <- gsub("NUR seed", "| seed:nurse", output_vector)
  return(output_vector)
}


# Loop through each matching object
for (obj in sorted_objects) {
  # Get the stats_outcomes$formatted element from the current object
  formatted_element <- get(obj)[['stats_outcomes']][['formatted']]
  formatted_name <- transform_vector(obj)
  
  # Construct the output string and print using cat()
  output <- paste("####",formatted_name,"####", "\n", paste(formatted_element, collapse = '\n'), "\n\n")
  vec <- c(vec,cat(output))
}


vec1 <- NULL
sim_RandObs_results <- list(dol_RandObs$formatted_outcomes,
                            net_RandObs$formatted_outcomes)
sim_RandObs_results <- c(sim_RandObs_results,sim_RandObs_formatted_outputs)

# For the simulation outcomes, different organisation
for (obj in sim_RandObs_results) {
  # Construct the output string and print using cat()
  output <- paste(paste(obj, collapse = '\n'), "\n\n")
  vec1 <- c(vec1,cat(output))
}

















########## EXTRAS

#----------------------------------
# preferential groomers

root_path <- paste(disk_path,"/main_experiment_grooming",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <-        c("duration_grooming_received_min") #GROOMING  ouside is negligibile as only 58/6016 events happen outside , "prop_duration_grooming_received_outside_min","duration_grooming_received_min_zone2"
names(variable_list) <- c("duration grooming received (min)") # , "Prop. duration grooming received outside min","duration grooming received outside min"
#transf_variable_list <- c("log"                           )   ######"none", "sqrt" "log","power2"

dur_groom_rec_data <- read_data(data_path,which_individuals=c("forager","nurse"))

dur_groom_rec_data$task_group <- factor(dur_groom_rec_data$task_group , levels=task_group_order[which(task_group_order%in%dur_groom_rec_data$task_group )])


# Define a function to calculate the proportion of top outliers
calculate_top_outliers <- function(x) {
  Q3 <- quantile(x, 0.75)
  IQR <- IQR(x)
  sum(x > Q3 + 1.5*IQR) / length(x)
}

# Calculate the proportion of top outliers per task_group and treatment in each colony
outlier_data <- dur_groom_rec_data %>%
  group_by(colony, task_group, treatment) %>%
  summarise(proportion_outliers = calculate_top_outliers(duration_grooming_given_to_treated_min),
            .groups = "drop")

# Calculate the standard error
outlier_data <- outlier_data %>%
  group_by(task_group, treatment) %>%
  summarise(mean_proportion = mean(proportion_outliers),
            se = sd(proportion_outliers) / sqrt(n()),
            .groups = "drop")

# Plot the proportion of top outliers with standard errors
ggplot(outlier_data, aes(x = treatment, y = mean_proportion, fill=treatment)) +
  geom_errorbar(aes(ymin = mean_proportion - se, ymax = mean_proportion + se), width = 0.2) +
  geom_bar(stat = "identity") +
  facet_wrap(~task_group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(#title = "Proportion of Top Outliers per Treatment and Task Group",
    x = "",
    y = "Proportion of top 25% Groomers")+
  colFill_treatment +
  STYLE +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#EXTRA SOLITAIRE
Entropy_size$formatted_outcomes
