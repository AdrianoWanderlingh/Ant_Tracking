library(dplyr)
library(caret)
library(MLeval)
library(randomForest)
library(datasets)
library(data.table)
library(ggplot2)
library(viridis)
library(stringr)


#VARIABLES
#named vector for all the variables in your dataset:

var_names_ind <- c("transfmean_speed_BLpersec" = "GM Speed (bodylength/sec)",
                   "transfSD_speed_BLpersec" = "SD Speed (bodylength/sec)",
                   "transfmean_abs_accel_BLpersec2" = "GM Abs. acceleration (bodylength/sec²)",
                   "transfSD_abs_accel_BLpersec2" = "SD Abs. acceleration (bodylength/sec²)",
                   "transfmean_abs_jerk_BLPerSec3" = "GM Abs. jerk (bodylength/sec³)",
                   "transfSD_abs_jerk_BLPerSec3" = "SD Abs. jerk (bodylength/sec³)",
                   "transfmean_abs_TurnAngle" = "RMS Abs. turn angle",
                   "transfSD_abs_TurnAngle" = "SD Abs. turn angle",
                   "stDev_turnAngle" = "AD Turn angle",
                   "transfmean_abs_ang_Velocity_Movement" = "RMS Abs. angular velocity movement",
                   "transfSD_abs_ang_Velocity_Movement" = "SD Abs. angular velocity movement",
                   "stDev_ang_Velocity_Movement" = "AD Angular velocity movement",
                   "transfmean_abs_Body_Rotation" = "RMS Abs. body rotation",
                   "transfSD_abs_Body_Rotation" = "SD Abs. body rotation",
                   "stDev_Body_Rotation" = "AD Body rotation",
                   "transfmean_abs_ang_Velocity_Body" = "RMS Abs. angular velocity body",
                   "transfSD_abs_ang_Velocity_Body" = "SD Abs. angular velocity body",
                   "stDev_ang_Velocity_Body" = "AD Angular velocity body",
                   "StDev_Body_angle" = "AD Body Angle",
                   "StDev_Movement_angle" = "AD Movement Angle",
                   "transfmean_abs_Movement_Body_angle_diff" = "RMS Abs. movement body angle difference",
                   "transfSD_abs_Movement_Body_angle_diff" = "SD Abs. movement body angle difference",
                   "stDev_Movement_Body_angle_diff" = "AD Movement body angle difference",
                   "transfmean_abs_Movement_Body_Inclination_angle" = "RMS Abs. movement body inclination angle",
                   "transfSD_abs_Movement_Body_Inclination_angle" = "SD Abs. movement body inclination angle",
                   "stDev_Movement_Body_Inclination_angle" = "AD Movement body inclination angle",
                   "root_mean_square_deviation_BL2" = "RMSD",
                   "prop_time_undetected" = "Proportion Time undetected",
                   "sum_moved_distance_BL" = "moved distance (bodylength)",
                   "chull_area_BL2" = "Convex Hull (bodylength²)"
)

var_names_pair <- c("duration_sec" = "Duration (sec)",
                    "mean_strghtline_pair_dist_BL" = "AM Straight-line Distance (bodylength)",
                    "mean_pair_body_orient_diff_abs" = "AM Absolute difference in body orientation",
                    "mean_pair_body_inclination_diff" = "AM Difference in body inclination")


#In this vector, GM stands for Geometric Mean, SD for Standard Deviation, RMS for Root Mean Square, AD for Angular Deviation, RMSD for Root Mean Square Deviation, and AM for Arithmetic Mean. You can use this vector to rename the variables in your plots as shown in the previous example.

# Create a named vector for actor variables
var_names_ind_ACT <- paste0(var_names_ind, " ACT")
names(var_names_ind_ACT) <- paste0(names(var_names_ind), "_ACT")

# Create a named vector for receiver variables
var_names_ind_REC <- paste0(var_names_ind, " REC")
names(var_names_ind_REC) <- paste0(names(var_names_ind), "_REC")

# Combine the 3 vectors
var_names_all <- c(var_names_ind_ACT, var_names_ind_REC,var_names_pair)




#if (USER=="Supercomputer1") {
WORKDIR <- "/home/cf19810/Ant_Tracking/Data/PhD-Ant_Colonies_Tracking_Analysis/Automated_Behavioural_Inference"
#WORKDIR <- "/home/cf19810/Desktop/Final_Classifier_&_Scripts"
DATADIR <- paste(WORKDIR,"MachineLearning_outcomes_FINAL",sep="/")
#SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")
#EXPDATADIR <- "/media/bzniks/DISK4/ADRIANO/EXPERIMENT_DATA" #"/home/cf19810/Documents/Ants_behaviour_analysis/Data"
#}

###############################################################################
#### READ CHOSEN METHOD #######################################################
###############################################################################
chosen <- read.table(paste(DATADIR,"/quality_scores_Fbeta_1_priorityoverall_CHOSEN.txt",sep=""),header=T,stringsAsFactors = F)


# Get the list of files in the directory
file_list <- list.files(DATADIR, full.names = TRUE,pattern = "\\.txt$")

# Initialize an empty data frame to store the combined data
all_methods <- data.frame()

# Loop through each file
for (file in file_list) {
  # Check if the file name contains the text "Fbeta"
  if (!grepl("Fbeta", file)) {
    #print(file)
    # Read the file and bind it to the combined data frame
    data <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
    all_methods <- rbind(all_methods, data)
  }
}



#all_methods <- read.table(paste(DATADIR,"/quality_scores.txt",sep=""),header=T,stringsAsFactors = F)

###############################################################################
###### EXTRACT CHOSEN PARAMETERS FROM CHOSEN ##################################
subDir                      <- paste0("Loop_ID_",chosen[,"Loop_ID"]) 

##Load_classifier
rf_model        <- readRDS (file.path(DATADIR,subDir,"fits",paste(chosen[,"classifier"],".rds",sep="")))
rf_model_LIST       <- list(readRDS (file.path(DATADIR,subDir,"fits",paste(chosen[,"classifier"],".rds",sep=""))))
names(rf_model_LIST) <- chosen[,"classifier"]

# #plot model using randomforest package tools
# ###Error plot
# rf_model.legend <- if (is.null(rf_model$test$err.rate)) {colnames(rf_model$err.rate)} else {colnames(rf_model$test$err.rate)}
# plot(rf_model, log="y", main = "log of error rate" )
# legend("top", cex =1, legend=rf_model.legend, lty=c(1,2,3), col=c(1,2,3), horiz=T)

#####error plot in ggplot
# Get OOB data from plot and coerce to data.table
oobData <- as.data.table(plot(rf_model))
# Define trees as 1:ntree
oobData[, trees := .I]
# Cast to long format
oobData2 <- melt(oobData, id.vars = "trees")
setnames(oobData2, "value", "error")
# Plot using ggplot
ggplot(data = oobData2, aes(x = trees, y = error, color = variable)) + 
  theme_bw() +
  geom_line() +
  #scale_x_discrete(name = "Class", labels = c("OOB", "Non-grooming", "Grooming"))
  scale_color_manual(name="class",
                     labels=c("Overall error", "Non-allogrooming", "Allogrooming"),
                     values=c("black","brown1","forestgreen"))+
  labs(x="N of trees", y= "OOB error")+
  scale_x_continuous( expand = c(0, 0)) 

### Variables importance plot
#measure of the Mean decrease Gini
imp <- varImpPlot(rf_model) # let's save the varImp object

# this part just creates the data.frame for the plot part
imp <- as.data.frame(imp)
imp$varnames <- rownames(imp) # row names to column
rownames(imp) <- NULL
#assign category
imp$var_categ  <- ifelse(grepl('ACT', imp$varnames), 'Actor',
                         ifelse(grepl('REC', imp$varnames), 'Receiver', "Pair"))
#remove transformation
imp$varnames <- sub("\\..*", "", imp$varnames)

# Add a new column to the summary table
imp <- imp %>%
  mutate(new_names = var_names_all[match(varnames, names(var_names_all))])


# First, sort the data frame by MeanDecreaseGini in descending order
imp_sorted <- imp[order(-imp$MeanDecreaseGini),]

# Then, select the top 10 rows
imp_top10 <- head(imp_sorted, 10)

# Now, plot the data
ggplot(imp_top10, aes(x=reorder(new_names, MeanDecreaseGini), y=MeanDecreaseGini, color=as.factor(var_categ))) + 
  geom_point() +
  geom_segment(aes(x=new_names,xend=new_names,y=0,yend=MeanDecreaseGini)) +
  scale_color_discrete(name="Individual measured") +
  theme_bw() +
  theme(legend.position="top") +
  ylab("Mean Decrease Gini") +
  xlab("top 10  most influential variables") +
  coord_flip()



varUsed(rf_model, by.tree=FALSE, count=TRUE)



##### precision vs sensitivity
#ggplot(all_methods$precision_test,all_methods$sensitivity_test)


min_x <- min(all_methods$precision_test * 100)
min_y <- min(all_methods$sensitivity_test * 100)

# Create a new column with modified classifier names
# Create a new column with modified classifier names
# Create new columns for classifier name and class balance technique
all_methods$CLASSIFNAME <- sapply(strsplit(as.character(all_methods$classifier), "_"), function(x) x[2])
all_methods$CLASSBALANCE <- sapply(strsplit(as.character(all_methods$classifier), "_"), function(x) {
  class_balance <- paste(gsub("_", " ", x[3:length(x)]), collapse = " ")
  class_balance <- gsub("\\bsmote\\b", "SMOTE", class_balance, ignore.case = TRUE)
  class_balance <- gsub("\\bros\\b", "ROS", class_balance, ignore.case = TRUE)
  class_balance <- gsub("\\brus\\b", "RUS", class_balance, ignore.case = TRUE)
  class_balance <- gsub("\\bsbc\\b", "SBC", class_balance, ignore.case = TRUE)
  class_balance
})

# Now use these new columns for the colour and shape aesthetics in your ggplot call
ggplot(all_methods, aes(x = precision_test * 100, y = sensitivity_test * 100)) + 
  geom_point(aes(shape = CLASSIFNAME, colour = CLASSBALANCE)) +
  scale_colour_viridis_d(option = "plasma", name = "Class balancing technique") + # Set legend title here and number of columns
  scale_shape_discrete(name = "Classifier") +
  geom_segment(data = chosen,
               aes(x = precision_test * 100, y = min_y, xend = precision_test * 100, yend = sensitivity_test * 100),
               linetype = "dashed", color = "green3") +
  geom_segment(data = chosen,
               aes(x = min_x, y = sensitivity_test * 100, xend = precision_test * 100, yend = sensitivity_test * 100),
               linetype = "dashed", color = "green3") +
  theme_bw() +
  geom_point(data = chosen,
             aes(x = precision_test * 100, y = sensitivity_test * 100), 
             fill = NA, 
             color = 'green3',
             size = 3,
             shape = 21) + # Shape 21 is a circle with a border +# Shape 21 is a circle with a border
  geom_text(data = chosen,
            aes(x = precision_test * 100, y = sensitivity_test * 100, 
                label = paste0(round(chosen$precision_test * 100), ", ", round(chosen$sensitivity_test * 100))),
            hjust = 0.5, vjust = -3, color = "green3") +
  labs(caption = paste0("N=", nrow(all_methods)),
       x ="Precision (%)",
       y ="Sensitivity (%)") +
  coord_fixed(ratio = 1)  # Set aspect ratio to 1:1
#coord_cartesian(expand = T)




# other option, split by class balancing technique


# Create a new data frame for the overall distribution
overall_distribution <- data.frame(
  precision_test = all_methods$precision_test,
  sensitivity_test = all_methods$sensitivity_test
)

# Now use these new columns for the colour and shape aesthetics in your ggplot call
ggplot() + 
  geom_point(data = overall_distribution, aes(x = precision_test * 100, y = sensitivity_test * 100), colour = "gray", alpha = 0.5, shape = 16) + # Add gray points at the back
  geom_point(data = all_methods, aes(x = precision_test * 100, y = sensitivity_test * 100, colour = CLASSIFNAME), alpha = 0.5, shape = 16) + # Add transparency and use filled circle
  scale_colour_viridis_d(option = "plasma", name = "Classifier") + # Set legend title here and number of columns
  facet_wrap(~CLASSBALANCE) +
  geom_segment(data = chosen,
               aes(x = precision_test * 100, y = min_y, xend = precision_test * 100, yend = sensitivity_test * 100),
               linetype = "dashed", color = "red") +
  geom_segment(data = chosen,
               aes(x = min_x, y = sensitivity_test * 100, xend = precision_test * 100, yend = sensitivity_test * 100),
               linetype = "dashed", color = "red") +
  theme_bw() +
  # Set the number of rows and columns in the legend
  theme(legend.position = c(0.7, 0.1), 
        legend.box = "horizontal", 
        legend.key.size = unit(1.5, "lines")) +
  guides(colour = guide_legend(nrow = 1, ncol = 2, byrow = TRUE,override.aes = list(size = 5))) +
  geom_point(data = chosen,
             aes(x = precision_test * 100, y = sensitivity_test * 100), 
             fill = NA, 
             color = 'red',
             size = 3,
             shape = 21) +
  labs(caption = paste0("N=", nrow(all_methods)),
       x ="Precision (%)",
       y ="Sensitivity (%)") +
  coord_fixed(ratio = 1)



#-------------------------------------

#Assuming you already have a vector of probabilities (called probs) computed with your model and the true class labels are in your data frame as df$label (0 and 1) this code should work:
# https://stats.stackexchange.com/questions/10501/calculating-aupr-in-r
fg <- probs[df$label == 1]
bg <- probs[df$label == 0]

# ROC Curve    
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(roc)

# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr)

#PS: The only disconcerting thing is you use scores.class0 = fg when fg is computed for label 1 and not 0.

#-------------------------------------



