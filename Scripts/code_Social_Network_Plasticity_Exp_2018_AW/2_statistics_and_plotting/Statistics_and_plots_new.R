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
exposure_order  <- c("control","pathogen")
size_order      <- c("small","big")
period_order    <- c("pre","post")
task_group_order <- c("queen","nurse","forager","untreated","treated")

####source programs
source(paste(source_path,"/libraries.R",sep=""))
source(paste(source_path,"/plotting_parameters.R",sep=""))
source(paste(source_path,"/functions_new.R",sep=""))
source(paste(source_path,"/analysis_parameters.R",sep=""))


###################################################################
###I - Example collective analysis - without rescaling ############
###################################################################
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="colony_data.txt"
variable_list <- c("modularity","clustering","task_assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
names(variable_list) <- c("modularity","clustering","assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
transf_variable_list <- c("log"      ,"none"      ,"none"         ,"log"      ,"log"       ,"log"          ,"none"   ,"log")   ######"none", "sqrt" "log","power2"

###1. read data
setwd(data_path)
file_list <- list.files(pattern=pattern)
print(file_list)
data <- NULL
for (file in file_list){
  data <- rbind(data,read.table(file,header=T,stringsAsFactors = F))  
}
##remove any duplicated line
data <- data[which(!duplicated(data)),]

##2. Extract exposure and size from treatment column
data$exposure <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[1]  ))
data$size     <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))


###2. Loop over variables
data_ori <- data
stats_outcomes    <- NULL
post_hoc_outcomes <- NULL

for (i in 1:length(variable_list)){
  # for (i in c(1)){
  print(paste0("######## ",variable_list[i]," ########"))
  data <- data_ori
  
  ###create a variable column
  data$untransformed_variable <- data[,variable_list[i]]
  
  ###transform variable
  if (transf_variable_list[i]=="log"){
    print("Logging variable...")
    data[!is.na(data$untransformed_variable),"variable"] <- log_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (grepl("power",transf_variable_list[i])){
    data[!is.na(data$untransformed_variable),"variable"]  <- (data[!is.na(data$untransformed_variable),"untransformed_variable"] )^as.numeric(gsub("power","",transf_variable_list[i]))
  }else if (transf_variable_list[i]=="sqrt"){
    data[!is.na(data$untransformed_variable),"variable"]  <- sqrt_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (transf_variable_list[i]=="none"){
    data$variable <- data$untransformed_variable
  }
  
  ###statistics
  ###make sure treatment, exposure, size and period are factors with levels in the desired order
  data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
  data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
  data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
  data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
  
  ###fit model - using treatment
  model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) ,data=data)
  anov  <- anova(model)
  p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
  
  stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),pval=anov["period:treatment","Pr(>F)"],stringsAsFactors = F))
  ###check interaction
  if (p_interaction_treatment<=0.05){
    test_norm(residuals(model))
    # print("Significant interaction between 4-way treatment and period.")
    # print(anov)
    
    contrast_matrix <- rbind(
      "Delta_control.small minus Delta_control.big"=c(0,0,0,0,0,-1,0,0)
      ,
      "Delta_control.small minus Delta_pathogen.small"=c(0,0,0,0,0,0,-1,0)
      ,
      "Delta_control.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,0,-1)
      ,
      "Delta_control.big minus Delta_pathogen.small"=c(0,0,0,0,0,1,-1,0)
      ,
      "Delta_control.big minus Delta_pathogen.big"=c(0,0,0,0,0,1,0,-1)
      ,
      "Delta_pathogen.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,1,-1)
    )
    
    posthoc_groups_treatments        <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix,which_levels="treatment_order",dataset=data))
    names(posthoc_groups_treatments) <- variable_list[i]
    post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_treatments)
    
    
  }else{
    # print( "Interaction between 4-way treatment and period is not significant:")
    # print(anov["period:treatment",])
    
    
    ###if interaction with treatment is not significant, repeat analysis, but separating size and exposure
    ###fit model 
    model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) ,data=data)
    anov  <- anova(model)
    p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
    p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
    
    if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
      stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
      stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
      
      
      test_norm(residuals(model))
      # print("Significant interaction between exposure (pathogen vs solvent) and period, AND significant interaction beween size (small vs big) and period.")
      # print(anov)
      
      contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,0,-1,0))
      contrast_matrix_size     <- rbind("Delta_small minus Delta_big"=c(0,0,0,0,0,-1))
      
      posthoc_groups_exposure <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data)
      posthoc_groups_size     <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size    ,which_levels="size_order",dataset=data)
      
    }else{
      if (p_interaction_size>p_interaction_exposure){ 
        ###if interaction with size is not significant and greater than interaction with exposure, repeat analysis, but exposure only
        # print( "Interaction between size and period is not significant:")
        # print(anov["period:size",])
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
        
        ###fit model 
        model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) ,data=data)
        anov  <- anova(model)
        test_norm(residuals(model))
        for (rowi in 1:nrow(anov)){
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
        }
        
        p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
        if (p_interaction_exposure>0.05){
          # print( "Interaction between exposure and period is not significant:")
          # print(anov["period:exposure",])
          
        }else{
          contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,-1))
          
          posthoc_groups_exposure <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data))
          names(posthoc_groups_exposure) <- variable_list[i]
          post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_exposure)
          
          # print( "Significant interaction between exposure and period:")
          # print(anov)
        }
        
      }else{
        ###if interaction with exposure is not significant and greater than interaction with size, repeat analysis, but size only
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
        
        # print( "Interaction between exposure and period is not significant:")
        # print(anov["period:exposure",])
        # 
        ###fit model 
        model <- lmer(   variable ~ period*size + (1|time_of_day) + (1|colony) ,data=data)
        anov  <- anova(model)
        test_norm(residuals(model))
        for (rowi in 1:nrow(anov)){
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
        }
        
        p_interaction_size <- anov["Pr(>F)"]["period:size","Pr(>F)"]
        if (p_interaction_size>0.05){
          # print( "Interaction between size and period is not significant:")
          # print(anov)
          
        }else{
          # print( "Significant interaction between size and period:")
          # print(anov)
          contrast_matrix_size <- rbind("Delta_small minus Delta_big"=c(0,0,0,-1))
          
          posthoc_groups_size <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size,which_levels="size_order",dataset=data))
          names(posthoc_groups_size) <- variable_list[i]
          post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_size)
          
          
        }
        
      }
    }
    
  }
  rm(list=ls()[which(grepl("p_interaction",ls()))])
  rm(list=ls()[which(grepl("posthoc_groups_",ls()))])
  
  #plot
  barplot_delta_collective <- barplot_delta(dataset=data,predictor="treatment",type="collective",collective=T,plot_untransformed=T,diff_type="absolute_difference") #form_stat=NULL,
  print(barplot_delta_collective)
  #ADD SAVING OF THE PLOT on single pdf file as done in plot_grooming file
  }
rownames(stats_outcomes) <- NULL


# Reshape the data
stats_outcomes_reshaped <- reshape(stats_outcomes[,which(names(stats_outcomes)!="df")], idvar = "predictor", timevar = "variable", direction = "wide")
#colnames(reshaped_data)[-1] <- gsub("pval.", "", colnames(reshaped_data)[-1])



#######################################################################################################
###II - Example collective analysis - with rescaling by pre-exposure mean for each size##############
#######################################################################################################
root_path <- paste(disk_path,"/main_experiment/processed_data",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming/processed_data",sep="")
data_path <- paste(root_path,"/network_properties/pre_vs_post_treatment/all_workers",sep="")
pattern="colony_data.txt"
variable_list <- c("modularity","clustering","task_assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
names(variable_list) <- c("modularity","clustering","assortativity","efficiency","degree_mean","degree_maximum","density","diameter")
transf_variable_list <- c("log"       ,"sqrt"      ,"none"          ,"log"       ,"log"        ,"log"           ,"log"   ,"power0.01")   ######"none", "sqrt" "log","power2"
# NOTE: task_assortativity is hard to normalise because of the bimodal distribution (small vs big cols)

###1. read data
setwd(data_path)
file_list <- list.files(pattern=pattern)
print(file_list)
data <- NULL
for (file in file_list){
  data <- rbind(data,read.table(file,header=T,stringsAsFactors = F))  
}
##remove any duplicated line
data <- data[which(!duplicated(data)),]

##2. Extract exposure and size from treatment column
data$exposure <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[1]  ))
data$size     <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))


###2. Loop over variables
data_ori <- data
stats_outcomes    <- NULL
post_hoc_outcomes <- NULL

for (i in 1:length(variable_list)){
  # for (i in c(1)){
  print(variable_list[i])
  data <- data_ori
  
  ###create a variable column
  data$untransformed_variable <- data[,variable_list[i]]

  ###rescale, for each treatment, using "pre" period as reference
  means <- aggregate(untransformed_variable~size,data=data[which(data$period=="pre"),],FUN=mean)
  names(means)[2] <- "pre_mean"
  data <- merge(data,means)
  data$untransformed_variable <- data$untransformed_variable/data$pre_mean
  
    
  ###transform variable
  if (transf_variable_list[i]=="log"){
    print("Logging variable...")
    data[!is.na(data$untransformed_variable),"variable"] <- log_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (grepl("power",transf_variable_list[i])){
    data[!is.na(data$untransformed_variable),"variable"]  <- (data[!is.na(data$untransformed_variable),"untransformed_variable"] )^as.numeric(gsub("power","",transf_variable_list[i]))
  }else if (transf_variable_list[i]=="sqrt"){
    data[!is.na(data$untransformed_variable),"variable"]  <- sqrt_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (transf_variable_list[i]=="none"){
    data$variable <- data$untransformed_variable
  }
  
  
  ###statistics
  ###make sure treatment, exposure, size and period are factors with levels in the desired order
  data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
  data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
  data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
  data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
  
  ###fit model - using treatment
  model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) ,data=data)
  anov  <- anova(model)
  p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
  
  stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),pval=anov["period:treatment","Pr(>F)"],stringsAsFactors = F))
  ###check interaction
  if (p_interaction_treatment<=0.05){
    test_norm(residuals(model))
    # print("Significant interaction between 4-way treatment and period.")
    # print(anov)
    
    contrast_matrix <- rbind(
      "Delta_control.small minus Delta_control.big"=c(0,0,0,0,0,-1,0,0)
      ,
      "Delta_control.small minus Delta_pathogen.small"=c(0,0,0,0,0,0,-1,0)
      ,
      "Delta_control.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,0,-1)
      ,
      "Delta_control.big minus Delta_pathogen.small"=c(0,0,0,0,0,1,-1,0)
      ,
      "Delta_control.big minus Delta_pathogen.big"=c(0,0,0,0,0,1,0,-1)
      ,
      "Delta_pathogen.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,1,-1)
    )
    
    posthoc_groups_treatments        <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix,which_levels="treatment_order",dataset=data))
    names(posthoc_groups_treatments) <- variable_list[i]
    post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_treatments)
    
    
  }else{
    # print( "Interaction between 4-way treatment and period is not significant:")
    # print(anov["period:treatment",])
    
    
    ###if interaction with treatment is not significant, repeat analysis, but separating size and exposure
    ###fit model 
    model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) ,data=data)
    anov  <- anova(model)
    p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
    p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
    
    if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
      stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
      stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
      
      
      test_norm(residuals(model))
      # print("Significant interaction between exposure (pathogen vs solvent) and period, AND significant interaction beween size (small vs big) and period.")
      # print(anov)
      
      contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,0,-1,0))
      contrast_matrix_size     <- rbind("Delta_small minus Delta_big"=c(0,0,0,0,0,-1))
      
      posthoc_groups_exposure <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data)
      posthoc_groups_size     <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size    ,which_levels="size_order",dataset=data)
      
    }else{
      if (p_interaction_size>p_interaction_exposure){ 
        ###if interaction with size is not significant and greater than interaction with exposure, repeat analysis, but exposure only
        # print( "Interaction between size and period is not significant:")
        # print(anov["period:size",])
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
        
        ###fit model 
        model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) ,data=data)
        anov  <- anova(model)
        test_norm(residuals(model))
        for (rowi in 1:nrow(anov)){
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
        }
        
        p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
        if (p_interaction_exposure>0.05){
          # print( "Interaction between exposure and period is not significant:")
          # print(anov["period:exposure",])
          
        }else{
          contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,-1))
          
          posthoc_groups_exposure <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data))
          names(posthoc_groups_exposure) <- variable_list[i]
          post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_exposure)
          
          # print( "Significant interaction between exposure and period:")
          # print(anov)
        }
        
      }else{
        ###if interaction with exposure is not significant and greater than interaction with size, repeat analysis, but size only
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
        
        # print( "Interaction between exposure and period is not significant:")
        # print(anov["period:exposure",])
        # 
        ###fit model 
        model <- lmer(   variable ~ period*size + (1|time_of_day) + (1|colony) ,data=data)
        anov  <- anova(model)
        test_norm(residuals(model))
        for (rowi in 1:nrow(anov)){
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
        }
        
        p_interaction_size <- anov["Pr(>F)"]["period:size","Pr(>F)"]
        if (p_interaction_size>0.05){
          # print( "Interaction between size and period is not significant:")
          # print(anov)
          
        }else{
          # print( "Significant interaction between size and period:")
          # print(anov)
          contrast_matrix_size <- rbind("Delta_small minus Delta_big"=c(0,0,0,-1))
          
          posthoc_groups_size <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size,which_levels="size_order",dataset=data))
          names(posthoc_groups_size) <- variable_list[i]
          post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_size)
          
          
        }
        
      }
    }
    
  }
  rm(list=ls()[which(grepl("p_interaction",ls()))])
  rm(list=ls()[which(grepl("posthoc_groups_",ls()))])
  
  #plot
  barplot_delta_collective <- barplot_delta(dataset=data,predictor="treatment",type="collective",collective=T,plot_untransformed=T,diff_type="") #form_stat=NULL,
  print(barplot_delta_collective)
}
rownames(stats_outcomes) <- NULL

# Reshape the data
stats_outcomes_reshaped <- reshape(stats_outcomes[,which(names(stats_outcomes)!="df")], idvar = "predictor", timevar = "variable", direction = "wide")

###################################################################################################################################
###III - Example individual analysis - for ONE type of individuals (e.g. treated workers only, or queens only)#####################
###################################################################################################################################
root_path <- paste(disk_path,"/main_experiment",sep="") # root_path <- paste(disk_path,"/main_experiment_grooming",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <- c("prop_time_outside","proportion_time_active","duration_of_contact_with_treated_min","inter_caste_contact_duration")
names(variable_list) <- c("prop_time_outside","proportion_time_active","duration_of_contact_with_treated_min","inter_caste_contact_duration")
transf_variable_list <- c("power0.01"        ,"none"                  ,"sqrt"                                ,"none"      )   ######"none", "sqrt" "log","power2"

which_individuals <- "treated" ## "treated","queen","nurse","forager"

###1. read data
setwd(data_path)
file_list <- list.files(pattern=pattern)
print(file_list)
data <- NULL
for (file in file_list){
  data <- rbind(data,read.table(file,header=T,stringsAsFactors = F))  
}
##remove any duplicated line
data <- data[which(!duplicated(data)),]

##2a. Extract exposure and size from treatment column
data$exposure <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[1]  ))
data$size     <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))

##2b. add information on task group
task_groups <- read.table(paste(root_path,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
task_groups <- task_groups[which(!duplicated(task_groups)),]

data        <- merge(data,task_groups[c("colony","tag","task_group")],all.x=T,all.y=F)
data[which(data$status=="treated"),"task_group"] <- "treated"

##2c. keep only target individuals
data <- data[which(data$task_group%in%which_individuals),]

##2d.  ###add a unique antid column
data <- within(data,antID <- paste(colony,tag,sep="_"))

###2. Loop over variables
data_ori <- data
stats_outcomes    <- NULL
post_hoc_outcomes <- NULL

for (i in 1:length(variable_list)){
  print(variable_list[i])
  data <- data_ori
  
  ###create a variable column
  data$untransformed_variable <- data[,variable_list[i]]
  
  ###transform variable
  if (transf_variable_list[i]=="log"){
    print("Logging variable...")
    data[!is.na(data$untransformed_variable),"variable"] <- log_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (grepl("power",transf_variable_list[i])){
    data[!is.na(data$untransformed_variable),"variable"]  <- (data[!is.na(data$untransformed_variable),"untransformed_variable"] )^as.numeric(gsub("power","",transf_variable_list[i]))
  }else if (transf_variable_list[i]=="sqrt"){
    data[!is.na(data$untransformed_variable),"variable"]  <- sqrt_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (transf_variable_list[i]=="none"){
    data$variable <- data$untransformed_variable
  }
  
  ###statistics
  ###make sure treatment, exposure, size and period are factors with levels in the desired order
  data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
  data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
  data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
  data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
  data$antID     <- factor( data$antID )
  
  ###fit model - using treatment
  model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) + (1|antID) ,data=data)
  anov  <- anova(model)
  p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
  
  stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),pval=anov["period:treatment","Pr(>F)"],stringsAsFactors = F))
  ###check interaction
  if (p_interaction_treatment<=0.05){
    test_norm(residuals(model))
    # print("Significant interaction between 4-way treatment and period.")
    # print(anov)
    
    contrast_matrix <- rbind(
      "Delta_control.small minus Delta_control.big"=c(0,0,0,0,0,-1,0,0)
      ,
      "Delta_control.small minus Delta_pathogen.small"=c(0,0,0,0,0,0,-1,0)
      ,
      "Delta_control.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,0,-1)
      ,
      "Delta_control.big minus Delta_pathogen.small"=c(0,0,0,0,0,1,-1,0)
      ,
      "Delta_control.big minus Delta_pathogen.big"=c(0,0,0,0,0,1,0,-1)
      ,
      "Delta_pathogen.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,1,-1)
    )
    
    posthoc_groups_treatments        <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix,which_levels="treatment_order",dataset=data))
    names(posthoc_groups_treatments) <- variable_list[i]
    post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_treatments)
    
    
  }else{
    # print( "Interaction between 4-way treatment and period is not significant:")
    # print(anov["period:treatment",])
    
    
    ###if interaction with treatment is not significant, repeat analysis, but separating size and exposure
    ###fit model 
    model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) + (1|antID),data=data)
    anov  <- anova(model)
    p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
    p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
    
    if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
      stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
      stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
      
      
      test_norm(residuals(model))
      # print("Significant interaction between exposure (pathogen vs solvent) and period, AND significant interaction beween size (small vs big) and period.")
      # print(anov)
      
      contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,0,-1,0))
      contrast_matrix_size     <- rbind("Delta_small minus Delta_big"=c(0,0,0,0,0,-1))
      
      posthoc_groups_exposure <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data)
      posthoc_groups_size     <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size    ,which_levels="size_order",dataset=data)
      
    }else{
      if (p_interaction_size>p_interaction_exposure){ 
        ###if interaction with size is not significant and greater than interaction with exposure, repeat analysis, but exposure only
        # print( "Interaction between size and period is not significant:")
        # print(anov["period:size",])
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
        
        ###fit model 
        model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) + (1|antID),data=data)
        anov  <- anova(model)
        test_norm(residuals(model))
        for (rowi in 1:nrow(anov)){
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
        }
        
        p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
        if (p_interaction_exposure>0.05){
          # print( "Interaction between exposure and period is not significant:")
          # print(anov["period:exposure",])
          
        }else{
          contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,-1))
          
          posthoc_groups_exposure <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data))
          names(posthoc_groups_exposure) <- variable_list[i]
          post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_exposure)
          
          # print( "Significant interaction between exposure and period:")
          # print(anov)
        }
        
      }else{
        ###if interaction with exposure is not significant and greater than interaction with size, repeat analysis, but size only
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
        
        # print( "Interaction between exposure and period is not significant:")
        # print(anov["period:exposure",])
        # 
        ###fit model 
        model <- lmer(   variable ~ period*size + (1|time_of_day) + (1|colony) + (1|antID),data=data)
        anov  <- anova(model)
        test_norm(residuals(model))
        for (rowi in 1:nrow(anov)){
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
        }
        
        p_interaction_size <- anov["Pr(>F)"]["period:size","Pr(>F)"]
        if (p_interaction_size>0.05){
          # print( "Interaction between size and period is not significant:")
          # print(anov)
          
        }else{
          # print( "Significant interaction between size and period:")
          # print(anov)
          contrast_matrix_size <- rbind("Delta_small minus Delta_big"=c(0,0,0,-1))
          
          posthoc_groups_size <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size,which_levels="size_order",dataset=data))
          names(posthoc_groups_size) <- variable_list[i]
          post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_size)
          
          
        }
        
      }
    }
    
  }
  rm(list=ls()[which(grepl("p_interaction",ls()))])
  rm(list=ls()[which(grepl("posthoc_groups_",ls()))])
  
  #plot
  barplot_delta_collective <- barplot_delta(dataset=data,predictor="treatment",type="individual",collective=F,plot_untransformed=T,diff_type="absolute_difference") #form_stat=NULL,
  print(barplot_delta_collective)
  #ADD SAVING OF THE PLOT on single pdf file as done in plot_grooming file
  
}
rownames(stats_outcomes) <- NULL


###################################################################################################################################
###IV - Example individual analysis - for TWO types of individuals (e.g. nurses and foragers; or "untreated" and "treated") #######
###################################################################################################################################
root_path <- paste(disk_path,"/main_experiment",sep="")
data_path=paste(root_path,"/processed_data/individual_behaviour/pre_vs_post_treatment",sep="")
pattern="individual_behavioural_data"
variable_list <- c("prop_time_outside")
names(variable_list) <- c("prop_time_outside")
transf_variable_list <- c("power0.01")   ######"none", "sqrt" "log","power2"

which_individuals <- c("nurse","forager") ## "treated","queen","nurse","forager"

###1. read data
setwd(data_path)
file_list <- list.files(pattern=pattern)
print(file_list)
data <- NULL
for (file in file_list){
  data <- rbind(data,read.table(file,header=T,stringsAsFactors = F))  
}
##remove any duplicated line
data <- data[which(!duplicated(data)),]

##2a. Extract exposure and size from treatment column
data$exposure <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[1]  ))
data$size     <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))

##2b. add information on task group
task_groups <- read.table(paste(root_path,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
task_groups <- task_groups[which(!duplicated(task_groups)),]

data        <- merge(data,task_groups[c("colony","tag","task_group")],all.x=T,all.y=F)
data[which(data$status=="treated"),"task_group"] <- "treated"

##2c. keep only target individuals
data <- data[which(data$task_group%in%which_individuals),]

##2d.  ###add a unique antid column
data <- within(data,antID <- paste(colony,tag,sep="_"))

###2. Loop over variables
data_ori <- data
stats_outcomes    <- NULL
post_hoc_outcomes <- NULL

for (i in 1:length(variable_list)){
  print(variable_list[i])
  data <- data_ori
  
  ###create a variable column
  data$untransformed_variable <- data[,variable_list[i]]
  
  ###transform variable
  if (transf_variable_list[i]=="log"){
    print("Logging variable...")
    data[!is.na(data$untransformed_variable),"variable"] <- log_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (grepl("power",transf_variable_list[i])){
    data[!is.na(data$untransformed_variable),"variable"]  <- (data[!is.na(data$untransformed_variable),"untransformed_variable"] )^as.numeric(gsub("power","",transf_variable_list[i]))
  }else if (transf_variable_list[i]=="sqrt"){
    data[!is.na(data$untransformed_variable),"variable"]  <- sqrt_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
  }else if (transf_variable_list[i]=="none"){
    data$variable <- data$untransformed_variable
  }
  
  ###statistics
  ###make sure treatment, exposure, size, period and task_group are factors with levels in the desired order
  data$treatment  <- factor(data$treatment , levels=treatment_order [which(treatment_order %in%data$treatment )])
  data$size       <- factor(data$size      , levels=size_order      [which(size_order      %in%data$size )])
  data$exposure   <- factor(data$exposure  , levels=exposure_order  [which(exposure_order  %in%data$exposure )])
  data$period     <- factor(data$period    , levels=period_order    [which(period_order    %in%data$period )])
  data$antID      <- factor(data$antID )
  data$task_group <- factor(data$task_group, levels=task_group_order[which(task_group_order%in%data$task_group )])
  
  ###fit model - using treatment and task_group in 3-way interaction
  model <- lmer(   variable ~ period*treatment*task_group + (1|time_of_day) + (1|colony) + (1|antID) ,data=data)
  anov  <- anova(model)
  p_interaction_treatment_group <- anov["period:treatment:task_group","Pr(>F)"]
  stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment:task_group",df=paste(round(anov["period:treatment:task_group","NumDF"]),round(anov["period:treatment:task_group","DenDF"]),sep=","),pval=anov["period:treatment:task_group","Pr(>F)"],stringsAsFactors = F))
  
  if (p_interaction_treatment_group<=0.05){
    test_norm(residuals(model))
    
    level_names <- expand.grid(levels(data$treatment),levels(data$task_group))
    level_names <- within(level_names,Var3<-paste(Var1,Var2,sep="."))$Var3
    names(level_names) <- level_names
    
    ####contrast matrix will depend on how many groups are included here!
    if (length(levels(data$task_group))==2){
      contrast_matrix <- rbind(
      "1 minus 2"=c(0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0)
      ,
      "1 minus 3"=c(0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0)
      ,
      "1 minus 4"=c(0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0)
      ,
      "1 minus 5"=c(0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0)
      ,
      "1 minus 6"=c(0,0,0,0,0,0,-1,0,0,-1,0,0,0,-1,0,0)
      ,
      "1 minus 7"=c(0,0,0,0,0,0,0,-1,0,-1,0,0,0,0,-1,0)
      ,
      "1 minus 8"=c(0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0,-1)
      ,
      "2 minus 3"=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0)
      ,
      "2 minus 4"=c(0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0)
      ,
      "2 minus 5"=c(0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0)
      ,
      "2 minus 6"=c(0,0,0,0,0,0,0,0,0,-1,0,0,0,-1,0,0)
      ,
      "2 minus 7"=c(0,0,0,0,0,0,1,-1,0,-1,0,0,0,0,-1,0)
      ,
      "2 minus 8"=c(0,0,0,0,0,0,1,0,-1,-1,0,0,0,0,0,-1)
      ,
      "3 minus 4"=c(0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0)
      ,
      "3 minus 5"=c(0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0)
      ,
      "3 minus 6"=c(0,0,0,0,0,0,-1,1,0,-1,0,0,0,-1,0,0)
      ,
      "3 minus 7"=c(0,0,0,0,0,0,0,0,0,-1,0,0,0,0,-1,0)
      ,
      "3 minus 8"=c(0,0,0,0,0,0,0,1,-1,-1,0,0,0,0,0,-1)
      ,
      "4 minus 5"=c(0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0)
      ,
      "4 minus 6"=c(0,0,0,0,0,0,-1,0,1,-1,0,0,0,-1,0,0)
      ,
      "4 minus 7"=c(0,0,0,0,0,0,0,-1,1,-1,0,0,0,0,-1,0)
      ,
      "4 minus 8"=c(0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1)
      ,
      "5 minus 6"=c(0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0)
      ,
      "5 minus 7"=c(0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0)
      ,
      "5 minus 8"=c(0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1)
      ,
      "6 minus 7"=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,1,-1,0)
      ,
      "6 minus 8"=c(0,0,0,0,0,0,1,0,-1,0,0,0,0,1,0,-1)
      ,
      "7 minus 8"=c(0,0,0,0,0,0,0,1,-1,0,0,0,0,0,1,-1)
      )
      for (i in 1:length(level_names)){
        rownames(contrast_matrix) <- gsub(i,level_names[i],rownames(contrast_matrix))
      }
      
    }else{print("TOO MANY TASK GROUP LEVELS _ CANNOT COPE")}
    
     
    posthoc_groups_treatment_task_groups        <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix,which_levels="level_names",dataset=data))
    names(posthoc_groups_treatment_task_groups) <- variable_list[i]
    post_hoc_outcomes                           <- c(post_hoc_outcomes,posthoc_groups_treatment_task_groups)
    
  }else{
    
    ###fit model - using treatment
    model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) + (1|antID) ,data=data)
    anov  <- anova(model)
    p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
    
    stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),pval=anov["period:treatment","Pr(>F)"],stringsAsFactors = F))
    ###check interaction
    if (p_interaction_treatment<=0.05){
      test_norm(residuals(model))
      # print("Significant interaction between 4-way treatment and period.")
      # print(anov)
      
      contrast_matrix <- rbind(
        "Delta_control.small minus Delta_control.big"=c(0,0,0,0,0,-1,0,0)
        ,
        "Delta_control.small minus Delta_pathogen.small"=c(0,0,0,0,0,0,-1,0)
        ,
        "Delta_control.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,0,-1)
        ,
        "Delta_control.big minus Delta_pathogen.small"=c(0,0,0,0,0,1,-1,0)
        ,
        "Delta_control.big minus Delta_pathogen.big"=c(0,0,0,0,0,1,0,-1)
        ,
        "Delta_pathogen.small minus Delta_pathogen.big"=c(0,0,0,0,0,0,1,-1)
      )
      
      posthoc_groups_treatments        <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix,which_levels="treatment_order",dataset=data))
      names(posthoc_groups_treatments) <- variable_list[i]
      post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_treatments)
      
      
    }else{
      # print( "Interaction between 4-way treatment and period is not significant:")
      # print(anov["period:treatment",])
      
      
      ###if interaction with treatment is not significant, repeat analysis, but separating size and exposure
      ###fit model 
      model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) + (1|antID),data=data)
      anov  <- anova(model)
      p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
      p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
      
      if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
        
        
        test_norm(residuals(model))
        # print("Significant interaction between exposure (pathogen vs solvent) and period, AND significant interaction beween size (small vs big) and period.")
        # print(anov)
        
        contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,0,-1,0))
        contrast_matrix_size     <- rbind("Delta_small minus Delta_big"=c(0,0,0,0,0,-1))
        
        posthoc_groups_exposure <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data)
        posthoc_groups_size     <- get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size    ,which_levels="size_order",dataset=data)
        
      }else{
        if (p_interaction_size>p_interaction_exposure){ 
          ###if interaction with size is not significant and greater than interaction with exposure, repeat analysis, but exposure only
          # print( "Interaction between size and period is not significant:")
          # print(anov["period:size",])
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
          
          ###fit model 
          model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) + (1|antID),data=data)
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
          }
          
          p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
          if (p_interaction_exposure>0.05){
            # print( "Interaction between exposure and period is not significant:")
            # print(anov["period:exposure",])
            
          }else{
            contrast_matrix_exposure <- rbind("Delta_control minus Delta_pathogen"=c(0,0,0,-1))
            
            posthoc_groups_exposure <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_exposure,which_levels="exposure_order",dataset=data))
            names(posthoc_groups_exposure) <- variable_list[i]
            post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_exposure)
            
            # print( "Significant interaction between exposure and period:")
            # print(anov)
          }
          
        }else{
          ###if interaction with exposure is not significant and greater than interaction with size, repeat analysis, but size only
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
          
          # print( "Interaction between exposure and period is not significant:")
          # print(anov["period:exposure",])
          # 
          ###fit model 
          model <- lmer(   variable ~ period*size + (1|time_of_day) + (1|colony) + (1|antID),data=data)
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
          }
          
          p_interaction_size <- anov["Pr(>F)"]["period:size","Pr(>F)"]
          if (p_interaction_size>0.05){
            # print( "Interaction between size and period is not significant:")
            # print(anov)
            
          }else{
            # print( "Significant interaction between size and period:")
            # print(anov)
            contrast_matrix_size <- rbind("Delta_small minus Delta_big"=c(0,0,0,-1))
            
            posthoc_groups_size <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix_size,which_levels="size_order",dataset=data))
            names(posthoc_groups_size) <- variable_list[i]
            post_hoc_outcomes                <- c(post_hoc_outcomes,posthoc_groups_size)
          }
        }
      }
    }
  }
  rm(list=ls()[which(grepl("p_interaction",ls()))])
  rm(list=ls()[which(grepl("posthoc_groups_",ls()))])
}
rownames(stats_outcomes) <- NULL

