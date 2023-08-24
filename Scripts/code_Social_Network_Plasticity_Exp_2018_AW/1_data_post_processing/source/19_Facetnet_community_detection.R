####Facetnet_community_detection.R#####

#### Takes an interaction list as an input, builds a network to provide a continuous measure of workers' social maturity

#### requires the following:
#### https://c4science.ch/source/facet_unil/
# TO Richardson, T Kay, R Braunschweig, OA Journeau, M Rüegg, ... Ant behavioral maturation is mediated by a stochastic transition between two fundamental states. Current Biology 31, 1-8

###Created by TO Richardson and Adriano Wanderlingh

## FACENET community partition parameters
m            <- 2   ## how many communities do we want to instruct FACETNET to search for...?
alpha        <- 0.5 ## used for modulating the `memory` - how much the community structure @ time t influences that at t+1 - only matters when we ask facetnet to pay attention to this...
t_step       <- 0   ## not important
N_ITERATIONS <- 50

#define main dirs (THESE ARE ALSO DEFINED IN A_MAIN EXPERIMENT)
if(!exists("DISK")){DISK <- "Seagate Portable Drive"}
data_path <- paste0("/media/cf19810/",DISK,"/Lasius-Bristol_pathogen_experiment/main_experiment")
code_path  <- "/home/cf19810/Ant_Tracking/Scripts/code_Social_Network_Plasticity_Exp_2018_AW/1_data_post_processing/source"
source(paste(code_path,"/functions_and_parameters.R",sep=""))

FACETNET_DIR <- "/home/cf19810/facetNet"  ## the Python script is here - can be anywhere, but this must point to it...

input_path           <-  paste(data_path,"/intermediary_analysis_steps/full_interaction_lists/PreTreatment/observed",sep="")
network_files <- list.files(path=paste(input_path),full.names=T)
## load the file containing the % time each ant spent outside
TaskStats_File <- paste(data_path, "/processed_data/individual_behaviour/pre_treatment/network_position_vs_time_outside.dat", sep="/")
TaskStats_all      <- read.table(TaskStats_File , header=T, stringsAsFactors = F)

#where scores are saved
WORKDIR      <- paste(data_path,"/Soft_community_scores",sep="")

to_keep <- c(ls(),"to_keep","network_files")

for (network_file in network_files){
  #### re-open the task_groups (where results are saved) as it gets updated at the end of the loop (the loaded version has to be the newest one)
  task_groups    <- read.table(paste(data_path,"original_data/task_groups.txt",sep="/"),header=T,stringsAsFactors = F)
  # Create new columns if they don't exist in 'task_groups'
  if(!"task_group_FACETNET_0.5" %in% names(task_groups)) {task_groups$task_group_FACETNET_0.5 <- NA}
  if(!"Forager_score" %in% names(task_groups)) {task_groups$Forager_score <- NA}
  
  
  ####get file metadata
  root_name          <- gsub("_interactions.txt","",unlist(strsplit(network_file,split="/"))[grepl("interactions",unlist(strsplit(network_file,split="/")))]) # LS: replace grepl("colony", ...) with grepl("interactions")
  Cassette           <- paste(root_name,"_interactions", sep="")
  components         <- unlist(strsplit(root_name,split="_"))
  # colony             <- components[grepl("colony",components)]
  colony             <- unlist(strsplit(root_name,split="_"))[1] # LS
  treatment          <- info[which(info$colony==colony),"treatment"] #AW: no need for as.numeric() 
  colony_size        <- info[which(info$colony==colony),"colony_size"]
  period             <- "Pre"
  
  print(basename(root_name)) #"\r",
  
  #MAKE A FOLDER PER COLONY
  FACETNET_REP_folder <- paste(WORKDIR, paste(root_name,"_FACETNET_iters",sep=""), sep="/"); if (!file.exists(FACETNET_REP_folder)){dir.create(FACETNET_REP_folder)}
  Module_File        <- paste(FACETNET_REP_folder, paste(Cassette,"Modularities.txt",sep=""), sep="/")
  
  ####get appropriate task_group list, treated list and tag
  tag <- read.tag(tag_list)
  alive <- tag$tag #AW
  colony_task_group  <- task_groups[which(task_groups$colony==colony), c("tag","task_group")];   colony_task_group <- subset(colony_task_group, tag %in% alive)
  queenid            <- as.character(colony_task_group[which(colony_task_group$task_group=="queen"),"tag"])
  TaskStats          <- TaskStats_all[which(TaskStats_all$colony==colony),]
  TaskStats          <- subset(TaskStats, tag %in% alive)
  
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### PART 1: re-format the aggregated interactions list for facetnet soft community labelling  #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

if (!file.exists("FACETNET_INPUT")){
  ####read interactions
  interactions       <- read.table(network_file,header=T,stringsAsFactors = F)
  #remove dead ants from interactions list  #AW
  interactions <- subset(interactions, Tag1 %in% alive)
  interactions <- subset(interactions, Tag2 %in% alive)
  ## time-aggregated edge list
  EdgeList              <- aggregate(Box ~ Tag1 + Tag2, FUN=length, data=interactions)  ; colnames(EdgeList)[ncol(EdgeList)] <- "Count"       ## pairwise contact count
  
  ## DEFINE INPUT FILE - needs to be formatted for facetnet python to understand
  FACETNET_INPUT <- paste(FACETNET_REP_folder,"/",Cassette,"interactions_Time-aggregated_network,FACETNET_format.txt",sep="")
  
  ## Export the aggregated edge list
  write.table(EdgeList, file=FACETNET_INPUT, row.names=F, col.names = F, quote=F)
  }

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### PART 2: Repeatedly apply facetnet to generate a community partition ### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

for (ITER in 1:N_ITERATIONS)  {  ## the modularity of the found solutions vary quite a bit... so repeat many times & select the highest-modularity solution 
  ## create new subdirectory for community results with m
  FACETNET_REP_OUTPUT_DIR_M <- paste(FACETNET_REP_folder, paste(Cassette,",m=",m,",Iteration=",ITER,sep=""), sep="/"); if (!file.exists(FACETNET_REP_OUTPUT_DIR_M)){dir.create(FACETNET_REP_OUTPUT_DIR_M)}
  ## check if the outputs already exist
  if ( file.exists( paste(FACETNET_REP_OUTPUT_DIR_M, "soft_comm_step_alpha0.5_nw0.csv", sep="/") ))
    {print(paste("Facetnet output exists for m=",m,"modules, iteration #",ITER))
    }else{cat(paste("Facetnet output missing for m=",m, "modules, iteration #",ITER),"\r")
    
      FACETNET_INPUT_brackers <- paste0("'",FACETNET_INPUT,"'") #make sure the command is not broken
      FACETNET_REP_OUTPUT_DIR_M_brackers <- paste0("'",FACETNET_REP_OUTPUT_DIR_M,"'")
      
    ## construct facetnet command 
    command <- paste("python3", paste(FACETNET_DIR, "facetnet_step.py", sep = "/"),
                     FACETNET_INPUT_brackers,
                     alpha,
                     m,
                     FACETNET_REP_OUTPUT_DIR_M_brackers,
                     t_step,
                     sep = " ")
    ## jump through hoops to extract the modularity which is just printed to the prompt
    OutPut <- system(command, intern=TRUE) 
    
    ## if error, escape the loop & stop
    if ("TRUE" %in% grepl("rror",OutPut)){print("ERROR in OutPut"); break; print(OutPut)
      }else{
      Item       <- grep("modularity",OutPut) ## 
      Modularity <- as.numeric(gsub(" ", "", gsub("\\(","", gsub("\\)","", strsplit(OutPut[Item], "=")[[1]][2]))))  ## extract the modularity from the output
      }
    ## stack modularity for each m
    Modules <- data.frame(colony, treatment,period, ITER, alpha, MODULARITY=Modularity )
    ## check if the modularity for m communities has already been calculated and recorded
    if (file.exists(Module_File)){ 
      Modules_precomputed <- read.table(Module_File, header = T)
      if (!m %in% Modules_precomputed$N_modules) {write.table (Modules, file=Module_File, row.names=F, col.names=F, quote=F, append=T)  } ## only add a line to the file if m is not already there
      }else{
        write.table (Modules, file=Module_File, row.names=F, col.names=T, quote=F, append=F)  }  ## if the file doesn't exist, create it
    }
  }##ITER


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### PART 3: Find the top-modularity solution & assign biological labels to both communities # ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

## reduce these files to the current colony, treatment,period
TaskStats_ColTreatPer <- TaskStats[which(TaskStats$colony==colony & TaskStats$treatment==treatment & tolower(TaskStats$period)==tolower(period) ), c("tag","time_hours","prop_time_outside")]

## **KEY STAGE** average the % time outside across the 3 hour pre-introduction bins
## IS THIS RIGHT?
TaskStats_ColTreatPer_MEAN <- aggregate(prop_time_outside ~ tag, FUN=mean, na.rm=TRUE, na.action=NULL, TaskStats_ColTreatPer)
if (nrow(TaskStats_ColTreatPer_MEAN) != nrow(colony_task_group)) {print("ERROR : the number of tags in TaskStats_ColTreatPer_MEAN is different from that in colony_task_group")}

if (file.exists(Module_File)){
  Modules <- read.table(Module_File, header = T)
  ITER    <- Modules$ITER [which.max(Modules$MODULARITY)]

  ## load community scores for the top solution 
  FACETNET_REP_OUTPUT_DIR_M <- paste(FACETNET_REP_folder, paste(Cassette,",m=",m,",Iteration=",ITER,sep=""), sep="/")
  ## load the community scores
  Ant_Modules <- read.csv( file = paste(FACETNET_REP_OUTPUT_DIR_M, "soft_comm_step_alpha0.5_nw0.csv", sep="/")); colnames(Ant_Modules)[match("id",colnames(Ant_Modules))] <- "tag"
 
  ## convert continuous -> binary communities
  Ant_Modules$cluster_binary <- apply( Ant_Modules[,c("cluster_0","cluster_1")], 1, function(x) {c(0,1)[which.max(x)]} ) 

  ## add the caste information 
  Combined <- merge(x = TaskStats_ColTreatPer_MEAN, y = colony_task_group, by="tag", all =TRUE)
  ## add the scores for the 2 as-yet unlabelled communities to the df containing % time outside
  Combined <- merge(x = Combined,                   y = Ant_Modules,    by="tag", all =TRUE)

  ##  Label the foragers using average time spent outside (& maybe also the distance to the queen)
  ClusterMeans <- aggregate(cbind(prop_time_outside) ~ cluster_binary, FUN=mean, na.rm=T, Combined)
  
  ## which community has the lower mean % time outside? --> THE NURSES
  Nurses_by_prop_time_outside <-  ClusterMeans$cluster_binary [which.min(ClusterMeans$prop_time_outside)]
  ## which community contains the queen? ---> THE NURSES
  Nurses_by_queen_affiliation <-  Combined$cluster_binary [ which(Combined$task_group=="queen") ]
  
  ##  check agreement between the 2 diagnosis methods - if they don't agree, we will need to discuss
  if (Nurses_by_prop_time_outside != Nurses_by_queen_affiliation){
    print("PROBLEM: community identification by % time outside DISAGREES with that by queen affiliation")}
  else{  
    ## the nurse community is the one that spends the least time outside and that includes the queen (assuming the 2 methods agree...)
    NURSE_COMMUNITY   <- ClusterMeans$cluster_binary [ which.min(ClusterMeans$prop_time_outside)]
    FORAGER_COMMUNITY <- ClusterMeans$cluster_binary [ which.max(ClusterMeans$prop_time_outside)]
    
    ## alter the column names in Ant_Modules to reflect task group labels - a continuous score in the range 0-1 for both modules, SUMMING TO 1 ACROSS MODULES!!
    colnames(Combined) [match( paste("cluster_",FORAGER_COMMUNITY,sep=""), colnames(Combined))] <- "Forager_score"
    colnames(Combined) [match( paste("cluster_",NURSE_COMMUNITY,sep=""),   colnames(Combined))] <- "Nurse_score"
    
    ## add categorical versions of the FACETNET scores - just to compare with the original threshold-based task_group:
    Combined$cluster_binary <- NULL
    Combined$task_group_FACETNET_0.5 <- NA
    Combined$task_group_FACETNET_0.5 [which(Combined$Forager_score > Combined$Nurse_score)] <- "forager"
    Combined$task_group_FACETNET_0.5 [which(Combined$Forager_score < Combined$Nurse_score)] <- "nurse"
    Combined$task_group_FACETNET_0.5 [which(Combined$tag == queenid)] <- "queen"
    }
  
  ## tally the agreement between the original and new (facetnet) task labels
  Agreement <- table(Combined$task_group, Combined$task_group_FACETNET_0.5); rownames(Agreement) <- paste("orig", rownames(Agreement),sep="_"); colnames(Agreement) <- paste("facet", colnames(Agreement), sep="_"); print(Agreement)
  
# Save output to the task_groups file
Combined$colony <- colony

# Merge by specified columns
merged_data <- merge(task_groups, Combined[ ,which(names(Combined) %in%  c("colony", "tag", "Forager_score", "task_group_FACETNET_0.5"))],
                     by=c("colony","tag"), all.x=TRUE, suffixes=c("", ".update"))

# Update FACETNET results based on the specific colony
# !is.na(merged_data$task_group_FACETNET_0.5.update) checks if there is an updated value available for the task_group_FACETNET_0.5 column
merged_data$task_group_FACETNET_0.5 <- ifelse(merged_data$colony == colony & !is.na(merged_data$task_group_FACETNET_0.5.update), merged_data$task_group_FACETNET_0.5.update, merged_data$task_group_FACETNET_0.5)
merged_data$Forager_score <- ifelse(merged_data$colony == colony & !is.na(merged_data$Forager_score.update), merged_data$Forager_score.update, merged_data$Forager_score)
merged_data$Forager_score <- round(merged_data$Forager_score,3)

merged_data <- merged_data[order(merged_data$colony),]

write.table(merged_data[, names(task_groups)],
            file=paste(data_path,"original_data/task_groups.txt",sep="/"),
            append = F, col.names = T, row.names = F, quote = F, sep = "\t")
}
clean()
}


#PLOT results
task_groups_A    <- read.table("/media/cf19810/Seagate Portable Drive/Lasius-Bristol_pathogen_experiment/main_experiment/original_data/task_groups.txt",header=T,stringsAsFactors = F)

task_groups_A$size     <- unlist(lapply( task_groups_A$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))

# Generate plots
ggplot(task_groups_A, aes(x = Forager_score, colour = size)) + 
  #geom_density(alpha = 0.6,size=1.5, adjust=1/1.2) + # 'adjust' changes the smoothing
  geom_line(aes(color=size,group = colony), stat="density", size=1.5, alpha=0.15, adjust=1/1.2) +
  geom_line(aes(color=size), stat="density", size=2.5, alpha=1, adjust=1/1.2) +
  geom_vline(aes(xintercept = 0.5), linetype = "dashed",colour="grey20") + 
  geom_vline(aes(xintercept = 1/4), linetype = "dashed",colour="grey20") + 
  #facet_wrap(~ colony, scales = "free") +
  theme_minimal() +
  xlab("Social maturity")




warning("rewrite table by calling the task_group file, faceting for colony")
## tally the agreement between the original and new (facetnet) task labels
Agreement <- table(Combined$task_group, Combined$task_group_FACETNET_0.5); rownames(Agreement) <- paste("orig", rownames(Agreement),sep="_"); colnames(Agreement) <- paste("facet", colnames(Agreement), sep="_"); print(Agreement)


# TO DOS
#test different thresholds?
# - rerun qpcr script with new labs (merge from where those are saved)
