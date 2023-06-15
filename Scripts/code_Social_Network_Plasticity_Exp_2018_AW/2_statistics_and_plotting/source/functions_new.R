clean <- function() {
  rm(list = ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv) %in% to_keep], envir = .GlobalEnv)
  no_print <- gc(verbose = F)
}

GetColorHex <- function(color) {
  clr <- col2rgb(color)
  hex_and_col <-
    sprintf("#%02X%02X%02X %3d %3d %3d", clr[1], clr[2], clr[3], clr[1], clr[2], clr[3])
  hex <- unlist(strsplit(hex_and_col, split = " "))[1]
  return(hex)
}

get_posthoc_groups <-
  function(model,
           contrast_matrix,
           which_levels,
           dataset) {
    levels <- get(which_levels)
    if (is.null(names(levels))) {
      names(levels) <- paste("Delta_", levels, sep = "")
    }
    
    post_hoc <-
      summary(glht(model, contrast_matrix), test = adjusted("BH"))
    # print("z value");print(post_hoc$test$tstat);print("Pr>|z|");print(post_hoc$test$pvalues);
    p_values <-
      as.numeric(post_hoc$test$pvalues)
    names(p_values) <- names(post_hoc$test$coefficients)
    
    post_hoc_levels <- names(levels)
    post_hoc_mat <-
      matrix(NA,
             nrow = length(post_hoc_levels) - 1,
             ncol = length(post_hoc_levels) - 1)
    rownames(post_hoc_mat) <-
      post_hoc_levels[2:length(post_hoc_levels)]
    colnames(post_hoc_mat) <-
      post_hoc_levels[1:(length(post_hoc_levels) - 1)]
    for (i in 1:nrow(post_hoc_mat)) {
      for (j in 1:i) {
        post_hoc_mat[i, j] <-
          as.logical(as.numeric(p_values[paste(colnames(post_hoc_mat)[j],
                                               " minus ",
                                               rownames(post_hoc_mat)[i],
                                               sep = "")]) > 0.05)
      }
    }
    g <- post_hoc_mat
    g <- cbind(rbind(NA, g), NA)
    g <- replace(g, is.na(g), FALSE)
    g <- g + t(g)
    diag(g) <- 1
    n <- length(post_hoc_levels)
    rownames(g) <- 1:n
    colnames(g) <- 1:n
    #g
    same <- which(g == 1)
    topology <-
      data.frame(N1 = ((same - 1) %% n) + 1, N2 = ((same - 1) %/% n) + 1)
    topology <-
      topology[order(topology[[1]]), ] # Get rid of loops and ensure right naming of vertices
    g3 <- simplify(graph.data.frame(topology, directed = FALSE))
    #get.data.frame(g3)
    #plot(g3)
    res <- maximal.cliques(g3)
    
    # Reorder given the smallest level
    clique_value <- NULL
    if (which_levels == "level_names") {
      form <-
        as.formula(paste("variable ~ ", paste(
          c("treatment", "task_group"), collapse = "+"
        )))
    } else{
      form <-
        as.formula(paste("variable ~ ", paste(
          gsub("_order", "", which_levels), collapse = "+"
        )))
    }
    
    means <- aggregate(form, FUN = mean, data = dataset)
    means$predictor <- names(levels)
    for (i in 1:length(res)) {
      clique_value <-
        c(clique_value, mean(means[as.numeric(unlist(res[[i]])), "variable"]))
    }
    res <- res[order(clique_value)]
    
    # Get group letters
    lab.txt <- vector(mode = "list", n)
    lab <- letters[seq(res)]
    for (i in seq(res)) {
      for (j in res[[i]]) {
        lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
      }
    }
    post_hoc_groups <-
      unlist(lab.txt)
    names(post_hoc_groups) <- levels[post_hoc_levels]
    return(post_hoc_groups)
  }

log_transf <- function(x) {
  if (all(x > 0)) {
    replac_val <- 0
  } else if (all(x >= 0)) {
    replac_val <- (min(x[x != 0], na.rm = T)) / sqrt(2)
  } else{
    replac_val_1 <- -min(x, na.rm = T)
    y <- x + replac_val_1
    replac_val <- replac_val_1 + (min(y[y != 0], na.rm = T)) / sqrt(2)
  }
  return(log10(x + replac_val))
}

sqrt_transf <- function(x) {
  if (all(x >= 0)) {
    replac_val <- 0
  } else{
    replac_val <- -min(x, na.rm = T)
  }
  return(sqrt(x + replac_val))
}

test_norm <- function(resids) {
  print("Testing normality")
  if (length(resids) <= 300) {
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
  } else{
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print(
      "Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero"
    )
    print(kurtosis(resids))
  }
}

# convert significance levels to stars
from_p_to_ptext <- function(pvalue) {
  if (is.na(pvalue)) {
    pvaluetext <- "p = NA"
  } else{
    if # (pvalue< 0.00001){
    #   pvaluetext <- "*****"
    # }else if (pvalue< 0.0001){
    #   pvaluetext <- "****"
    # }else if
    (pvalue < 0.0001) {
      pvaluetext <- "***"
    } else if (pvalue < 0.005) {
      pvaluetext <- "**"
    } else if (pvalue < 0.05) {
      pvaluetext <- "*"
    } else if (pvalue < 0.1) {
      pvaluetext <- paste("p=", sprintf("%.3f", pvalue), sep = "")
    } else if (pvalue < 1) {
      pvaluetext <- paste("p=", sprintf("%.2f", pvalue), sep = "")
    } else{
      pvaluetext <- "p = 1"
    }#if (pvalue< 10^-6)
  }
  return(pvaluetext)
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#AW adapted Functions

collective_analysis_no_rescal <- function(data_path=data_path){
  
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
  barplot_delta_period_list <- list()
  
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
    barplot_delta_period <- barplot_delta(dataset=data,predictor="treatment",post_hoc_outcomes=post_hoc_outcomes,stats_outcomes=stats_outcomes,i=i,type="collective",collective=T,plot_untransformed=T,diff_type="") #form_stat=NULL,
    print(barplot_delta_period) 
    
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period
    #names(barplot_delta_period_list[[i]]) <- variable_list[i]
  }
  rownames(stats_outcomes) <- NULL
  
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
}

collective_analysis_rescal <- function(data_path=data_path){
  
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
  barplot_delta_period_list <- list()
  
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
    barplot_delta_period <- barplot_delta(dataset=data,predictor="treatment",post_hoc_outcomes=post_hoc_outcomes,stats_outcomes=stats_outcomes,i=i,type="collective",collective=T,plot_untransformed=T,diff_type="") #form_stat=NULL,
    print(barplot_delta_period)
    
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period
  }
  rownames(stats_outcomes) <- NULL
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
}

individual_ONE_analysis_rescal <- function(data_path=data_path,which_individuals){
  # for ONE type of individuals (e.g. treated workers only, or queens only)
  
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
  if (pattern=="individual_behavioural_data") {
    task_groups <- read.table(paste(root_path,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
    task_groups <- task_groups[which(!duplicated(task_groups)),]
    data        <- merge(data,task_groups[c("colony","tag","task_group")],all.x=T,all.y=F)
    
  } else if (pattern=="individual_data") {
    #make sure that the task_group is named correctly
    data$task_group <- data$status 
    data$status     <- NULL
    treated_worker_list <- read.table(paste(root_path,"original_data","treated_worker_list.txt",sep="/"),header=T,stringsAsFactors = F)
    treated_worker_list <- treated_worker_list[which(!duplicated(treated_worker_list)),]
    treated_worker_list$status <- "treated"
    data        <- merge(data,treated_worker_list[c("colony","tag","status")],all.x=T,all.y=F) 
    data[which(is.na(data$status)),"status"] <- "untreated"
  } 
  
  data[which(data$status=="treated"),"task_group"] <- "treated"

  ##2c. keep only target individuals
  data <- data[which(data$task_group%in%which_individuals),]
  
  ##2d.  ###add a unique antid column
  data <- within(data,antID <- paste(colony,tag,sep="_"))
  
  ###2. Loop over variables
  data_ori <- data
  stats_outcomes    <- NULL
  post_hoc_outcomes <- NULL
  barplot_delta_period_list <- list()
  
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
    barplot_delta_period <- barplot_delta(dataset=data,predictor="treatment",post_hoc_outcomes=post_hoc_outcomes,stats_outcomes=stats_outcomes,i=i,type="individual",collective=F,plot_untransformed=T,diff_type="absolute_difference") #form_stat=NULL,
    print(barplot_delta_period)
    
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period 
    #ADD SAVING OF THE PLOT on single pdf file as done in plot_grooming file
    
  }
  rownames(stats_outcomes) <- NULL
  
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
  
}

individual_TWO_analysis_rescal <- function(data_path=data_path,which_individuals){
  # for ONE type of individuals (e.g. treated workers only, or queens only)
  
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
  if (pattern=="individual_behavioural_data") {
    task_groups <- read.table(paste(root_path,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
    task_groups <- task_groups[which(!duplicated(task_groups)),]
    data        <- merge(data,task_groups[c("colony","tag","task_group")],all.x=T,all.y=F)
    
  } else if (pattern=="individual_data") {
    #make sure that the task_group is named correctly
    data$task_group <- data$status 
    data$status     <- NULL
    treated_worker_list <- read.table(paste(root_path,"original_data","treated_worker_list.txt",sep="/"),header=T,stringsAsFactors = F)
    treated_worker_list <- treated_worker_list[which(!duplicated(treated_worker_list)),]
    treated_worker_list$status <- "treated"
    data        <- merge(data,treated_worker_list[c("colony","tag","status")],all.x=T,all.y=F) 
    data[which(is.na(data$status)),"status"] <- "untreated"
  } 
  
  data[which(data$status=="treated"),"task_group"] <- "treated"
  
  
  
  ##2c. keep only target individuals
  data <- data[which(data$task_group%in%which_individuals),]
  
  ##2d.  ###add a unique antid column
  data <- within(data,antID <- paste(colony,tag,sep="_"))
  
  ###2. Loop over variables
  data_ori <- data
  stats_outcomes    <- NULL
  post_hoc_outcomes <- NULL
  barplot_delta_period_list <- list()
  
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
    
    #plot
    barplot_delta_period <- barplot_delta(dataset=data,predictor="treatment",post_hoc_outcomes=post_hoc_outcomes,stats_outcomes=stats_outcomes,i=i,type="individual",collective=F,plot_untransformed=T,diff_type="absolute_difference") #form_stat=NULL,
    print(barplot_delta_period)
    
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period 
    
  }
  rownames(stats_outcomes) <- NULL
  
  
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
  
}


barplot_delta <-
  function(dataset = data,
           predictor,
           post_hoc_outcomes,
           stats_outcomes,
           i, #variable number
           type,
           collective,
           plot_untransformed, #is always TRUE as the the variable fed to the plotting is transformed beforehand (see section: transform variable)
           diff_type) {
    #form_stat=NULL,
    
    #get difference on RAW data
    diff <-
      create_diff(dataset,
                  predictor,
                  type,
                  collective,
                  plot_untransformed,
                  diff_type) #form_stat=NULL,
    diff["time"] <- "Delta"
    
    if (!collective) {
      # I SUM BY INDIVIDUAL, THEN MEAN AFTER
      diff <-
        aggregate(
          na.rm = T,
          na.action = "na.pass",
          variable ~ predictor + colony_size + colony + time + tag,
          FUN = mean,
          data = diff
        )
    }
    
    
    #mean by colony (drop time_of_day)
    diff <-
      aggregate(
        na.rm = T,
        na.action = "na.pass",
        variable ~ predictor + colony_size + colony + time,
        FUN = mean,
        data = diff
      )
    
    means <-
      aggregate(
        na.rm = T,
        na.action = "na.pass",
        variable ~ time + predictor,
        FUN = "mean",
        data = diff
      )
    ses <-
      aggregate(
        na.rm = T,
        na.action = "na.pass",
        variable ~ time + predictor,
        FUN = "std.error",
        data = diff
      )
    
    names(means)[names(means) == "variable"] <-
      "mean"
    names(ses)[names(ses) == "variable"] <-
      "se"
    means <- merge(means, ses)
    means <-
      means[order(match(means$predictor, status_order),
                  match(means$time, status_order)), ]
    #to_plot <- unique(means[c("time","predictor")])
    
    
    
    
    ##### PLOT
    # significant letters/stars position
    max_mean_se <- 2 * max(means$mean + means$se)
    # rename vars
    
    
    
    # if interaction period*treatment
    if (all(names(post_hoc_outcomes[[variable_list[i]]]) == levels(dataset$treatment))) {
      print("period*treatment interaction significant")
      
      # Create the base ggplot object
      plot_var <-
        ggplot(means,
               aes_string(x = "predictor", y = "mean", fill = "predictor")) +
        geom_errorbar(
          aes_string(ymin = paste0("mean - se"), ymax = paste0("mean + se")),
          position = position_dodge2(width = 0.8, preserve = "single")
        ) +
        geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
        STYLE +
        colFill_treatment +
        labs(y = paste("\u0394", names(variable_list[i]), sep = " "), x = predictor)
      
      if (variable_list[i] %in% names(post_hoc_outcomes)) {
        additional_geoms <- list(geom_text(
          aes(label = post_hoc_outcomes[[variable_list[i]]][predictor], y = max_mean_se),
          position = position_dodge2(width = 0.8, preserve = "single"),
          vjust = -0.4
        ))
        plot_var <- plot_var + additional_geoms
      }
      
      if (predictor == "treatment") {
        scale_geom <- list(scale_x_discrete(
          labels = function(x)
            str_wrap(gsub("big", "large", gsub("\\.", " ", x)), width = 4)
        ))
        plot_var <- plot_var + scale_geom
      } else{
        scale_geom <- list(scale_x_continuous(
          labels = function(x)
            str_wrap(x, width = 4),
          breaks = 1:length(levels(means$predictor)) + c(0, 0, 1, 1)
        ))
        plot_var <- plot_var + scale_geom
      }
      
      
      
    } else if (all(names(post_hoc_outcomes[[variable_list[i]]]) == levels(dataset$size))) {
      print("Only period*size interaction significant")
      
      # Reorder the vector by size
      treatment_order_size <-
        treatment_order[order(grepl("small", treatment_order, fixed = TRUE), decreasing = TRUE)]
      means$predictor <-
        factor(means$predictor , levels = treatment_order_size[which(treatment_order_size %in% means$predictor)])
      # Add a new variable for grouping and x-axis
      means$group <-
        ifelse(grepl("small", means$predictor), "small", "big")
      #create a gap in the layout add a +1 to the big colonies
      means$x_axis <-
        as.numeric(means$predictor) + ifelse(means$group == "big", 1, 0)
      
      # Create the base ggplot object
      plot_var <-
        ggplot(means, aes(x = x_axis, y = mean, fill = predictor)) +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                      position = position_dodge2(width = 0.8, preserve = "single")) +
        geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
        STYLE +
        colFill_treatment +
        labs(y = paste("\u0394", names(variable_list[i]), sep = " "), x = predictor) +
        theme(axis.text.x = element_text())
      
      if (predictor == "treatment") {
        scale_geom <- list(scale_x_continuous(
          labels = function(x)
            str_wrap(gsub("big", "large", gsub("\\.", " ", x)), width = 4),
          breaks = 1:length(levels(means$predictor)) + c(0, 0, 1, 1)
        ))
        plot_var <- plot_var + scale_geom
      } else{
        scale_geom <- list(scale_x_continuous(
          labels = function(x)
            str_wrap(x, width = 4),
          breaks = 1:length(levels(means$predictor)) + c(0, 0, 1, 1)
        ))
        plot_var <- plot_var + scale_geom
      }
      
      
      # Add geom_segment and geom_text if the condition is true (for more than 1 element in the ggplot if statement, there is the error: Error in `+.gg`:! Cannot add ggproto objects together )
      if (variable_list[i] %in% names(post_hoc_outcomes)) {
        additional_geoms <- list(geom_segment(
          aes(
            x = sum(means$x_axis[means$group == "small"]) / 2,
            xend = sum(means$x_axis[means$group == "big"]) / 2,
            y = max_mean_se,
            yend = max_mean_se
          )
        ),
        geom_text(
          aes(
            x = mean(means$x_axis),
            y = max_mean_se + 1 / 3 * (max_mean_se),
            label = from_p_to_ptext(stats_outcomes[which(
              stats_outcomes$variable == variable_list[i] &
                stats_outcomes$predictor == "period:size"
            ), "pval"])
          )
        ))
        plot_var <- plot_var + additional_geoms
      }
      
    } else if (all(names(post_hoc_outcomes[[variable_list[i]]]) == levels(dataset$exposure))) {
      print("Only period*exposure interaction significant")
      
      # Reorder the vector by size
      treatment_order_exp <-
        treatment_order[order(grepl("control", treatment_order, fixed = TRUE),
                              decreasing = TRUE)]
      means$predictor <-
        factor(means$predictor , levels = treatment_order_exp[which(treatment_order_exp %in% means$predictor)])
      # Add a new variable for grouping and x-axis
      means$group <-
        ifelse(grepl("control", means$predictor), "control", "pathogen")
      #create a gap in the layout add a +1 to the pathogen colonies
      means$x_axis <-
        as.numeric(means$predictor) + ifelse(means$group == "pathogen", 1, 0)
      
      # Create the base ggplot object
      plot_var <-
        ggplot(means, aes(x = x_axis, y = mean, fill = predictor)) +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                      position = position_dodge2(width = 0.8, preserve = "single")) +
        geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
        STYLE +
        colFill_treatment +
        labs(y = paste("\u0394", names(variable_list[i]), sep = " "), x = predictor) +
        theme(axis.text.x = element_text())
      
      if (predictor == "treatment") {
        scale_geom <- list(scale_x_continuous(
          labels = function(x)
            str_wrap(gsub("big", "large", gsub("\\.", " ", x)), width = 4),
          breaks = 1:length(levels(means$predictor)) + c(0, 0, 1, 1)
        ))
        plot_var <- plot_var + scale_geom
      } else{
        scale_geom <- list(scale_x_continuous(
          labels = function(x)
            str_wrap(x, width = 4),
          breaks = 1:length(levels(means$predictor)) + c(0, 0, 1, 1)
        ))
        plot_var <- plot_var + scale_geom
      }
      
      # Add geom_segment and geom_text if the condition is true (for more than 1 element in the ggplot if statement, there is the error: Error in `+.gg`:! Cannot add ggproto objects together )
      if (variable_list[i] %in% names(post_hoc_outcomes)) {
        additional_geoms <- list(geom_segment(
          aes(
            x = sum(means$x_axis[means$group == "control"]) / 2,
            xend = sum(means$x_axis[means$group == "pathogen"]) / 2,
            y = max_mean_se,
            yend = max_mean_se
          )
        ),
        geom_text(
          aes(
            x = mean(means$x_axis),
            y = max_mean_se + 1 / 3 * (max_mean_se),
            label = from_p_to_ptext(stats_outcomes[which(
              stats_outcomes$variable == variable_list[i] &
                stats_outcomes$predictor == "period:size"
            ), "pval"])
          )
        ))
        plot_var <- plot_var + additional_geoms
      }
      
      
      
      
      
      
      
      
      
    }
    
    #remove legend?
    plot_var <- plot_var + theme(legend.position = "none")   
    
    return(plot_var)
  }



#AW simplified Stroeymeyt FUNCTIONS
create_diff <-
  function(dataset,
           predictor,
           type,
           collective,
           plot_untransformed,
           diff_type) {
    #form_stat=NULL,
    if (is.null(dataset)) {
      return(NULL)
    } else{
      #################################################################################
      ####plot
      #################################################################################
      dataset$predictor <- dataset[, predictor]
      befores <-
        dataset[dataset$period == "pre", ]
      afters <- dataset[dataset$period == "post", ]
      # print(plot_untransformed)
      if (!plot_untransformed) {
        names(befores)[names(befores) == "variable"] <-
          "variable_before"
        names(afters)[names(afters) == "variable"] <- "variable_after"
      } else{
        names(befores)[names(befores) == "untransformed_variable"] <-
          "variable_before"
        names(afters)[names(afters) == "untransformed_variable"] <-
          "variable_after"
      }
      # if ((!all(!grepl(":treatment:predictor",as.character(form_stat))))|(!all(!grepl(" \\* treatment \\* predictor",as.character(form_stat))))){
      #   befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+treatment+antid+time_of_day,FUN=mean,data=befores)
      #   afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+treatment+antid+time_of_day,FUN=mean,data=afters)
      # }else{
      # if (!grepl("age",root_path)){
      if (collective | type != "individual") {
        befores <-
          aggregate(
            na.rm = T,
            na.action = "na.pass",
            variable_before ~ colony_size + colony + predictor + size + time_of_day,
            FUN = mean,
            data = befores
          )
        afters <-
          aggregate(
            na.rm = T,
            na.action = "na.pass",
            variable_after ~ colony_size + colony + predictor + size + time_of_day,
            FUN = mean,
            data = afters
          )
      } else{
        befores <-
          aggregate(
            na.rm = T,
            na.action = "na.pass",
            variable_before ~ colony_size + colony + predictor + size + tag + time_of_day,
            FUN = mean,
            data = befores
          )
        afters <-
          aggregate(
            na.rm = T,
            na.action = "na.pass",
            variable_after ~ colony_size + colony + predictor + size + tag + time_of_day,
            FUN = mean,
            data = afters
          )
      }
      # } #else{
      #   if (collective){
      #     befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+treatment+time_of_day,FUN=mean,data=befores)
      #     afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+treatment+time_of_day,FUN=mean,data=afters)
      #   }else{
      #     if(type!="individual"){
      #       befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+treatment+time_of_day,FUN=mean,data=befores)
      #       afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+treatment+time_of_day,FUN=mean,data=afters)
      #     }else{
      #       befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+treatment+antid+time_of_day,FUN=mean,data=befores)
      #       afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+treatment+antid+time_of_day,FUN=mean,data=afters)
      #     }
      #   }
      # }
      # }
      # befores["average_before"] <- mean(befores$variable_before,na.rm=T)
      
      
      
      diff <- merge(befores, afters, all = T)
      # if (diff_type=="absolute_difference"){
      diff["variable"] <-
        diff$variable_after - diff$variable_before
      # diff[""] <- "diff"
      # }
      # if (diff_type=="normalised_by_average_before"){
      #   standardisation <- aggregate(variable_before ~ size , FUN=mean,data=befores)
      #   names(standardisation)[which(names(standardisation)=="variable_before")] <- "average_before"
      #   befores <- merge(standardisation,befores,all.x=T)
      #   diff["variable"] <- (diff$variable_after-diff$variable_before)/diff$average_before; # diff[""] <- "diff"
      # }
      # if (diff_type=="relative_difference"){
      #   diff["variable"] <- diff$variable_after-diff$variable_before; diff[""] <- "diff"
      #   if ("treatment"%in%names(diff)){
      #     diff_1 <- diff[which(as.character(diff$treatment)==levels(diff$treatment)[1]),]
      #     diff <- diff[which(!as.character(diff$treatment)==levels(diff$treatment)[1]),]
      #     diff_1 <- aggregate(na.rm=T,na.action="na.pass",variable~predictor,FUN=get(relative_function),data=diff_1)
      #     names(diff_1)[names(diff_1)=="variable"] <- "mean_diff_1"
      #     diff <- merge(diff,diff_1,all.x=T)
      #     diff$variable <- diff$variable-diff$mean_diff_1
      #     diff <- diff[,which(names(diff)!="mean_diff_1")]
      #   }else{
      #     diff_1 <- diff[which(as.character(diff$predictor)==levels(diff$treatment)[1]),]
      #     diff <- diff[which(!as.character(diff$predictor)==levels(diff$treatment)[1]),]
      #
      #     diff_1 <- aggregate(na.rm=T,na.action="na.pass",variable~1,FUN=get(relative_function),data=diff_1)
      #     names(diff_1)[names(diff_1)=="variable"] <- "mean_diff_1"
      #
      #     diff <- merge(diff,diff_1,all.x=T)
      #     diff$variable <- diff$variable-diff$mean_diff_1
      #     diff <- diff[,which(names(diff)!="mean_diff_1")]
      #
      #   }
      #
      #
      # }
      #if (!grepl("age",root_path)){diff$predictor <- factor(diff$predictor)}else{diff$treatment <- factor(diff$treatment)}
      
      return(diff)
    }
  }



#@@@@@@@ STYLING FUNCTIONS @@@@@@@#


# Import system fonts
#font_import()
# font_import(paths = "/usr/share/texmf/fonts/opentype/public/tex-gyre", pattern = "texgyretermes")
# # Set the default font family
# loadfonts()

#font_import(pattern = "Liberation", prompt= FALSE)
#loadfonts(device = "pdf", quiet = TRUE)

# #Create a custom color scale FOR COLONIES + treatments
# FullPal <- scales::viridis_pal(option = "D")(20)
# myColorsSmall  <- tail(FullPal,5)
# myColorsLarge  <- head(FullPal,5)
# Cols_treatments <- c("#440154FF","#FDE725FF") #"#31688EFF"
# myColors      <- c(myColorsLarge,myColorsSmall, Cols_treatments)
# names(myColors) <- c("R3BP","R5BP","R7BP","R8BP","R12BP","R1SP", "R2SP", "R5SP", "R7SP","R11SP","pathogen.big","pathogen.small")
# colScale <- scale_colour_manual(name = "Colony",values = myColors,drop=TRUE)
# fillScale <- scale_fill_manual(name = "Colony",values = myColors,drop=TRUE)

#### CREATE CONSISTENT COLORING FOR ALL THE PLOTTING
# GENERATE 4 MAIN COLORS FOR THE 4 TREATMENTS BS,BP,SS,SP + SHADES FOR THE RESPECTIVE COLONIES

# Create a color palette with 4 colors as distant as possible
colors_full <- scales::viridis_pal(option = "D")(100)
# Create a list to store the subsets of colors
color_subsets <- list()
Shades <- list()
# Loop over the 4 colors to get shades of the colour in a +5, -5 range
for (i in c(10, 30, 70, 90)) {
  # BE CAREFUL: IF CHANGING THIS, ALSO CHANGE Treat_colors
  color_ramp <- colorRampPalette(colors_full[(i - 5):(i + 5)])
  #color_ramp <- colorRampPalette(colors_full[(i-1):(i+1)])
  color_subsets[[i]] <- color_ramp(12)
}


Treat_colors <-
  structure(list(
    Shades = c(colors_full[10], colors_full[30], colors_full[70], colors_full[90]),
    Treatment = c(
      "pathogen.big",
      "control.big",
      "pathogen.small",
      "control.small"
    )
  ),
  row.names = c(NA, 4L),
  class = "data.frame")

#show_col(c(colors_full[10],colors_full[30],colors_full[70],colors_full[90]))


#clean the output
color_subsets <- color_subsets[lapply(color_subsets, length) > 0]

#Darken the colours progressively to increase contrast among the colonies of the same treatment
for (i in 1:4) {
  # Define the color gradient
  color_shades <- color_subsets[[i]]
  
  # Convert the colors from hexadecimal to HSL (hue, saturation, lightness)
  
  colors_lightGrad <- c()
  # Decrease the lightness of each color by an increasing amount
  lightness_decrease <-
    rev(seq(
      from = 0,
      to = 0.2,
      length.out = length(color_shades)
    ))
  lightness_increase <-
    seq(from = 0,
        to = 0.2,
        length.out = length(color_shades))
  
  for (j in 1:length(color_shades)) {
    hsl_colors <- col2hsl(color_shades[j])
    hsl_colors[3] <- hsl_colors[3] - lightness_decrease[j]
    hsl_colors[3] <- hsl_colors[3] + lightness_increase[j]
    colors_lightGrad <- c(colors_lightGrad, hsl2col(hsl_colors))
  }
  
  Shades[[i]] <- colors_lightGrad
}

# #inspect output
# par(mfrow = c(2, 2))
# for (i in 1:4) {
#   plot(1:12, 1:12,
#        col = Shades[[i]], # color_subsets
#        pch = 19,
#        cex = 5,
#        xaxt = "n",
#        yaxt = "n",
#        xlab = "",
#        ylab = "")
# }

### ADD THE METADATA
#meta.data <- read.table(paste("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
# HARD CODE THEM
meta.data_size_treat <- c("BS", "SP", "SS", "BP")
meta.data_REP_treat <-
  c(
    "R1BS",
    "R1SP",
    "R1SS",
    "R10BP",
    "R11BP",
    "R11BS",
    "R11SP",
    "R11SS",
    "R12BP",
    "R12BS",
    "R12SP",
    "R13BS",
    "R13SP",
    "R13SS",
    "R14BP",
    "R14BS",
    "R14SP",
    "R14SS",
    "R2BP",
    "R2BS",
    "R2SP",
    "R2SS",
    "R3BP",
    "R3BS",
    "R3SP",
    "R3SS",
    "R4BP",
    "R4SS",
    "R5BP",
    "R5SP",
    "R5SS",
    "R6SP",
    "R6SS",
    "R7BP",
    "R7BS",
    "R7SP",
    "R7SS",
    "R8BP",
    "R8BS",
    "R8SP",
    "R8SS",
    "R9BS",
    "R9SP",
    "R9SS"
  )

# create groups to assign the colours
Cols <- list()
# divide each size_treat into a list element with its colonies inside
for (i in 1:length(meta.data_size_treat)) {
  treatment_labs <- meta.data_size_treat[i]
  treatment_vector <-
    unique(meta.data_REP_treat[grepl(treatment_labs, meta.data_REP_treat)])
  Cols[[i]] <- treatment_vector
  names(Cols)[i] <- treatment_labs
}
#name list elements according to the favoured pattern (the colour order I have been using since the first plots)
names(Shades)[1] <- "BP"
names(Shades)[2] <- "BS"
names(Shades)[3] <- "SP"
names(Shades)[4] <- "SS"

# dput((list_of_vectors))
# dput((Shades))

# Create an empty dataframe to store the results
colour_palette <- data.frame()

# bind together colours, REP_treats and treatments
# Loop over the list "Cols" to create the dataframe
for (group in names(Cols)) {
  group_cols <- Cols[[group]]
  group_shades <- Shades[[group]]
  # Take a random subset of the shades that matches the length of the cols
  rand_shades <- sample(group_shades, length(group_cols))
  # Create a dataframe with two columns: "Cols" and "Shades"
  group_colour_palette <-
    data.frame(Cols = group_cols,
               Shades = rand_shades,
               Treatment = group)
  # Append the current group dataframe to the overall dataframe
  colour_palette <- rbind(colour_palette, group_colour_palette)
}

# #visualise output
# ggplot(colour_palette, aes(x=Treatment, y=Shades, fill=Shades)) +
#   geom_tile(colour="white", size=0.5) +
#   scale_fill_identity() +
#   theme_void() +
#   labs(title="Colours by Treatment", x="Treatment", y="Shades") +
#   theme(axis.text.y=element_text(angle=0, hjust=1)) +
#   facet_wrap(~Treatment, scales = "free_y")




myColors_Colony <- colour_palette$Shades
names(myColors_Colony) <- colour_palette$Cols
myColors_Treatment <- Treat_colors$Shades
names(myColors_Treatment) <- Treat_colors$Treatment
# PERIOD
myColors_period <- rev(hue_pal()(2))
names(myColors_period) <- c("pre", "post")
#### DEFINE THE FILL AND COLOR AS STANDALONE

#COLOR
# geom_point(aes(color = Sample)) +
colScale_Colony <-
  scale_colour_manual(name = "Colony",
                      values = myColors_Colony,
                      drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colScale_Treatment <-
  scale_color_manual(name = "Treatment", values = myColors_Treatment) #for lines
colScale_treatment <-
  scale_color_manual(name = "treatment", values = myColors_Treatment) #for lines
####FILL
# geom_point(aes(color = Sample)) +
colFill_Colony <-
  scale_fill_manual(name = "Colony",
                    values = myColors_Colony,
                    drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colFill_Treatment <-
  scale_fill_manual(name = "Treatment", values = myColors_Treatment) #for lines
colFill_treatment <-
  scale_fill_manual(name = "treatment", values = myColors_Treatment) #for lines

colScale_period <-
  scale_colour_manual(name = "period",
                      values = myColors_period,
                      drop = TRUE)
colFill_period  <-
  scale_fill_manual(name = "period", values = myColors_period) #for lines


#### DEFINE REMAINING PLOT STYLE
# ggplot PLOT STYLE
STYLE <- list(#colScale, fillScale,
  theme_bw(),
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(family = "Liberation Serif")
  ))#,
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long)
  
  STYLE_CONT <- list(#colScale, fillScale,
    theme_bw(),
    theme(
      panel.grid.minor = element_blank(),
      text = element_text(family = "Liberation Serif")
    ))
  
  
  remove <-
    c(
      "color_ramp",
      "color_shades",
      "color_subsets",
      "colors_full",
      "colors_lightGrad",
      "colour_palette",
      "Cols",
      "group",
      "group_colour_palette",
      "group_cols",
      "group_shades",
      "hsl_colors",
      "i",
      "j",
      "lightness_decrease",
      "lightness_increase",
      "meta.data",
      "myColors_Colony",
      "rand_shades",
      "treatment_labs",
      "Shades"
    )
  # cleaning
  rm(list = ls()[which(ls() %in% remove)])
  gc()
  