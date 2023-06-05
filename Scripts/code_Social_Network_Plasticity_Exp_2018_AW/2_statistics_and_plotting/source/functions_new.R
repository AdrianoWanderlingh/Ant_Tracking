clean <- function(){
  rm(list=ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv)%in%to_keep], envir = .GlobalEnv)
  no_print <- gc(verbose=F)
}

GetColorHex <- function(color){
  clr <- col2rgb(color)
  hex_and_col <- sprintf("#%02X%02X%02X %3d %3d %3d", clr[1],clr[2],clr[3], clr[1], clr[2], clr[3])
  hex <- unlist(strsplit(hex_and_col,split=" "))[1]
  return(hex)
}

get_posthoc_groups <- function(model,contrast_matrix,which_levels,dataset){
  levels <- get(which_levels)
  if (is.null(names(levels))){
    names(levels) <- paste("Delta_",levels,sep="")
  }
  
  post_hoc <- summary(glht(model,contrast_matrix),test=adjusted("BH"))
  # print("z value");print(post_hoc$test$tstat);print("Pr>|z|");print(post_hoc$test$pvalues);
  p_values <- as.numeric(post_hoc$test$pvalues);names(p_values) <- names(post_hoc$test$coefficients)
  
  post_hoc_levels <- names(levels)
  post_hoc_mat <- matrix(NA,nrow=length(post_hoc_levels)-1,ncol=length(post_hoc_levels)-1)
  rownames(post_hoc_mat) <- post_hoc_levels[2:length(post_hoc_levels)]
  colnames(post_hoc_mat) <- post_hoc_levels[1:(length(post_hoc_levels)-1)]
  for (i in 1:nrow(post_hoc_mat)){
    for (j in 1:i){
        post_hoc_mat[i,j] <- as.logical(as.numeric(p_values[paste(colnames(post_hoc_mat)[j]," minus ",rownames(post_hoc_mat)[i],sep="")])>0.05)
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
  same <- which(g==1)
  topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
  topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
  g3 <- simplify(graph.data.frame(topology,directed = FALSE))
  #get.data.frame(g3)
  #plot(g3)
  res <- maximal.cliques(g3)
  
  # Reorder given the smallest level
  clique_value <- NULL
  if (which_levels=="level_names"){
    form <- as.formula(paste("variable ~ ", paste(c("treatment","task_group"), collapse= "+")))
  }else{
    form <- as.formula(paste("variable ~ ", paste(gsub("_order","",which_levels), collapse= "+")))
  }
  
  means <- aggregate(form,FUN=mean,data=dataset)
  means$predictor <- names(levels)
  for (i in 1:length(res)){
    clique_value <- c(clique_value,mean(means[as.numeric(unlist(res[[i]])),"variable"]))
  }
  res <- res[order(clique_value)]
  
  # Get group letters
  lab.txt <- vector(mode="list", n)
  lab <- letters[seq(res)]
  for(i in seq(res)){
    for(j in res[[i]]){
      lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
    }
  }
  post_hoc_groups <- unlist(lab.txt); names(post_hoc_groups) <- levels[post_hoc_levels]
  return(post_hoc_groups)
}

log_transf <- function(x){
  if (all(x>0)){
    replac_val <- 0
  }else if (all(x>=0)){
    replac_val <- (min(x[x!=0],na.rm=T))/sqrt(2)
  }else{
    replac_val_1 <- -min(x,na.rm=T)
    y <- x+replac_val_1
    replac_val <- replac_val_1 + (min(y[y!=0],na.rm=T))/sqrt(2)
  }
  return(log10(x+replac_val))
}

sqrt_transf <- function(x){
  if(all(x>=0)){
    replac_val <- 0
  }else{
    replac_val <- -min(x,na.rm=T)
  }
  return(sqrt(x+replac_val))
}

test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
  }else{
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}


##############################################################################
#AW Functions

barplot_delta <- function(dataset=data,predictor,type,collective,plot_untransformed,diff_type){ #form_stat=NULL,
  
  #get difference on RAW data
  diff <- create_diff(dataset,predictor,type,collective,plot_untransformed,diff_type) #form_stat=NULL,
  diff["time"] <- "Delta"

  if (!collective){
    # I SUM BY INDIVIDUAL, THEN MEAN AFTER
    diff <- aggregate(na.rm=T,na.action="na.pass",variable~predictor+colony_size+colony+time+tag,FUN=mean,data=diff)
  }
    
  
  #mean by colony (drop time_of_day)
  diff <- aggregate(na.rm=T,na.action="na.pass",variable~predictor+colony_size+colony+time,FUN=mean,data=diff)
  
    means <- aggregate(na.rm=T,na.action="na.pass",variable~time+predictor,FUN="mean",data=diff);ses <- aggregate(na.rm=T,na.action="na.pass",variable~time+predictor,FUN="std.error",data=diff);
    names(means)[names(means)=="variable"] <- "mean";names(ses)[names(ses)=="variable"] <- "se";means <- merge(means,ses)
    means <- means[order(match(means$predictor,status_order),match(means$time,status_order)),]
    #to_plot <- unique(means[c("time","predictor")])

    

  
  ##### PLOT
  max_mean_se <- max(means$mean + means$se)
  
  plot_var <- ggplot(means, aes_string(x = "predictor", y = "mean", fill = "predictor")) +
    geom_errorbar(aes_string(ymin = paste0("mean - se"), ymax =paste0("mean + se")),
                  position = position_dodge2(width = 0.8, preserve = "single")
    ) +
    geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
    STYLE +
    colFill_treatment +
    labs(
      y = paste("Delta", names(variable_list[i]), sep = " ")
    ) +
    if (variable_list[i] %in% names(post_hoc_outcomes)) {
      geom_text(aes(label = post_hoc_outcomes[[variable_list[i]]][predictor], y= max_mean_se), position = position_dodge2(width = 0.8, preserve = "single"), vjust = -0.4)
      #geom_text(data = post_hoc_outcomes[[variable_list[i]]], aes(x = period, y = 70, group = treatment, label = letters, fontface = "bold"), position = position_dodge2(width = 0.9, preserve = "single"))
    }
  return(plot_var)
}



#AW simplified Stroeymeyt FUNCTIONS
create_diff <- function(dataset,predictor,type,collective,plot_untransformed,diff_type){ #form_stat=NULL,
  if (is.null(dataset)){
    return(NULL)
  }else{
    #################################################################################
    ####plot
    #################################################################################
    dataset$predictor <- dataset[,predictor]
    befores <- dataset[dataset$period=="pre",];afters <-dataset[dataset$period=="post",]
    # print(plot_untransformed)
    if (!plot_untransformed){
      names(befores)[names(befores)=="variable"] <- "variable_before";names(afters)[names(afters)=="variable"] <- "variable_after"
    }else{
      names(befores)[names(befores)=="untransformed_variable"] <- "variable_before";names(afters)[names(afters)=="untransformed_variable"] <- "variable_after"
    }
    # if ((!all(!grepl(":treatment:predictor",as.character(form_stat))))|(!all(!grepl(" \\* treatment \\* predictor",as.character(form_stat))))){
    #   befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+treatment+antid+time_of_day,FUN=mean,data=befores) 
    #   afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+treatment+antid+time_of_day,FUN=mean,data=afters) 
    # }else{
    # if (!grepl("age",root_path)){
    if (collective|type!="individual"){
      befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+size+time_of_day,FUN=mean,data=befores)
      afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+size+time_of_day,FUN=mean,data=afters)
    }else{
        befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+size+tag+time_of_day,FUN=mean,data=befores)
        afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+size+tag+time_of_day,FUN=mean,data=afters)
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
    

    
    diff <- merge(befores,afters,all=T)
    # if (diff_type=="absolute_difference"){
      diff["variable"] <- diff$variable_after-diff$variable_before; # diff[""] <- "diff"
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



########## STYLING FUNCTIONS ###############


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
for (i in c(10, 30, 70, 90)) { # BE CAREFUL: IF CHANGING THIS, ALSO CHANGE Treat_colors
  color_ramp <- colorRampPalette(colors_full[(i-5):(i+5)])
  #color_ramp <- colorRampPalette(colors_full[(i-1):(i+1)])
  color_subsets[[i]] <- color_ramp(12)
}


Treat_colors <- structure(list(Shades = c(colors_full[10],colors_full[30],colors_full[70],colors_full[90]), 
                               Treatment = c("pathogen.big", "control.big", "pathogen.small", "control.small")),
                          row.names = c(NA, 4L), class = "data.frame")

#show_col(c(colors_full[10],colors_full[30],colors_full[70],colors_full[90]))


#clean the output
color_subsets <- color_subsets[lapply(color_subsets,length)>0]

#Darken the colours progressively to increase contrast among the colonies of the same treatment
for (i in 1:4) {
  # Define the color gradient
  color_shades <- color_subsets[[i]]
  
  # Convert the colors from hexadecimal to HSL (hue, saturation, lightness)
  
  colors_lightGrad <- c()
  # Decrease the lightness of each color by an increasing amount
  lightness_decrease <- rev(seq(from = 0, to = 0.2, length.out = length(color_shades)))
  lightness_increase <- seq(from = 0, to = 0.2, length.out = length(color_shades))
  
  for (j in 1:length(color_shades)) {
    hsl_colors <- col2hsl(color_shades[j])
    hsl_colors[3] <- hsl_colors[3] - lightness_decrease[j]
    hsl_colors[3] <- hsl_colors[3] + lightness_increase[j]
    colors_lightGrad <- c(colors_lightGrad,hsl2col(hsl_colors))
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
meta.data <- read.table(paste("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")

# create groups to assign the colours
Cols <- list()
# divide each size_treat into a list element with its colonies inside
for (i in 1:length(unique(meta.data$size_treat))) {
  treatment_labs <- unique(meta.data$size_treat)[i]
  treatment_vector <- unique(meta.data$REP_treat[grepl(treatment_labs, meta.data$REP_treat)])
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
  group_colour_palette <- data.frame(Cols = group_cols, Shades = rand_shades, Treatment = group)
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
names(myColors_period) <- c("pre","post")
#### DEFINE THE FILL AND COLOR AS STANDALONE

#COLOR
# geom_point(aes(color = Sample)) +
colScale_Colony <- scale_colour_manual(name = "Colony", values = myColors_Colony, drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colScale_Treatment <- scale_color_manual(name = "Treatment", values = myColors_Treatment) #for lines
colScale_treatment <- scale_color_manual(name = "treatment", values = myColors_Treatment) #for lines
####FILL
# geom_point(aes(color = Sample)) +
colFill_Colony <- scale_fill_manual(name = "Colony", values = myColors_Colony, drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colFill_Treatment <- scale_fill_manual(name = "Treatment", values = myColors_Treatment) #for lines
colFill_treatment <- scale_fill_manual(name = "treatment", values = myColors_Treatment) #for lines

colScale_period <- scale_colour_manual(name = "period", values = myColors_period, drop = TRUE)
colFill_period  <- scale_fill_manual(name = "period", values = myColors_period) #for lines


#### DEFINE REMAINING PLOT STYLE
# ggplot PLOT STYLE
STYLE <- list(
  #colScale, fillScale,
  theme_bw(),
  theme( panel.grid.minor = element_blank(),text=element_text(family="Liberation Serif")),
  scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

STYLE_CONT <- list(
  #colScale, fillScale,
  theme_bw(),
  theme(panel.grid.minor = element_blank(),text=element_text(family="Liberation Serif"))
)


remove <- c("color_ramp", "color_shades", "color_subsets", "colors_full", 
            "colors_lightGrad", "colour_palette", "Cols",
            "group", "group_colour_palette", 
            "group_cols", "group_shades", "hsl_colors", "i", "j", "lightness_decrease", 
            "lightness_increase", "meta.data", "myColors_Colony", 
            "rand_shades", "treatment_labs", 
            "Shades")
# cleaning
rm(list = ls()[which(ls() %in% remove)])
gc()















