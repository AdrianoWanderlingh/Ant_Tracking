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
           dataset,
           level_names=NULL) {
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


#Box-cox transformation
Box_Cox <- function(x) {
  if (all(x > 0)) {
    replac_val <- 0
  } else if (all(x >= 0)) {
    replac_val <- (min(x[x != 0], na.rm = T)) / sqrt(2)
  } else{
    replac_val_1 <- -min(x, na.rm = T)
    y <- x + replac_val_1
    replac_val <- replac_val_1 + (min(y[y != 0], na.rm = T)) / sqrt(2)
  }
  X <- x + replac_val
  
  bc <- MASS::boxcox(X ~ 1,plotit = FALSE)
  #computes the log-likelihood for a range of lambda values and returns the lambda that maximizes the log-likelihood
  lambda <- bc$x[which.max(bc$y)]
  (X^lambda - 1) / lambda

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

# convert significance levels to stars # AW modded
from_p_to_ptext <- function(pvalue) {
  sapply(pvalue, function(p) {
    if (length(p) == 0) {
      pvaluetext <- "if pval not here, the model failed"
    } else if (is.na(p)) {
      pvaluetext <- "p = NA"
    } else {
      if   (p < 0.0001) {
        pvaluetext  <- "***"
      } else if (p < 0.005) {
        pvaluetext <- "**"
      } else if (p < 0.05) {
        pvaluetext <- "*"
      } else if (p < 0.1) {
        pvaluetext <- paste("p=", sprintf("%.3f", p), sep = "")
      } else if (p < 1) {
        pvaluetext <- paste("p=", sprintf("%.2f", p), sep = "")
      } else {
        pvaluetext <- "p = 1"
      }
    }
    return(pvaluetext)
  })
}

# convert significance levels to rounded vals for text # AW modded
from_p_to_prounded <- function(pvalue) {
  sapply(pvalue, function(p) {
    if (length(p) == 0) {
      pvaluetext <- "if pval not here, the model failed"
    } else if (is.na(p)) {
      pvaluetext <- "p = NA"
    } else {
      if   (p < 0.0001) {
        pvaluetext <- 0.0001
      } else if (p < 0.005) {
        pvaluetext <- 0.005
      } else if (p < 0.05) {
        pvaluetext <- 0.05
      } else if (p < 0.1) {
        pvaluetext <- paste( sprintf("%.3f", p), sep = "")
      } else if (p < 1) {
        pvaluetext <- paste( sprintf("%.2f", p), sep = "")
      } else {
        pvaluetext <- "1"
      }
    }
    return(pvaluetext)
  })
}

closest_match <- function(x,y){
  return(min(which(abs(x-y)==min(abs(x-y))),na.rm=T))
}

plot_arrows <- function(means,plotx,plot_type,LWD,LENGTH,colz=NULL,direction="normal"){
  options(warn=-1)
  if (grepl("points",plot_type)){
    LENGTH <- LENGTH*2
  }
  if (grepl("bars",plot_type)){
    if (direction=="normal"){
      arrows_low <- means$mean-1*as.numeric(sign(means$mean)<=0)*means$se
      arrows_high <- means$mean+1*as.numeric(sign(means$mean)>=0)*means$se
    }else{
      arrows_low <- means$mean-1*as.numeric(sign(means$mean)>=0)*means$se
      arrows_high <- means$mean+1*as.numeric(sign(means$mean)<=0)*means$se
    }
    code1 <- which(arrows_high==means$mean&arrows_low<means$mean)
    code2 <- which(arrows_low==means$mean&arrows_high>means$mean)
    code3 <- which(arrows_low==means$mean&arrows_high==means$mean)
    
    arrows (plotx[code1],arrows_low[code1],plotx[code1],arrows_high[code1],code=1,angle=90,col="black",lwd=LWD,length=LENGTH)
    arrows (plotx[code2],arrows_low[code2],plotx[code2],arrows_high[code2],code=2,angle=90,col="black",lwd=LWD,length=LENGTH)
    arrows (plotx[code3],arrows_low[code3],plotx[code3],arrows_high[code3],code=3,angle=90,col="black",lwd=LWD,length=LENGTH)
  }else{
    arrows_low <- means$mean-means$se
    arrows_high <- means$mean+means$se
    arrows (plotx,arrows_low,plotx,arrows_high,code=3,angle=90,col=colz,lwd=1.5*LWD,length=1.5*LENGTH)
  }
  options(warn=0)
}


VioPlot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                     horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                     lwd = 1, rectCol = "black", colMed = "black", pchMed = 21, bgMed = "white",
                     at, add = FALSE, wex = 1, drawRect = TRUE, mode="Median",cexMed=min_cex){
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)){ 
    at <- 1:n
  }
  std_low <- vector(mode = "numeric", length = n)
  std_high <- vector(mode = "numeric", length = n)
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  meaN <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))){
    args <- c(args, h = h)
  } 
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    meaN[i] <- mean(data)
    std <- std.error(data)
    std_low[i] <- meaN[i]-std
    std_high[i] <- meaN[i]+std
    
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        if (mode=="Median"){
          lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                lty = lty)
          rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
               q3[i], col = rectCol)
          points(at[i], med[i], pch = pchMed, col = colMed, bg=bgMed, cex=cexMed)
        }else if (mode=="Mean"){
          # rect(at[i] - boxwidth/2, std_low[i], at[i] + boxwidth/2, 
          # std_high[i], col = rectCol)
          plot_arrows(means=data.frame(mean=meaN[i],se=std),plotx=at[i],plot_type="violinplot",LWD=line_max,LENGTH=0.05,colz="black")
          points(at[i], meaN[i], pch = pchMed, col = colMed, bg=bgMed, cex=cexMed)
        }
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

scatterplot_violin_forpaper <- function(formula_stat,formula_plot,ylabel,xlabel,title,dat,ymin=NULL,ymax=NULL,xmin=NULL,xmax=NULL,sorting="status",time_point,IC=NULL,output=F,means=T,input_color=NULL,violin_params,point_cex=NULL,predict=NULL){
  violin_params <- as.numeric(unlist(violin_params))
  ##read violin param
  range <- violin_params[1]
  ylim_fac1 <- violin_params[2]
  ylim_fac2 <- violin_params[3]
  wex <- violin_params[4]
  h <- violin_params[5]
  
  # if (all(dat$predictor==dat$colony_size)){
  formula_stat <- update(formula_stat,~.-colony_size)
  # }
  dep_var <- row.names(attr(terms(formula_plot),"factors"))[1]
  pf <- parent.frame()
  dat["variable"] <- eval(parse(text=dep_var),dat,pf)
  
  ##########plotting########
  if (is.numeric(dat$predictor)){
    categories <- sort(unique(dat$predictor_plot))
  }else{
    categories <- unique(c(dat$predictor_plot))
    categories <- categories[order(match(categories,status_order))]
  }
  if (is.null(input_color)){
    if (sorting=="treatment"){
      colour_pal <- NULL
      for (category in categories){
        colour_pal <- c(colour_pal,get(paste(category,"_colour",sep="")))
      }
    }else{
      colour_pal <- rev(colorRampPalette(colour_palette_workers)(length(categories)))
    }
  }else{
    colour_pal <- rev(colorRampPalette(input_color)(length(categories)))
  }
  
  names(colour_pal) <- categories
  dat["colour"] <-  colour_pal[match(dat$predictor_plot,names(colour_pal))]
  forplot <- aggregate(na.rm=T,na.action="na.pass",colour~predictor_plot,FUN=unique,data=dat)
  forplot_med <- aggregate(na.rm=T,na.action="na.pass",variable~predictor_plot,FUN=median,data=dat);names(forplot_med)[names(forplot_med)=="variable"] <- "median"
  forplot_mean <- aggregate(na.rm=T,na.action="na.pass",variable~predictor_plot,FUN=mean,data=dat);names(forplot_mean)[names(forplot_mean)=="variable"] <- "mean"
  forplot <- merge(merge(forplot,forplot_med),forplot_mean)
  if (is.null(ymin)){
    ymin <- min(dat$variable,na.rm=T) - ylim_fac1*(max(dat$variable,na.rm=T)-min(dat$variable,na.rm=T))
    # ymin <- min(dat$variable,na.rm=T)
  }
  if (is.null(ymax)){
    ymax <- max(dat$variable,na.rm=T) + ylim_fac2*(max(dat$variable,na.rm=T)-min(dat$variable,na.rm=T))
    # ymax <- max(dat$variable,na.rm=T) 
  }
  
  ###prepare plot shell
  par_mar_ori <- par()$mar
  if(is.character(dat$predictor)){
    forplot ["pred"] <- as.numeric(factor(forplot$predictor_plot,levels=categories))/2
    dat["pred"] <- as.numeric(factor(dat$predictor_plot,levels=categories))/2
    par(bty="n",xaxt = "n")
  }else{
    forplot ["pred"] <- forplot$predictor_plot
    dat["pred"] <- dat$predictor_plot
    par(bty='l')
  }
  par(mar=par_mar_ori+c(0,0,0,0.5))
  values <- sort(unique(forplot$pred))
  if (is.null(xmin)){
    xlim <- c(min(forplot$pred),max(forplot$pred))+mean(diff(values,lag=1))*c(-0.5,0.5)
  }else{
    xlim <- c(xmin,xmax)
  }
  
  plot(dat$pred,dat$variable,xlab="",ylab="",xaxt="n",yaxt="n",cex.main=inter_cex,cex.lab=inter_cex,font.lab=1,cex.axis=min_cex,xlim=xlim,ylim=c(ymin,ymax),pch=21,type="n")
  if (!all(!grepl("Log",xlabel))){
    where <- axisTicks(c(par("usr")[1],par("usr")[2]),log=F)
    where <- where[which(where==round(where))]
    axis(1,at=where,labels=format(10^(where),scientific=T),cex.lab=inter_cex,cex.axis=min_cex)
    xlab <- substr(gsub("Log\\(","",xlabel),1,-1+nchar(gsub("Log\\(","",xlabel)))
    if (xlab=="Measured pathogen load"){
      title(xlab=expression(paste("Measured pathogen load (ng/", mu, "L)")),cex.lab=inter_cex,mgp=par()$mgp+c(0.2,0,0))
    }else{
      title(xlab=xlab,cex.lab=inter_cex,mgp=par()$mgp+c(0.2,0,0))
    }
  }else{
    axis(1,cex.lab=inter_cex,cex.axis=min_cex)
    title(xlab=xlabel,cex.lab=inter_cex)
  }
  
  if(all(!grepl("log",ylabel))){
    axis(2,cex.lab=inter_cex,cex.axis=min_cex)
    title(ylab=ylabel,cex.lab=inter_cex)
  }else{
    where <- axisTicks(c(par("usr")[3],par("usr")[4]),log=F)
    where <- where[which(where==round(where))]
    axis(2,at=where,labels=format(10^(where),scientific=T),cex.lab=inter_cex,cex.axis=min_cex)
    ylab <- as.character(ylabel[2])
    if (ylab=="Measured pathogen load"){
      title(ylab=expression(paste("Measured pathogen load (ng/", mu, "L)")),cex.lab=inter_cex,mgp=par()$mgp+c(0.1,0,0))
    }else{
      title(ylab=ylab,cex.lab=inter_cex,mgp=par()$mgp+c(0.1,0,0))
    }
    
    
  }
  
  if(is.character(dat$predictor)){
    par(xaxt = "s")
    axis(side=1,at=sort(unique(forplot$pred)),labels=full_statuses_names[categories],tick=F,lty=0,cex.axis=inter_cex)
    ###plus add grid lines
    par(xpd=F)
    # for (idx in 1:nrow(forplot)){
    #   abline(h=forplot[idx,"median"],lwd=line_min,lty=3,col="black")
    # }
    abline(h=0,lwd=line_min,lty=3,col="black")
  }
  if(!is.null(title)){title(main=title,cex.main=inter_cex,line=2.5)}
  ###if input, add 95%IC
  if(!is.null(IC)){
    polygon(c(par("usr")[c(1,2)],par("usr")[c(2,1)]),c(IC[1],IC[1],IC[2],IC[2]),border=NA,col=alpha("orange",0.1))
  }
  par_cex_ori <- par()$cex
  par(cex=0.3)
  if (is.null(point_cex)){
    cexMed <- min_cex
  }else{
    cexMed <- point_cex
  }
  ###add violins
  for (i in 1:nrow(forplot)){
    subset <- dat[which(dat$pred==forplot[i,"pred"]),"variable"]
    if (is.na(range)){
      VioPlot(na.omit(subset),col=alpha(forplot[i,"colour"],0.7), horizontal=F, at=forplot[i,"pred"], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Median",cexMed=cexMed)
    }else{
      VioPlot(na.omit(subset),range=range, h=h,col=alpha(forplot[i,"colour"],0.7), horizontal=F, at=forplot[i,"pred"], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Median",cexMed=cexMed)
    }
  }
  par(cex=par_cex_ori)
  
  
  ##########stats############
  ###initialise
  to_plot <- data.frame(intercept=as.numeric(as.character()),slope=as.numeric(as.character()),colour=as.character(),stringsAsFactors=F)
  to_mtext <- data.frame(effect=as.character(),pvalue=as.numeric(as.character()),stringsAsFactors=F)
  pred<- attr(terms(formula_stat),"term.labels")
  pred <- pred[!grepl("\\|",pred)]
  
  ###get final stats formula (step by step in case of a multiple terms model); i.e., reduce model
  try(model_temp <- lmer(formula_stat,data=dat),silent=T)
  if (exists("model_temp")){
    coeff <- Anova(model_temp,type="III")
    if ("Pr(>Chisq)"%in%colnames(coeff)){
      ###first test if colony size is significant; if not remove it
      if ("colony_size" %in%pred){
        if (coeff["colony_size","Pr(>Chisq)"]>0.05){
          formula_stat <- update(formula_stat,~.-colony_size)
          model_temp <- lmer(formula_stat,data=dat)
          coeff <- Anova(model_temp,type="III")
        }#if (coeff["colony_size","Pr(>Chisq)"]>0.05)
        ###and now that it is dealt with, remove it from pred
        pred <- pred[pred!="colony_size"]
      }#("colony_size" %in%pred)
      
      if (time_point=="comparison"){
        for (effect in pred){
          to_mtext <- rbind(to_mtext,data.frame(effect=effect,pvalue=coeff[effect,"Pr(>Chisq)"],stringsAsFactors=F))
        }
        interaction_effect <- pred[grepl("\\:",pred)]
        p_interaction <- coeff[interaction_effect,"Pr(>Chisq)"]
        if (p_interaction>0.05){
          if (sorting=="status"){
            formula_stat <- update(formula_stat,~.-predictor:status)
          }
          if (sorting=="period"){
            formula_stat <- update(formula_stat,~.-predictor:period)
          }
          pred<- attr(terms(formula_stat),"term.labels")
          pred <- pred[!grepl("\\|",pred)&!grepl("colony_size",pred)]
          try(model_bis <- lmer(formula_stat,data=dat),silent=T)
          if (exists("model_bis")){
            coeff <- Anova(model_bis,type="III")
            pvalue <-   coeff[pred[!grepl("predictor",pred)],"Pr(>Chisq)"]
            to_mtext[to_mtext$effect==pred[!grepl("predictor",pred)],"pvalue"] <- pvalue
            if (pvalue>0.05){
              if (sorting=="status"){
                formula_stat <- update(formula_stat,~.-status)
              }
              if (sorting=="period"){
                formula_stat <- update(formula_stat,~.-period)
              }
              pred<- attr(terms(formula_stat),"term.labels")
              pred <- pred[!grepl("\\|",pred)&!grepl("colony_size",pred)]
            }
          }
        }#if (p_interaction>0.05)
      }#if (time_point=="comparison")
    }#if ("Pr(>Chisq)"%in%colnames(coeff))
    
    rm(list=c("model_temp"))
  }#(exists("model_temp"))
  
  
  ####get final pvalues based on final model
  ########get names of the predictors in the table for extracting coefficients
  formula_simple <- update(formula_stat,~.-(1|colony)-(1|antid_1)-(1|antid_2)-(1|antid))
  pred2 <- Names(  formula_simple,dat);pred2 <- pred2[!grepl("Intercept",pred2)&!grepl("colony_size",pred2)];pred2 <- pred2[!grepl("\\|",pred2)]
  try(model_final <- lmer(formula_stat,data=dat),silent=T)
  #########testing normality of residuals
  resids <- residuals(model_final)
  test_norm(resids)
  
  if (exists("model_final")){
    coeff <- Anova(model_final,type="III")
    print(coeff)
    coeff2 <- summary(model_final)$coefficients
    if ("Pr(>Chisq)"%in%colnames(coeff)){
      if (length(pred)==1){
        pvalue <- coeff[pred,"Pr(>Chisq)"]
        to_output <- list(coeff2[pred2,"Estimate"])
        names(to_output) <- time_point
        if (nrow(to_mtext)==0){
          to_mtext <- rbind(to_mtext,data.frame(effect=pred,pvalue=pvalue,stringsAsFactors=F))
        }else{
          to_mtext[to_mtext$effect==pred,"pvalue"] <- pvalue
        }
        if ( is.numeric(dat$predictor)){
          if (pvalue < 0.05){
            if ("colony_size"%in%Names(formula_simple,dat)){
              to_plot <- rbind(to_plot,data.frame(
                intercept = coeff2["(Intercept)","Estimate"]+coeff2["colony_size","Estimate"]*mean(dat$colony_size),
                slope = coeff2[pred2,"Estimate"],
                #colour = statuses_colours[time_point],stringsAsFactors=F))
                colour = "black",stringsAsFactors=F))
              
            }else{
              to_plot <- rbind(to_plot,data.frame(intercept = coeff2["(Intercept)","Estimate"],slope = coeff2[pred2,"Estimate"],
                                                  #colour = statuses_colours[time_point],stringsAsFactors=F))
                                                  colour = "black",stringsAsFactors=F))
            }
          }#if (pvalue < 0.05)         
        }
      }#if (length(pred)==1)
    }#if ("Pr(>Chisq)"%in%colnames(coeff))
  }#if (exists("model_final"))
  
  ###plot ablines
  if (!is.null(predict)){
    predicted_value <- to_plot[1,"intercept"] + to_plot[1,"slope"]*predict
    segments(x0=predict,y0=ymin-0.1*(ymax-ymin),y1=predicted_value,col="springgreen2",xpd=F,lty=2)
    segments(x0=xlim[1]-0.1*(xlim[2]-xlim[1]),y0=predicted_value,x1=predict,col="red",xpd=F,lty=2)
  }
  par(xpd=F)
  if (nrow(to_plot)>=1){
    for (i in 1:nrow(to_plot)){
      abline(a=to_plot[i,"intercept"],b=to_plot[i,"slope"],col=to_plot[i,"colour"],lwd=line_inter)
    }# i
  }#(nrow(to_plot)>=1)
  ###plot mtext
  if (nrow(to_mtext)>=1){
    for (i in 1:nrow(to_mtext)){
      pvalue <- to_mtext[i,"pvalue"];effect <- to_mtext[i,"effect"]
      if(grepl("\\:",effect)){effect_ <- "Interaction: "}
      if(!grepl(sorting,effect)){
        effect_ <- paste(ylabel[1],": ",sep="")
        if (nchar(effect_)>30){
          effect_ <- "main: "
        }
      }
      if (sorting=="status"){
        if(!grepl("predictor",effect)){effect_ <- "treated vs. nestmates: "}
      }
      if (sorting=="period"){
        if(!grepl("predictor",effect)){effect_ <- "before vs. after: "}
      }
      p_value <- from_p_to_ptext(pvalue)
      if (pvalue>0.05){p_cex <- inter_cex}else{p_cex <- max_cex}
      par(xpd=T)
      if (nrow(to_mtext)==1){
        title(main=p_value,cex.main=p_cex,font.main=2,line=stat_line-((i-1)*0.75),xpd=T)
      }else{
        title(main=paste(effect_,p_value,sep=""),cex.main=p_cex,font.main=2,line=stat_line-((i-1)*0.75),xpd=T)
      }
      par(xpd=T)
    } # i
  }#if (nrow(to_mtext)>=1)
  if(output){return(to_output)}
  
  par(mar=par_mar_ori)
  if(!is.null(predict)){return(predicted_value)}else{return(predict)}
}#scatterplot_bubbles_qpcr

plot_regression <- function(data,time_point,analysis,n_cat_horiz,n_cat_vertic,pool=F,prepare=F,status=NULL,collective=NULL,pool_plot=F,input_color=NULL,plot_untransformed=F,boldy=F,aligned=F,ymin=NULL,ymax=NULL,xmin=NULL,xmax=NULL,point_cex=NULL,adjust_title_line=0,predict=NULL){
  adjust_title_line_ori <- adjust_title_line
  data_ori <- data
  #### plot regression for each desired combination of variable and predictor
  for (i in 1:length(analysis[["variable_list"]])){
    data <- data_ori
    adjust_title_line <- adjust_title_line_ori
    ####get variable and predictor
    variable <- analysis[["variable_list"]][i]
    
    ####if necessary: convert pixels to mm
    if (grepl("changetomm",names(variable))){
      if (grepl("changetomm2",names(variable))){
        data[,variable] <- data[,variable]*pix_to_mm_ratio*pix_to_mm_ratio
        names(variable) <- gsub("_changetomm2","",names(variable))
      }else{
        data[,variable] <- data[,variable]*pix_to_mm_ratio
        names(variable) <- gsub("_changetomm","",names(variable))
      }
    }
    
    print(variable)
    
    predictor <- analysis[["predictor_list"]][i]
    transf_variable <- analysis[["transf_variable_list"]][i]
    transf_predictor <- analysis[["transf_predictor_list"]][i]
    pooli <- pool[i]
    
    ####if necessary: include queen and treated into predictor function
    if (!is.null(predictor)){
      if ((!collective&refine!=""&predictor!="")&("tag"%in%names(data))){
        ####first change the content of predictor column 
        data["predictor"] <- data[,predictor]
        ###second add treated
        data[which(data$status=="treated"),"predictor"] <- "treated"
        ####fourth copy back into predictor column
        data[,predictor] <- data[,"predictor"]
        ####fifth if necessary remove queen
        if (length(unique(data$predictor))>1){
          if (!queen){
            data <- data[which(data$task_group!="queen"),]
          }
          if (!treated){
            data <- data[which(data$predictor!="treated"),]
          }
          if (!nurses){
            data <- data[which(data$predictor!="nurse"),]
          }
          if (!foragers){
            data <- data[which(data$predictor!="forager"),]
          }
        }
      }else if (predictor!=""){
        data["predictor"] <- data[,predictor]
      }  
    }
    
    ####if necessary: apply prepare dataset function
    if (prepare){
      data <- prepare_dataset(data,variable)
    }else{
      ####process variable
      data["variable"] <- data[,variable]
      data[which(!is.finite(data$variable)),"variable"] <- NA
    }
    data["untransformed_variable"] <- data$variable
    ####transform variable
    ylabel <- names(variable)
    ylabel <- capitalize(ylabel)
    if (transf_variable=="log"){
      print("Logging variable...")
      data[!is.na(data$variable),"variable"] <- log_transf(data[!is.na(data$variable),"variable"] )
      if (!plot_untransformed){
        
        if(boldy){
          ylabel <- substitute(paste(bold(ylabel),bolditalic(" (log)")),list(ylabel=ylabel))
          adjust_title_line <- 0.17
        }else{
          ylabel <- substitute(paste(ylabel,italic(" (log)")),list(ylabel=ylabel))
          adjust_title_line <- 0.17
        }
        
      }
    }else if (grepl("power",transf_variable)){
      data[!is.na(data$variable),"variable"]  <- (data[!is.na(data$variable),"variable"] )^as.numeric(gsub("power","",transf_variable))
      if (!plot_untransformed){
        
        if(boldy){
          ylabel <- substitute(paste(bold(ylabel),bolditalic(" (") ^pow,bolditalic(")")),list(ylabel=ylabel,pow=as.numeric(gsub("power","",transf_variable))))
          adjust_title_line <- 0
        }else{
          ylabel <- substitute(paste(ylabel,italic(" (") ^pow,italic(")")),list(ylabel=ylabel,pow=as.numeric(gsub("power","",transf_variable))))
          adjust_title_line <- 0
        }
        
      }
    }else if (transf_variable=="sqrt"){
      data[!is.na(data$variable),"variable"]  <- sqrt_transf(data[!is.na(data$variable),"variable"] )
      if (!plot_untransformed){
        
        if (boldy){
          ylabel <- substitute(paste(bold(ylabel),bolditalic(" ("),sqrt(bolditalic(")"))),list(ylabel=ylabel))
          adjust_title_line <- 0.24
        }else{
          ylabel <- substitute(paste(ylabel,italic(" ("),sqrt(italic(")"))),list(ylabel=ylabel))
          adjust_title_line <- 0.24
        }
      }
      
    }else if (transf_variable=="Box_Cox"){
      data[!is.na(data$variable),"variable"]  <- Box_Cox(data[!is.na(data$variable),"variable"] )
      if (!plot_untransformed){
        
        if (boldy){
          ylabel <- substitute(paste(bold(ylabel),bolditalic(" ("),sqrt(bolditalic(")"))),list(ylabel=ylabel))
          adjust_title_line <- 0.24
        }else{
          ylabel <- substitute(paste(ylabel,italic(" ("),sqrt(italic(")"))),list(ylabel=ylabel))
          adjust_title_line <- 0.24
        }
      }
      
    }
    ####for comparison: use perform_analysis_combined_function
    if (time_point=="comparison"){
      if ("randy"%in%names(data)){
        case <- "case3"
      }else{
        
        if (is.null(predictor)){
          case <- "case1"
        }else if (length(unique(data$predictor))==1){
          case <- "case1"
        }else if (!((!collective&refine!=""&predictor!=""))){
          case <- "case1"
        }else{
          case <- "case2"
        }
        
      }
      if (case=="case1"){
        perform_barplot_analysis(root_path=root_path,collective=collective,dataset=data,lab_title=ylabel,excluded=NULL,pool=pool,violin_params=analysis[["violin_plot_param"]][i],pool_plot=pool_plot,adjust_title_line=adjust_title_line)
      }else if (case=="case2"){
        perform_barplot_analysis_refined(root_path=root_path,collective=collective,dataset=data,lab_title=ylabel,excluded=NULL,pool=pool,violin_params=analysis[["violin_plot_param"]][i],pool_plot=pool_plot,plot_untransformed=plot_untransformed,aligned=aligned,adjust_title_line=adjust_title_line)
      }else{
        perform_barplot_analysis_simple(root_path=root_path,collective=collective,dataset=data,lab_title=ylabel,excluded=NULL,pool=pool,violin_params=analysis[["violin_plot_param"]][i],pool_plot=pool_plot,plot_untransformed=plot_untransformed,aligned=aligned,adjust_title_line=adjust_title_line)
      }
    }else{
      ####process predictor
      data["predictor"] <- data[,predictor]
      if (is.numeric(data$predictor)){
        data[which(!is.finite(data$predictor)),"predictor"] <- NA
        if (transf_predictor=="log"){
          if (!is.null(predict)){
            transfi <- log_transf(c(data[!is.na(data$predictor),"predictor"] , predict))
            predict <- transfi[length(transfi)]
            data[!is.na(data$predictor),"predictor"]  <- transfi[1:(length(transfi)-1)]
          }else{
            data[!is.na(data$predictor),"predictor"]  <- log_transf(data[!is.na(data$predictor),"predictor"] )
          }
          xlabel<- paste("Log(", names(predictor),")",sep="")
        }else if (transf_predictor=="power2"){
          data[!is.na(data$predictor),"predictor"]  <- (data[!is.na(data$predictor),"predictor"] )^2
          if (!is.null(predict)){predict <- predict^2}
          xlabel <- substitute( xlabely ^2,list(xlabely=names(predictor)))
        }else if (transf_predictor=="sqrt"){
          if (!is.null(predict)){
            transfi <- sqrt_transf(c(data[!is.na(data$predictor),"predictor"] , predict))
            predict <- transfi[length(transfi)]
            data[!is.na(data$predictor),"predictor"]  <- transfi[1:(length(transfi)-1)]
          }else{
            data[!is.na(data$predictor),"predictor"]  <- sqrt_transf(data[!is.na(data$predictor),"predictor"] )
          }
          xlabel <- substitute(sqrt ( xlabely),list(xlabely=names(predictor)))
        }else{
          xlabel <- names(predictor)
        }
      }else{
        xlabel <- ""
      }
      
      title <- ""
      if(!"period"%in%names(data)){data["period"] <- data$time}
      
      ###process pool argument
      if (pooli){
        dat <- aggregate(na.rm=T,na.action="na.pass",variable~predictor+colony+antid+status+colony_size+period,FUN="mean",data=data)
      }else{
        dat <- data
      }
      
      ###process horizontal bin argument
      if (length(unique(dat[,"predictor"]))>n_cat_horiz){
        dat[which(!is.na(dat[,"predictor"])),"predictor_plot"] <- as.numeric(gsub("\\(","",unlist(strsplit(as.character(cut(dat$predictor,n_cat_horiz,include_lowest=T)),split=","))[grepl("\\(",unlist(strsplit(as.character(cut(dat$predictor,n_cat_horiz,include_lowest=T)),split=",")))]))
      }else{#if (length(unique(dat[,"predictor"]))>n_cat_horiz)
        dat[,"predictor_plot"] <- dat[,"predictor"]
      }##else
      
      ###define formula
      formula_stat <- as.formula(paste("variable"," ~ ", paste(c("predictor","colony_size","(1|colony)","(1|antid)"), collapse= "+")))
      if (length(unique(na.omit(aggregate(variable~antid,FUN=length,data=dat)$variable)))==1&unique(na.omit(aggregate(variable~antid,FUN=length,data=dat)$variable))[1]==1){
        formula_stat <- update(formula_stat,~.-(1|antid))
      }
      formula_plot <- as.formula(paste("variable"," ~ ", paste(c("predictor_plot","period"), collapse= "+")))
      ###plot
      predicted_value <- scatterplot_violin_forpaper(formula_stat=formula_stat,formula_plot=formula_plot,ylabel=ylabel,xlabel=xlabel,title=title,dat=dat,sorting="period",time_point=time_point,output=F,violin_params = analysis[["violin_plot_param"]][i],input_color=input_color,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,point_cex=point_cex,predict=predict)
      if (!is.null(predicted_value)){
        predicted_value <- data[,variable][closest_match(predicted_value,data$variable) ]
        return(predicted_value)
      } 
    }
  }
}



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#AW adapted Functions

model_n_interations <- lmerControl(optCtrl = list(maxfun = 2e5))

#this function performs no rescaling but has if-statements to skip time_of_day for simulation data
collective_analysis_no_rescal <- function(data_path=data_path,showPlot=T){
  
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
  # stats_outcomes1    <- NULL
  # post_hoc_outcomes1 <- NULL
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
    }else if (transf_variable_list[i]=="Box_Cox"){
      data[!is.na(data$untransformed_variable),"variable"]  <- Box_Cox(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
    }else if (transf_variable_list[i]=="none"){
      data$variable <- data$untransformed_variable
    }
    
    # Check if 'before' or 'after' is present in the 'period' column
    data$period <- ifelse(grepl("before", data$period), "pre", 
                          ifelse(grepl("after", data$period), "post", data$period))
    
    
    ###statistics
    ###make sure treatment, exposure, size and period are factors with levels in the desired order
    data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
    data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
    data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
    data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
    
    ###fit model - using treatment
    if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
      model <- lmer(   variable ~ period*treatment + (1|colony) ,data=data)
    }else{
      model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) ,data=data)
    }
    anov  <- anova(model)
    p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
    
    # Check if the model converged
    if(length(summary(model)$optinfo$conv$lme4$messages) != 0) {
      # If the model fails to converge, calculate VIF
      # The Variance Inflation Factor (VIF) is a measure of multicollinearity. A rule of thumb is that if the VIF is greater than 5, then the explanatory variables are highly correlated, which can affect the stability and interpretability of the model
      vif_values <- car::vif(model)
      
      # If any VIF is > 10, assign p_interaction_treatment = 1, and try simpler model 
      if(any(vif_values[, "GVIF"] > 10)) {
        p_interaction_treatment <- 1
        warning("model multicollinearity detected, try simpler model")
      }
    }
    
    stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),Fvalue=anov["period:treatment","F value"],pval=p_interaction_treatment,stringsAsFactors = F))
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
      if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
        model <- lmer(   variable ~ period*exposure + period*size + (1|colony) ,data=data,control = model_n_interations)
      }else{
        model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) ,data=data,control = model_n_interations)
      }
      
      anov  <- anova(model)
      p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
      p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
      
      if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
        
        
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
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
          
          ###fit model
          if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
            model <- lmer(   variable ~ period*exposure + (1|colony) ,data=data,control = model_n_interations)
          }else{
            model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) ,data=data,control = model_n_interations)
          }
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
          
          # print( "Interaction between exposure and period is not significant:")
          # print(anov["period:exposure",])
          # 
          ###fit model
          if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
            model <- lmer(   variable ~ period*size + (1|colony) ,data=data ,control = model_n_interations)
          }else{
            model <- lmer(   variable ~ period*size + (1|time_of_day) + (1|colony) ,data=data ,control = model_n_interations)
          }
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
    if (showPlot) {print(barplot_delta_period)}
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period
    #names(barplot_delta_period_list[[i]]) <- variable_list[i]
    
    # # add cleaning step
    # post_hoc_outcomes1 <- c(post_hoc_outcomes1,post_hoc_outcomes)
    # stats_outcomes1    <- c(stats_outcomes1,stats_outcomes)
    # 
    # rm(list = ls()[which(names(ls()) == "post_hoc_outcomes")])
    # rm(list = ls()[which(names(ls()) == "stats_outcomes")])
  }
  rownames(stats_outcomes) <- NULL
  
  #add  formatted output
  stats_outcomes$formatted <- paste0("(GLMM,treatment-induced changes (", stats_outcomes$variable,", ",stats_outcomes$predictor,"): F(",stats_outcomes$df,") = ", round(stats_outcomes$Fvalue,3), ", P ≤ ",format(from_p_to_prounded(stats_outcomes$pval),scientific=F))
  
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
}

#this function handles rescaling by pre-exposure mean for each size
collective_analysis_rescal <- function(data_path=data_path,showPlot=T){
  
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
  
  print(names(data))
  
  ###2. Loop over variables
  data_ori <- data
  stats_outcomes    <- NULL
  post_hoc_outcomes <- NULL
  # stats_outcomes1    <- NULL
  # post_hoc_outcomes1 <- NULL
  barplot_delta_period_list <- list()
  
  for (i in 1:length(variable_list)){
    # for (i in c(1)){
    print(paste0("######## ",variable_list[i]," ########"))
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
    }else if (transf_variable_list[i]=="Box_Cox"){
      data[!is.na(data$untransformed_variable),"variable"]  <- Box_Cox(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
    }else if (transf_variable_list[i]=="none"){
      data$variable <- data$untransformed_variable
    }
    
    # Check if 'before' or 'after' is present in the 'period' column
    data$period <- ifelse(grepl("before", data$period), "pre", 
                          ifelse(grepl("after", data$period), "post", data$period))
    
    ###statistics
    ###make sure treatment, exposure, size and period are factors with levels in the desired order
    data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
    data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
    data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
    data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
    
    ###fit model - using treatment
    model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) ,data=data,control = model_n_interations)
    anov  <- anova(model)
    p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
    
    # Check if the model converged
    if(length(summary(model)$optinfo$conv$lme4$messages) != 0) {
      # If the model fails to converge, calculate VIF
      # The Variance Inflation Factor (VIF) is a measure of multicollinearity. A rule of thumb is that if the VIF is greater than 5, then the explanatory variables are highly correlated, which can affect the stability and interpretability of the model
      vif_values <- car::vif(model)
      
      # If any VIF is > 10, assign p_interaction_treatment = 1, and try simpler model 
      if(any(vif_values[, "GVIF"] > 10)) {
        p_interaction_treatment <- 1
        warning("model multicollinearity detected, try simpler model")
      }
    }
    
    stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),Fvalue=anov["period:treatment","F value"],pval=p_interaction_treatment,stringsAsFactors = F))
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
      model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) ,data=data ,control = model_n_interations)
      anov  <- anova(model)
      p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
      p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
      
      if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
        
        
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
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
          
          ###fit model 
          model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) ,data=data ,control = model_n_interations)
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
          
          # print( "Interaction between exposure and period is not significant:")
          # print(anov["period:exposure",])
          # 
          ###fit model 
          model <- lmer(   variable ~ period*size + (1|time_of_day) + (1|colony) ,data=data ,control = model_n_interations)
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
    if (showPlot) {print(barplot_delta_period)}
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period
  
    # # add cleaning step
    # post_hoc_outcomes1 <- c(post_hoc_outcomes1,post_hoc_outcomes)
    # stats_outcomes1    <- c(stats_outcomes1,stats_outcomes)
    # 
    # rm(list = ls()[which(names(ls()) == "post_hoc_outcomes")])
    # rm(list = ls()[which(names(ls()) == "stats_outcomes")])
    
    }
  rownames(stats_outcomes) <- NULL
  
  #add  formatted output
  stats_outcomes$formatted <- paste0("(GLMM,treatment-induced changes (", stats_outcomes$variable,", ",stats_outcomes$predictor,"): F(",stats_outcomes$df,") = ", round(stats_outcomes$Fvalue,3), ", P ≤ ",format(from_p_to_prounded(stats_outcomes$pval),scientific=F))
  
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
}

individual_ONE_analysis <- function(data_path=data_path,which_individuals,showPlot=T){
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
  
  print(names(data))
  
  
  ##2a. Extract exposure and size from treatment column
  data$exposure <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[1]  ))
  data$size     <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))
  
  ##2b. add information on task group
  if (pattern %in% c("individual_behavioural_data", "individual_simulation_results_observed")) {
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
  
  print(paste("",which_individuals,"",sep=" xxxxxxxxxx "))
  
  ##2d.  ###add a unique antid column
  data <- within(data,antID <- paste(colony,tag,sep="_"))
  
  ###2. Loop over variables
  data_ori <- data
  stats_outcomes    <- NULL
  post_hoc_outcomes <- NULL
  # stats_outcomes1    <- NULL
  # post_hoc_outcomes1 <- NULL
  barplot_delta_period_list <- list()
  
  for (i in 1:length(variable_list)){
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
    }else if (transf_variable_list[i]=="Box_Cox"){
      data[!is.na(data$untransformed_variable),"variable"]  <- Box_Cox(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
    }else if (transf_variable_list[i]=="none"){
      data$variable <- data$untransformed_variable
    }
    
    
    # Check if 'before' or 'after' is present in the 'period' column
    data$period <- ifelse(grepl("before", data$period), "pre", 
                          ifelse(grepl("after", data$period), "post", data$period))
    
    #hist(data$variable, breaks = 100)
    ###statistics
    ###make sure treatment, exposure, size and period are factors with levels in the desired order
    data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
    data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
    data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
    data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
    data$antID     <- factor( data$antID )
    
    ###fit model - using treatment
    if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
      model <- lmer(   variable ~ period*treatment + (1|colony) + (1|antID) ,data=data ,control = model_n_interations)
    }else{
      model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) + (1|antID) ,data=data ,control = model_n_interations)
    }
    
    anov  <- anova(model)
    p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
    
    # Check if the model converged
    if(length(summary(model)$optinfo$conv$lme4$messages) != 0) {
      # If the model fails to converge, calculate VIF
      # The Variance Inflation Factor (VIF) is a measure of multicollinearity. A rule of thumb is that if the VIF is greater than 5, then the explanatory variables are highly correlated, which can affect the stability and interpretability of the model
      vif_values <- car::vif(model)
      
      # If any VIF is > 10, assign p_interaction_treatment = 1, and try simpler model 
      if(any(vif_values[, "GVIF"] > 10)) {
        p_interaction_treatment <- 1
      warning("model multicollinearity detected, try simpler model")
      }
    }

    stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),Fvalue=anov["period:treatment","F value"],pval=p_interaction_treatment,stringsAsFactors = F))
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
      if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
        model <- lmer(   variable ~ period*exposure + period*size + (1|colony) + (1|antID),data=data ,control = model_n_interations)
      }else{
        model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) + (1|antID),data=data ,control = model_n_interations)
      }
      
      anov  <- anova(model)
      p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
      p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
      
      if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
        stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
        
        
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
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
          
          ###fit model
          if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
            model <- lmer(   variable ~ period*exposure + (1|colony) + (1|antID),data=data ,control = model_n_interations)
          }else{
            model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) + (1|antID),data=data ,control = model_n_interations)
          }
          
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
          
          # print( "Interaction between exposure and period is not significant:")
          # print(anov["period:exposure",])
          # 
          ###fit model 
          if (grepl("simulation_results",pattern)) { #simulation data won't have time_of_day
            model <- lmer(   variable ~ period*size + (1|colony) + (1|antID),data=data ,control = model_n_interations)
          }else{
            model <- lmer(   variable ~ period*size + (1|colony) + (1|antID),data=data ,control = model_n_interations)
          }
          
          anov  <- anova(model)
          test_norm(residuals(model))
          for (rowi in 1:nrow(anov)){
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
    if (showPlot) {print(barplot_delta_period)}
    
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period 
    #ADD SAVING OF THE PLOT on single pdf file as done in plot_grooming file
    
    # # add cleaning step
    # post_hoc_outcomes1 <- c(post_hoc_outcomes1,post_hoc_outcomes)
    # stats_outcomes1    <- c(stats_outcomes1,stats_outcomes)
    # 
    # rm(list = ls()[which(names(ls()) == "post_hoc_outcomes")])
    # rm(list = ls()[which(names(ls()) == "stats_outcomes")])
  }
  rownames(stats_outcomes) <- NULL
  
  #add  formatted output
  stats_outcomes$formatted <- paste0("(GLMM,treatment-induced changes (", stats_outcomes$variable,", ",stats_outcomes$predictor,"): F(",stats_outcomes$df,") = ", round(stats_outcomes$Fvalue,3), ", P ≤ ",format(from_p_to_prounded(stats_outcomes$pval),scientific=F))
  
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
  
}

# this function currently does not handle grepl("simulation_results",pattern) for models without time_of_day
individual_TWO_analysis <- function(data_path=data_path,which_individuals,showPlot=T){
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
  
  print(names(data))
  
  ##2a. Extract exposure and size from treatment column
  data$exposure <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[1]  ))
  data$size     <- unlist(lapply( data$treatment, function(x)  unlist(strsplit(x,split="\\.") )[2]  ))
  
  ##2b. add information on task group
  if (pattern %in% c("individual_behavioural_data", "individual_simulation_results_observed")) {
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
  # stats_outcomes1    <- NULL
  # post_hoc_outcomes1 <- NULL
  barplot_delta_period_list <- list()
  
  for (i in 1:length(variable_list)){
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
    }else if (transf_variable_list[i]=="Box_Cox"){
      data[!is.na(data$untransformed_variable),"variable"]  <- Box_Cox(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
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
    model <- lmer(   variable ~ period*treatment*task_group + (1|time_of_day) + (1|colony) + (1|antID) ,data=data ,control = model_n_interations)
    anov  <- anova(model)
    p_interaction_treatment_group <- anov["period:treatment:task_group","Pr(>F)"]
    
    # Check if the model converged
    if(length(summary(model)$optinfo$conv$lme4$messages) != 0) {
      # If the model fails to converge, calculate VIF
      # The Variance Inflation Factor (VIF) is a measure of multicollinearity. A rule of thumb is that if the VIF is greater than 5, then the explanatory variables are highly correlated, which can affect the stability and interpretability of the model
      vif_values <- car::vif(model)
      
      # If any VIF is > 10, assign p_interaction_treatment = 1, and try simpler model 
      if(any(vif_values[, "GVIF"] > 10)) {
        p_interaction_treatment <- 1
        warning("model multicollinearity detected, try simpler model")
      }
    }
    
    stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment:task_group",df=paste(round(anov["period:treatment:task_group","NumDF"]),round(anov["period:treatment:task_group","DenDF"]),sep=","),Fvalue=anov["period:treatment:task_group","F value"],pval=anov["period:treatment:task_group","Pr(>F)"],stringsAsFactors = F))
    
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
        for (j in 1:length(level_names)){
          rownames(contrast_matrix) <- gsub(j,level_names[j],rownames(contrast_matrix))
        }
        
      }else{print("TOO MANY TASK GROUP LEVELS _ CANNOT COPE")}
      
      
      posthoc_groups_treatment_task_groups        <- list(get_posthoc_groups(model=model,contrast_matrix=contrast_matrix,which_levels="level_names",dataset=data, level_names=level_names))
      names(posthoc_groups_treatment_task_groups) <- variable_list[i]
      post_hoc_outcomes                           <- c(post_hoc_outcomes,posthoc_groups_treatment_task_groups)
      
    }else{
      
      ###fit model - using treatment
      model <- lmer(   variable ~ period*treatment + (1|time_of_day) + (1|colony) + (1|antID) ,data=data ,control = model_n_interations)
      anov  <- anova(model)
      p_interaction_treatment <- anov["period:treatment","Pr(>F)"]
      
      stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:treatment",df=paste(round(anov["period:treatment","NumDF"]),round(anov["period:treatment","DenDF"]),sep=","),Fvalue=anov["period:treatment","F value"],pval=p_interaction_treatment,stringsAsFactors = F))
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
        model <- lmer(   variable ~ period*exposure + period*size + (1|time_of_day) + (1|colony) + (1|antID),data=data ,control = model_n_interations)
        anov  <- anova(model)
        p_interaction_exposure <- anov["Pr(>F)"]["period:exposure","Pr(>F)"]
        p_interaction_size     <- anov["Pr(>F)"]["period:size","Pr(>F)"] 
        
        if (p_interaction_size<=0.05&p_interaction_exposure<=0.05){
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
          stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
          
          
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
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:size",df=paste(round(anov["period:size","NumDF"]),round(anov["period:size","DenDF"]),sep=","),Fvalue=anov["period:size","F value"],pval=anov["period:size","Pr(>F)"],stringsAsFactors = F))
            
            ###fit model 
            model <- lmer(   variable ~ period*exposure + (1|time_of_day) + (1|colony) + (1|antID),data=data ,control = model_n_interations)
            anov  <- anova(model)
            test_norm(residuals(model))
            for (rowi in 1:nrow(anov)){
              stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
            stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor="period:exposure",df=paste(round(anov["period:exposure","NumDF"]),round(anov["period:exposure","DenDF"]),sep=","),Fvalue=anov["period:exposure","F value"],pval=anov["period:exposure","Pr(>F)"],stringsAsFactors = F))
            
            # print( "Interaction between exposure and period is not significant:")
            # print(anov["period:exposure",])
            # 
            ###fit model 
            model <- lmer(   variable ~ period*size + (1|time_of_day) + (1|colony) + (1|antID),data=data ,control = model_n_interations)
            anov  <- anova(model)
            test_norm(residuals(model))
            for (rowi in 1:nrow(anov)){
              stats_outcomes <- rbind(stats_outcomes,data.frame(variable=variable_list[i],predictor=rownames(anov[rowi,]),df=paste(round(anov[rowi,"NumDF"]),round(anov[rowi,"DenDF"]),sep=","),Fvalue=anov[rowi,"F value"],pval=anov[rowi,"Pr(>F)"],stringsAsFactors = F))
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
    
  
    barplot_delta_period <- barplot_delta(dataset=data,predictor="treatment",post_hoc_outcomes=post_hoc_outcomes,stats_outcomes=stats_outcomes,i=i,type="individual",collective=F,plot_untransformed=T,diff_type="absolute_difference") #form_stat=NULL,
    if (showPlot) {print(barplot_delta_period)}
    
    barplot_delta_period_list[[variable_list[i]]]        <- barplot_delta_period

    
    ## add cleaning step
    # post_hoc_outcomes1 <- c(post_hoc_outcomes1,post_hoc_outcomes)
    # stats_outcomes1    <- c(stats_outcomes1,stats_outcomes)
    # 
    # rm(list = ls()[which(names(ls()) == "post_hoc_outcomes")])
    # rm(list = ls()[which(names(ls()) == "stats_outcomes")])
  }
  
  
  rownames(stats_outcomes) <- NULL
  
  #add  formatted output
  stats_outcomes$formatted <- paste0("(GLMM,treatment-induced changes (", stats_outcomes$variable,", ",stats_outcomes$predictor,"): F(",stats_outcomes$df,") = ", round(stats_outcomes$Fvalue,3), ", P ≤ ",format(from_p_to_prounded(stats_outcomes$pval),scientific=F))
  
  return(list(stats_outcomes=stats_outcomes,    post_hoc_outcomes=post_hoc_outcomes,  barplot_delta_period_list=barplot_delta_period_list))
  
}


line_plot <- function(data_path, which_individuals,showPlot=T){
  
  ###1. read data
  setwd(data_path)
  file_list <- list.files(pattern=pattern)
  print(file_list)
  warning("this function is not well generalised")
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
  
  line_plot_obj_list <- list()
  
  for (i in 1:length(variable_list)){
    print(paste0("######## ",variable_list[i]," ########"))
    data <- data_ori
    
    
    # ###create a variable column
    data$untransformed_variable <- data[,variable_list[i]]
    # 
    # ###transform variable
    # if (transf_variable_list[i]=="log"){
    #   print("Logging variable...")
    #   data[!is.na(data$untransformed_variable),"variable"] <- log_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
    # }else if (grepl("power",transf_variable_list[i])){
    #   data[!is.na(data$untransformed_variable),"variable"]  <- (data[!is.na(data$untransformed_variable),"untransformed_variable"] )^as.numeric(gsub("power","",transf_variable_list[i]))
    # }else if (transf_variable_list[i]=="sqrt"){
    #   data[!is.na(data$untransformed_variable),"variable"]  <- sqrt_transf(data[!is.na(data$untransformed_variable),"untransformed_variable"] )
    # }else if (transf_variable_list[i]=="none"){
    #   data$variable <- data$untransformed_variable
    # }
    #hist(data$variable, breaks = 100)
    ###statistics
    ###make sure treatment, exposure, size and period are factors with levels in the desired order
    data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
    data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
    data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
    data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
    data$antID     <- factor( data$antID )
    
    
  
  # 1. Calculate the mean of the untransformed_variable by tag
  mean_data <- aggregate(untransformed_variable ~ period + time_hours + treatment + colony, 
                         FUN = mean, na.rm = T, na.action = na.pass, data)
  
  # 2. Calculate the grand mean and standard error dropping the colony factor
  grand_mean_data <- mean_data %>%
    group_by(period, time_hours, treatment) %>%
    summarise(grand_mean = mean(untransformed_variable),
              standard_error = sd(untransformed_variable) / sqrt(n()))
  
  # Add NA values at time_hours == -3
  unique_treatments <- unique(grand_mean_data$treatment)
  unique_periods <- unique(grand_mean_data$period)
  na_rows <- expand.grid(period = unique_periods,
                         time_hours = -3,
                         treatment = unique_treatments,
                         grand_mean = NA,
                         standard_error = NA)
  
  grand_mean_data <- rbind(grand_mean_data, na_rows) %>%
    arrange(treatment, period, time_hours)
  
  # 3. Create a ggplot geom_line with geom_ribbon for the untransformed_variable
  line_plot_obj <- ggplot(grand_mean_data, aes(x = time_hours, y = grand_mean, color = treatment, group = treatment)) +
    geom_line() +
    geom_ribbon(aes(ymin = grand_mean - standard_error, ymax = grand_mean + standard_error, fill = treatment), alpha = 0.2) +
    scale_x_continuous(limits = c(min(grand_mean_data$time_hours), max(grand_mean_data$time_hours)), expand = c(0, 0)) +
    labs(#title = "Prop Time Outside by Time Hours and Treatment",
         x = "Time Hours since treatment exposure",
         y = names(variable_list[i])
         ) +
    STYLE +
    colFill_treatment +
    colScale_treatment
    
  if (showPlot) {print(line_plot_obj)}
  line_plot_obj_list[[variable_list[i]]]        <- line_plot_obj
  }
  return(line_plot_obj_list)
}


barplot_delta <-
  function(dataset,
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
          variable ~ predictor + colony_size + colony + time + tag + task_group,
          FUN = mean,
          data = diff
        )
    }


    #mean by colony (drop time_of_day)
    diff <-
      aggregate(
        na.rm = T,
        na.action = "na.pass",
        as.formula(paste("variable ~ predictor + colony_size + colony + time",
                         if (!collective) "+ task_group")),
        FUN = mean,
        data = diff
      )
    
    means <-
      aggregate(
        na.rm = T,
        na.action = "na.pass",
        as.formula(paste("variable ~ time + predictor",
                         if (!collective) "+ task_group")),
        FUN = "mean",
        data = diff
      )
    ses <-
      aggregate(
        na.rm = T,
        na.action = "na.pass",
        as.formula(paste("variable ~ time + predictor",
                         if (!collective) "+ task_group")),
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
    max_mean_se <- 1.2 * max(means$mean + means$se)
    # rename vars
    
    
    
    # if interaction period*treatment
    if (all(grepl(paste0(levels(dataset$treatment), collapse = "|"), names(post_hoc_outcomes[[variable_list[i]]])))
      #all(names(post_hoc_outcomes[[variable_list[i]]]) == levels(dataset$treatment))
      ) {
      print("period*treatment interaction significant")
      
      # Create the base ggplot object
      plot_var <-
        ggplot(means,
               aes_string(x = "predictor", y = "mean", fill = "predictor")) +
        STYLE +
        colFill_treatment +
        labs(y = paste("\u0394", names(variable_list[i]),
                       ifelse(transf_variable_list[i] == "none", "", paste0("(", transf_variable_list[i], ")")), #transformation
                       sep = " "), x = predictor)
      
      
      # if there are two ants groups
      if (any(grepl("task_group", stats_outcomes[which(stats_outcomes$variable == variable_list[i]),"predictor"])) ) {

        col_geom <-         list( 
          geom_errorbar(
          aes_string(ymin = paste0("mean - se"), ymax = paste0("mean + se"),group="task_group"),
          position = position_dodge2(width = 0.8, preserve = "single"),
          color="black"
        ) ,
          geom_col(position = position_dodge2(width = 0.8, preserve = "single"), aes_string(col = "task_group"), size = 0.8
                   ),
        scale_color_manual(values = c("#00FFFF", "#FF00FF"))
        )

        plot_var <- plot_var + col_geom
      } else {
        
        col_geom <- list(
          geom_errorbar(
          aes_string(ymin = paste0("mean - se"), ymax = paste0("mean + se")),
          position = position_dodge2(width = 0.8, preserve = "single")
        ) ,
          geom_col(position = position_dodge2(width = 0.8, preserve = "single"))
        )
        plot_var <- plot_var + col_geom
      }
      
      
      # add posthoc letters
      if (variable_list[i] %in% names(post_hoc_outcomes)) {
        # if there are two ants groups
        # first if: check if task_group appears in the predictors
        # second if: return true if stats_outcomes$pval is lower than 0.05 for row that grep "task" in stats_outcomes$predictor
        if (any(grepl("task_group", stats_outcomes[which(stats_outcomes$variable == variable_list[i]),"predictor"])
            &  ifelse(grepl("task_group", stats_outcomes$predictor) & stats_outcomes$pval < 0.05, TRUE, FALSE))
            ) {

         reshaped_posthocs <- data.frame(predictor = names(post_hoc_outcomes[[variable_list[i]]]), letters = post_hoc_outcomes[[variable_list[i]]], stringsAsFactors = FALSE)
          
         reshaped_posthocs$task_group <-  sub(".*\\.(.*)", "\\1", reshaped_posthocs$predictor)
         reshaped_posthocs$predictor <- sub("\\.[^.]*$", "", reshaped_posthocs$predictor)
         
         warning("geom_text faulty, as it shows letters inverted in each treatment. temporary fix: coloured the letters to help readability")
          additional_geoms <- list(geom_text(data = reshaped_posthocs,
            aes(label = letters, x = predictor, y = max_mean_se, col  = task_group, group= task_group),
            position = position_dodge2(width = 0.8, preserve = "single"),
            vjust = -0.4
          ))
          plot_var <- plot_var + additional_geoms

          #normal condition (4 groups)
        }else{
          
          reshaped_posthocs <- data.frame(predictor = names(post_hoc_outcomes[[variable_list[i]]]), letters = post_hoc_outcomes[[variable_list[i]]], stringsAsFactors = FALSE)
          reshaped_posthocs <- merge(reshaped_posthocs,means)
          
          # adjust position of letters if vals are negative
          reshaped_posthocs$se <- ifelse(reshaped_posthocs$mean+reshaped_posthocs$se<0,
                                         0,
                                         reshaped_posthocs$se)
          reshaped_posthocs$mean <- ifelse(reshaped_posthocs$mean<0,0,reshaped_posthocs$mean)
          
        additional_geoms <- list(geom_text(data = reshaped_posthocs,
                                           aes(label = letters, x = predictor, y = mean+se),
                                           position = position_dodge2(width = 0.8, preserve = "single"),
                                           vjust = -0.4
        ))
          
        #   list(geom_text(
        #   aes(label = post_hoc_outcomes[[variable_list[i]]][predictor], y = max_mean_se),
        #   position = position_dodge2(width = 0.8, preserve = "single"),
        #   vjust = -0.4
        # ))
        plot_var <- plot_var + additional_geoms
        }
         }
      
      # manage x axis labels
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
      
      
    } else if (all(names(post_hoc_outcomes[[variable_list[i]]]) == levels(dataset$size))
               ) {
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
        labs(y = paste("\u0394", names(variable_list[i]),
                       ifelse(transf_variable_list[i] == "none", "", paste0("(", transf_variable_list[i], ")")), #transformation
                       sep = " "), x = predictor) +
        theme(axis.text.x = element_text())
      
      # manage x axis labels
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
      
      
      # Add bar with significant stars 
      # (for more than 1 element in the ggplot if statement, there is the error: Error in `+.gg`:! Cannot add ggproto objects together )
      if (variable_list[i] %in% names(post_hoc_outcomes)) {
        additional_geoms <- list(geom_segment(
          aes(
            x = sum(x_axis[group == "small"]) / 2,
            xend = sum(x_axis[group == "big"]) / 2,
            y = abs(max_mean_se),
            yend = abs(max_mean_se)
          )
        ),
        geom_text(
          aes(
            x = mean(x_axis),
            y = abs(max_mean_se + 1 / 2 * (max_mean_se)),
            label = paste(from_p_to_ptext(stats_outcomes[which(
              stats_outcomes$variable == variable_list[i] &
                stats_outcomes$predictor == "period:size"
            ), "pval"]),collapse="")
          )
        ))
        plot_var <- plot_var + additional_geoms
      }
      
    } else if (
      all(names(post_hoc_outcomes[[variable_list[i]]]) %in% levels(dataset$exposure)) #&&
      #all(levels(dataset$exposure) %in% names(post_hoc_outcomes[[variable_list[i]]])) #extra recursivness
      ) {
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
        labs(y = paste("\u0394", names(variable_list[i]),
                       ifelse(transf_variable_list[i] == "none", "", paste0("(", transf_variable_list[i], ")")), #transformation
                       sep = " "), x = predictor) +
        theme(axis.text.x = element_text())
      
      # manage x axis labels
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
      # 
      # Add bar with significant stars
      # for more than 1 element in the ggplot if statement, there is the error: Error in `+.gg`:! Cannot add ggproto objects together )
      if (variable_list[i] %in% names(post_hoc_outcomes)) {
        
        
         additional_geoms <- list(geom_segment(
          aes(
            x = sum(x_axis[group == "control"]) / 2,
            xend = sum(x_axis[group == "pathogen"]) / 2,
            y = abs(max_mean_se),
            yend = abs(max_mean_se)
          )
        ),
        
        geom_text(
          aes(
            x = mean(x_axis),
            y = abs(max_mean_se + 1 / 2 * (max_mean_se)),
            label = paste(from_p_to_ptext(stats_outcomes[which(
              stats_outcomes$variable == variable_list[i] &
                stats_outcomes$predictor == "period:exposure"
            ), "pval"]),collapse="")
          )
        ))
        plot_var <- plot_var + additional_geoms #+ geom_text(aes(label = ''))
        
        # print(paste0("*****************************************
        #              \nFAULTY GEOM!!!\n",stats_outcomes[which(
        #       stats_outcomes$variable == variable_list[i] &
        #         stats_outcomes$predictor == "period:exposure"
        #     ), "pval"]))
      }

    }else{
      plot_var <- plot_var + geom_text(aes(label = '')) #add empty label
    }
    
    if(!all(grepl(paste0(levels(dataset$treatment), collapse = "|"), names(post_hoc_outcomes[[variable_list[i]]])))){#not main int sig
    if (all(names(post_hoc_outcomes[[variable_list[i]]]) == levels(dataset$size)) && # size is sign
        all(names(post_hoc_outcomes[[variable_list[i]]]) %in% levels(dataset$exposure)) #exp is sign
        ) {
      warning("both period*exposure and period*size are significant. There is currently no plotting contidion for that")
    }}
    
    
    #legend modifiers (2 rows and changed labs)
    plot_var <- plot_var +
    guides(fill = guide_legend(nrow = 2))
    # +guides(fill = "none")   # remove treatment legend

    
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
            variable_before ~ colony_size + colony + predictor + size + tag + time_of_day + task_group,
            FUN = mean,
            data = befores
          )
        afters <-
          aggregate(
            na.rm = T,
            na.action = "na.pass",
            variable_after ~ colony_size + colony + predictor + size + tag + time_of_day + task_group,
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


read_data <- function(data_path, which_individuals){
  
  ###1. read data
  setwd(data_path)
  file_list <- list.files(pattern=pattern)
  print(file_list)
  warning("this function is not well generalised")
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
  
  data_list <- list()
  
  # for (i in 1:length(variable_list)){
  #   print(paste0("######## ",variable_list[i]," ########"))
  #   data <- data_ori
  #   
  #   
  #   # ###create a variable column
  #   data$untransformed_variable <- data[,variable_list[i]]
  #   
  ###statistics
  ###make sure treatment, exposure, size and period are factors with levels in the desired order
  data$treatment <- factor(data$treatment , levels=treatment_order[which(treatment_order%in%data$treatment )])
  data$size      <- factor(data$size    , levels=size_order   [which(size_order%in%data$size )])
  data$exposure  <- factor(data$exposure , levels=exposure_order[which(exposure_order%in%data$exposure )])
  data$period    <- factor(data$period    , levels=period_order   [which(period_order%in%data$period )])
  data$antID     <- factor( data$antID )
  
  
  # data_list[[variable_list[i]]] <- data
  
  #}
  return(data)
}

plot_qpcr <- function(experiments){
  plot1 <- function(data,experiment,ylim=ylim){
    par_mari <- par("mar")
    par(mar=par_mari-c(0,0.35,0,par_mari[4]))
    
    ###Plot (load)=f(caste)
    data["variable"] <- data$measured_load_ng_per_uL;data["predictor"] <- data$task_group
    replac_val <- min(data$variable[data$variable!=0],na.rm=T)/sqrt(2)
    
    ####transform 
    data <- aggregate(na.action=na.pass,na.rm=T,log10(variable+replac_val)~predictor+colony,FUN=mean,data=data)
    names(data)[grepl("variable",names(data))] <-"variable"
    
    means <- aggregate(na.rm=T,na.action="na.pass",variable~predictor,FUN="mean",data=data);ses <- aggregate(na.rm=T,na.action="na.pass",variable~predictor,FUN="std.error",data=data);
    
    names(means)[grepl("variable",names(means))] <- "mean";names(ses)[grepl("variable",names(ses))] <- "se";means <- merge(means,ses)
    means <- means[order(match(means$predictor,status_order)),]
    
    means[is.na(means$se),"se"] <- 0
    if (is.null(ylim)){
      ymin <- min(c(means$mean-means$se),na.rm=T);ymax<- max(c(means$mean+means$se),na.rm=T)
      ymin <- floor(ymin)
      ymax <- ceiling(ymax)
      rangy <- ymax-ymin
      ymax <- ymax+0.1*rangy
      yrange <- c(ymin,ymax)
    }else{
      ymin <- ylim[1]
      ymax <- ylim[2]
      yrange <- ymax-ymin
      yrange <- c(ymin-0.04*yrange,ymax+0.04*yrange)
    }
    
    barwidth <- 0.5; barwidth_fac_within <- 0.5; barwidth_fac_between <- 2
    barspace <- c(barwidth_fac_between,barwidth_fac_within,barwidth_fac_within)
    
    plotx <- barplot(means$mean,plot=F,width=barwidth,space=barspace)
    ####empty plot
    plot(plotx,means$mean,ylim=yrange,xlim=c(min(plotx)- 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ 0.6*(barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab=expression(paste("Mean measured pathogen load (ng/", mu, "L)")),bty="n",xaxs="i",yaxs="i",type="n",cex.axis=min_cex,cex.lab=inter_cex,lwd=line_min,xaxt="n",yaxt="n")
    ####arrows
    plot_arrows(means=means,plotx=plotx,plot_type="means",LWD=line_max,LENGTH=0.025,colz=statuses_colours[as.character(means$predictor)])
    ####points
    points(plotx,means$mean,col=statuses_colours[as.character(means$predictor)],pch=16,cex=max_cex*0.8,lwd=line_min)
    ####Y-axis
    axis(2,at=ymin:floor(ymax),labels=format(10^(ymin:floor(ymax)),scientific=T),cex.axis=min_cex,cex.lab=inter_cex,lwd=0,lwd.ticks=1)
    
    par(xpd=T)
    labses <- full_statuses_names[as.character(means$predictor)]
    labses <- c(labses[1],rep("",length(labses)-1))
    axis(1,at=plotx,labels=labses,tick=F,cex.axis=inter_cex,las=1,mgp=c(0.8,0.8,0))
    
    for (labse in 2:length(labses)){
      print(par("mgp")[2] + (labse-1)*1)
      if (labse/2==round(labse/2)){
        mtext(full_statuses_names[as.character(means$predictor)][labse],side=1,at=plotx[labse],line=0+par("mgp")[1] + 1,cex=par("cex")*inter_cex)
      }else{
        mtext(full_statuses_names[as.character(means$predictor)][labse],side=1,at=plotx[labse],line=0+par("mgp")[1] ,cex=par("cex")*inter_cex)
      }
    }
    
    segments(x0=min(plotx)- 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),y0=yrange[1],x1=max(plotx)+ 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),lwd=line_max)
    segments(x0=min(plotx)- 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),y0=yrange[1],y1=yrange[2],lwd=line_max)
    par(xpd=F)
    
    
    data$predictor <- factor(data$predictor,levels=means$predictor)
    
    model_lmer <- lmer(variable~predictor+(1|colony),data=data ,control = model_n_interations)
    test_norm(residuals(model_lmer))
    pvalue <- Anova(model_lmer)["predictor","Pr(>Chisq)"]
    print(Anova(model_lmer))
    # title(ylab=lab_title,cex.lab=1.5,mgp=c(5,1,0))
    par(xpd=F)
    abline(h=0)
    if (pvalue>0.05){p_cex <- inter_cex}else{p_cex <- max_cex}
    title(main=from_p_to_ptext(pvalue),cex.main=p_cex,font.main=2,line=stat_line,xpd=T)
    
    post_hoc <- summary(glht(model_lmer,linfct = mcp(predictor="Tukey")),test=adjusted("BH"))
    print(post_hoc)
    post_hoc_groups <- cld(post_hoc)$mcletters$Letters
    for (idx in 1:length(post_hoc_groups)){
      group <- post_hoc_groups[as.character(means[idx,"predictor"])]
      text(x=plotx[idx],y=ymax-0.1*(ymax-ymin),adj=c(0.5,0),labels=as.character(group),xpd=T,cex=inter_cex)
    }
    par(mar=par_mari)
  }
  
  all_qpcr_data <- NULL
  data_for_plot <- NULL
  for (experiment in experiments){
    ###read qpcr data
    warning("rename file qPCR_results.txt instead of qPCR_file.txt")
    data <- read.table(paste(disk_path,"/",experiment,"/original_data/qPCR/qPCR_results.txt",sep=""),header=T,stringsAsFactors = F)
    ###keep only ants
    data <- data[which(!is.na(as.numeric(data$tag))),]
    ###remove workers that died before the end
    data <- data[which(data$alive_at_sampling_time),]
    
    ###read task group
    task_groups <-  read.table(paste(disk_path,"/",experiment,"/original_data/",task_group_file,sep=""),header=T,stringsAsFactors = F)
    ### add task groups info to data
    data <- merge(data,task_groups,all.x=T,all.y=F)
    
    ###keep only pathogen treatments
    data <- data[grep("pathogen", data$treatment),] #AW
    
    ###add metadata
    data <- data.frame(experiment=experiment,data,stringsAsFactors = F)
    data$period <- "after"
    ####list desired variables and transformations
    data$antid <- as.character(data$colony,data$tag)
    all_qpcr_data <- rbind(all_qpcr_data,data)
    data$colony <- as.character(interaction(data$experiment,data$colony))
    
    ###remove treated workers
    data <- data[which(data$status!="treated"),]
    data_for_plot <- rbind(data_for_plot,data)
  }
  
  all_sim_data <- NULL
  for (experiment in experiments){
    data <-read.table(paste(disk_path,experiment,"transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds","individual_simulation_results_observed.txt",sep="/"),header=T,stringsAsFactors=F) #AW # in stroeymeyt script there is no trace on how the "calibration/individual_simulation_results.dat" is generated
    all_sim_data <- rbind(all_sim_data,data.frame(experiment=experiment,data[names(data)!="age"],stringsAsFactors = F))
  
    #AW select (as reported in paper figure 2A) load received by untreated ants in post-treatment simulation ("Simulationswere run over the posttreatment networks using experimentally treated foragers as disease origin")
    all_sim_data <- all_sim_data[grep("pathogen", all_sim_data$treatment),]
    all_sim_data <- all_sim_data[grep("untreated", all_sim_data$status),]
    all_sim_data <- all_sim_data[grep("after", all_sim_data$period),]
    }
  
  ###Now join qPCR and simulation data into single data frame
  all_qpcr_data <- all_qpcr_data[which(!is.na(as.numeric(all_qpcr_data$tag))),]
  all_qpcr_data <- all_qpcr_data[which(all_qpcr_data$alive_at_sampling_time),]
  
  all_sim_qpcr_data <- merge(all_qpcr_data[c("experiment","colony","treatment","tag","status","task_group","age","measured_load_ng_per_uL")],all_sim_data[c("experiment","colony","tag","simulated_load")])
  
  ###remove treated individuals
  all_sim_qpcr_data <- all_sim_qpcr_data[which(all_sim_qpcr_data$status!="treated"),]
  all_sim_qpcr_data["antid"] <- as.character(interaction(all_sim_qpcr_data$experiment,all_sim_qpcr_data$colony,all_sim_qpcr_data$tag))
  
  varb <- "simulated_load"
  variable_list <- c("measured_load_ng_per_uL")
  names(variable_list) <-  c("Measured pathogen load")
  predictor_list <- c( "simulated_load")
  names(predictor_list) <- c("Simulated pathogen load")
  transf_variable_list <- c("log") #AW c("log")
  transf_predictor_list <- c("log") #AW c("log")
  
  ymin <- floor(log10(min(all_sim_qpcr_data$measured_load_ng_per_uL[all_sim_qpcr_data$measured_load_ng_per_uL!=0])/sqrt(2)))
  ymax <- ceiling(log10(max(all_sim_qpcr_data$measured_load_ng_per_uL)))
  xmin <- floor(log10(min(all_sim_qpcr_data[all_sim_qpcr_data[,varb]!=0,varb])/sqrt(2)))
  xmax <- ceiling(log10(max(all_sim_qpcr_data[,varb])))
  
  analysis <- list(variable_list=variable_list,
                   predictor_list=predictor_list,
                   transf_variable_list=transf_variable_list,
                   transf_predictor_list=transf_predictor_list,
                   violin_plot_param = list(c(1,0,-0.025,0.2,0.2)))
  
  ###
  # all_sim_qpcr_data_big <- all_sim_qpcr_data[grep("big", all_sim_qpcr_data$treatment),]
  # all_sim_qpcr_data_small <- all_sim_qpcr_data[grep("small", all_sim_qpcr_data$treatment),]
  # 
  # hist(all_sim_qpcr_data_small[,c("measured_load_ng_per_uL","simulated_load")],breaks = 100)
  # hist(all_sim_qpcr_data_big[,c("measured_load_ng_per_uL","simulated_load")],breaks = 100)
  
 #  # Plotting the variables by treatment
 # print( ggplot(all_sim_qpcr_data, aes(x = measured_load_ng_per_uL, fill = treatment)) +
 #          geom_histogram(bins = 30, alpha = 0.7) +
 #          labs(x = "Measured Load (ng/uL)", y = "Frequency") +
 #          ggtitle("Frequency Histogram of Measured Load by Treatment")
 # )
 # 
 # print(
 #   ggplot(all_sim_qpcr_data, aes(x = simulated_load, fill = treatment)) +
 #     geom_histogram(bins = 30, alpha = 0.7) +
 #     labs(x = "Simulated Load", y = "Frequency") +
 #     ggtitle("Frequency Histogram of Simulated Load by Treatment")
 # )
  
  predicted_value <- plot_regression(data=all_sim_qpcr_data,time_point="after",analysis=analysis,n_cat_horiz=20,n_cat_vertic=11,pool=c(F,F),collective=T,input_color=colour_palette_age,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,point_cex=1.5,predict=NULL) #high_threshold
  par(xpd=NA)
  if (varb =="simulated_load"){
    title(sub=expression(italic("(Prop. exposure dose)")),cex.sub=min_cex,font.sub=1,mgp=c(1,0.1,0))
  }else if (varb =="predicted_measured_load_ng_per_uL_SI"){
    title(sub=expression(paste("(ng/", mu, "L)")),cex.sub=min_cex,font.sub=1,mgp=c(1,0.1,0))
  }
  par(xpd=F)
  
  if (!exists("ymin")){
    ylim <- NULL
  }else{
    ylim <- c(ymin,ymax)
  }
  plot1(data=data_for_plot,experiment="both",ylim=NULL)
  if (exists("predicted_value")){return(predicted_value)}else{return(NULL)} 
}

plot_distribution <- function(experiments,desired_treatments){
  par_mar_ori <- par()$mar
  par(mar=par_mar_ori+c(0,0,0,0.5))
  if (experiments=="all"){
    experiments <- c("age_experiment","survival_experiment","main_experiment")
  }
  if (experiments =="both"){
    experiments <- c("age_experiment","main_experiment")
  }
  transf <- function(x){
    return(x^(1/2))
  }
  rev_transf <- function(x){
    return(x^2)
  }
  xlabel <- substitute(sqrt ( xlabely),list(xlabely="Simulated load"))
  ####read data
  infection_data <- NULL
  for (experiment in experiments){
    setwd(paste(disk_path,experiment,"transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="/"))
    si_outcome <- read.table("individual_simulation_results_observed.txt",header=T,stringsAsFactors = T)
    si_outcome["ant_id"] <- as.character(interaction(experiment,si_outcome$colony,si_outcome$tag))
    
    for (time_point in c("before","after")){
      si_outcome_temp <- si_outcome[which(si_outcome$period==time_point),c("colony","treatment","tag","status","antid","simulated_load","transmission_rank")]
      names(si_outcome_temp)[names(si_outcome_temp)%in%c("simulated_load","transmission_rank")] <- paste(names(si_outcome_temp)[names(si_outcome_temp)%in%c("simulated_load","transmission_rank")],time_point,sep="_")
      assign(paste("si_outcome_",time_point,sep=""),si_outcome_temp)
    }
    si_outcome <- merge(si_outcome_before,si_outcome_after)
    infection_data <- rbind(infection_data,data.frame(experiment=experiment,si_outcome,stringsAsFactors = F))
  }
  
  ####fill in missing untreated and queen info
  # infection_data[which(infection_data$tag==queenid),"tag"] <- "queen"
  #infection_data <- infection_data[which(infection_data$task_group=="queen"),"tag"] <- "queen" # AW this line does not make sense
  
  ###read task group #AW
  task_groups <-  read.table(paste(disk_path,"/",experiment,"/original_data/",task_group_file,sep=""),header=T,stringsAsFactors = F)
  ### add task groups info to data #AW
  infection_data <- merge(infection_data,task_groups,all.x=T,all.y=F)

  ###modify data
  infection_data[,"colony"] <- as.character(interaction(infection_data$experiment,infection_data$colony))
  
  ####make new datasets for further analyses######
  infection_data$status <- as.character(infection_data$status)
  infection_data <- infection_data[infection_data$status!="treated",]###contains queens and untreated workers
  infection_data <- infection_data[infection_data$task_group!="queen",]###contains untreated workers

  #####1. Make bins
  xmin <- min(c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after)),na.rm=T)
  xmin_bis <- min(
    c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after))
    [
      c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after))
      >0
    ]
    ,
    na.rm=T
  )
  xmax <- max(c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after)),na.rm=T)
  xmin <- min(c(floor((xmin)*100)/100))
  xmax <- ceiling((xmax)*100)/100
  breaks <- seq(from=xmin,to=xmax,length.out=23)
  bw <- 0.09
  xmax <- ceiling((xmax)*10)/10
  
  ###2. plot histograms plus densities
  #infection_data <- infection_data[which(infection_data$treatment%in%desired_treatments),]
  infection_data <- infection_data[grep(paste(desired_treatments, collapse = "|"), infection_data$treatment),] #AW
  #infection_data <- infection_data[grep("small", infection_data$treatment),] #TEMP
  
  
  afters_dens <- density(transf(infection_data$simulated_load_after),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
  befores_dens <- density(transf(infection_data$simulated_load_before),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
  
  percolony_density <- NULL
  for (colony in sort(unique(infection_data$colony[!is.na(infection_data$colony)]))){
    if (colony!="age_experiment.colony021"){
      subsett <- infection_data[which(infection_data$colony==colony),]
      
      afters <- density(transf(subsett$simulated_load_after),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
      afters$y <- afters$y/sum(afters$y)
      befores <- density(transf(subsett$simulated_load_before),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
      befores$y <- befores$y/sum(befores$y)
      percolony_density <- rbind(percolony_density,data.frame(colony=colony,xcoor=afters$x,density_after=afters$y,density_before=befores$y,density_diff=afters$y-befores$y))
    }
  }
  
  ####plot actual distributions
  forplot <- data.frame(as.matrix(aggregate(cbind(density_after,density_before,density_diff)~xcoor,function(x)cbind(mean(x),std.error(x)),data=percolony_density)))
  names(forplot) <- c("x","mean_after","std.error_after","mean_before","std.error_before","mean","std.error")
  forplot["lower_y"] <- forplot$mean-forplot$std.error
  forplot["top_y"] <- forplot$mean+forplot$std.error
  forplot <- forplot[order(forplot$x),]
  
  xshade <- c(forplot$x,rev(forplot$x))
  yshade <- c(forplot$lower_y,rev(forplot$top_y))
  
  ##get ylims for plots ####
  ######first get an idea of how the density plots will be distributed
  ymin_dens <- 2*min(forplot$lower_y)
  ymax_dens <- max(c(forplot$mean_after,forplot$mean_before))+0.1*(max(c(forplot$mean_after,forplot$mean_before))-min(forplot$mean))
  neg <- abs(ymin_dens)/abs(ymax_dens)
  
  
  ######second get an idea of how the histogram will be distributed
  afters <- hist(transf(infection_data$simulated_load_after),breaks=breaks,plot=F)
  befores <- hist(transf(infection_data$simulated_load_before),breaks=breaks,plot=F)
  ymax <- max(c(afters$density,befores$density))+0.1*max(c(afters$density,befores$density))
  #####...and deduce ymin
  ymin <- -neg*ymax
  
  ###then plot histograms; in frequency
  befores <- hist(transf(infection_data$simulated_load_before),breaks=breaks,plot=T,col=alpha("blue",0.4),prob=T,ylim=c(ymin,ymax),xlim=c(min(0,xmin),xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",main="",bty="l",cex.axis=min_cex,cex.lab=inter_cex,xlab="",lwd=line_min/10,border="black")
  afters <- hist(transf(infection_data$simulated_load_after),breaks=breaks,col=alpha("red",0.3),add=T,plot=T,prob=T,lwd=line_min/10,border="black")
  prospected_ats <- axisTicks(c(min(0,xmin),xmax),log=F)
  prospected_ats <- c(prospected_ats,prospected_ats[length(prospected_ats)]+prospected_ats[length(prospected_ats)]-prospected_ats[-1+length(prospected_ats)])
  axis(1,at=prospected_ats,cex.axis=min_cex,cex.lab=inter_cex)
  axis(2,cex.axis=min_cex,cex.lab=inter_cex)
  title(xlab=xlabel,cex.axis=min_cex,cex.lab=inter_cex,mgp=par()$mgp+c(0.2,0,0))
  ###then add densities; with new axes
  par(new=T)
  #####get new ymin, ymax
  plot(forplot$x,forplot$mean, type="n", axes=F, xlab=NA, ylab=NA, cex=1.2,col=statuses_colours[desired_treatments],ylim=c(ymin_dens,ymax_dens),xaxs="i",yaxs="i",main="",bty="l",cex.axis=min_cex,cex.lab=inter_cex,xlim=c(min(0,xmin),xmax))
  polygon(xshade,yshade, border = alpha("yellow",0),col=alpha("yellow",0.6))
  lines(forplot$x,forplot$mean,lwd=line_max,col="black")
  
  ####write legend
  legend(x=par("usr")[2]+0.05*(par("usr")[2]-par("usr")[1]),y=par("usr")[4]-0.025*(par("usr")[4]-par("usr")[3]),xjust=1,yjust=1,legend=c("Pre-treatment","Post-treatment","Difference"),pt.bg=c(alpha("blue",0.4),alpha("red",0.3),alpha("yellow",0.6)),col=c("black","black",alpha("yellow",0)),bty='n',pch=22,lty=0,lwd=0,pt.lwd=1,pt.cex=1.5,text.col="white",cex=min_cex)
  legend(x=par("usr")[2]+0.05*(par("usr")[2]-par("usr")[1]),y=par("usr")[4]-0.025*(par("usr")[4]-par("usr")[3]),xjust=1,yjust=1,legend=c("Pre-treatment","Post-treatment","Difference"),col=c(alpha("blue",0),alpha("red",0),"black"),bty='n',lty=c(0,0,1),lwd=c(0,0,1),cex=min_cex)
  
  
  print("KS-test:")
  ks_test <- ks.test(transf(infection_data$simulated_load_after),transf(infection_data$simulated_load_before))
  print(ks_test)
  p_value <- ks_test$p.value
  
  where_to_print_stat <- median(c(transf(infection_data$simulated_load_after),transf(infection_data$simulated_load_before)))
  
  par(xpd=T) 
  mtext(full_statuses_names[desired_treatments],side=3,line=stat_line,adj=0.5,cex=par("cex") *inter_cex,font=2)
  mtext(from_p_to_ptext(p_value),side=3,line=stat_line-1,adj=0.5,cex=par("cex") *max_cex,font=2,at=where_to_print_stat)
  
  # print("Thresholds at which after becomes lower than before")
  forplot["positive"] <- forplot$mean>=0
  change <- min(which(diff(forplot$positive)==-1))
  threshold1 <- rev_transf(forplot[change,"x"])
  threshold2 <- rev_transf(forplot[change+1,"x"])
  if(!exists("threshold")){threshold <-round(((threshold1+threshold2)/2)*10000)/10000}
  
  par(xpd=F)
  lines(x=c(transf(threshold),transf(threshold)),y=c(0,ymin_dens),col="springgreen3",lty=1,lwd=2*line_max)
  ###Now write down "high level", "low_level"
  arrows(x0=transf(threshold),y0=1.7*min(forplot$lower_y),x1=par("usr")[2],y1=1.7*min(forplot$lower_y),col="springgreen4",code=3,length=0.025)
  text(x=mean(c(transf(threshold),par("usr")[2])),y=1.5*min(forplot$lower_y),labels="high load",adj=c(0.5,0),cex=min_cex,col="springgreen4",font=3)
  xmin <-min(c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after)),na.rm=T)
  arrows(x1=transf(threshold),y0=1.7*min(forplot$lower_y),x0=xmin_bis,y1=1.7*min(forplot$lower_y),col="springgreen2",code=3,length=0.025)
  text(x=mean(c(transf(threshold),xmin)),y=1.5*min(forplot$lower_y),labels="low load",adj=c(0.5,0),cex=min_cex,col="springgreen2",font=3)
  
  par(mar=par_mar_ori)
}


########################################################################
##########  ADRIANO EXP1: ANALYSIS AND STYLING FUNCTIONS ###############

# inherited from FUNCTIONS_Analysis_and_Styling.R

########## STATS FUNCTIONS ###############

# function to test normality of residuals
test_norm <- function(resids) {
  print("Testing normality")
  if (length(resids) <= 300) {
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
    print("below 0.05, the data significantly deviate from a normal distribution")
  } else {
    print("More than 300 data points so using the skewness and kurtosis
approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}

# function to report a model output
output_lmer <- function(model,show_report = F) {
  print(paste("############### MODEL :",deparse(substitute(m1)),"###############",sep = " "))
  print("------------RESIDUALS NORMALITY------------")
  test_norm(residuals(model))
  print("------------SUMMARY------------")
  print(summary(model))
  print("------------ANOVA------------")
  print(anova(model)) 
  warning("this uses Type III Analysis of variance, for when there are interactions")
  print("------------RSQUARED------------")
  print(r.squaredGLMM(model))
  print("------------REPORT------------")
  if (show_report) {
    print(report(model))
  }else{ cat("Set show_report to TRUE to obtain a human readable model report")}
  #tab_model(model)
}

# function to simplify a model if the interaction is not significant
simplify_model <- function(model) {
  # Extract the model formula
  model_formula <- formula(model)
  # Check the significance of the interaction using anova()
  anova_m1 <- as.data.frame(Anova(model))
  print(anova_m1)
  # Find the interaction term in the model formula
  # define a regular expression pattern to match the desired substring
  pattern <- "\\b\\w+\\s*\\*\\s*\\w+\\b"
  # use the sub() function to extract the first match of the pattern
  interaction_term <- sub(paste0(".*(", pattern, ").*"), "\\1", as.character(model_formula)[3])
  interaction_vars <- unlist(strsplit(interaction_term, " * "))
  anova_term <- gsub("\\s*\\*\\s*", ":", interaction_term)
  # If the Anova of the interaction is not significant, simplify the model by removing the interaction
  if (anova_m1[which(rownames(anova_m1)== anova_term),"Pr(>Chisq)"] > 0.05) {
    cat("\n#\nModel interaction NOT significant, simplify\n#\n")
    model_no_interaction_formula <-  as.formula(gsub("\\*", "+", deparse(model_formula)))
    model_no_interaction <- update(model, formula = model_no_interaction_formula)
    #print(summary(m1_no_interaction))
    print(Anova(model_no_interaction))
    return(model_no_interaction)
  } else {
    cat("\n#\nModel interaction significant, don't simplify\n#\n")
    return(model)
  }
}


#check if model has interaction
has_interaction <- function(model) {
  formula_str <- as.character(formula(model))
  return(any(grepl(":", formula_str) | grepl("\\*", formula_str)))
}


# function to perform posthocs
posthoc_list <- list()
interactions_to_explore <- list()
compute_posthocs <- function(model) {
  #warning("this function has only been tested with lmer()")
  warning("How to use: \nA.Run on models without interactions. If interaction is present, run 'simplify_model' first. 
          \nB.If has_interaction= T, paste your variables naming the new var in this format 'VAR1_VAR2'
          \nC. assign the output of this function to 'posthoc_list'
          \nD.Optional: provide an ID_model [i.e. formatted as paste(GROUP,VAR,sep='-')] to later recall the posthoc.\n\n")
  print(paste("model predictor:", paste0(row.names(Anova(model))), sep = " "))
  # check that there are significant outputs
  if (length(row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) == 0) {
    print("there are no significant vars.")
  } else {
    for (SIG.VAR in row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) {
      if (grepl(":", SIG.VAR)) {
        warning(paste0(SIG.VAR, "is an interaction, currently this situation is not handled by the function. Re-run model using pasted variables"))
        #interactions_to_explore <- c(interactions_to_explore, list(paste(GENE,GROUP,SIG.VAR,deparse(substitute(model)), sep = "-") ))
      } else {
        # check if the variable is not numeric . to do so, we need to access the dataframe from the model
        if (!is.numeric(get(gsub("\\[.*", "", as.character(model@call)[3]))[, SIG.VAR])) {
          print(paste0("Performing posthocs for the significant var: ", SIG.VAR))
          arg <- list("Tukey")
          names(arg) <- SIG.VAR
          # glht (mcp) : General linear hypotheses Testing (glht) and multiple comparisons for parametric models (mcp)
          cmp <- do.call(mcp, arg)
          posthoc_SIG.VAR <- summary(glht(model, linfct = cmp), test = adjusted("BH"))
          # Set up a compact letter display of all pair-wise comparisons
          model_means_cld <- cld(posthoc_SIG.VAR)
          # create dataframe usable with ggplot geom_text
          model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters)))
          # add column name
          model_means_cld$newcol <- NA
          colnames(model_means_cld)[which(names(model_means_cld) == "newcol")] <- SIG.VAR
          colnames(model_means_cld)[which(names(model_means_cld) == "V1")] <- "letters"
          model_means_cld[, SIG.VAR] <- row.names(model_means_cld)
          rownames(model_means_cld) <- NULL
          # if interaction term, split columns
          if (grepl("_", SIG.VAR)) {
            # Split the column into two new columns using "_"
            #SIG.VAR.STRIP <- names(model_means_cld[grep( "_",model_means_cld)])
            model_means_cld[, strsplit(SIG.VAR, "_")[[1]]] <- t(sapply(model_means_cld[,SIG.VAR], function(x) strsplit(x, "_")[[1]]))
          }
          # add to list
          posthoc_list <- c(posthoc_list, list(model_means_cld))
          if (exists("ID_model")) {
            names(posthoc_list)[length(posthoc_list)] <- paste(ID_model,SIG.VAR,  deparse(substitute(model)),sep = "-")
          } else {
            warning("if you provide an ID_model [i.e. formatted as paste(GROUP,VAR,sep='-')] before running the function, it will be easier to later call the posthoc output for plotting")
            names(posthoc_list)[length(posthoc_list)] <- paste(SIG.VAR,deparse(substitute(model)), sep = "-")
          }
          print(paste(deparse(substitute(model)), SIG.VAR, sep = "_"))
          print(model_means_cld)
        } # SIG.VAR LOOP
      } # check if is an interaction
    } # check if numeric
    warning("call 'posthoc_list' to get posthocs")
  } # if significant vars exist
  return(posthoc_list)
}


# # The calculate_weights function takes as input a categorical variable group and an optional argument ratio, which specifies the desired ratio between the weights of the largest and smallest groups
# calculate_weights <- function(group) {
#   print("function to calculate_weights for imbalanced datasets")
#   #make sure that group is a factor
#   group <- as.factor(group)
#   # Calculate the proportions of each group
#   proportions <- table(group) / length(group)
#   print(proportions)
#   # Calculate the weights for each group
#   weights <- numeric(length(group))
#   
#   for (g in levels(group)) {
#     weights[group == g] <- 1 / proportions[g]
#   }
#   # Normalize the weights so that they sum up to 1
#   weights <- weights / sum(unique(weights))
#   return(weights)
# }

## pretty print the model selection output
nice_print_model_sel <- function(model_output) {
  # clean output
  sel.table <- round(as.data.frame(model_output)[-c(1:5)], 3)
  # number of parameters (df) should be K
  names(sel.table)[1] <- "K"
  sel.table$Model <- rownames(sel.table)
  rownames(sel.table) <- NULL
  # replace Model name with formulas
  for (i in 1:nrow(sel.table)) sel.table$formula[i] <- as.character(formula(get(sel.table$Model[i])))[3]
  return(sel.table)
}

# convert significance levels to stars
add_star <- function(p) {
  if (p<0.001) {
    return('***')
  } else if (p<0.01) {
    return('**')
  } else if (p<0.05) {
    return('*')
  } else {
    return('ns')
  }
}




#### individual immunity ###

#### MAXIMUM ACCEPTED DIFF IN CT between duplicates WITH 15% PIPETTING ERROR# #Efficiency is between 1.95 and 2.05, differing CT values.
# FUNCTION MODIFIED FROM  # https://rnajournal.cshlp.org/content/23/5/811.full.pdf, BASED ON THE ABOVE
assign_CTDiff_15PipErr <- function(mean_Ct) {
  if (is.na(mean_Ct)) {
    return(0)
  } else if (mean_Ct > 34.5) {return(1.9)
  } else if (mean_Ct > 33.5) {return(1.3)
  } else if (mean_Ct > 32.5) {return(0.9)
  } else if (mean_Ct > 31.5) {return(0.7)  
  } else if (mean_Ct >    0) {return(0.5)
  } else                     {return(0)
  }
}

# # Loop through transformations
# for (trans_name in names(transformations)) {
#   # Apply transformation
#   trans_func <- transformations[[trans_name]]
#   transformed_data <- trans_func(sample_data)
#   
# 
#   # Test normality using Shapiro-Wilk test
#   shapiro_test <- shapiro.test(transformed_data)
#   cat(trans_name, "Transformation:\n")
#   cat("Shapiro-Wilk Test p-value:", shapiro_test$p.value, "\n\n")
# }

####################################################
####################################################
####################################################

###########    PLOT SAVING    ###############
SavePrint_plot <- function(plot_obj, plot_name, dataset_name, save_dir, plot_size = c(7, 4), dpi = 300, font_size = 30) {
  # Create the directory if it doesn't exist
  if (!dir.exists(save_dir_plots)) {
    dir.create(save_dir_plots, recursive = TRUE)
  }
  # Modify the plot object to adjust the font size for jpg
  plot_obj_jpg <- plot_obj + theme(text = element_text(size = font_size, lineheight = .3))
  # Check if the directory is writable
  if (!file.access(save_dir_plots, 2)) {
    # Save plot as png
    ggsave(paste0(save_dir_plots, dataset_name, "_", plot_name, "_", Sys.Date(), ".png"), plot = plot_obj_jpg, width = plot_size[1], height = plot_size[2], dpi = dpi)
    # Save plot as pdf
    ggsave(paste0(save_dir_plots, dataset_name, "_", plot_name, "_", Sys.Date(), ".pdf"), plot = plot_obj, width = plot_size[1], height = plot_size[2])
    # Print the plot to the currently open device (the cumulative PDF file)
    print(plot_obj)
  } else {
    cat("Error: The directory is not writable.")
  }
}


######### STYLING FUNCTIONS ###########

## general
remove_y_labs <-  list(theme(axis.title.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.ticks.y = element_blank()))

remove_x_labs <-  list(theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank()))


# Function to split text into two lines
split_title <- function(title) {
  words <- str_split(title, " ")[[1]]
  half <- length(words) %/% 2
  paste(paste(words[1:half], collapse = " "), "\n", paste(words[(half+1):length(words)], collapse = " "), sep = "")
}


multi_plot_comps <- function(plot_list,ylim_extra){
  
  # Set the same y-axis limits for all plots
  y_limits <- c(min(sapply(plot_list, function(x) ggplot_build(x)$data[[1]]$ymin)), 
                max(sapply(plot_list, function(x) ggplot_build(x)$data[[1]]$ymax)) + ylim_extra)
  
  leg <-cowplot::get_legend(plot_list[[1]] + theme(legend.direction="horizontal"))
  
  multi_plot_comps_list <- list(y_limits=y_limits,leg=leg)
  
  return(multi_plot_comps_list)
}




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

##OVERWRITE THE COLOUR SHADING AS I NEVER USE THE INDIVIDUAL DATA POINTS ANYMORE
myColors_Treatment <- scales::viridis_pal(option = "D")(4)

myColors_Colony <- colour_palette$Shades
names(myColors_Colony) <- colour_palette$Cols
#myColors_Treatment <- Treat_colors$Shades
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
  scale_fill_manual(name = "Treatment", values = myColors_Treatment,labels = function(x) str_to_title(gsub("big", "large", gsub("\\.", " ", x)))) #for lines
colFill_treatment <-
  scale_fill_manual(name = "treatment", values = myColors_Treatment,labels = function(x) str_to_title(gsub("big", "large", gsub("\\.", " ", x)))) #for lines

colScale_period <-
  scale_colour_manual(name = "period",
                      values = myColors_period,
                      drop = TRUE,
                      labels = function(x) str_to_title(x))
colFill_period  <-
  scale_fill_manual(name = "period", values = myColors_period,labels = function(x) str_to_title(x)) #for lines


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
  