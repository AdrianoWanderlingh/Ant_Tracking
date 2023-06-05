####
#  PLOTS, TO DO:
#  - CALC THE DELTA OF THE DATA
#  - ASSIGN POST-HOC LETTERS ONLY IF PRESENT
#  - 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

###################################################################
##########             I    collective analysis        ############
###################################################################

collective = T
type <-  "collective"
predictor <- "treatment"
diff_type <- "absolute_difference"




diff <- create_diff(dataset=data,predictor="treatment",type="collective",form_stat=NULL,collective=F,plot_untransformed=F,diff_type="absolute_difference")


























data1 <- data.frame(lapply(data, function(x) if(is.numeric(x)) round(x, 1) else x))

dput(head(data1))
