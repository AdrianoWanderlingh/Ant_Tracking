library(BSDA)
library(car)
library(coxme)
library(e1071)
library(extrafont)
loadfonts()
library(Hmisc)
library(gtools)
library(kSamples) 
library(igraph)
library(lme4)
library(lmerTest)
library(MASS)
library(measurements)
library(metap)
library(multcomp)
library(nlme)
library(plotrix)
library(pls)
library(RColorBrewer)
library(scales)
library(spatstat)
library(stringr)
library(survival)
library(survcomp)
library(vioplot)
library(viridisLite)
# AW
require(report)
require(sjPlot)
require(lme4)
require(car)
require(MuMIn)
require(diptest) #dip test for bimodularity

# in case of installation issues, use mamba install r-PACKAGENAME

### AW # for plotting
require(RColorBrewer)
require(shades)
require(colorspace)
require(plotwidgets)
require(ggplot2)
require(ggnewscale)
require(gridExtra)
require(extrafont)
require(cowplot)
library(Cairo) # to ensure working export to pdf with non-standard fonts
#font_add_google("Crimson Text", "crimson", db_cache = FALSE)
require(dplyr)
library(forcats)  # for fct_rev() reverses order of factors for plotting needs
#font_import()
#fonts()
library(grid)
library(gridBase)
library(gridGraphics)

######library(showtext)
######showtext_auto()#must be called to indicate that showtext is going to be automatically invoked to draw text whenever a plot is created.


## to install a specifi version of a package when "Warning in install.packages :package ‘ XXX’ is not available for this version of R":
# library(remotes)
# remotes::install_version("MuMIn", version = "1.46.0")
# remotes::install_version("sjPlot", version = "1.46.0")


# To add the "Liberation Serif" font to the `pdfFonts` in Linux, you can use the `extrafont` package in R. Here are the steps:
#   
#   1. **Install and load the `extrafont` package**
#   
# install.packages("extrafont")
# library(extrafont)
# 2. **Import system fonts**
#   
# font_import()
# 
# This function will import all fonts on your system into R. It may take a while to complete.
# 
# 3. **Check if "Liberation Serif" has been imported**
#   
# fonts()
# 
# This function will list all the fonts that have been imported into R. Check if "Liberation Serif" is in the list.
# 
# 4. **Load the fonts for the PDF device**
#   
# loadfonts(device = "pdf")
# Now, "Liberation Serif" should be available to the PDF device in R.
# CHECK THAT > loadfonts(device = "pdf") DOESN'T RETURN:
# More than one version of regular/bold/italic found for Liberation Serif. Skipping setup for this font.
