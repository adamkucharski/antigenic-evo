# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)

library(plot3D)
library(colorspace)


# LOAD Functions

setwd("~/Documents/antigenic-landscapes/immunity_model/")
source("model_functions.R")  #Set up functions



dataload="Fluscape_08"
#dataload="Australia_98"

landscape.build(dataload,0.5) # Generate maps

landscape.plot(dataload,radius1=4,2009)
titre.plot(dataload)
proct.plot(dataload)
