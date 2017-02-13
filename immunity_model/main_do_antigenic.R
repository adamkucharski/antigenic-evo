# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016-)

library(plot3D)
library(colorspace)
library(plotly)
library(magrittr)

# LOAD Functions

setwd("~/Documents/antigenic-landscapes/immunity_model/")
#setwd("/Users/marcbaguelin/GitHub/antigenic-landscapes/immunity_model/")
source("model_functions.R")  #Set up functions

# Define dataset to load
dataload="Fluscape_08"

# - - - - - - - - - - - - 
# Generate maps
landscape.build(dataload,d.step=0.25,extendD=5,bandW=22) 


# Plot titre and reproduction number landscapes
landscape.plot(dataload,radius1=5,yearload=2009,groupN=2,borderA=F)

reproduction.number.plot(dataload,rR=2,borderA=F)
