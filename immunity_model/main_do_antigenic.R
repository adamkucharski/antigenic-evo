# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)

library(plot3D)
library(colorspace)
library(plotly)
library(magrittr)

# LOAD Functions

setwd("~/Documents/antigenic-landscapes/immunity_model/")
#setwd("/Users/marcbaguelin/GitHub/antigenic-landscapes/immunity_model/")
source("model_functions.R")  #Set up functions



dataload="Fluscape_08"

landscape.build(dataload,d.step=0.25,extendD=5,bandW=22) # Generate maps

# Plot data
landscape.plot(dataload,radius1=5,yearload=2009,groupN=2,borderA=F)

reproduction.number.plot(dataload,rR=2,borderA=F)
