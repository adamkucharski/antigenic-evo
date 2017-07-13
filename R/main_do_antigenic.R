# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Main code to accompany
# Kucharski AJ, Baguelin M. The role of human immunity and social behavior in shaping influenza evolution. PLOS Pathogens
#

library(plot3D)
library(colorspace)
library(plotly)
library(magrittr)

rm(list=ls())

# LOAD Functions

setwd("~/Documents/antigenic-evo/R/")
source("model_functions.R")  #Set up functions

# Define dataset to load
dataload="China_08"

# - - - - - - - - - - - - 
# Generate antigenic maps
landscape.build(dataload,d.step=0.25,extendD=5,bandW=22) 


# Plot titre and reproduction number landscapes
landscape.plot(dataload,radius1=5,yearload=2009,groupN=2,borderA=F)

reproduction.number.plot(dataload,rR=2,borderA=F)
