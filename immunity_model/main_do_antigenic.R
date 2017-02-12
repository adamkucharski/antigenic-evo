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
#dataload="Australia_98"

# - - - - - - - - - - - - 
# Generate maps
landscape.build(dataload,d.step=0.25,extendD=5,bandW=22) 

# - - - - - - - - - - - - 
# Plot titre landscapes
landscape.plot(dataload,radius1=5,yearload=2009,groupN=2)

# Plot reproduction number landscapes
reproduction.number.plot(dataload,rR=2)

#titre.plot(dataload)
#proct.plot(dataload)


# Cross validate to explore bandwidth

store.crossVal = NULL
for(kk in seq(5,70,5)){
  rmsq=cross.validation(Data.load=dataload,d.step=0.5,extendD=0,bandW=kk, Nsamp = 10, bootstrap = 100)
  store.crossVal=cbind(store.crossVal,c(kk,rmsq))
}

write.csv(store.crossVal,"output_data/storeCross.csv")

# Plot cross validation

data <- read.csv("output_data/storeCross.csv")
data <- data[,-1]
plot(as.numeric(data[1,]),as.numeric(data[2,]))

