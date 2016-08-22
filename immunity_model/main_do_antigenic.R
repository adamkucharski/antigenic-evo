# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)

library(plot3D)
library(colorspace)


# LOAD Functions

setwd("~/Documents/antigenic-landscapes/immunity_model/")
#setwd("/Users/marcbaguelin/GitHub/antigenic-landscapes/immunity_model/")
source("model_functions.R")  #Set up functions



dataload="Fluscape_08"
#dataload="Australia_98"

landscape.build(dataload,0.5,extendD=7,bandW=25) # Generate maps

# Plot data
landscape.plot(dataload,radius1=4,yearload=2009,groupN=2)

reproduction.number.plot(dataload)

#titre.plot(dataload)
#proct.plot(dataload)


# Cross validate to find bandwidth

store.crossVal = NULL
for(kk in seq(26,30,1)){
  rmsq=cross.validation(Data.load=dataload,d.step=0.5,extendD=0,bandW=kk, Nsamp = 30, bootstrap = 200)
  store.crossVal=cbind(store.crossVal,c(kk,rmsq))
}

write.csv(store.crossVal,"output_data/storeCross.csv")