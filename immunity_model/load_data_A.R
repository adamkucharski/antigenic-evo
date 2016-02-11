# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)
# Load data (Fonville et al.)

options(StringsAsFactors=F)
data0=read.csv("datasets/Australia_98.csv", as.is=T,head=F)

# List test strains
data1=data0[,!data0[2,]=="POST"] # Focus on prevaccination
names(data1)=data1[1,] # Add names
data1=data1[-c(1,2),] # Remove duplicate row

nstrains=length(data1)-2 # remove subject and sample year
strain_names=names(data1)[3:(nstrains+2)]
test.index=c(1:nstrains)

# Convert to log titres and set missing data = NA
data1[data1=="*"]=NA
data1[,strain_names]=apply(data1[,strain_names],2,function(x){log2(as.numeric(x)/10)+1}) 


# Gather participants and infection years

n_part=max(data1$Subject.number) # number of participants

test_years=unique(data1$Sample.year) # year of testing
test.n=length(test_years) # number of test years

inf_years=seq(min(strain_years),max(c(test_years,strain_years))) #annual infection model
inf.n=length(inf_years) # number of possible infecting strains

# Set up list of test data for quick access

data.Test=data1[,strain_names]
test.list=list()

for(ii in 1:n_part){

subjectn=ii
i.list=list()

for(jj in 1:test.n){

testyr=test_years[jj]
dataI=data.Test[data1$Subject.number==subjectn & data1$Sample.year==testyr,]
i.list[[jj]]=rbind(rep(testyr,nstrains),
      dataI[,!is.na(dataI)],
      strain_years[!is.na(dataI)],
      strain_years[test.index[!is.na(dataI)]]-min(strain_years)+1
      )

}

test.list[[ii]]=i.list

save(test_years,inf_years,strain_years,n_part,test.list,file=paste("R_datasets/HaNam_data_V.RData",sep=""))

}