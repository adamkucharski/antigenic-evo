# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)
# Load data (Fonville et al.)

yr.load=98

options(StringsAsFactors=F)
data0=read.csv(paste("datasets/Australia_",yr.load,".csv",sep=""), as.is=T,head=F)


# - - - - - - - - - - - - - - - - 
# List test strains
data1=data0[,!data0[2,]=="POST"] # Focus on prevaccination
names(data1)=data1[1,] # Add names
data1=data1[-c(1,2),] # Remove duplicate row

nstrains=length(data1)-2 # remove subject and sample year
strain_names=names(data1)[3:(nstrains+2)]
test.index=c(1:nstrains)
npart=length(data1[,1])

# Convert to log titres and set missing data = NA
data1[data1=="*"]=NA
data1[data1=="<10"]=5
data1[,strain_names]=apply(data1[,strain_names],2,function(x){log2(as.numeric(x)/10)+1}) 


# - - - - - - - - - - - - - - - - 
# Set up antigenic coords

ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
ag.coord=ag.coord[match(strain_names,ag.coord$viruses),]



# Convert antigenic coords into cluster centroids
strain_years=as.numeric(sapply(strain_names,function(x){
  a1=max(which(strsplit(x, "")[[1]]=="/"))
  lstr=nchar(x)
  yr1=substr(x, a1+1, lstr)
  
  if(nchar(yr1)>4){yr1=substr(yr1, 1, 4)}
  year=yr1
  if(nchar(yr1)==2 & as.numeric(yr1)>15){year=paste("19",yr1,sep="")}
  if(nchar(yr1)==2 & as.numeric(yr1)<15){year=paste("20",yr1,sep="")}
  year
}
))

strain_unique=unique(strain_years)

strain_centre=data.frame(t(sapply(strain_unique,function(x){
  aa=ag.coord[!is.na(match(strain_years,x)),c("AG_x","AG_y")]
  c(x,colSums(aa)/length(aa[,1]))
})))
names(strain_centre)=c("year","AG_x","AG_y")

save(ag.coord,data1,npart,nstrains,strain_years,yr.load,strain_centre,file=paste("R_datasets/Australia_",yr.load,"_V.RData",sep=""))



# Load antigenic coordinate data



#ag.coord0=ag.coord[strain_years=="1997",]
#plot(ag.coord$AG_y,ag.coord$AG_x,col='white')
#text(ag.coord0$AG_y,ag.coord0$AG_x,ag.coord0$viruses,cex=1)

#dev.copy(pdf,paste("plots/antigenic_map.pdf",sep=""),width=20,height=10)
#dev.off()
