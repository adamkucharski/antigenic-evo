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
npart=length(data1[,1])

# Convert to log titres and set missing data = NA
data1[data1=="*"]=NA
data1[data1=="<10"]=5
data1[,strain_names]=apply(data1[,strain_names],2,function(x){log2(as.numeric(x)/10)+1}) 

ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
ag.coord=ag.coord[match(strain_names,ag.coord$viruses),]

save(ag.coord,data1,npart,nstrains,file=paste("R_datasets/Australia_98_V.RData",sep=""))


# Load antigenic coordinate data




#plot(ag.coord$AG_y,ag.coord$AG_x,col='white')
#text(ag.coord$AG_y,ag.coord$AG_x,ag.coord$viruses,cex=1)

#dev.copy(pdf,paste("plots/antigenic_map.pdf",sep=""),width=20,height=10)
#dev.off()
