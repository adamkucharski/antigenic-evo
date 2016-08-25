# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Build antigenic landscape from data 

landscape.build<-function(Data.load,d.step=0.5,extendD=3,bandW=20){
  
  # Follows Section 1.2.3 of Supplement of Fonville et al (2015) Science
  # Data.load = dataload ; d.step=0.5 ; extendD=3 ; bandW=20
  
  load(paste("R_datasets/",Data.load,"_V.RData",sep=""))
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Construct matrix of map coords
  x.range=seq(floor(min(ag.coord$AG_x)-1),ceiling(max(ag.coord$AG_x)),d.step)
  y.range=seq(floor(min(ag.coord$AG_y)),ceiling(max(ag.coord$AG_y)+extendD),d.step)
  points.j=expand.grid(x.range,y.range) # Define list of points to evaluate
  names(points.j)=c("agx","agy")
  npointsj=length(points.j[,1])
  aA=bandW # Define local bandwidth (set as 11 in paper)
  
  ag.weights=matrix(NA,nrow=npointsj,ncol=nstrains)
  
  # How much weighting to give to strains in the fitting
  for(ii in 1:npointsj){
    for(jj in 1:nstrains){
      a_ij=sqrt((points.j[ii,"agx"]-ag.coord[jj,"AG_x"])^2+(points.j[ii,"agy"]-ag.coord[jj,"AG_y"])^2)
      ag.weights[ii,jj]=ifelse(a_ij<=aA,(1-(a_ij/aA)^3)^3,0) #a_ij #
    }
  }
  
  image(ag.weights)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Set up fitting data set for group of participants
  age.list=as.numeric(data1[,2])
  #age.gp1=c(1:npart)[age.list<=15] #define young age
  #age.gp2=c(1:npart)[age.list>15 & age.list<=30] #define young age
  #age.gp3=c(1:npart)[age.list>=30] #define old age
  #group.list=list(age.gp1,age.gp2,age.gp3)
  #group.names=c("<15","15-30","30+")
  
  age.gp1=c(1:npart)[age.list<=19] #define young age
  age.gp2=c(1:npart)[age.list>19] #define young age
  group.list=list(age.gp1,age.gp2)
  group.names=c("<20","20+")
  n.groups=length(group.list)
  
  lm.data.group=NULL
  
  for(kk in 1:n.groups){ #Loop of groups of interest
    
    lm.d0=NULL
   
    for(ii in group.list[[kk]]){ # Gather data for all relevant participants
  
    p.data=data1[ii,3:(nstrains+2)]
    lm.data=data.frame(matrix(NA,nrow=nstrains,ncol=3))
    names(lm.data)=c("agx","agy","titre")
    
    lm.data$titre=as.numeric(p.data)
    #lm.data$titre=probability.protection(2^as.numeric(p.data)*10)
    lm.data$agx=ag.coord$AG_x
    lm.data$agy=ag.coord$AG_y
    
    lm.d0=rbind(lm.d0,lm.data)
  }
  
  lm.data.group[[kk]]=lm.d0
    
  }
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Estimate antigenic surface with linear model
  
  for(kk in 1:n.groups){ #Loop over groups of interest
  
    #pick a point from points.j
    pred.table=data.frame(matrix(NA,nrow=npointsj,ncol=3))
    names(pred.table)=c("agx","agy","pred.titre")
    lm.data=lm.data.group[[kk]]
    
    for(pp in 1:npointsj){
    
      agW=rep(ag.weights[pp,],length(group.list[[kk]]))
      
      fit.model<-lm(titre ~ agx*agy,data=lm.data,weights=agW)
      #summary(fit.model); predict(fit.model,newdata=points.j[pp,])
      
      # predict titre for that point
      pred.table[pp,]=c(points.j[pp,],predict(fit.model,newdata=points.j[pp,]))
    
    }
    
    pred.tableP=sapply(pred.table$pred.titre,function(x){min(max(x,0),8)})
    pred.matrix=matrix(pred.tableP,byrow=F,nrow=length(x.range))
    
    save(pred.matrix,lm.data,x.range,y.range,group.names,file=paste("R_datasets/maps/AGmap_",Data.load,"_gp",kk,".RData",sep=""))
  }
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Cross-validation of bandwidth of map

cross.validation <- function(Data.load,d.step=0.5,extendD=3,bandW=20, Nsamp = 10, bootstrap = 5){
  
  # Follows Section 1.2.5 of Supplement of Fonville et al (2015) Science
  # Data.load = dataload ; d.step=1 ; extendD=3 ; bandW=20; bootstrap = 5; Nsamp = 10
  
  load(paste("R_datasets/",Data.load,"_V.RData",sep=""))
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Construct matrix of map coords
  x.range=seq(floor(min(ag.coord$AG_x)-1),ceiling(max(ag.coord$AG_x)),d.step)
  y.range=seq(floor(min(ag.coord$AG_y)),ceiling(max(ag.coord$AG_y)+extendD),d.step)
  points.j=expand.grid(x.range,y.range) # Define list of points to evaluate
  names(points.j)=c("agx","agy")
  npointsj=length(points.j[,1])
  aA=bandW # Define local bandwidth (set as 11 in paper)
  
  ag.weights=matrix(NA,nrow=npointsj,ncol=nstrains)
  
  
  for(ii in 1:npointsj){
    for(jj in 1:nstrains){
      a_ij=sqrt((points.j[ii,"agx"]-ag.coord[jj,"AG_x"])^2+(points.j[ii,"agy"]-ag.coord[jj,"AG_y"])^2)
      ag.weights[ii,jj]=ifelse(a_ij<=aA,(1-(a_ij/aA)^3)^3,0) #a_ij #
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Sample 10 random individuals and calculate titres
  
  lm.data.group=NULL
  partN = length(data1$subject)
  
  for(kk in 1:bootstrap){
    picks = sort(sample(1:partN,Nsamp))
    lm.d0=NULL
  
    for(ii in picks){ # Gather data for all relevant participants
        
      p.data=data1[ii,3:(nstrains+2)]
      lm.data=data.frame(matrix(NA,nrow=nstrains,ncol=3))
      names(lm.data)=c("agx","agy","titre")
        
      lm.data$titre=as.numeric(p.data)
      lm.data$agx=ag.coord$AG_x; lm.data$agy=ag.coord$AG_y
        
      lm.d0=rbind(lm.d0,lm.data)
      lm.data.group[[kk]]=lm.d0
      
    }
  
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Estimate antigenic surface with linear model and calculate error
  
  rmsq.error.tab = NULL
  
  for(kk in 1:bootstrap){ #Loop over boostrap groups
    
    #pick a point from points.j
    pred.table=data.frame(matrix(NA,nrow=npointsj,ncol=3))
    names(pred.table)=c("agx","agy","pred.titre")
    lm.data=lm.data.group[[kk]]
    
    # pick training data
    option.L = c(1:length(lm.data$titre))
    train.picks = sapply(sort(sample(1:9,6)),function(x){ seq(x,9*Nsamp,9) }) # Need to pick groups of viruses
    dim(train.picks)=NULL; train.picks = sort(train.picks)
    test.picks = option.L[-train.picks]
    subset.data = lm.data.group
    
    for(pp in 1:npointsj){
      
      agW=rep(ag.weights[pp,],length(picks))
      agW = agW[train.picks] # Only pick weights for training strains
      fit.model<-lm(titre ~ agx*agy,data=lm.data[train.picks,],weights=agW)
      #summary(fit.model); predict(fit.model,newdata=points.j[pp,])
      
      # predict titre for that point
      pred.table[pp,]=c(points.j[pp,],predict(fit.model,newdata=points.j[pp,]))
      
    }
    
    pred.tableP=sapply(pred.table$pred.titre,function(x){min(max(x,0),8)})
    pred.matrix=matrix(pred.tableP,byrow=F,nrow=length(x.range))
    
    # calculate error in test data by matching location in map
    test.data = lm.data[test.picks,]
    xMatch = sapply(test.data$agx,function(x){ c(1:length(x.range))[abs(x-x.range)==min(abs(x-x.range)) ] })
    yMatch = sapply(test.data$agy,function(x){ c(1:length(y.range))[abs(x-y.range)==min(abs(x-y.range)) ] })

    p.test = NULL
    for(ii in 1:length(test.picks)){
      
      p.test = c(p.test, (test.data[ii,"titre"]-pred.matrix[xMatch[ii],yMatch[ii]])^2)# Pick out nearest point in matrix
    }
    
    rmsq.error = sqrt(sum(p.test))
    rmsq.error.tab = c(rmsq.error.tab,rmsq.error)

    #save(pred.matrix,lm.data,x.range,y.range,group.names,file=paste("R_datasets/maps/ValMap_",Data.load,"_gp",kk,".RData",sep=""))
  }
  
  mean(rmsq.error.tab)
  
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot titre landscape (load from above function)

landscape.plot<-function(Data.load,radius1,yearload,groupN=3,circleShow=F){
  
  # Data.load=dataload; radius1=4; yearload=2009; groupN=2; circleShow=F
  
  load(paste("R_datasets/",Data.load,"_V.RData",sep=""))
  
  par(mfrow=c(groupN,1))
  par(mgp=c(1.8,0.6,0))
  
  for(kk in 1:groupN){
  
    load(paste("R_datasets/maps/AGmap_",Data.load,"_gp",kk,".RData",sep=""))  
  
    # Calculate mean titre for each sample strain
    group.size=length(lm.data$titre)/nstrains
    group_m1=matrix(lm.data$titre,nrow=group.size,byrow=T)
    mean.titre=apply(group_m1,2,function(x){mean(x[!is.na(x)])})
    
    # Plot antigenic surface

    if(kk<groupN){
      par(mar = c(2,4,2,2))
    }else{
      par(mar = c(3,4,2,2))
    }
    MTx=c(332,372)
    MTy=c(245,262)
    
    image2D(z = t(pred.matrix), x = y.range, y = x.range, xlab=ifelse(kk<groupN,"","antigenic dimension 1"), ylab="antigenic dimension 2", zlim = c(0, 8),
            main=paste("Age ", group.names[kk],sep=""),col=rev(ramp.col (col = c("blue",rgb(0.4,0.6,1),"white"), n = 100, alpha = 1))) #paste("Landscape ",Data.load, ". Age ", group.names[kk],sep="")
    
    #title(main=LETTERS[kk],adj=0)
    
    ag.coordALL=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
    #vals1=predict(am.spl,scalemap(inf_years))
    ag.coordALL=ag.coordALL[ag.coordALL$AG_y<=max(y.range) & ag.coordALL$AG_y>=min(y.range)&
                ag.coordALL$AG_x<=max(x.range) & ag.coordALL$AG_x>=min(x.range),]
    
    # All points
    #points(ag.coordALL$AG_y,ag.coordALL$AG_x,xlim=c(min(y.range),max(y.range)),col="black",xlab="strain dimension 1", ylab="strain dimension 2",pch=19)
    
    # Points of interest
    #points(ag.coord$AG_y,ag.coord$AG_x,cex=1.2*(mean.titre+1),col=rgb(0.1,0.1,0.1),lwd=2)
    ofs=0.03
    tx1=strain_centre$AG_y
    tx1[1]=tx1[1]+1
    text(tx1+ofs,strain_centre$AG_x-ofs,labels=strain_centre$year,col="white",cex=1.2)
    text(tx1+2*ofs,strain_centre$AG_x-2*ofs,labels=strain_centre$year,col="white",cex=1.2)
    text(tx1,strain_centre$AG_x,labels=strain_centre$year,col="black",cex=1.2)
    
    
    # Plot circle of test points centred on (yearload-1)
    if(circleShow==T){
      c.pts=c(1:100); c.tot=length(c.pts)
      yr.test=as.numeric(yearload)
      yr.testC=strain_centre[strain_centre==(yr.test-1),]
      c.test=c(AG_x=yr.testC$AG_x,AG_y=yr.testC$AG_y,rad=radius1)
      
      circle.coord=t(sapply(c.pts, function(x){
        c(c.test[["rad"]] *cos(2*pi*x/c.tot)+c.test[["AG_y"]],c.test[["rad"]] *sin(2*pi*x/c.tot)+c.test[["AG_x"]])
      }))
      
      lines(circle.coord,pch=19,col='white',lwd=2,lty=2)
      
      # Match points and calculate titres
      pred1=apply(circle.coord,1,function(zz){
        ydist=abs(y.range-zz[1])
        xdist=abs(x.range-zz[2])
        pred.matrix[xdist==min(xdist),ydist==min(ydist)]
      })
  
      save(circle.coord,pred1,file=paste("R_datasets/maps/predC_",Data.load,"_gp",kk,".RData",sep=""))
    }
  }
  
  dev.copy(png,paste("plots/antigenic_map",Data.load,".png",sep=""),width=1500,height=1800,res=190)
  dev.off()

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot reproduction number landscape (load from above function)

build.china.matrix<-function(r0=2){
  
  contacts <- data.frame(read.csv("datasets/fluscape_participant_contact_aggr_v1b.csv",stringsAsFactors = F))
  npart <- length(contacts$age)
  
  age.list=as.numeric(contacts$age)
  
  # Define mixing matrix
  n.u20 = length(contacts[age.list<=19,"age"])
  n.o20 = length(contacts[age.list>19,"age"])
  contact.u20 = c(sum(contacts[age.list<=19,c("c.age.under5","c.age.6to19")]),sum(contacts[age.list<=19,c("c.age.20to64","c.age.over65")]))/n.u20
  contact.o20 = c(sum(contacts[age.list>19,c("c.age.under5","c.age.6to19")]),sum(contacts[age.list>19,c("c.age.20to64","c.age.over65")]))/n.o20
  
  r.matrix <- matrix(c(contact.u20,contact.o20),nrow=2)
  rage.matrix <- r0*r.matrix/max(eigen(r.matrix)$values)
  
  rage.matrix
}

proct <- function(x){probability.protection(2^x*10)}

strain_years <- function(){
  
  ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)

  # Convert antigenic coords into cluster centroids
  strain_years=as.numeric(sapply(ag.coord$viruses,function(x){
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
  
}

reproduction.number.plot<-function(Data.load){
  
  r.matrix <- build.china.matrix(2)

  # Data.load=dataload; radius1=4; yearload=2009; groupN=2; circleShow=F
  
  load(paste("R_datasets/",Data.load,"_V.RData",sep=""))
  
  
  load(paste("R_datasets/maps/AGmap_",Data.load,"_gp",1,".RData",sep=""))  
  pred.matrix.U20 <- proct(pred.matrix)
  load(paste("R_datasets/maps/AGmap_",Data.load,"_gp",2,".RData",sep=""))  
  pred.matrix.O20 <- proct(pred.matrix)
  
  # Compile R matrix
  pred.matrixR <- NA+pred.matrix*0
  contrib.mat <- NA+pred.matrix*0
  for(xx in 1:length(x.range)){
    for(yy in 1:length(y.range)){
      
      s.u20 <- pred.matrix.U20[xx,yy]
      s.o20 <- pred.matrix.O20[xx,yy]
      rmat <- r.matrix*matrix(1-c(s.u20,s.u20,s.o20,s.o20),ncol=2,byrow=T)
      pred.matrixR[xx,yy] <- max(eigen(rmat)$values)
      contrib.mat[xx,yy] <- sum(rmat[,1])/sum(rmat[,2])
      
    }
  }

  # Calculate mean titre for each sample strain
  group.size=length(lm.data$titre)/nstrains
  group_m1=matrix(lm.data$titre,nrow=group.size,byrow=T)

  # - - - - - 
  # Plot antigenic surface
  par(mfrow=c(1,1))
  par(mar = c(4,4,2,2))
  #par(mgp=c(1.8,0.6,0))
  MTx=c(332,372)
  MTy=c(245,262)
    
  image2D(z = t(pred.matrixR), x = y.range, y = x.range, xlab="antigenic dimension 1", ylab="antigenic dimension 2", zlim = c(0, 2),
            main=paste("",sep=""),labels="") #paste("Landscape ",Data.load, ". Age ", group.names[kk],sep="")

 # image2D(z = t(contrib.mat), x = y.range, y = x.range, contour = list(levels = 1, col = "black", lwd = 2) , xlab="antigenic dimension 1", ylab="antigenic dimension 2", zlim = c(0, 3),
  #        main=paste("",sep="")) #paste("Landscape ",Data.load, ". Age ", group.names[kk],sep="")
  
  
  ag.coordALL=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)

  #vals1=predict(am.spl,scalemap(inf_years))
  #ag.coordALL=ag.coordALL[ag.coordALL$AG_y<=max(y.range) & ag.coordALL$AG_y>=min(y.range)&
  #                            ag.coordALL$AG_x<=max(x.range) & ag.coordALL$AG_x>=min(x.range),]
    
  # All locations
  #
  #plot(ag.coordALL$AG_y,ag.coordALL$AG_x,xlim=c(min(y.range),max(y.range)),col="black",xlab="strain dimension 1", ylab="strain dimension 2",pch=19)
  
  # - - - - - 
  # Post 2008 strain locations

  # Select post XX year strains
  strainY <- strain_years()
  years.plot=c(1968:2011)
  ag.coordPick = NULL;ag.coordPickY=NULL
  for(ii in 1:length(ag.coordALL$viruses)){
    yy=if(sum(strainY[ii]==years.plot)>0){
      ag.coordPick=rbind(ag.coordPick,ag.coordALL[ii,])
      ag.coordPickY=rbind(ag.coordPickY,strainY[ii])
    }
  }

  points(ag.coordPick[,c("AG_y","AG_x")],pch=19,col="black")
  #text(ag.coordP2008$AG_y,ag.coordP2008$AG_x,labels=ag.coordP2008$viruses)
  
  # Plot vaccine selection?   
  #vaccine2008 = ag.coordALL[match(c("A/BRISBANE/10/2007","A/URUGUAY/716/2007"),ag.coordALL$viruses),c("AG_y","AG_x")]
  #vaccine2009 = ag.coordALL[match("A/PERTH/16/2009",ag.coordALL$viruses),c("AG_y","AG_x")]
  vaccine2010 = ag.coordALL[match("A/VICTORIA/361/2011",ag.coordALL$viruses),c("AG_y","AG_x")]
  #points(vaccine2008[,1],vaccine2008[,2],col="white",pch=15,cex=2)
  #points(vaccine2009[1],vaccine2009[2],col="white",pch=19,cex=2)
  points(vaccine2010[1],vaccine2010[2],col="white",pch=17,cex=2)
  
  # Add year text labels
  #points(ag.coord$AG_y,ag.coord$AG_x,cex=1.2*(mean.titre+1),col=rgb(0.1,0.1,0.1),lwd=2)
  ofs=0.03
  tx1=strain_centre$AG_y
  tx1[1]=tx1[1]+1
  text(tx1+ofs,strain_centre$AG_x-ofs,labels=strain_centre$year,col=rgb(0,0,0),cex=1.2)
  text(tx1+2*ofs,strain_centre$AG_x-2*ofs,labels=strain_centre$year,col=rgb(0,0,0),cex=1.2)
  text(tx1,strain_centre$AG_x,labels=strain_centre$year,col=rgb(1,1,1),cex=1.2)
  
  #title(main=LETTERS[3],adj=0)
  dev.copy(png,paste("plots/reproduction_number_map",Data.load,".png",sep=""),width=1500,height=1200,res=190)
  dev.off()
  
  #
  pred1=apply(ag.coordPick,1,function(zz){
    ydist=abs(y.range-as.numeric(zz[3]))
    xdist=abs(x.range-as.numeric(zz[2]))
    mean(pred.matrixR[xdist==min(xdist),ydist==min(ydist)])
  })
  
  plot(ag.coordPickY,pred1,col="blue",xlim=c(1968,2012))
  v2008 = match(c("A/BRISBANE/10/2007","A/URUGUAY/716/2007"),ag.coordPick$viruses)
  v2010 = match("A/PERTH/16/2009",ag.coordPick$viruses)
  points(ag.coordPickY[v2010],pred1[v2010],col="black",pch=19,xlim=c(1968,2015))
  points(ag.coordPickY[v2008],pred1[v2008],col="black",pch=17,xlim=c(1968,2015))
  
  #text(ag.coordPickY,pred1,labels=ag.coordPick$viruses,col="black",cex=0.8)

  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Plot protection by age

proct.plot<-function(Data.load){
  
  col.palette=rainbow_hcl(3, c = 100, l = 60)
  
  par(mfrow=c(1,1))
  par(mar = c(5,5,2,2))
  
  for(kk in 1:3){
    
    load(paste("R_datasets/maps/predC_",Data.load,"_gp",kk,".RData",sep=""))  
    load(paste("R_datasets/maps/AGmap_",Data.load,"_gp",kk,".RData",sep=""))  
    
    proct <- probability.protection(2^pred1*10)
    if(kk==1){plot(proct,type="l",ylim=c(0,1),lwd=2,col=col.palette[kk],ylab="AB mediated protection")}
    else{lines(proct,lwd=2,col=col.palette[kk])}
    
    
    text(x=round(0.85*length(proct)),y=4-0.4*kk,group.names[kk],adj=0,col=col.palette[kk])
    
  }
  
  dev.copy(pdf,paste("plots/vaccine",Data.load,".pdf",sep=""),width=10,height=6)
  dev.off()
  
}

#antibody protection function for an individual with a given titre T
#a (alpha) and b (beta) are the log(50% PT) and steepness of the function
#Defaults values are taken from Coudeville et al. (model ALL) 
probability.protection<-function(T,a=2.844,b=1.299)
{
  return(1-1/(1+exp(b*(log(T)-a))))
}
