# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)


landscape.build<-function(Data.load,d.step=0.5){
  
  # Follows Section 1.2.3 of Supplement of Fonville et al (2015) Science
  
  load(paste("R_datasets/",Data.load,"_V.RData",sep=""))
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Construct matrix of map coords
  x.range=seq(floor(min(ag.coord$AG_x)),ceiling(max(ag.coord$AG_x)),d.step)
  y.range=seq(floor(min(ag.coord$AG_y)),ceiling(max(ag.coord$AG_y)+5),d.step)
  points.j=expand.grid(x.range,y.range) # Define list of points to evaluate
  names(points.j)=c("agx","agy")
  npointsj=length(points.j[,1])
  aA=20 # Define local bandwidth (set as 11 in paper)
  
  ag.weights=matrix(NA,nrow=npointsj,ncol=nstrains)
  
  
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
  age.gp1=c(1:npart)[age.list<=15] #define young age
  age.gp2=c(1:npart)[age.list>15 & age.list<=30] #define young age
  age.gp3=c(1:npart)[age.list>=30] #define old age
  
  group.list=list(age.gp1,age.gp2,age.gp3)
  group.names=c("<15","15-30","30+")
  n.groups=length(group.list)
  
  lm.data.group=NULL
  
  for(kk in 1:n.groups){ #Loop of groups of interest
    
    lm.d0=NULL
   
    for(ii in group.list[[kk]]){ # Gather data for all relevant participants
  
    p.data=data1[ii,3:(nstrains+2)]
    lm.data=data.frame(matrix(NA,nrow=nstrains,ncol=3))
    names(lm.data)=c("agx","agy","titre")
    
    lm.data$titre=as.numeric(p.data)
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
      
      summary(fit.model)
      predict(fit.model,newdata=points.j[pp,])
      
      # predict titre for that point
      pred.table[pp,]=c(points.j[pp,],predict(fit.model,newdata=points.j[pp,]))
    
    }
    
    pred.tableP=sapply(pred.table$pred.titre,function(x){min(max(x,0),8)})
    pred.matrix=matrix(pred.tableP,byrow=F,nrow=length(x.range))
    
    save(pred.matrix,lm.data,x.range,y.range,group.names,file=paste("R_datasets/maps/AGmap_",Data.load,"_gp",kk,".RData",sep=""))
  }
  
}



landscape.plot<-function(Data.load,radius1,yearload){
  
  load(paste("R_datasets/",Data.load,"_V.RData",sep=""))
  
  par(mfrow=c(3,1))
  
  for(kk in 1:3){
  
    load(paste("R_datasets/maps/AGmap_",Data.load,"_gp",kk,".RData",sep=""))  
  
    # Calculate mean titre for each sample strain
    group.size=length(lm.data$titre)/nstrains
    group_m1=matrix(lm.data$titre,nrow=group.size,byrow=T)
    mean.titre=apply(group_m1,2,function(x){mean(x[!is.na(x)])})
    
    # Plot antigenic surface
    par(mar = c(1,3,2,2))
    
    image2D(z = t(pred.matrix), x = y.range, y = x.range, zlim = c(0, 8),
            main=paste("Landscape 19",yr.load, ". Age ", group.names[kk],sep=""))
    
    points(ag.coord$AG_y,ag.coord$AG_x,cex=1.2*(mean.titre+1),col=rgb(0.1,0.1,0.1),lwd=2)
    ofs=0.05
    text(strain_centre$AG_y+ofs,strain_centre$AG_x-ofs,labels=strain_centre$year,col=rgb(0,0,0),cex=0.8)
    text(strain_centre$AG_y,strain_centre$AG_x,labels=strain_centre$year,col=rgb(1,1,1),cex=0.8)
    
    # Plot circle of test points centred on (yearload-1)
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
  
  dev.copy(png,paste("plots/antigenic_map",Data.load,".png",sep=""),width=1500,height=1800,res=150)
  dev.off()

}


titre.plot<-function(Data.load){
  
  col.palette=rainbow_hcl(3, c = 100, l = 60)
  
  par(mfrow=c(1,1))
  par(mar = c(5,5,2,2))
  
  for(kk in 1:3){
    
    load(paste("R_datasets/maps/predC_",Data.load,"_gp",kk,".RData",sep=""))  
    load(paste("R_datasets/maps/AGmap_",Data.load,"_gp",kk,".RData",sep=""))  
    
    if(kk==1){plot(pred1,type="l",ylim=c(0,8),lwd=2,col=col.palette[kk],ylab="titre")}
    else{lines(pred1,lwd=2,col=col.palette[kk])}
    
    
    text(x=round(0.85*length(pred1)),y=4-0.4*kk,group.names[kk],adj=0,col=col.palette[kk])
    
  }
  
  dev.copy(pdf,paste("plots/vaccine",Data.load,".pdf",sep=""),width=10,height=6)
  dev.off()
  
}

proct.plot<-function(Data.load){
  
  col.palette=rainbow_hcl(3, c = 100, l = 60)
  
  par(mfrow=c(1,1))
  par(mar = c(5,5,2,2))
  
  for(kk in 1:3){
    
    load(paste("R_datasets/maps/predC_",Data.load,"_gp",kk,".RData",sep=""))  
    load(paste("R_datasets/maps/AGmap_",Data.load,"_gp",kk,".RData",sep=""))  
    
    proct<-probability.protection(2^pred1*10)
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
