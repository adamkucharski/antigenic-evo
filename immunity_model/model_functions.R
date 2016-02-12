# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Code by Adam Kucharski (2016)


landscape.plot<-function(){
  
  # Follows Section 1.2.3 of Supplement of Fonville et al (2015) Science
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Construct matrix of map coords
  x.range=seq(floor(min(ag.coord$AG_x)),ceiling(max(ag.coord$AG_x)),0.5)
  y.range=seq(floor(min(ag.coord$AG_y)),ceiling(max(ag.coord$AG_y)),0.5)
  points.j=expand.grid(x.range,y.range) # Define list of points to evaluate
  names(points.j)=c("agx","agy")
  npointsj=length(points.j[,1])
  aA=11 # Define local bandwidth (set as 11 in paper)
  
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
  age.gp1=c(1:npart)[age.list<=20] #define young age
  age.gp2=c(1:npart)[age.list>20 & age.list<=50] #define young age
  age.gp3=c(1:npart)[age.list>=50] #define old age
  
  group.list=list(age.gp1,age.gp2,age.gp3)
  group.names=c("18-20","20-50","50+")
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
    
    # Plot antigenic surface
    par(mar = c(4,4,2,2))
    
    image2D(z = t(pred.matrix), x = y.range, y = x.range, zlim = c(0, 8),
            main=paste("Landscape #",kk, ". Age ", group.names[kk],sep=""))
    
    points(lm.data$agy,lm.data$agx,cex=(lm.data$titre+1),col=rgb(0,0,0,0.1))
    
    dev.copy(png,paste("plots/antigenic_map",kk,".png",sep=""),width=1500,height=800,res=150)
    dev.off()
    
  }
  
  #npart
  #nstrains
  



}