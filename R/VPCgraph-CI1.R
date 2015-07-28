#' draw the original visual predictive check plot
#' 
#' @import ggplot2 plyr
#' @title Visual predictive checks
#' @param orig.data NONMEM data 
#' @param sim.data simulated data from NONMEM
#' @param N.timebin number of time bin
#' @param N.sim number of simulation
#' @param q.list list of quantiles for VPC plot
#' @param alpha significance level of CI for each quantile
#' @param X.name x label in VPC plot
#' @param Y.name y label in VPC plot
#' @param main.title title of plot
#' @param opt.DV.point option for drawing data points
#' @param opt.DV.quantile.line option for drawing quantiles of the original data
#' @param opt.SIM.quantile.line option for drawing quantiles of simulated data
#' @param opt.SIM.quantile.CI.area option for drawing confidence area of 
#'                                  quantiles for simulated data
#' @param Y.min minimum of y axis in VPC plot
#' @param Y.max maximum of y axis in VPC plot
#' @param plot.flag draw plot if TRUE; generate data for drawing plot if FALSE
#' @return plot or the values to draw plot
#' @export
#' @examples
#' data(origdata)
#' data(simdata)
#' VPC.graph.CI(origdata,simdata,10,100)

VPC.graph.CI1<-function(orig.data,sim.data,N.timebin=NULL,N.sim,
                       q.list=c(0.05,0.5,0.95),
                       alpha=0.05,
                       X.name="TIME",Y.name="DV", ID.name="ID",
                       main.title=NULL,
                       opt.DV.point=FALSE,
                       opt.DV.quantile.line=TRUE,
                       opt.SIM.quantile.line=FALSE,
                       opt.SIM.quantile.CI.area=TRUE,
                       Y.min=NULL,Y.max=NULL,
                       X.min=NULL,X.max=NULL,                   
                       plot.flag=TRUE,
                       Y0.adjust=TRUE,
                       maxK=NULL,beta=0.2,lambda=0.3,
                       CIarea.box=FALSE,
                       bin.grid=TRUE,orig.flag=TRUE,...){
  SIM.CIarea.1<-NULL
  SIM.CIarea.2<-NULL
  SIM.CIarea.3<-NULL
  DV.point<-NULL
  DV.quant<-NULL
  SIM.quant<-NULL 
  ID<-NULL;G<-NULL
  
  optK<-function(orig.data,X.name="TIME",Y.name="DV",ID.name="ID",
                 maxK=NULL,beta=0.2,lambda=0.3){
    plot.data<-data.frame(X=orig.data[,X.name],
                          Y=orig.data[,Y.name],ID=orig.data[,ID.name])
    if(is.null(maxK)) 
      maxK <- max(c(round(nrow(orig.data)/10),
                    tapply(plot.data$X,plot.data$ID,length)))
    keep.tot<-NULL
    for(k in 2:maxK){
      time.bin<-makeCOVbin(plot.data$X,k)  
      J.temp<-sum((tapply(plot.data$X,time.bin$COV.bin,var)^beta)*
                    table(time.bin$COV.bin))
      temp.U<-ifelse(J.temp!=0,log(J.temp),0)+lambda*beta*k
      keep.tot<-rbind(keep.tot,c(k,temp.U))
    }
    colnames(keep.tot)<-c("k","U")
    keep.tot<-data.frame(keep.tot)
    return(keep.tot$k[which(keep.tot$U==min(keep.tot$U))])
  }
  
  findQuantile<-function(DV.data,TIME.data,time.bin.temp,q){
    TIME<-NULL
    DV.quantile<-NULL
    for(i in 1:nrow(time.bin.temp$COV.bin.summary)){
      if(time.bin.temp$COV.bin.summary$n.bin[i]!=0){
        sel.id<-which(time.bin.temp$COV.bin==
                        time.bin.temp$COV.bin.summary[i,1])
        DV.temp<-c(as.matrix(DV.data)[sel.id,])
        DV.quantile<-rbind(DV.quantile,quantile(DV.temp,prob=q,na.rm=TRUE))        
      } else{
        DV.quantile<-rbind(DV.quantile,NA) 
      }
    }
    DV.q<-data.frame(time.bin.temp$COV.bin.summary,DV.quantile)
    colnames(DV.q)[-(1:(ncol(DV.q)-length(q)))]<-
      paste("Q",round(q*100),"th",sep="")
    return(DV.q)
  }
  
  findQuantileCI<-function(sim.data,TIME.data,time.bin.temp,q,alpha,N.sim){
    orig.sim.Q.temp<-apply(sim.data,2,function(x) 
      findQuantile(x,TIME.data,time.bin.temp,q)[,
                                                -(1:ncol(time.bin.temp$COV.bin.summary))])
    keep.name<-NULL#names(orig.sim.Q.temp[[1]])
    orig.sim.Q.temp<-array(unlist(orig.sim.Q.temp),
                           dim=c(nrow(time.bin.temp$COV.bin.summary),
                                 length(q),N.sim))
    Q.CI<-list()
    for(i in 1:length(q)){
      Q.CI[[i]]<-cbind(findQuantile(sim.data,TIME.data,
                                    time.bin.temp,q[i])[,4],
                       t(apply(orig.sim.Q.temp[,i,],1,
                               function(x) quantile(x,
                                                    probs=c(alpha/2,1/2,1-alpha/2),
                                                    na.rm=TRUE))))
      rownames(Q.CI[[i]])<-time.bin.temp$COV.bin.summary[,1]    
      keep.name<-c(keep.name,paste("Q",round(q[i]*100),"th",sep=""))
      colnames(Q.CI[[i]])[1]<-keep.name[i]
    }  
    names(Q.CI)<-keep.name
    return(list(TIME.bin.summary=time.bin.temp$COV.bin.summary,Q.CI=Q.CI))
  } 
  
  if(is.null(N.timebin)) 
    N.timebin<-optK(orig.data,X.name,Y.name,ID.name,maxK,beta,lambda)
  if(is.null(main.title))
    main.title<-paste("N.timebin=",N.timebin,sep="")
  
  
  plot.data<-data.frame(orig.data,X=orig.data[,X.name],Y=orig.data[,Y.name])
  time.bin<-makeCOVbin(plot.data$X,N.timebin)   
  if(is.null(Y.min)) Y.min<-min(c(unlist(sim.data),plot.data$Y),na.rm=T)
  if(is.null(Y.max)) Y.max<-max(c(unlist(sim.data),plot.data$Y),na.rm=T)
  
  P.temp<-ggplot(plot.data,aes(x=X,y=Y))+ylim(Y.min,Y.max)+
    labs(x=X.name,y=Y.name,title=main.title)
  if(bin.grid){
    temp.tick<-unique(unlist(time.bin$COV.bin.summary[,4:5]))
    P.temp<-P.temp+geom_vline(xintercept=temp.tick,linetype=3,col="grey50")
  }
  
  if(opt.SIM.quantile.CI.area){
    
    temp.simCI<-findQuantileCI(sim.data,plot.data$X,
                               time.bin,q=q.list,alpha=alpha,N.sim)
    if(is.null(X.min)) X.min<-min(c(plot.data$X,temp.simCI$Q.CI[[1]][1,2]),
                                  na.rm=T)-range(plot.data$X)[2]/20
    if(is.null(X.max)) X.max<-max(c(plot.data$X,
                                    temp.simCI$Q.CI[[1]][nrow(
                                      temp.simCI$Q.CI[[1]]),3]),
                                  na.rm=T)+range(plot.data$X)[2]/20
    P.temp<-P.temp+xlim(X.min,X.max)
    
    if(!CIarea.box){
      NN<-nrow(temp.simCI$TIME.bin.summary)
      test.LU<-c(temp.simCI$TIME.bin.summary[1,4],
                 temp.simCI$TIME.bin.summary[,6],         
                 # test.LU<-c(temp.simCI$TIME.bin.summary[1,2],
                 #            temp.simCI$TIME.bin.summary[,2],
                 #apply(temp.simCI$TIME.bin.summary[,4:5],1,mean),
                 temp.simCI$TIME.bin.summary[NN,5])
      test.data.tot<-temp.simCI$Q.CI
      n.temp<-length(test.LU)
      X<-c(test.LU,test.LU[n.temp:1])
      test.data<-test.data.tot[[1]]
      NN<-nrow(test.data)
      Y<-c(test.data[1,2],test.data[,2],test.data[NN,2],
           test.data[NN,4],test.data[NN:1,4],test.data[1,4])
    } else{
      test.LU<-temp.simCI$TIME.bin.summary[,4:5]
      test.data.tot<-temp.simCI$Q.CI
      X.temp<-c(test.LU[,1],test.LU[nrow(test.LU),2])
      n.temp<-nrow(test.LU)
      X<-c(test.LU[1,1],rep(test.LU[2:n.temp,1],each=2),test.LU[n.temp,2])
      X<-c(X,X[length(X):1])
      test.data<-test.data.tot[[1]]
      Y<-c(rep(test.data[,2],each=2),rep(test.data[(n.temp:1),4],each=2))
    }
    if(Y0.adjust) Y[Y<0]<-0
    SIM.CIarea.1<-data.frame(X=X,Y=Y,ID=1)
    
    P.temp<-P.temp+geom_polygon(data= SIM.CIarea.1,
                                aes(x=X,y=Y,group=ID,fill=ID),
                                fill="lightblue",colour="lightblue")
    #fill="gray80",colour="gray80")
    test.data<-test.data.tot[[3]]
    if(!CIarea.box){
      NN<-nrow(test.data)
      Y<-c(test.data[1,2],test.data[,2],test.data[NN,2],
           test.data[NN,4],test.data[NN:1,4],test.data[1,4])
    } else{
      Y<-c(rep(test.data[,2],each=2),
           rep(test.data[(n.temp:1),4],each=2))
    }
    if(Y0.adjust) Y[Y<0]<-0      
    SIM.CIarea.3<-data.frame(X=X,Y=Y,ID=1)
    P.temp<-P.temp+geom_polygon(data= SIM.CIarea.3,
                                aes(x=X,y=Y,group=ID,fill=ID),
                                fill="lightblue",colour="lightblue")
    #fill="gray80",colour="gray80")
    test.data<-test.data.tot[[2]]
    if(!CIarea.box){
      NN<-nrow(test.data)      
      Y<-c(test.data[1,2],test.data[,2],test.data[NN,2],
           test.data[NN,4],test.data[NN:1,4],test.data[1,4])
    } else{
      Y<-c(rep(test.data[,2],each=2),
           rep(test.data[(n.temp:1),4],each=2))
    }
    if(Y0.adjust) Y[Y<0]<-0      
    SIM.CIarea.2<-data.frame(X=X,Y=Y,ID=1)
    P.temp<-P.temp+geom_polygon(data=SIM.CIarea.2,
                                aes(x=X,y=Y,group=ID,fill=ID),
                                fill="pink",colour="pink")  
    #fill="gray50",colour="gray50")  
  } 
  if(opt.DV.point){
    P.temp<-P.temp+geom_point(color="grey30",size=2,alpha=0.5)
    DV.point<-data.frame(X=orig.data[,X.name],Y=orig.data[,Y.name])
  } 
  
  if(opt.DV.quantile.line){
    temp<-findQuantile(plot.data$Y,plot.data$X,time.bin,q=q.list)
    if(Y0.adjust) 
      temp[,6:8]<-apply(temp[,6:8],2,function(x) {x[x<0]<-0;x})
    temp.group<-factor(c(temp[,7]>=test.data.tot[[1]][,2] & temp[,7]<=test.data.tot[[1]][,4],
                         temp[,8]>=test.data.tot[[2]][,2] & temp[,8]<=test.data.tot[[2]][,4],
                         temp[,9]>=test.data.tot[[3]][,2] & temp[,9]<=test.data.tot[[3]][,4]),
                       levels=c(TRUE,FALSE),labels=c("In","Out"))
    DV.quant<-data.frame(X=rep(apply(temp[,4:5],1,mean),length(q.list)),
                         G=factor(rep(paste("Q",round(q.list*100),
                                            "th",sep=""),each=nrow(temp))),
                         Y=unlist(temp[,-(1:(ncol(temp)-length(q.list)))]),
                         groupS=temp.group)
    
    P.temp<-P.temp+geom_line(data = DV.quant[DV.quant$G!="Q50th",],
                             aes(x=X,y=Y,group=G),
                             linetype=1,size=1,color="blue")+
      #linetype=1,size=1,color="black")+
      geom_line(data = DV.quant[DV.quant$G=="Q50th",],
                aes(x=X,y=Y,group=G),
                linetype=1,size=1,color="red")
    #linetype=1,size=1,color="black")
    if(orig.flag){
      P.temp<-P.temp+geom_point(data=DV.quant[DV.quant$G!="Q50th",],
                                aes(x=X,y=Y,group=G,shape=groupS),
                                size=4,color="blue")+
        #linetype=1,size=1,color="black")+
        geom_point(data=DV.quant[DV.quant$G=="Q50th",],
                   aes(x=X,y=Y,group=G,shape=groupS),
                   size=4,color="red")
      
    }
  } 
  
  if(opt.SIM.quantile.line){
    temp.sim.Q<-findQuantile(sim.data,plot.data$X,time.bin,q=q.list)
    if(Y0.adjust) 
      temp.sim.Q[,6:8]<-apply(temp.sim.Q[,6:8],2,function(x) {x[x<0]<-0;x})      
    SIM.quant<-data.frame(X=rep(temp.sim.Q$med.COV,length(q)),
                          G=factor(rep(paste("Q",round(q.list*100),
                                             "th",sep=""),each=nrow(temp.sim.Q))),
                          Y=unlist(temp.sim.Q[,-(1:(ncol(temp.sim.Q)-
                                                      length(q.list)))]))
    
    P.temp<-P.temp+geom_line(data = SIM.quant[SIM.quant$G!="Q50th",],
                             aes(x=X,y=Y,group=G),linetype=2,
                             size=1,color="blue")+
      geom_line(data = SIM.quant[SIM.quant$G=="Q50th",],
                aes(x=X,y=Y,group=G),linetype=2,
                size=1,color="red")
  }
  if(plot.flag){
    P.temp+theme_bw()+theme(panel.grid.major=element_line(colour = "white"))+
      theme(panel.grid.minor=element_line(colour="white"))
  } else{
    return(list(SIM.CIarea.1=SIM.CIarea.1,SIM.CIarea.2=SIM.CIarea.2, 
                SIM.CIarea.3=SIM.CIarea.3,DV.point=DV.point,DV.quant=DV.quant,
                SIM.quant=SIM.quant))
  }
}