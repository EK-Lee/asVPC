#' calculate percentiles of original data using bin-related weight 
#' percentiles of simulated data with corresponding confidence interval 
#' and draw the asVPC plot
#' test
#' 
#' @import ggplot2  plyr
#' @title asVPC with bin-related weight
#' @usage asVPC.binW(orig.data,sim.data,N.timebin,n.sim,n.hist,
#'                   q.list=c(0.05,0.5,0.95),conf.level=0.95,
#'                   X.name="TIME",Y.name="DV",opt.DV.point=FALSE,
#'                   weight.flag=FALSE,Y.min=NULL,Y.max=NULL,
#'                   only.med=FALSE,plot.flag=TRUE)
#' @param orig.data the original data for model fitting
#' @param sim.data the simulated data 
#' @param N.timebin the number of bin of X axis
#' @param n.sim the number of simulation in the simulated data
#' @param n.hist the number of shifted 
#' @param q.list numeric vector of probabilities with values in [0,1]
#' @param conf.level confidence level of the interval
#' @param X.name the name of X variable in the original scatter plot
#' @param Y.name the name of Y variable in the original scatter plot
#' @param opt.DV.point option to put data point in the plot
#' @param weight.flag option to use weight in average shifted calculation 
#' @param Y.min minimum of Y axis in the plot
#' @param Y.max maximum of Y axis in the plot
#' @param only.med use only median if TRUE
#' @param plot.flag draw plot if TRUE; generate data for drawing plot if FALSE
#' @return plot or the values to draw plot
#' @export
#' @seealso \code{\link{asVPC.distanceW}}
#' @examples
#' data(origdata)
#' data(simdata)
#' asVPC.binW.CI(origdata,simdata,N.timebin=10, n.sim=100,n.hist=3)

asVPC.binW.CI<-function(orig.data,sim.data,n.sim,
                        N.timebin=NULL,n.hist=NULL,
                        q.list=c(0.05,0.5,0.95),
                        conf.level=0.95,
                        X.name="TIME",Y.name="DV",ID.name="ID",
                        main.title=NULL,
                        opt.DV.point=FALSE,
                        weight.flag=FALSE,
                        Y.min=NULL,
                        Y.max=NULL,
                        X.min=NULL,X.max=NULL,
                        only.med=FALSE,
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
   
   if(is.null(N.timebin)) 
     N.timebin<-optK(orig.data,X.name,Y.name,ID.name,maxK,beta,lambda)
   if(is.null(n.hist))
     n.hist<-ifelse(N.timebin==length(table(orig.data[,X.name])),
                    2,round(N.timebin/2))
   bintot.N<-N.timebin*n.hist
   if(is.null(main.title))
     main.title<-paste("BIN: N.timebin=",N.timebin,", n.hist=",n.hist,sep="")
   
   time.bin<-makeCOVbin(orig.data[,X.name],N.covbin=bintot.N)
   alpha<-1-conf.level   
   Q.CI<-vector("list",3)
   orig.Q<-NULL
   bintot.N<-nrow(time.bin$COV.bin.summary)
   for(i in 1:bintot.N){
      if(i<n.hist){
         sel.id<-which(as.numeric(time.bin$COV.bin)<=i+n.hist-1)
         sel.id1<-which(as.numeric(time.bin$COV.bin)==i)
         mid.point<-median(orig.data[sel.id1,X.name])
         low.point<-time.bin$COV.bin.summary$lower.COV[i]
         upper.point<-time.bin$COV.bin.summary$upper.COV[i]
      } else if(i>(bintot.N-n.hist+1)){
         sel.id<-which(as.numeric(time.bin$COV.bin)>=i-(n.hist-1))
         sel.id1<-which(as.numeric(time.bin$COV.bin)==i)
         mid.point<-median(orig.data[sel.id1,X.name])
         low.point<-time.bin$COV.bin.summary$lower.COV[i]
         upper.point<-time.bin$COV.bin.summary$upper.COV[i]
      } else{
         sel.id <-which(as.numeric(time.bin$COV.bin)>i-n.hist & 
                        as.numeric(time.bin$COV.bin)<i+n.hist)
         sel.id1<-which(as.numeric(time.bin$COV.bin)==i)
         mid.point<-median(orig.data[sel.id1,X.name])
         low.point<-time.bin$COV.bin.summary$lower.COV[i]
         upper.point<-time.bin$COV.bin.summary$upper.COV[i]
      }   

      A<-as.numeric(time.bin$COV.bin[sel.id])
      temp<-abs(A-as.numeric(time.bin$COV.bin[sel.id1[1]]))
      temp.weight<-(max(temp)+1)-temp
      temp.weight<-temp.weight/max(temp.weight)
      
      if(weight.flag){
         temp.quantile<-t(apply(sim.data[sel.id,],2,function(x) 
           Hmisc::wtd.quantile(x,weight=temp.weight,
                                               prob=q.list,na.rm=TRUE)))
         temp.orig.q<-Hmisc::wtd.quantile(orig.data[,Y.name][sel.id],
                                   weight=temp.weight,prob=q.list,na.rm=TRUE)
      } else {
         temp.quantile<-t(apply(sim.data[sel.id,],2,function(x) 
                                         quantile(x,prob=q.list,na.rm=TRUE)))
         temp.orig.q<-quantile(orig.data[,Y.name][sel.id],
                               prob=q.list,na.rm=TRUE)
      } 
      orig.Q<-rbind(orig.Q,c(mid.point,temp.orig.q))
      temp<-t(apply(temp.quantile,2,
                    function(x) quantile(x,prob=c(alpha/2,0.5,1-alpha/2),
                                         na.rm=TRUE)))
      for(j in 1:length(q.list))
         Q.CI[[j]]<-rbind(Q.CI[[j]],c(mid.point,low.point,upper.point,temp[j,]))
   } 

   keep.name<-NULL
   for(j in 1:length(q.list)){
      keep.name<-c(keep.name,paste("Q",round(q.list[j]*100),"th",sep=""))
      colnames(Q.CI[[j]])   <-c("mid","Lower","upper",colnames(Q.CI[[j]])[4:6])
   }    
   names(Q.CI)<-keep.name
   colnames(orig.Q)<-c("mid","Y1","Y2","Y3")
   orig.Q<-data.frame(orig.Q)
 
   plot.data<-data.frame(orig.data,X=orig.data[,X.name],Y=orig.data[,Y.name])
   if(is.null(Y.min)) Y.min<-min(c(plot.data$Y,Q.CI[[1]][,4]),na.rm=T)
   if(is.null(Y.max)) Y.max<-max(c(plot.data$Y,Q.CI[[length(Q.CI)]][,6]),na.rm=T)
   if(is.null(X.min)) X.min<-min(c(plot.data$X,Q.CI[[1]][1,2]),
                                 na.rm=T)-range(plot.data$X)[2]/20
   if(is.null(X.max)) X.max<-max(c(plot.data$X,Q.CI[[1]][bintot.N,3]),
                                 na.rm=T)+range(plot.data$X)[2]/20
   
   P.temp<-ggplot(plot.data,aes(x=X,y=Y))+ylim(Y.min,Y.max)+xlim(X.min,X.max)+
               labs(x=X.name,y=Y.name,title=main.title)+
               theme_bw()+theme(panel.grid.major=element_line(colour="white"))+
               theme(panel.grid.minor=element_line(colour="white"))
   
   if(bin.grid){
     temp.tick<-unique(unlist(time.bin$COV.bin.summary[,4:5]))
     P.temp<-P.temp+geom_vline(xintercept=temp.tick,linetype=3,col="grey50")
   }   
   
   NN<-nrow(Q.CI[[1]])   
   test.LU<-c(Q.CI[[1]][1,2],Q.CI[[1]][,1],Q.CI[[1]][NN,3],
              Q.CI[[1]][NN,3],Q.CI[[1]][NN:1,1],Q.CI[[1]][1,2])
   test.data.tot<-Q.CI
   n.temp<-nrow(test.LU)
   if(CIarea.box){
     test.LU<-Q.CI[[1]][,2:3]
     test.data.tot<-Q.CI
     X.temp<-c(test.LU[,1],test.LU[nrow(test.LU),2])
     n.temp<-nrow(test.LU)
     X<-c(test.LU[1,1],rep(test.LU[2:n.temp,1],each=2),test.LU[n.temp,2])
     X<-c(X,X[length(X):1])
   } else{
      X<-c(Q.CI[[1]][1,2],Q.CI[[1]][,1],Q.CI[[1]][NN,3],
           Q.CI[[1]][NN,3],Q.CI[[1]][NN:1,1],Q.CI[[1]][1,2]) 
   }
   if(!only.med){ 
      test.data<-test.data.tot[[1]]
      NN<-nrow(test.data)
      if(!CIarea.box){
        Y<-c(test.data[1,4],test.data[,4],test.data[NN,4],
             test.data[NN,6],test.data[NN:1,6],test.data[1,6])
      } else{ 
        Y<-c(rep(test.data[,4],each=2),rep(test.data[(n.temp:1),6],each=2))
      }
      SIM.CIarea.1<-data.frame(X=X,Y=Y,ID=1)
      P.temp<-P.temp+geom_polygon(data= SIM.CIarea.1,
                                  aes(x=X,y=Y,group=ID,fill=ID),
                                  fill="lightblue",colour="lightblue")
                                 #fill="gray80",colour="gray80")
      test.data<-test.data.tot[[3]]
      NN<-nrow(test.data)  
      if(!CIarea.box){
        Y<-c(test.data[1,4],test.data[,4],test.data[NN,4],
             test.data[NN,6],test.data[NN:1,6],test.data[1,6])
      } else{
        Y<-c(rep(test.data[,4],each=2),rep(test.data[(NN:1),6],each=2))
      }
      SIM.CIarea.3<-data.frame(X=X,Y=Y,ID=1)
      P.temp<-P.temp+geom_polygon(data= SIM.CIarea.3,
                                  aes(x=X,y=Y,group=ID,fill=ID),
                                  fill="lightblue",colour="lightblue")
                                  #fill="gray80",colour="gray80")
   }
   test.data<-test.data.tot[[2]]
   NN<-nrow(test.data)   
   if(!CIarea.box){
     Y<-c(test.data[1,4],test.data[,4],test.data[NN,4],
          test.data[NN,6],test.data[NN:1,6],test.data[1,6])    
   } else{
     Y<-c(rep(test.data[,4],each=2),rep(test.data[(n.temp:1),6],each=2))
   }
   SIM.CIarea.2<-data.frame(X=X,Y=Y,ID=1)
   P.temp<-P.temp+geom_polygon(data= SIM.CIarea.2,
                               aes(x=X,y=Y,group=ID,fill=ID),
                               fill="pink",colour="pink")  
                              #fill="gray50",colour="gray50")     
   if(opt.DV.point==TRUE){
      P.temp<-P.temp+geom_point(color="grey30",size=2,alpha=0.5) 
      DV.point <- data.frame(X=orig.data[,X.name],Y=orig.data[,Y.name])
   }
   DV.quant<-data.frame(X=rep(orig.Q$mid,length(q.list)),
                        G=factor(rep(paste("Q",round(q.list*100),"th",sep=""),
                                     each=nrow(orig.Q))),
                        Y=unlist(orig.Q[,-1]))
#   
      P.temp<-P.temp+geom_line(data=DV.quant[DV.quant$G!="Q50th",],
                               aes(x=X,y=Y,group=G),
                               linetype=1,size=0.8,color="blue")+
                               #linetype=1,size=1,color="black")+
                     geom_line(data=DV.quant[DV.quant$G=="Q50th",],
                               aes(x=X,y=Y,group=G),
                               linetype=1,size=0.8,color="red")
                              #linetype=1,size=1,color="black")
     if(orig.flag){
        P.temp<-P.temp+geom_point(data=DV.quant[DV.quant$G!="Q50th",],
                                 aes(x=X,y=Y,group=G),
                                 size=4,color="blue")+
       #linetype=1,size=1,color="black")+
          geom_point(data=DV.quant[DV.quant$G=="Q50th",],
                    aes(x=X,y=Y,group=G),
                    size=4,color="red")
     #linetype=1,size=1,color="black")     
   }
   
 
   colnames(orig.Q)<-c("X.mid",paste("Q",round(q.list*100),"th",sep="")) 
   if(plot.flag){
      P.temp
   } else{
      return(list(SIM.CIarea.1=SIM.CIarea.1,SIM.CIarea.2=SIM.CIarea.2, 
                  SIM.CIarea.3=SIM.CIarea.3,DV.point=DV.point,
                  DV.quant=DV.quant,SIM.quant=SIM.quant))
   }
}