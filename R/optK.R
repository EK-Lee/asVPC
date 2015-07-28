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


optK<-function(orig.data,X.name="TIME",Y.name="DV",ID.name="ID",maxK=NULL,beta=0.2,lambda=0.3){
  
   plot.data<-data.frame(X=orig.data[,X.name],Y=orig.data[,Y.name],ID=orig.data[,ID.name])
   if(is.null(maxK)) 
     maxK <- max(c(round(nrow(orig.data)/10),tapply(plot.data$X,plot.data$ID,length)))
   keep.tot<-NULL
   for(k in 2:maxK){
     time.bin<-makeCOVbin(plot.data$X,k)  
     J.temp<-sum((tapply(plot.data$X,time.bin$COV.bin,var)^beta)*
                   table(time.bin$COV.bin))
     temp.U<-log(J.temp)+lambda*beta*k
     keep.tot<-rbind(keep.tot,c(k,temp.U))
   }
   colnames(keep.tot)<-c("k","U")
   keep.tot<-data.frame(keep.tot)
   return(keep.tot$k[which(keep.tot$U==min(keep.tot$U))])
    
}

