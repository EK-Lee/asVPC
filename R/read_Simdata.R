#' read the simulated data file generated from NONMEM
#'
#' @title Read NONMEM simulation data
#' @param sim.file.name file name of simulation, generated from NONMEM with 'NOAPPEND ONEHEAD' options in TABLE statement
#' @param data.n number of observations in the original data 
#' @param sim.n number of simulation
#' @param name.DV name of dependent variable in simulated data file
#' @return data.n * sim.n matrix with simulated data
#' @export

read_Simdata<-function(sim.file.name,data.n,sim.n,name.DV){ 
   data.n<-data.n+2
   sim.data<-NULL
   for(i in 1:sim.n){
      sel.id<-(i-1)*data.n+(3:data.n)
      temp.data<-read.table(sim.file.name,skip=(i-1)*data.n+1,
                            header=TRUE,nrows=data.n-2)
      sim.data<-cbind(sim.data,temp.data[,name.DV])
   }
   return(sim.data)
}
