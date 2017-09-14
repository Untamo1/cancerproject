library(grDevices)
library(xlsx)
library(JADE)
library(ggplot2)
library(reshape2)

source("https://raw.githubusercontent.com/Untamo1/fractalseparation/master/complexica_fun.R")
# setwd("C:/Users/Niko/Desktop/cancerproject/syopasaatiossa")
setwd("C:/Users/lietz/OneDrive/cancerproject/syopasaatiossa")

reset <- function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="Time", ylab="", axes=FALSE)
}




cervicalcancerprin <- function(D,U,UM,col,lab,str=1953,inter=3,nam ){
  hei <- 10
  wei <- 12
  
  oricent <- sweep(D,2,colMeans(D),FUN="-")
  
  n <- dim(D)[1]
  p <- dim(D)[2]
  
  
  meancurve <- ts(as.matrix(rowMeans(D)),start=str)
  
  pdf("cervicaloriginal.pdf",paper="special",height=16,width=23)
  par(mar = c(1.5, 1.5, 1.5, 1.5), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=(c(1,1)),cex=2.5)
  ts.plot(D,xlab="Time", ylab="Incidence",main="Cervical Cancer Incidence in Finland between 1953-2014",col=var,type="l",lwd=2)
  lines(meancurve,col="black",lwd=3)
  mtext("Incidence", side = 2, outer = TRUE, cex = 2.5, line = 2.2,col = "black")
  reset()
  # add mean to lab
  legend(x="bottom",legend = c(lab,"mean"), col=c(col,"black"), pch=16,horiz = TRUE,bty="n", cex = 2.5)
  dev.off()
  
  
  for(i in 1:inter){
    nimi <- paste(nam,"comp",i,".pdf",sep="")
    pdf(nimi,paper="special",height=hei,width=wei)
    par(mar = c(1.5, 1.5, 1.5, 1.5), oma = c(0.25, 0.25, 0.25, 0.25))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    par(mfrow=(c(1,1)),cex=2.5)
    man <- paste("tICS Component ",i,sep="")
    plot(ts(U[,i],start=str),xlab="Time", ylab="Incidence",main=man,lwd=5,type="l")
    reset()
    dev.off()
  }
  
  nimi <- paste(nam,"_uninteresting",sep="")
  pdf("cervicaluninteresting.pdf",paper="special",height=10,width=16)
  par(mar = c(1.5, 1.5, 1.5, 1.5), oma = c(0.25, 0.25, 0.25, 0.25))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=(c(2,4)),cex=1.5)
  for(i in (inter+1):p){
    man <- paste("tICS Component ",i,sep="")
    plot(ts(U[,i],start=str),xlab="Time", ylab="Incidence",main=man,lwd=2.5,type="l")
    
  }
  mtext("Incidence", side = 2, outer = TRUE, cex = 2.5, line = 2.2,col = "black")
  reset()
  dev.off()

  
  pdf("dimticsfinal.pdf",paper="special",height=24,width=20)
  par(mfrow = c(4, 3))
  par(cex = 1.5)
  par(mar = c(3, 1.5, 1.5, 1.5), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  for(i in 1:p){
    inv <- solve(UM)
    varisyopa = matrix(data=NA, nrow=n, ncol=(inter+1))
    tmp <- 0
    
    for(k in 1:(inter)){
      tmp <- tmp + inv[i,k] * U[,k]
      varisyopa[,k] = tmp
    }
    
    varisyopa[,(inter+1)] = U %*% inv[i,]
    varisyopa1 <- ts(as.matrix(varisyopa),start=str)
    tex <- paste("Age ",lab[i],sep="")
    ts.plot(varisyopa1,main=tex,col=c("blue","darkgreen","red","black"),lty=c(3,4,2,1),lwd=2.5)
    
    # print(varisyopa[,(inter+1)] - oricent[,i])
    
  }
  mtext("Incidence", side = 2, outer = TRUE, cex = 2.5, line = 2.2,col = "black")
  reset()
  legend(x=0.387,y=0.1375, legend = c(1:3,10),col=c("blue","darkgreen","red","black"),lwd=12, lty=c(3,4,2,1),horiz = TRUE, 
         cex = 3,title="Number of tICS Components",box.lty=1, box.lwd=2, box.col="black")
  dev.off()
  
  # if(p >6){
  #   pdf("dimtics2.pdf",paper="special",height=12,width=20)
  #   par(mfrow = c(2, 3))
  #   par(cex = 1.5)
  #   par(mar = c(3, 1.5, 1.5, 1.5), oma = c(4, 4, 0.5, 0.5))
  #   par(tcl = -0.25)
  #   par(mgp = c(2, 0.6, 0))
  #   for(i in 6:p){
  #     inv <- solve(UM)
  #     varisyopa = matrix(data=NA, nrow=n, ncol=(p-5))
  #     tmp <- 0
  #     
  #     for(k in 1:(inter)){
  #       tmp <- tmp + inv[i,k] * U[,k]
  #       varisyopa[,k] = tmp
  #     }
  #     
  #     varisyopa[,(inter+1)] = U %*% inv[i,]
  #     varisyopa1 <- ts(as.matrix(varisyopa),start=str)
  #     tex <- paste("Age ",lab[i],sep="")
  #     ts.plot(varisyopa1,main=tex,col=c("blue","darkgreen","red","black"),lty=c(3,4,2,1),lwd=2.5)
  #     
  #     # print(varisyopa[,(inter+1)] - oricent[,i])
  #   }
  #   mtext("Incidence", side = 2, outer = TRUE, cex = 2.5, line = 2.2,col = "black")
  #   reset()
  #   legend(x="bottom", legend = c(1:3,10),col=c("blue","darkgreen","red","black"),lwd=8, lty=c(3,4,2,1),horiz = TRUE,bty="n", cex = 2.5)
  #   dev.off()
  # }
  

  pdf("clusteri.pdf",paper="special",height=10,width=25)
  par(mfrow = c(1, 3))
  par(cex = 1.5)
  par(mar = c(3, 1.5, 1.5, 1.5), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  for(k in 1:inter){
    # nam <- paste("cluster",i,".pdf",sep="")
    # pdf(nam,paper="special",height=hei,width=wei)

    
    # clusterplot(cervicalAMUSE,"ccluster1",cervicalAMUSEmatrix,i,"")
    tex <- paste("tICS Component",k)
    varisyopa = matrix(data=NA, nrow=n, ncol=p) #Transformoitu matriisi
    
    
    for(i in 1:p){
      
      varisyopa[,i] = as.matrix(solve(UM)[i,k]*U[,k])
      
    }
    
    varisyopa1 <- ts(varisyopa,start=str)
    ts.plot(varisyopa1,col=var,main=tex,lwd=2.5)
  
    
  }
  
  reset()
  legend(x="bottom",legend = lab, col=col, pch=16,horiz = TRUE,bty="n", cex = 2.5)
  dev.off()
  
  # 
  # pdf("cluster_uninteresting.pdf",paper="special",height=16,width=20)
  # par(mfrow = c(2, 4))
  # par(cex = 1.5)
  # par(mar = c(3, 1.5, 1.5, 1.5), oma = c(4, 4, 0.5, 0.5))
  # par(tcl = -0.25)
  # par(mgp = c(2, 0.6, 0))
  # for(i in (inter+1):p){
  #   # clusterplot(cervicalAMUSE,"ccluster4",cervicalAMUSEmatrix,i,"")
  #   tex <- paste("tICS Component",i)
  #   vvarisyopa = matrix(data=NA, nrow=n, ncol=p) #Transformoitu matriisi
  #   
  #   
  #   for(k in 1:p){
  #     
  #     varisyopa[,k] = as.matrix(solve(UM)[p,k]*U[,i])
  #     
  #   }
  #   
  #   varisyopa1 <- ts(varisyopa,start=str)
  #   ts.plot(varisyopa1,col=var,main=tex,lwd=2.5)
  #  
  # }
  # mtext("Incidence", side = 2, outer = TRUE, cex = 2.5, line = 2.2,col = "black")
  # reset()
  # legend(x="bottom",legend = lab, col=col, pch=16,horiz = TRUE,bty="n", cex = 2.5)
  # dev.off()

  
  
}
dat.cer<-  read.xlsx("cervical2014.xlsx",1,header=TRUE,row.names=1)
# dat.cer <-  read.csv("cervical2014.csv",sep=",",header=TRUE,row.names=1)

cd <- cbind(dat.cer[,1:9], (dat.cer[,10]+ dat.cer[,11]))
lab <- c("0-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69",
         "70-74","75+")

colnames(cd) <- lab

v1 <- "deeppink"
v2 <- "magenta4"
v3 <- "blue4"
v4 <- "blue" 
v5 <- "mediumaquamarine"
v6 <- "green"
v7 <- "yellow"
v8 <- "orange"
v9 <- "red"
v10 <- "azure4"
var <- c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)

matrixori <- ts(as.matrix(cd),start=1953)
# oricent <- sweep(matrixori,2,colMeans(matrixori),FUN="-")

Nestimate <- NAMUSE(matrixori,1)

Dat <- Nestimate$Data

Dat[,1] <- -Dat[,1]

Gam <- Nestimate$Gamma

Gam[1,] <- -Gam[1,] 

cervicalcancerprin(matrixori,Dat,Gam,var,lab,str=1953,inter=3,"c" )
