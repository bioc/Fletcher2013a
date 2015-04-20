###############################################################################
### Wrapper functions to extract and plot retuls from intermediate objects      		    
###############################################################################
#Internal function for Fletcher2013a: a simple function to extract results
extractResults<-function(myContrasts,contrastsFit,resultsFit){
  #---report deg list
  degAll<-list()
  labs<-colnames(myContrasts)
  for(i in 1:length(colnames(myContrasts))){
    nm<-names(resultsFit[resultsFit[,i]!=0,i])
    if(is.null(nm))nm<-NA
    degAll[[labs[i]]]<-nm
  }
  degPos<-list()
  for(i in 1:length(colnames(myContrasts))){
    nm<-names(resultsFit[resultsFit[,i]==1,i])
    if(is.null(nm))nm<-NA
    degPos[[labs[i]]]<-nm
  }
  degNeg<-list()
  for(i in 1:length(colnames(myContrasts))){
    nm<-names(resultsFit[resultsFit[,i]==-1,i])
    if(is.null(nm))nm<-NA
    degNeg[[labs[i]]]<-nm
  }
  deglist<-list(all=degAll,pos=degPos, neg=degNeg)
  #report differential expression results
  resdeg<-list()
  resdeg$PROBEID  <- contrastsFit$ids$PROBEID
  resdeg$SYMBOL   <- contrastsFit$ids$SYMBOL
  resdeg$ENTREZ   <- contrastsFit$ids$ENTREZ
  resdeg$coef     <- contrastsFit$coefficients
  resdeg$p.value  <- contrastsFit$p.value
  resdeg$degenes  <- unclass(resultsFit)
  resdeg <- data.frame(resdeg, check.names=FALSE, stringsAsFactors=FALSE)
  return(list(deg=resdeg,deglist=deglist))
}
#-------------------------------------------------------------------------------
#Internal function for Fletcher2013a: a simple function plot a summary
plotExpDeg<-function(results){
  opt<-c('info.Exp1limma','info.Exp2limma','info.Exp3limma')%in%names(results)
  elab<-which(opt)
  #------------------------------------
  #plot summary from Exp1
  plot.deg.exp1<-function(results){
    deglist<-results$deglist
    legs<-c("E2 vs. VEH","E2+FGF10 vs. E2","E2+FGF10+PD vs E2","DT VEH")
    labs<-c("t6.E2-t6.UT","t24.E2-t24.UT","t6.E2FGF10-t6.E2","t24.E2FGF10-t24.E2",
            "t6.E2FGF10PD-t6.E2","t24.E2FGF10PD-t24.E2",
            "t6.UT-t0.UT","t24.UT-t6.UT")
    plotendsum<-function(xx){
      maior<-1200
      pdf('fig1a.pdf',5,2)
      layout(matrix(c(1,3,5,7,2,4,6,8), 2, 4, byrow=TRUE))
      cols<-c("black","grey")
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[1]]][1]
      x[2]<-xx[[labs[2]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2), axisnames=FALSE)
      mtext("EXP1", col="black", side = 3, line = -1.0, outer = T,cex=0.7,adj = 0,font=2)
      mtext(legs[1],side=3,line=0,cex=0.5,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[1]]][2]
      x[2]<-xx[[labs[2]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[3]]][1]
      x[2]<-xx[[labs[4]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[2],side=3,line=0,cex=0.5,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[3]]][2]
      x[2]<-xx[[labs[4]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))  
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[5]]][1]
      x[2]<-xx[[labs[6]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[3],side=3,line=0,cex=0.5,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[5]]][2]
      x[2]<-xx[[labs[6]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))  
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[7]]][1]
      x[2]<-xx[[labs[8]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[4],side=3,line=0,cex=0.5,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[7]]][2]
      x[2]<-xx[[labs[8]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))
      legend("bottomleft",legend=c("induced","inhibited"),bty="n",cex=0.7,fill=cols)
      dev.off()
    }
    pos<-unlist(lapply(deglist$pos,length))
    neg<-unlist(lapply(deglist$neg,length))
    mc<-list()
    for(i in names(deglist$all)){
      mc[[i]]<-c(pos[[i]],neg[[i]])
    }
    plotendsum(mc)
  }
  
  #------------------------------------
  #plot summary from Exp2
  plot.deg.exp2<-function(results){
    legs<-c("E2 vs. VEH","E2+AP20187 vs. E2","E2+AP20187+PD vs. E2", "E2+FGF10 vs. E2", "DT VEH")
    labs<-c("t6.E2-t6.UT","t24.E2-t24.UT",
            "t6.E2.AP20187-t6.E2","t24.E2.AP20187-t24.E2",
            "t6.E2.AP20187.PD-t6.E2","t24.E2.AP20187.PD-t24.E2",
            "t6.E2.FGF10-t6.E2","t24.E2.FGF10-t24.E2",
            "t6.UT-t0.UT","t24.UT-t6.UT")
    deglist<-results$deglist
    plotendsum<-function(xx){
      maior<-max(as.data.frame(lapply(mc, sum)))
      maior<-4400
      pdf('fig2a.pdf',7,2.2)
      layout(matrix(c(1,3,5,7,9,2,4,6,8,10), 2, 5, byrow=TRUE))
      cols<-c("black","grey")
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[1]]][1]
      x[2]<-xx[[labs[2]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2), axisnames=FALSE)
      mtext("EXP2", col="black", side = 3, line = -1.0, outer = T,cex=0.7,adj = 0,font=2)
      mtext(legs[1],side=3,line=0,cex=0.6,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[1]]][2]
      x[2]<-xx[[labs[2]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[3]]][1]
      x[2]<-xx[[labs[4]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[2],side=3,line=0,cex=0.6,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[3]]][2]
      x[2]<-xx[[labs[4]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[5]]][1]
      x[2]<-xx[[labs[6]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[3],side=3,line=0,cex=0.6,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[5]]][2]
      x[2]<-xx[[labs[6]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))  
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[7]]][1]
      x[2]<-xx[[labs[8]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[4],side=3,line=0,cex=0.6,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[7]]][2]
      x[2]<-xx[[labs[8]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))  
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[9]]][1]
      x[2]<-xx[[labs[10]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0, maior), cex.names=1.0, cex.axis=0.9, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[5],side=3,line=0,cex=0.6,adj=0)
      x<-cbind(t6=NA, t24=NA)
      x[1]<-xx[[labs[9]]][2]
      x[2]<-xx[[labs[10]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.9, yaxp=c(-maior,0,2))
      legend("bottomleft",legend=c("induced","inhibited"),bty="n",cex=0.9,fill=cols)
      dev.off()
    }
    pos<-unlist(lapply(deglist$pos,length))
    neg<-unlist(lapply(deglist$neg,length))
    mc<-list()
    for(i in names(deglist$all)){
      mc[[i]]<-c(pos[[i]],neg[[i]])
    }
    plotendsum(mc)
  }
  #------------------------------------
  #plot summary from Exp3
  plot.deg.exp3<-function(results){
    deglist<-results$deglist
    plotendsum<-function(expname,xx,labs,legs=1){
      if(legs==1){
        legs<-c("E2 vs. VEH","E2+FGF10 vs. E2", "DT VEH")
      } else {
        legs<-c("E2 vs. E2","E2+FGF10 vs. E2+FGF10", "VEH vs. VEH")
      }
      maior<-2000
      cols<-c("black","grey")
      x<-cbind(t3=NA, t6=NA, t12=NA, t24=NA)
      x[1]<-xx[[labs[1]]][1]
      x[2]<-xx[[labs[2]]][1]
      x[3]<-xx[[labs[3]]][1]
      x[4]<-xx[[labs[4]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0,maior), cex.names=1.0, cex.axis=0.85, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext("EXP3", col="black", side = 3, line = -1.0, outer = T,cex=0.7,adj = 0,font=2)
      mtext(legs[1],side=3,line=0,cex=0.6,adj=0)
      x<-cbind(t3=NA, t6=NA, t12=NA, t24=NA)
      x[1]<-xx[[labs[1]]][2]
      x[2]<-xx[[labs[2]]][2]
      x[3]<-xx[[labs[3]]][2]
      x[4]<-xx[[labs[4]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.85, yaxp=c(-maior,0,2))
      x<-cbind(t3=NA, t6=NA, t12=NA, t24=NA)
      x[1]<-xx[[labs[5]]][1]
      x[2]<-xx[[labs[6]]][1]
      x[3]<-xx[[labs[7]]][1]
      x[4]<-xx[[labs[8]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0,maior), cex.names=1.0, cex.axis=0.85, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[2],side=3,line=0,cex=0.6,adj=0)
      x<-cbind(t3=NA, t6=NA, t12=NA, t24=NA)
      x[1]<-xx[[labs[5]]][2]
      x[2]<-xx[[labs[6]]][2]
      x[3]<-xx[[labs[7]]][2]
      x[4]<-xx[[labs[8]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.85, yaxp=c(-maior,0,2))
      x<-cbind(t3=NA, t6=NA, t12=NA, t24=NA)
      x[1]<-xx[[labs[9]]][1]
      x[2]<-xx[[labs[10]]][1]
      x[3]<-xx[[labs[11]]][1]
      x[4]<-xx[[labs[12]]][1]
      par(mar=c(0.2, 2.1, 2.1, 2.1))
      barplot(x,col=cols[1], horiz=FALSE, ylim=c(0,maior), cex.names=1.0, cex.axis=0.85, yaxp=c(0,maior,2),axisnames=FALSE)
      mtext(legs[3],side=3,line=0,cex=0.6,adj=0)
      mtext(expname, col="black", side = 4, line = 0.5, outer = F,cex=0.7,adj = 1,font=2)
      x<-cbind(t3=NA, t6=NA, t12=NA, t24=NA)
      x[1]<-xx[[labs[9]]][2]
      x[2]<-xx[[labs[10]]][2]
      x[3]<-xx[[labs[11]]][2]
      x[4]<-xx[[labs[12]]][2]
      par(mar=c(2.1, 2.1, 0.2, 2.1))
      barplot(-x,col=cols[2], horiz=FALSE, ylim=c(-maior,0), cex.names=1.0, cex.axis=0.85, yaxp=c(-maior,0,2))
      legend("bottomleft",legend=c("induced","inhibited"),bty="n",cex=0.9,fill=cols)
    }
    pos<-unlist(lapply(deglist$pos,length))
    neg<-unlist(lapply(deglist$neg,length))
    mc<-list()
    for(i in names(deglist$all)){
      mc[[i]]<-c(pos[[i]],neg[[i]])
    }
    pdf('fig3a.pdf',width=5,height=6.2)
    tp<-c(1,3,5,2,4,6)
    layout(matrix(c(tp,tp+6,tp+12), 6, 3, byrow=TRUE))
    labs<-names(mc)
    plotendsum("MinusTet",mc,labs[c(9,11,13,15,10,12,14,16,1,2,3,4)])
    plotendsum("PlusTet",mc,labs[c(17,19,21,23,18,20,22,24,5,6,7,8)])
    plotendsum("PlusTet vs. MinusTet", mc,labs[c(26,29,32,35,27,30,33,36,25,28,31,34)],legs=2)
    dev.off()
  }
  if(elab==1){
    plot.deg.exp1(results)
  } else if(elab==2){
    plot.deg.exp2(results)
  } else if(elab==3){
    plot.deg.exp3(results)
  }
}
#-------------------------------------------------------------------------------
#Internal function for Fletcher2013a: a simple function to plot the PCA analysis
plotExpPca<-function(rpca,targets,deg,elab){
  plot.pca.exp1<-function(rpca,targets,deg){
    ff=c(1,16,2,3)
    pdf(file='fig1b.pdf',width=6,height=6)
    plot(
      rpca$x[,1],rpca$x[,2],xlab="PCA 1",ylab="PCA 2",type="p",
      pch=ff[as.factor(targets$Treatment)],cex=1.5,
      col=c('red','green','blue')[as.factor(targets$Time)],
      xlim=c(min(rpca$x[,1])*1.1,max(rpca$x[,1])*1.1),
      ylim=c(min(rpca$x[,2]*1.1),max(rpca$x[,2])*1.1)
    )
    legend("bottomleft",c("E2","E2.FGF10","E2.FGF10.PD","VEH"),pch=ff,bty="n")
    legend("bottomright",c("red=0h","green=6h","blue=24h"),col=c('red','green','blue'),bty="n")
    mtext(paste("DE genes (n=",length(deg),")",sep=""),side=3,adj=0,cex=1.2)
    dev.off()
  }
  plot.pca.exp2<-function(rpca,targets,deg){
    ff=c(1,16,2,4,3)
    pdf(file='fig2b.pdf',width=6, height=6)
    plot(rpca$x[,1], rpca$x[,2], xlab="PCA 1", ylab="PCA 2",type="p",
         pch=ff[as.factor(targets$Treatment)],
         col=c('red','green','blue')[as.factor(targets$Time)],
         cex=1.5, 
         xlim=c(min(rpca$x[,1])*1.1, max(rpca$x[,1])*1.1),
         ylim=c(min(rpca$x[,2]*1.1), max(rpca$x[,2])*1.1) )
    text(rpca$x[!is.na(targets$TecRep),1]-5,rpca$x[!is.na(targets$TecRep),2]-2, 
         c(targets$TecRep[!is.na(targets$TecRep)]),cex=0.9)
    legend("topleft",c("E2","E2.AP20187","E2.AP20187.PD","E2.FGF10","VEH"),pch=ff,bty="n", cex=1.0)
    legend("topright",c("red=0h","green=6h","blue=24h"),col=c('red','green','blue'),bty="n", cex=1.0)
    leg<-paste("DE genes (n=",length(deg),")",sep="")
    mtext(leg,side = 3,adj = 0, cex = 1.2)
    dev.off()
  }
  plot.pca.exp3<-function(rpca,targets,deg){
    ff=c(1,19,10,0,15,12)
    pdf(file='fig3b.pdf',width=6, height=6)
    plot(rpca$x[,1], rpca$x[,2], xlab="PCA 1", ylab="PCA 2",type="p",
         pch=ff[as.factor(targets$Treatment)],
         col=c('red','orange','green','cyan','blue')[as.factor(targets$Time)],
         cex=1.5, 
         xlim=c(min(rpca$x[,1])*1.1, max(rpca$x[,1])*1.1),
         ylim=c(min(rpca$x[,2]*1.1), max(rpca$x[,2])*1.1) )
    legend("bottomleft",levels(as.factor(targets$Treatment)),pch=ff,bty="n", cex=1.0)
    legend("bottomright",c("red=0h","orange=3h","green=6h","cyan=12h","blue=24h"),
           col=c('red','orange','green','cyan','blue'),bty="n", cex=1.0)
    leg<-paste("DE genes (n=",length(deg),")",sep="")
    mtext(leg,side = 3,adj = 0, cex = 1.2)
    dev.off()
  }
  if(elab==1){
    plot.pca.exp1(rpca,targets,deg)
  } else if(elab==2){
    plot.pca.exp2(rpca,targets,deg)
  } else if(elab==3){
    plot.pca.exp3(rpca,targets,deg)
  }
}
##-----------------------------------------------------------------------------
##Internal function for Fletcher2013a: plot overlap among signatures
plotOverlap<-function(){
  deExp1 <- Fletcher2013pipeline.deg(what="Exp1")
  deExp2 <- Fletcher2013pipeline.deg(what="Exp2")
  deExp3 <- Fletcher2013pipeline.deg(what="Exp3")
  vv <- list(Exp1=deExp1$E2FGF10, Exp2=deExp2$E2AP20187, Exp3=deExp3$TetE2FGF10)
  pdf(file="fig4a.pdf",5,5)
  grid.draw(venn.diagram(vv,cat.col=c("red","blue","green"),filename=NULL))
  dev.off()
}
##-----------------------------------------------------------------------------
##Internal function for Fletcher2013a: plot IL8 gene expression validation
plotFollowup<-function(){
  #---load gx data toguether with qPRC and ELISA results
  data(IL8vld)
  data(Exp1)
  IL8vld<-get("IL8vld")
  Exp1<-get("Exp1")
  #---extract data from eSets
  Exp1geneids<-fData(Exp1)
  Exp1targets<-pData(Exp1)
  Exp1gexp<-exprs(Exp1)
  #---get IL8 probe ids
  idx<-Exp1geneids$SYMBOL=="IL8" 
  idx[is.na(idx)]<-FALSE
  check_genes<-Exp1geneids[idx,]
  #---start plot
  pdf(file="fig4b.pdf",5,5)
  layout(matrix(c(1,3,2,4), 2, 2, byrow=TRUE))
  palette(c("grey","black","blue","red","green"))
  #---get and plot gx of the 1st IL8 probe
  checkg<-check_genes[1,1]
  tempVector<-is.element(rownames(Exp1gexp), checkg) 
  tempExpdata <- Exp1gexp[tempVector,]
  PlotGroups(
    tempExpdata, time = Exp1targets$Time, repvect = Exp1targets$Replicates, groups = Exp1targets[,4:7],
    main = "", sub=checkg, legpos="topright", inset=c(0.1,0), legsz=0.6, ylab="IL8 gene expression",
    ylim=c(7.8,13),lwd=1.4
  )
  #---get and plot gx of the 2nd IL8 probe
  checkg<-check_genes[2,1]
  tempVector<-is.element(rownames(Exp1gexp), checkg) 
  tempExpdata <- Exp1gexp[tempVector,]
  PlotGroups(
    tempExpdata, time = Exp1targets$Time, repvect = Exp1targets$Replicates, groups = Exp1targets[,4:7],
    main = "", sub=checkg, legpos="topright", inset=c(0.1,0), legsz=0.6,ylab="IL8 gene expression",
    ylim=c(7.8,13),lwd=1.4
  )
  #---plot qPCR results
  plot(
    x = IL8vld$qPCR$Time, y = IL8vld$qPCR$fold.increase, pch = 21, 
    ylab = "IL8 mRNA fold increase", xlab="time", type="n",
    xaxt = "n", main="", sub = NULL, ylim = c(0,160), xlim = NULL, cex.axis = 0.9,
    cex = 0.9, col = c("blue","red"), cex.sub=0.6, mgp = c(2.0, 0.6, 0), bty='n'
  )   
  axis(1,at=IL8vld$qPCR$Time, labels=IL8vld$qPCR$Time, cex.axis = 0.9, mgp=c(2.0, 0.6, 0))
  lines(x=IL8vld$qPCR$Time[1:4], y=IL8vld$qPCR$fold.increase[1:4],col="blue",lwd=1.4)
  lines(x=IL8vld$qPCR$Time[5:8], y=IL8vld$qPCR$fold.increase[5:8],col="red",lwd=1.4)
  plotCI(
    x=IL8vld$qPCR$Time, y=IL8vld$qPCR$fold.increase, uiw=pmax(12,IL8vld$qPCR$SE), 
    col=c("blue","red")[as.numeric(as.factor(IL8vld$qPCR$Treatment))],gap=0.6,
    add=TRUE, pch=21, cex=0.9, pt.bg="white"
  )
  legend("topright", legend=levels(as.factor(IL8vld$qPCR$Treatment)), 
         text.col = c("blue","red"), col=c("blue","red"), cex=0.6, 
         lty=1, yjust = 0, bty="n")
  #---plot ELISA results
  bp<-t(IL8vld$ELISA)
  barplot2(
    bp[1:3,], beside = TRUE,col = c("blue", "red","green"), ylim = c(0, 150), ylab="IL8 (ng/mL per 500k cells)",
    main = "", xlab="time",font.main = 4, sub = "", col.sub = "gray20",axes=FALSE,axisnames=FALSE,
    cex.names = 0.9, plot.ci = TRUE, ci.l = bp[1:3,]-bp[4:6,], ci.u = bp[1:3,]+bp[4:6,],
    plot.grid = TRUE,bty="n", grid.lwd=0.3,cex=1,mgp = c(2.0, 0.6, 0),ci.lwd=1.4
  )
  axis(2,at=c(20,60,100,140),labels=c(20,60,100,140),cex.axis = 0.9, mgp=c(2.0, 0.6, 0))
  axis(1,at=c(2.5,6.5,10.5), labels=c(24,48,72), cex.axis = 0.9, mgp=c(2.0, 0.6, 0))
  dev.off()
}
##-----------------------------------------------------------------------------
## Adapted "PlotGroups" function for Fletcher2013a
PlotGroups<-function (data, edesign = NULL, time = edesign[, 1], 
                      groups = edesign[,c(3:ncol(edesign))], repvect = edesign[, 2], 
                      dis = NULL, step.method = "backward", min.obs = 2, alfa = 0.05, 
                      nvar.correction = FALSE, summary.mode = "median", show.lines = TRUE, 
                      groups.vector = NULL, xlab = "time", cex.xaxis = 0.9, ylim = NULL, 
                      main = NULL, cexlab = 0.9, legend = TRUE, sub = NULL, legpos="topright", inset=0,
                      legsz=0.5, startcol=NULL, plotdata=FALSE, ylab=NULL, plotSEM=FALSE,lwd=1.0)
{
  #library(gplots)
  if (!is.vector(data)) {
    if (summary.mode == "representative") {
      distances <- apply(as.matrix(dist(data, diag = TRUE, 
                                        upper = TRUE)), 1, sum)
      representative <- names(distances)[distances == min(distances)]
      yy <- as.numeric(data[rownames(data) == representative,])
      if(is.null(ylab))ylab = "expression value"
      sub <- paste("Representative:", representative)
    } else if (summary.mode == "median") {
      yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
      if(is.null(ylab))ylab = "gene expression (median)"
      sub <- paste("Median profile of ", nrow(data), " genes")
    } else stop("not valid summary.mode")
    if (dim(data)[1] == 1) {
      sub <- rownames(data)
    }
  } else if (length(data) != 0) {
    yy <- as.numeric(data)
    if(is.null(ylab))ylab = "expression value"
  } else stop("empty data")
  if (is.null(ncol(groups))) {
    ncol = 1
    legend = FALSE
    codeg = "group"
  } else {
    ncol = ncol(groups)
    codeg <- as.character(colnames(groups))
  }
  i.rank<-function (x){
    xx <- x
    for (i in 1:length(xx)) {
      xx[i] <- c(1:length(unique(x)))[x[i] == unique(x)]
    }
    xx
  }
  reps <- i.rank(repvect)
  y <- vector(mode = "numeric", length = length(unique(reps)))
  x <- vector(mode = "numeric", length = length(unique(reps)))
  g <- matrix(nrow = length(unique(reps)), ncol = ncol)
  for (k in 1:length(y)) {
    y[k] <- mean(yy[reps == k], na.rm = TRUE)
    x[k] <- mean(time[reps == k])
    for (j in 1:ncol) {
      g[k, j] <- mean(groups[reps == k, j])
    }
  }
  if(is.vector(data) || nrow(data)<2) plotSEM=FALSE
  if(plotSEM){
    dtm=as.matrix(data)
    errb<-apply(dtm, 2, sd, na.rm = TRUE)
    errb<-errb/sqrt(nrow(data)) 
    errb<-pmax(0.2,errb)
  }
  if (is.null(ylim)){ 
    if(plotSEM){
      ymn<-min(as.numeric(yy-errb), na.rm = TRUE)
      ymx<-max(as.numeric(yy+errb), na.rm = TRUE)      	
    } else {
      ymn<-min(as.numeric(yy), na.rm = TRUE)
      ymx<-max(as.numeric(yy), na.rm = TRUE)
    }
    ylim = c(ymn, ymx+((ymx-ymn)*0.1))
  }
  abcissa <- x
  xlim = c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE)*1.1)
  color1 <- as.numeric(sort(factor(colnames(groups))))
  color2 <- groups
  for (j in 1:ncol) {
    color2[, j] <- color2[, j] * j
  }
  color2 <- as.vector(apply(color2, 1, sum) + 1)
  #order to default colors!
  unicol<-unique(color2)
  color2<-as.factor(color2)
  ordcol<-c(1:length(unicol))
  color2<-order(unicol)[color2]
  if(is.null(startcol)){
    color1=color1+1
    color2[color2==1]=1
  } else {
    if(!is.numeric(startcol))startcol=1
    color1=color1+1
    color2[color2==1]=startcol    
  }
  #
  par(mar = c(3.1, 3.1, 3.1, 1.1))
  plot(x = time, y = yy, pch = 21, ylab = ylab, xlab=xlab, type="n",
       xaxt = "n", main="", sub = NULL, ylim = ylim, xlim = xlim, cex.axis = cex.xaxis,
       cex = cexlab, col = color2, cex.sub=0.6, mgp = c(2.0, 0.6, 0), bty='n')    
  axis(1,at=unique(abcissa), labels=unique(abcissa), cex.axis = cex.xaxis, mgp=c(2.0, 0.6, 0))
  mtext(main, col="black", side = 3, line = 1,cex=0.7,adj = 0,font=2)
  mtext(sub, col="black", side = 3, line = 0,cex=0.6,adj = 0,font=2)
  #
  if(plotdata){
    ally=t(as.matrix(data))
    points(x=matrix(rep(time,ncol(ally)),ncol=ncol(ally)), y=ally, pch=19,col=color2,cex=cex.xaxis*0.6)
  }
  for (i in 1:ncol(groups)) {
    group <- g[, i]
    if (show.lines) {
      lx <- abcissa[group != 0]
      ly <- y[group != 0]
      ord <- order(lx)
      lxo <- lx[ord]
      lyo <- ly[ord]
      lines(lxo, lyo, col = color1[i],lwd=lwd)
    }
  }
  if(plotSEM){
    plotCI(x=time, y=yy, uiw=errb, add=TRUE, barcol=color2, pch=21, col=color2,cex=1.4, pt.bg="white")		
  } else {    
    points(x=time, y=yy, pch=21, col=color2, cex=cex.xaxis*1.2, bg="white")
  }
  op <- par(bg = "white")
  if (legend) 
    legend(legpos, legend=codeg, text.col = color1, col=color1, cex=legsz, lty=1, 
           yjust = 0, bty="n", inset=inset)
  par(op)
}
