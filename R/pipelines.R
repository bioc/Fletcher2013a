###############################################################################
### Pipeline function to run limma analysis
###############################################################################
Fletcher2013pipeline.limma<-function(exprSet, par=list(method="global", adjust.method="BH", p.value=1e-2, lfc=0)){
  cat("Running analysis pipeline ... \n\n")
  if(class(exprSet)!="ExpressionSet")
    stop("Not a 'ExpressionSet' object!")
  opt<-c('Exp1dataset','Exp2dataset','Exp3dataset','siOTHERSdataset','siESR1dataset')%in%names(notes(exprSet))
  if(sum(opt)!=1)stop("'ExpressionSet' object not valid for this pipeline!")
  elab<-which(opt)
  #Extract gene expression data
  geneids <- fData(exprSet)
  targets <- pData(exprSet)
  gexp <- exprs(exprSet)
  myContrasts <- notes(exprSet)$contrasts
  #prepare design model for limma
  f <- factor(targets$Target)
  design <- model.matrix(~0+f)
  colnames(design) <- levels(f)
  #Fit model, set contrasts and run eBayes correction
  fit <- lmFit(gexp,design)
  myContrasts <- makeContrasts(contrasts=myContrasts,levels=design)
  contrastsFit <- eBayes(contrasts.fit(fit, myContrasts))
  contrastsFit$ids <- geneids #add aligned ids!
  #Run decision
  resultsFit <- decideTests(contrastsFit,method=par$method,adjust.method=par$adjust.method,
                            p.value=par$p.value,lfc=par$lfc)
  #Extract results set exp. names!
  reslimma <- extractResults(myContrasts,contrastsFit,resultsFit)
  if(elab==1){
    reslimma$Exp1limma<-"Endogenous FGFRs experiments"
    plotExpDeg(reslimma)
    cat("-Fig1(a) generated!\n")
    save(reslimma,file="Exp1limma.rda")
    cat("-Exp1limma data file generated!\n")
  }
  else if(elab==2){
    reslimma$Exp2limma<-"iF2 construct experiments"
    plotExpDeg(reslimma)
    cat("-Fig2(a) generated!\n")
    save(reslimma,file="Exp2limma.rda")
    cat("-Exp2limma data file generated!\n")
  }
  else if(elab==3){
    reslimma$Exp3limma<-"FGFR2b over-expression experiments" 
    plotExpDeg(reslimma)
    cat("-Fig3(a) generated!\n")
    save(reslimma,file="Exp3limma.rda")
    cat("-Exp3limma data file generated!\n")
  } else if(elab==4){
    reslimma$siRNAlimma<-"siRNA experiments: siPTTG1, siSPDEF, siE2F2 and siELF3"
    save(reslimma,file="siOTHERSlimma.rda")
    cat("-siOTHERSlimma file generated!\n")
  } else if(elab==5){
    reslimma$siRNAlimma<-"siESR1 dataset"
    colnames(reslimma$deg)<-c("PROBEID","SYMBOL","ENTREZ","coef","p.value","degenes")
    save(reslimma,file="siESR1limma.rda")
    cat("-siESR1limma file generated!\n")
  }
  cat("-data preprocessing finished!\n\n")
}
###############################################################################
### Pipeline function to run PCA analysis              
###############################################################################
Fletcher2013pipeline.pca<-function(exprSet){
  opt<-c('Exp1dataset','Exp2dataset','Exp3dataset')%in%names(notes(exprSet))
  if(sum(opt)!=1)stop("'exprSet' not valid for this pipeline!")
  elab<-which(opt)
  #Extract gene expression data
  targets <- pData(exprSet)
  gexp <- exprs(exprSet)
  #Get pre-processed limma analysis
  if(elab==1){
    data("Exp1limma",package="Fletcher2013a")
    resLimma<-get("Exp1limma")
    #Run PCA analysis
    deg <- unique(unlist(resLimma$deglist$all))
    tempgexp <- gexp[is.element(rownames(gexp),deg),]
    tempgexp <- as.data.frame(t(na.omit(tempgexp)))
    rpca <- prcomp(tempgexp,retx=TRUE,center=TRUE,scale=TRUE)
    #Plot the first two principal components
    plotExpPca(rpca,targets,deg,elab)    
    cat("-Fig1(b) generated!\n")
  } else if(elab==2){
    data("Exp2limma",package="Fletcher2013a")
    resLimma<-get("Exp2limma")
    #Run PCA analysis
    deg <- unique(unlist(resLimma$deglist$all))
    tempgexp <- gexp[is.element(rownames(gexp),deg),]
    tempgexp <- as.data.frame(t(na.omit(tempgexp)))
    rpca <- prcomp(tempgexp,retx=TRUE,center=TRUE,scale=TRUE)
    #Plot the first two principal components
    plotExpPca(rpca,targets,deg,elab)
    cat("-Fig2(b) generated!\n")
  } else if(elab==3){
    data("Exp3limma",package="Fletcher2013a")
    resLimma<-get("Exp3limma")
    #Run PCA analysis
    deg <- unique(unlist(resLimma$deglist$all))
    tempgexp <- gexp[is.element(rownames(gexp),deg),]
    tempgexp <- as.data.frame(t(na.omit(tempgexp)))
    rpca <- prcomp(tempgexp,retx=TRUE,center=TRUE,scale=TRUE)
    #Plot the first two principal components
    plotExpPca(rpca,targets,deg,elab)
    cat("-Fig3(b) generated!\n")
  }
}

###############################################################################
### Pipeline function to extrac DE gene lists from limma analysis          
###############################################################################
Fletcher2013pipeline.deg<-function(what="Exp1", idtype="probeid",response="all"){
  opt<-c('Exp1','Exp2','Exp3')%in%what
  if(sum(opt)!=1)
    stop("'what' should be any of 'Exp1', 'Exp2' and 'Exp3'!")
  if(sum(idtype%in%c("probeid","entrez"))!=1)
    stop("'idtype' should be 'probeid' or 'entrez'!")
  if(sum(response%in%c("all","early","late"))!=1)
    stop("'response' should be 'all', 'early' or 'late'!")
  if(which(opt)==1){
    sig<-Exp1signatures(idtype=idtype,response=response)
  }
  else if(which(opt)==2){
    sig<-Exp2signatures(idtype=idtype,response=response)
  }
  else if(which(opt)==3){
    sig<-Exp3signatures(idtype=idtype,response=response)
  }
  return(sig)
}
###############################################################################
### Pipeline function for additional plots          
###############################################################################
Fletcher2013pipeline.supp<-function(){
  #Plot overlap among experiments
  plotOverlap()
  cat("-Fig4(a) generated!\n")
  #Plot followup on IL8
  plotFollowup()
  cat("-Fig4(b-d) generated!\n")
}

