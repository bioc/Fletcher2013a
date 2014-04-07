###############################################################################
### Wrapper functions to extract DE genes from pre-processed data objects  
###############################################################################
## Extract FGFR2 signatures from Exp1limma object
Exp1signatures<-function(idtype="probeid",response="all", mode="all"){
  data("Exp1limma",package="Fletcher2013a")
  dt<-get("Exp1limma")
  #---select relevant DE gene lists from Exp1
  if(response=="early"){
    DT<-unique(unlist(dt$deglist[[mode]][c("t6.UT-t0.UT")]))
    E2<-unique(unlist(dt$deglist[[mode]][c("t6.E2-t6.UT")]))
    E2FGF10<-unique(unlist(dt$deglist[[mode]][c("t6.E2FGF10-t6.E2")]))
    random<-rownames(dt$deg)[sample( 1:nrow(dt$deg) )[1:length(E2FGF10)] ]
  } else if(response=="late"){
    DT<-unique(unlist(dt$deglist[[mode]][c("t24.UT-t6.UT")]))
    E2<-unique(unlist(dt$deglist[[mode]][c("t24.E2-t24.UT")]))
    E2FGF10<-unique(unlist(dt$deglist[[mode]][c("t24.E2FGF10-t24.E2")]))
    random<-rownames(dt$deg)[sample( 1:nrow(dt$deg) )[1:length(E2FGF10)] ]
  } else {
    DT<-unique(unlist(dt$deglist[[mode]][c("t6.UT-t0.UT","t24.UT-t6.UT")]))
    E2<-unique(unlist(dt$deglist[[mode]][c("t6.E2-t6.UT","t24.E2-t24.UT")]))
    E2FGF10<-unique(unlist(dt$deglist[[mode]][c("t6.E2FGF10-t6.E2","t24.E2FGF10-t24.E2")]))
    random<-rownames(dt$deg)[sample( 1:nrow(dt$deg) )[1:length(E2FGF10)] ]
  }
  if(idtype=="entrez"){
    DT<-unique(dt$deg[DT,"ENTREZ"]);DT<-DT[!is.na(DT)]
    E2<-unique(dt$deg[E2,"ENTREZ"]);E2<-E2[!is.na(E2)]
    E2FGF10<-unique(dt$deg[E2FGF10,"ENTREZ"]);E2FGF10<-E2FGF10[!is.na(E2FGF10)]
    random<-unique(dt$deg[,"ENTREZ"]);random<-random[!is.na(random)]
    random<-random[ sample( 1:length(random) )[1:length(E2FGF10)] ]
  }
  #---combine signatures
  sig<-list()
  sig$DT<-DT 
  sig$E2<-E2
  sig$E2FGF10<-E2FGF10
  sig$random<-random
  return(sig)
}
##-----------------------------------------------------------------------------
## Extract FGFR2 signatures from Exp2limma object
Exp2signatures<-function(idtype="probeid",response="all", mode="all"){
  data("Exp2limma",package="Fletcher2013a")
  dt<-get("Exp2limma")
  #---select relevant DE gene lists from Exp2
  if(response=="early"){
    DT<-unique(unlist(dt$deglist[[mode]][c("t6.UT-t0.UT")]))
    E2<-unique(unlist(dt$deglist[[mode]][c("t6.E2-t6.UT")]))
    E2AP20187<-unique(unlist(dt$deglist[[mode]][c("t6.E2.AP20187-t6.E2")]))
  } else if(response=="late"){
    DT<-unique(unlist(dt$deglist[[mode]][c("t24.UT-t6.UT")]))
    E2<-unique(unlist(dt$deglist[[mode]][c("t24.E2-t24.UT")]))
    E2AP20187<-unique(unlist(dt$deglist[[mode]][c("t24.E2.AP20187-t24.E2")]))
  } else {
    DT<-unique(unlist(dt$deglist[[mode]][c("t6.UT-t0.UT","t24.UT-t6.UT")]))
    E2<-unique(unlist(dt$deglist[[mode]][c("t6.E2-t6.UT","t24.E2-t24.UT")]))
    E2AP20187<-unique(unlist(dt$deglist[[mode]][c("t6.E2.AP20187-t6.E2","t24.E2.AP20187-t24.E2")]))
  }
  random<-rownames(dt$deg)[sample( 1:nrow(dt$deg) )[1:length(E2AP20187)] ]
  if(idtype=="entrez"){
    DT<-unique(dt$deg[DT,"ENTREZ"]);DT<-DT[!is.na(DT)]
    E2<-unique(dt$deg[E2,"ENTREZ"]);E2<-E2[!is.na(E2)]
    E2AP20187<-unique(dt$deg[E2AP20187,"ENTREZ"]);E2AP20187<-E2AP20187[!is.na(E2AP20187)]
    random<-unique(dt$deg[,"ENTREZ"]);random<-random[!is.na(random)]
    random<-random[ sample( 1:length(random) )[1:length(E2AP20187)] ]
  }
  #---combine signatures
  sig<-list()
  sig$DT<-DT 
  sig$E2<-E2
  sig$E2AP20187<-E2AP20187
  sig$random<-random
  return(sig)
}
##-----------------------------------------------------------------------------
## Extract FGFR2 signatures from Exp3limma object
Exp3signatures<-function(idtype="probeid",response="all", mode="all"){
  data("Exp3limma",package="Fletcher2013a")
  dt<-get("Exp3limma")
  #---select relevant DE gene lists from Exp3
  if(response=="early"){
    PlusTet<-unique(
      unlist(
        dt$deglist[[mode]]["PlusTet.VEH.6h-MinusTet.VEH.6h"]
      )
    )    
    PlusTet.DT<-unique(
      unlist(
        dt$deglist[[mode]][c("PlusTet.VEH.3h-PlusTet.VEH.0h",
                          "PlusTet.VEH.6h-PlusTet.VEH.3h")]
      )
    )
    PlusTet.E2<-unique(
      unlist(
        dt$deglist[[mode]]["PlusTet.PlusE2.6h-PlusTet.VEH.6h"]
      )
    )
    PlusTet.E2.FGF10<-unique(
      unlist(
        dt$deglist[[mode]]["PlusTet.PlusE2PlusFGF10.6h-PlusTet.PlusE2.6h"]
      )
    )
    random<-rownames(dt$deg)[sample( 1:nrow(dt$deg) )[1:length(PlusTet.E2.FGF10)] ]
  } else if(response=="late"){
    PlusTet<-unique(
      unlist(
        dt$deglist[[mode]]["PlusTet.VEH.24h-MinusTet.VEH.24h"]
      )
    )
    PlusTet.DT<-unique(
      unlist(
        dt$deglist[[mode]][c("PlusTet.VEH.12h-PlusTet.VEH.6h",
                          "PlusTet.VEH.24h-PlusTet.VEH.12h")]
      )
    )
    PlusTet.E2<-unique(
      unlist(
        dt$deglist[[mode]]["PlusTet.PlusE2.24h-PlusTet.VEH.24h"]
      )
    )
    PlusTet.E2.FGF10<-unique(
      unlist(
        dt$deglist[[mode]]["PlusTet.PlusE2PlusFGF10.24h-PlusTet.PlusE2.24h"]
      )
    )
    random<-rownames(dt$deg)[sample( 1:nrow(dt$deg) )[1:length(PlusTet.E2.FGF10)] ]
  } else {
    PlusTet<-unique(
      unlist(
        dt$deglist[[mode]][c("PlusTet.VEH.3h-MinusTet.VEH.3h",
                         "PlusTet.VEH.6h-MinusTet.VEH.6h",
                         "PlusTet.VEH.12h-MinusTet.VEH.12h",
                         "PlusTet.VEH.24h-MinusTet.VEH.24h")]
      )
    )
    PlusTet.DT<-unique(
      unlist(
        dt$deglist[[mode]][c("PlusTet.VEH.3h-PlusTet.VEH.0h",
                          "PlusTet.VEH.6h-PlusTet.VEH.3h",
                          "PlusTet.VEH.12h-PlusTet.VEH.6h",
                          "PlusTet.VEH.24h-PlusTet.VEH.12h")]
      )
    )
    PlusTet.E2<-unique(
      unlist(
        dt$deglist[[mode]][c("PlusTet.PlusE2.3h-PlusTet.VEH.3h",
                          "PlusTet.PlusE2.6h-PlusTet.VEH.6h",
                          "PlusTet.PlusE2.12h-PlusTet.VEH.12h",
                          "PlusTet.PlusE2.24h-PlusTet.VEH.24h")]
      )
    )
    PlusTet.E2.FGF10<-unique(
      unlist(
        dt$deglist[[mode]][c("PlusTet.PlusE2PlusFGF10.3h-PlusTet.PlusE2.3h",
                          "PlusTet.PlusE2PlusFGF10.6h-PlusTet.PlusE2.6h",
                          "PlusTet.PlusE2PlusFGF10.12h-PlusTet.PlusE2.12h",
                          "PlusTet.PlusE2PlusFGF10.24h-PlusTet.PlusE2.24h")]
      )
    )
    random<-rownames(dt$deg)[sample( 1:nrow(dt$deg) )[1:length(PlusTet.E2.FGF10)] ]
  }
  if(idtype=="entrez"){
    PlusTet<-unique(dt$deg[PlusTet,"ENTREZ"])
    PlusTet.DT<-unique(dt$deg[PlusTet.DT,"ENTREZ"])
    PlusTet.DT<-PlusTet.DT[!is.na(PlusTet.DT)]
    PlusTet.E2<-unique(dt$deg[PlusTet.E2,"ENTREZ"])
    PlusTet.E2<-PlusTet.E2[!is.na(PlusTet.E2)]
    PlusTet.E2.FGF10<-unique(dt$deg[PlusTet.E2.FGF10,"ENTREZ"])
    PlusTet.E2.FGF10<-PlusTet.E2.FGF10[!is.na(PlusTet.E2.FGF10)]
    random<-unique(dt$deg[,"ENTREZ"]);random<-random[!is.na(random)]
    random<-random[ sample( 1:length(random) )[1:length(PlusTet.E2.FGF10)] ]
  }
  #---combine signatures
  sig<-list()
  sig$Tet<-PlusTet
  sig$TetDT<-PlusTet.DT
  sig$TetE2<-PlusTet.E2
  sig$TetE2FGF10<-PlusTet.E2.FGF10
  sig$random<-random
  return(sig)
}

