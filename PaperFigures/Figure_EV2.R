  source('./Conf.R')
  source("./Main.R")

  
  give.n <- function(x){
    return(c(y = median(x)+0.15, label = length(x))) 
  }
  
  ########### EV2 A
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  RNA_Protein_Combined$UTR5L <- sapply(RNA_Protein_Combined$UTR5Seq, nchar)
  atgCount <- unlist(vcountPattern("AUG", RNA_Protein_Combined$UTR5Seq))
  atgPos <- start(vmatchPattern("AUG", RNA_Protein_Combined$UTR5Seq))
  
  RNA_Protein_Combined$Motif_UTR5_AUG  <- as.factor(sapply(1:length(atgPos), function(i){
    if(length(atgPos[[i]])==0){
      return(0)
    }else{
      atgDist <- RNA_Protein_Combined[i, "UTR5L"] - atgPos[[i]]+1
      if(any((atgDist %% 3) != 0)){
        return(2)
      }else{
        return(1)
      }
    }
  }))
  
  
  allfeatureCols <- getAllFeatureColumns(RNA_Protein_Combined)
  allfeatureCols[c("UTR5_ATG","UTR5_Motifs","kozakSequence")] <- NULL
  
  RNA_Protein_Combined$rowMedianPTR <- getYPrime(inputDF=RNA_Protein_Combined, 
                                                 feature="",
                                                 responseCol = "rowMedianPTR",
                                                 allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatureCols)]))
  
  RNA_Protein_Combined$rowMedianPTR_temp <- 10^RNA_Protein_Combined$rowMedianPTR
  
  Fig_EV2A <- ggplot(RNA_Protein_Combined, aes(x=factor(Motif_UTR5_AUG),y=rowMedianPTR_temp))+
    geom_boxplot(fill=c("grey", "lightyellow","lightblue"),outlier.shape = NA, coef = 0) + 
    scale_y_continuous(trans='log10', 
                       breaks = c(0.1, 1,10), 
                       labels = trans_format("log10", 
                                             math_format(10^.x)))+
    stat_summary(fun.data = give.n, 
                 geom = "text", 
                 fun.y = "median",
                 position = position_dodge(width = 0.5),
                 size=3)+
    coord_cartesian(ylim=c(0.1,10))+
    labs(x="",y="Median PTR ratio \n corrected for other sequence features")+
    scale_x_discrete(labels=c("No uAUG", "Only in-frame \n uAUG(s)", ">=1 out-frame \n uAUG(s)"))
  
  Fig_EV2A <- add_pval(Fig_EV2A, 
                       pairs=list(c(1,2), c(1,3)), 
                       log=F, 
                       heights = c(4,6), 
                       fold_change = T, 
                       textsize = 3)
  
  ########### EV2 B
  
  findFirstCodonOccurancePosition <- function(codon, Sequence){
    
    a <- start(vmatchPattern(codon, RNAStringSet(Sequence), fixed=T))
    
    a<- lapply(a, function(x){
      if(length(x)==0)
      {return(0)} 
      else
      {return(x)}})
    
    return(a)
  }
  
  findInframeStopCodonDist <- function(upATGStarts, stopCodonStartsList, UTR5Length=NULL){
    
    inFrameStopCodDist = c()
    
    for(i in 1:length(upATGStarts)){
      if(upATGStarts[[i]] != 0 & stopCodonStartsList[[i]] != 0){
        
        upStartCodStart <- upATGStarts[[i]]
        stopCodonStarts <- stopCodonStartsList[[i]]
        
        if(!is.null(UTR5Length)){
          UTR5len <- UTR5Length[[i]]
          stopCodonStarts <- (stopCodonStarts+UTR5len)
        }
        
        stopCodonStarts <- stopCodonStarts[((stopCodonStarts-upStartCodStart) > 0)]
        stopCodonStarts <- stopCodonStarts[((stopCodonStarts-upStartCodStart) %% 3) == 0]
        if(length(stopCodonStarts) > 0){
          inFrameStopCodDist <- c(inFrameStopCodDist, min(stopCodonStarts-upStartCodStart-3)) 
        }else{
          inFrameStopCodDist <- c(inFrameStopCodDist, 0)
        }
        
      }else{
        inFrameStopCodDist <- c(inFrameStopCodDist, -1)
      }
      
    }
    return(inFrameStopCodDist)
  }
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatureCols <- getAllFeatureColumns(RNA_Protein_Combined)
  allfeatureCols[c("UTR5_ATG","UTR5_Motifs","kozakSequence")] <- NULL
  
  RNA_Protein_Combined$rowMedianPTR <- getYPrime(inputDF=RNA_Protein_Combined, 
                                                 feature="",
                                                 responseCol = "rowMedianPTR",
                                                 allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatureCols)]))
  
  RNA_Protein_Combined$rowMedianPTR <- 10^RNA_Protein_Combined$rowMedianPTR
  
  RNA_Protein_Combined$UTR5L <- sapply(RNA_Protein_Combined$UTR5Seq, nchar)
  RNA_Protein_Combined$atgCount <- unlist(vcountPattern("AUG", RNA_Protein_Combined$UTR5Seq))
  RNA_Protein_Combined_temp <- RNA_Protein_Combined[RNA_Protein_Combined$atgCount< 2,]
  
  RNA_Protein_Combined_temp$upATGStart <- unlist(findFirstCodonOccurancePosition("AUG", 
                                                                                 RNAStringSet(RNA_Protein_Combined_temp$UTR5Seq)))
  RNA_Protein_Combined_temp$upATGdistanceToCanStart <- RNA_Protein_Combined_temp$UTR5L - (RNA_Protein_Combined_temp$upATGStart -1)
  RNA_Protein_Combined_temp[RNA_Protein_Combined_temp$atgCount==0,"upATGdistanceToCanStart"] <- 0
  
  RNA_Protein_Combined_temp$upATGframeWrtCanStart <- as.integer(RNA_Protein_Combined_temp$upATGdistanceToCanStart) %% 3
  RNA_Protein_Combined_temp$upATGInFrameWithStart <- F
  RNA_Protein_Combined_temp[RNA_Protein_Combined_temp$atgCount == 1 & RNA_Protein_Combined_temp$upATGframeWrtCanStart == 0, 
                            "upATGInFrameWithStart"] <- T
  
  
  RNA_Protein_Combined_temp$upATGframeWrtCanStart <- factor(RNA_Protein_Combined_temp$upATGframeWrtCanStart, levels = c("NoATG",0,1,2))
  RNA_Protein_Combined_temp[RNA_Protein_Combined_temp$atgCount==0,"upATGframeWrtCanStart"] <- "NoATG"
  
  
  b <- findFirstCodonOccurancePosition("UAA", RNA_Protein_Combined_temp$UTR5Seq)
  RNA_Protein_Combined_temp$inFrameTAADists <- findInframeStopCodonDist(upATGStarts=RNA_Protein_Combined_temp$upATGStart,
                                                                        stopCodonStartsList=b)
  
  b <- findFirstCodonOccurancePosition("UGA", RNA_Protein_Combined_temp$UTR5Seq)
  RNA_Protein_Combined_temp$inFrameTGADists <- findInframeStopCodonDist(upATGStarts=RNA_Protein_Combined_temp$upATGStart, 
                                                                        stopCodonStartsList=b)
  
  b <- findFirstCodonOccurancePosition("UAG", RNA_Protein_Combined_temp$UTR5Seq)
  RNA_Protein_Combined_temp$inFrameTAGDists <- findInframeStopCodonDist(upATGStarts=RNA_Protein_Combined_temp$upATGStart, 
                                                                        stopCodonStartsList=b)
  
  
  RNA_Protein_Combined_temp$inFrameStopExistsIn5UTR <- (RNA_Protein_Combined_temp$inFrameTAADists > 0 |
                                                          RNA_Protein_Combined_temp$inFrameTGADists > 0 | 
                                                          RNA_Protein_Combined_temp$inFrameTAGDists > 0)
  
  RNA_Protein_Combined_temp <- RNA_Protein_Combined_temp[, c("rowMedianPTR", "atgCount", "upATGInFrameWithStart", "inFrameStopExistsIn5UTR")]
  RNA_Protein_Combined_temp$xAxis <- paste0(RNA_Protein_Combined_temp$atgCount,"_",
                                            RNA_Protein_Combined_temp$upATGInFrameWithStart, "_", 
                                            RNA_Protein_Combined_temp$inFrameStopExistsIn5UTR )
  
  RNA_Protein_Combined_temp$xAxis <- factor(RNA_Protein_Combined_temp$xAxis, 
                                            levels = c("0_FALSE_FALSE", "1_TRUE_FALSE", "1_TRUE_TRUE", "1_FALSE_FALSE", "1_FALSE_TRUE"), 
                                            labels = c("No uAUG", "1 in-frame uAUG", "1 in-frame uORF", "1 out-frame uAUG", "1 out-frame uORF"))
  
  
  Fig_EV2B <- ggplot(data=RNA_Protein_Combined_temp,
                     aes(x=xAxis, y=rowMedianPTR)) + 
    geom_boxplot(fill=c("grey", "lightyellow","lightyellow", "lightblue", "lightblue"), 
                 outlier.shape = NA,
                 coef = 0) +
    stat_summary(fun.data = give.n, 
                 geom = "text", 
                 fun.y = "median",
                 position = position_dodge(width = 0.5),
                 size=3)+
    labs(x="",y="Median PTR ratio \n corrected for other sequence features")+
    scale_y_continuous(trans='log10', 
                       breaks = c(0.1, 1,10), 
                       labels = trans_format("log10", 
                                             math_format(10^.x)))+
    coord_cartesian(ylim=c(0.1,500))+
    theme(axis.text.x=element_text(angle = 90, hjust = 1 ))
  
  Fig_EV2B <- add_pval(Fig_EV2B, 
                       pairs=list(c(1,2), c(1,3), c(1,4), c(1,5), c(4,5)), 
                       log=T, 
                       heights = c(4,10,50,150, 10), 
                       fold_change = T, 
                       textsize = 3)
  
  
  
  ######################################
  ### EV2 D
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  
  
  allFeaturesFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                               allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatures)]),
                                               responseStr=c("_PTR"))
  
  ############################################################################
  ######### Kozak sequence effect 
  
  kozak_TV <- getBetasOfFeature(allFeaturesFit, featureSubStr="Kozak_")
  posVec <- unlist(str_split(colnames(kozak_TV$betas), "Kozak_"))
  posVec <- posVec[posVec!=""]
  
  PosNames <- sapply(posVec, function(x){return(strsplit(x,"\\.")[[1]][1])})
  PosNames <- unique(PosNames)
  
  FigEV2_D<- plot_lm(lm_out=NULL,
                    posVec = posVec,
                    numPos=15,
                    PosNames=PosNames,
                    ylabStr="", 
                    mainStr=mainStr,
                    effectLog=T, 
                    allCoefs= as.data.frame(kozak_TV$betas), 
                    xlabStr="Position around start codon", 
                    breaksSeq=seq(0.6, 1.3, 0.1), 
                    limitSeq=c(0.6,1.3))
  

  FigEV2_D_2 <-myGGSeqLogo(sequences=RNA_Protein_Combined$kozakSeq, 
                        xlabStr="Position around start codon", 
                        motifLength=3, 
                        lrOffset=6)
  FigEV2_D_2 <- FigEV2_D_2 + theme(axis.title.x=element_text(size=7))
  

  ###############################################################################
  ####### stop codon context EV2E
  stop_TV <- getBetasOfFeature(allFeaturesFit, featureSubStr="stop_")
  posVec <- unlist(str_split(colnames(stop_TV$betas), "stop_"))
  posVec <- posVec[posVec!=""]
  
  PosNames <- sapply(posVec, function(x){return(strsplit(x,"\\.")[[1]][1])})
  PosNames <- unique(PosNames)
  
  FigEV2_E<- plot_lm(lm_out=NULL,
                    posVec = posVec,
                    numPos=15,
                    PosNames=PosNames,
                    ylabStr="", 
                    mainStr=mainStr,
                    effectLog=T, 
                    allCoefs=as.data.frame(stop_TV$betas), 
                    xlabStr="Position around stop codon", 
                    breaksSeq=seq(0.6, 1.3, 0.1), 
                    limitSeq=c(0.6,1.3))
  
  
  FigEV2_E_2 <-myGGSeqLogo(sequences=RNA_Protein_Combined$sCodContext, 
                        xlabStr="Position around stop codon", 
                        motifLength=3, 
                        lrOffset=6)
  
  FigEV2_E_2 <- FigEV2_E_2 + theme(axis.title.x=element_text(size=7))
  
  
  #################################################
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatureCols <- getAllFeatureColumns(RNA_Protein_Combined)
  allfeatureCols$stopCodonContext<- NULL
  RNA_Protein_Combined$rowMedianPTR <- getYPrime(inputDF=RNA_Protein_Combined, 
                                                 feature=NULL,
                                                 responseCol = "rowMedianPTR", 
                                                 allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatureCols)]))
  
  
  SopPlusOne <- apply(RNA_Protein_Combined, 1, function(x){ paste0(substr(x["CDSSeq"], nchar(x["CDSSeq"])-2,nchar(x["CDSSeq"])), 
                                                                   substr(x["UTR3Seq"], 1, 1 )  )})
  k<-which(substr(SopPlusOne,1,1) == "U")
  SopPlusOne <- SopPlusOne[k]
  RNA_Protein_Combined <- RNA_Protein_Combined[k,]
  
  stopCodon <- substr(SopPlusOne,1,3)
  stopPlusOneDF <- data.frame(SopPlusOne=SopPlusOne, 
                              response=(RNA_Protein_Combined$rowMedianPTR), 
                              stopCod=stopCodon)
  
  stopPlusOneDF$SopPlusOne <- factor(stopPlusOneDF$SopPlusOne, 
                                     levels = c("UAAG","UAAA","UAAU","UAAC","UAGG","UAGA","UAGU","UAGC","UGAU","UGAG","UGAA","UGAC") , 
                                     ordered = TRUE)
  
  
  give.n <- function(x){
    return(c(y = median(x)+0.2, label = length(x))) 
  }
  
  stopPlusOneDF$response <- 10^stopPlusOneDF$response
  
  
  Fig_EV2F <- ggplot(data=stopPlusOneDF, aes(x=stopCod, y=response, group=stopCod, fill=stopCodon)) + 
    geom_boxplot(outlier.shape = NA)+
    stat_summary(fun.data = give.n,
                 geom = "text", 
                 fun.y = "median",
                 position = position_dodge(width = 0.75),
                 size=rel(3)) +
    ylab(" Median PTR ratio \n corrected for other sequence features") +
    xlab("")+
    theme( legend.position = "None")+
    scale_y_continuous(trans='log10')
  
  Fig_EV2F <- add_pval(Fig_EV2F, 
                       fold_change = T,
                       pairs = list(c(1,2), c(1,3), c(2,3)), 
                       test='wilcox.test', 
                       heights = c(100,400,1200),
                       textsize = 3)
  
  
  ########### EV4 D
  Fig_EV2G <- ggplot(data=stopPlusOneDF, 
                     aes(x=reorder(SopPlusOne,stopCodon), 
                         y=response, 
                         group=SopPlusOne, 
                         fill=stopCod)) + 
    geom_boxplot(outlier.shape = NA)+
    stat_summary(fun.data = give.n, 
                 geom = "text", 
                 fun.y = "median",
                 position = position_dodge(width = 0.75),
                 size=rel(2)) +
    ylab(" Median PTR ratio \n corrected for other sequence features") +
    xlab("")+
    theme( legend.position = "None",
           axis.text.x=element_text(angle = 90, hjust = 1 ))+
    scale_y_continuous(trans='log10') 
  