  save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  

  give.n <- function(x){
    return(c(y = median(x)+0.4, label = length(x))) 
  }
  
  #############################################################
  ########## This function plots the boxplots of the motifs 
  
  boxplot_kmer <- function(motif, regionString, inputDF, titleStr, ylabStr, binary, response, xlabStr, colM="lightyellow"){
    
    motifWithRegion <- paste0(regionString,motif)
    
    if(binary){
      minCutOff= 0
    }else{
      minCutOff = min(as.integer(names(table(inputDF[[motifWithRegion]])[table(inputDF[[motifWithRegion]]) < 150])))-1
      if(minCutOff == Inf){
        minCutOff= 1
      }else{
        inputDF[ inputDF[[motifWithRegion]] > minCutOff,motifWithRegion] <- 100000
      }
    }
    
    foldChangeNatScale <- formattable(10^( median(inputDF[inputDF[,motifWithRegion] == 1, response]) - median(inputDF[inputDF[,motifWithRegion]==0, response])), digits = 2, format = "f")
    
    pvalue <- wilcox.test(2^(inputDF[[response]]) ~ cut(inputDF[[motifWithRegion]], c(-Inf,0,Inf)))$p.value
    
    if(pvalue < 0.00001){
      pvalue=0
    }else{
      pvalue = round(pvalue, digits=2)
    }
    
    inputDF[[response]] <- 10^inputDF[[response]]
    
    noMotifMedian <- median(inputDF[ inputDF[[motifWithRegion]] ==0 ,response])
    
    k <- data.frame(response=(inputDF[[response]]), group=inputDF[[motifWithRegion]] )
    p <- ggplot(k, aes(x=factor(group),y=response))+
      geom_boxplot(fill=c("grey", rep(colM,length(unique(factor(k$group)))-1)),outlier.shape = NA, varwidth = T, coef = 0) + 
      labs(x=xlabStr,y=ylabStr)+
      theme_minimal() +
      scale_x_discrete(labels=c(0:minCutOff,paste('≥',minCutOff+1))) +
      stat_summary(fun.data = give.n, geom = "text", fun.y = "median",position = position_dodge(width = 0.5),size=rel(6))+
      theme(legend.position = "none",
            axis.title = element_text(size=rel(2)), 
            axis.text = element_text(size=rel(2)),
            plot.title = element_text(size=rel(1.5), hjust = 0.5))+
      coord_trans(limy = quantile(inputDF[[response]], probs=c(0.05,0.95)), y="log10") +
      scale_y_continuous(trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))
    #+geom_hline(color="black", yintercept=noMotifMedian)
    
    return(p)
  }
  
  ############################################################
  ########## This function plots the beta values of a feature across tissues
  
  plotMotifEffectAcrossTissues <- function(betas, tissueNames, stdErrors,  ylabStr, breaksSeq, labelsSeq=NULL){
    
    if(is.null(labelsSeq)){
      labelsSeq <- breaksSeq
    }
    temp <- data.table(betas=betas, names= tissueNames, stdErr=stdErrors)
    colnames(temp) <- c("betas","names","stdErr")
    temp <- within(temp, names <- factor(names, levels=unique(names)))
    maxStd <- max(temp$stdErr)
    
    temp[,yMin:=10^(betas-1.96*stdErr),][,yMax:=10^(betas+1.96*stdErr)]
    temp[,betas:=10^betas]
    
    if (median(temp$betas) < 1){
      ylimC <- c(min(temp$yMin), max(1.02,max(temp$yMax)))
    }else{
      ylimC <- c(min(0.97,min(temp$yMin)), max(temp$yMax))
    }
    
    p<- ggplot(temp, aes(x = names, y = betas)) +
      geom_errorbar(aes(ymin = yMin, ymax = yMax), size=0.3) + 
      geom_point(size=0.5) +
      ylab(ylabStr)+
      xlab("")+
      scale_y_continuous(trans='log2', 
                         breaks = breaksSeq, 
                         labels = labelsSeq, 
                         limits = ylimC)+
      geom_hline(color="darkgrey", yintercept=1)+
      theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
    return(p)
  }
  
  plotMotifEffectAcrossTissuesForMultipleM <- function(betas, 
                                                       TissueNamesPrint, 
                                                       stdErrs,  
                                                       ylabStr="Effect size in multivariate model",
                                                       inputDF, breaksSeq, labelsSeq=NULL, changeNames=T){
    if(is.null(labelsSeq)){
      labelsSeq <- breaksSeq
    }
    motifOrder <- colnames(betas)
    
    betas$tissueNames <- TissueNamesPrint
    betas <- melt(betas)
    betas <- rename.vars(betas, from = c("variable", "value"), 
                         to = c("Motif", "betas"))
    
    nOfGenes <- getNumberOfGenesHavingTheMotifs(as.character(unique(betas$Motif)), inputDF)
    betas <- merge(betas, nOfGenes, by="Motif")
    betas$Motif <- paste0(betas$tissueNames,"_JJ_",betas$Motif)
    
    
    stdErrs$tissueNames <- TissueNamesPrint
    stdErrs <- melt(stdErrs, variable_name="Tissues")
    stdErrs <- rename.vars(stdErrs, from = c("variable", "value"), 
                           to = c("Motif", "stdErr"))
    stdErrs$Motif <- paste0(stdErrs$tissueNames,"_JJ_",stdErrs$Motif)
    
    
    
    allBetas <- merge(betas,stdErrs, by="Motif")
    if(changeNames){
      allBetas$consMotif <- apply(allBetas, 1, function(x){strsplit(x["Motif"], "_")[[1]][5]})
      allBetas$consMotif <- gsub("T", "U", allBetas$consMotif)
    }else{
      allBetas$consMotif <- apply(allBetas, 1, function(x){strsplit(x["Motif"], "_JJ_")[[1]][2]})
    }
    
    allBetas$consMotif <- paste0(allBetas$consMotif," (", allBetas$value, ")")
    allBetas$consMotif <- factor(allBetas$consMotif)
    allBetas$tissueNames.y <-  factor(allBetas$tissueNames.y, levels=unique(allBetas$tissueNames.y))
    
    allBetas$yMin <- allBetas$betas-2*allBetas$stdErr
    allBetas$yMax <- allBetas$betas+2*allBetas$stdErr
    allBetas[,c("betas", "yMin", "yMax")] <- 10^allBetas[,c("betas", "yMin", "yMax")]
    
    p<- ggplot(allBetas, aes(x = tissueNames.y, y = betas), size=0.3) +
      geom_errorbar(aes(ymin = yMin, ymax = yMax), size=0.3) + 
      facet_wrap(~ consMotif, nrow = 1)+
      geom_point(size=0.3) +
      geom_hline(yintercept = 1, col="darkgrey")+
      xlab("")+
      ylab(ylabStr)+
      scale_y_continuous(trans='log2', 
                         breaks = breaksSeq , 
                         labels = labelsSeq)+
      theme(panel.spacing = unit(panelSpacingUnit, "lines"), 
            #axis.text.x = element_text(angle = 90, vjust=0.5, size=rel(1.5), hjust=1)
            axis.text.x = element_blank())
  
    return(p)
  }
  
  plotCDSMotifEffectAcrossTissues <- function(betas, stdErrors, numberOfMotifsPerFrame,motif){
    
    tp <- data.frame(variable=names(numberOfMotifsPerFrame), 
                     noOfTrans=numberOfMotifsPerFrame)
    
    betas$tissueNames <- TissueNamesPrint
    betas <- melt(betas)
    betas <- rename.vars(betas, from = c("value"), to = c("betas"))
    
    betas$Motif <- motif
    betas$Frame <- unlist(lapply(as.character(betas$variable),
                                 function(x){strsplit(x,split = "_")[[1]][3]}))
    
    betas <- merge(betas, tp, by="variable")
    betas$variable <- paste0(betas$tissueNames,"_",betas$variable)
    
    
    stdErrors$tissueNames <- TissueNames
    stdErrors <- melt(stdErrors, variable_name="Tissues")
    stdErrors <- rename.vars(stdErrors, from = c("value"),  to = c("stdErr"))
    stdErrors$variable <- paste0(stdErrors$tissueNames,"_",stdErrors$variable)
    
    
    allBetas <- merge(betas,stdErrors, vy="variable")
    allBetas$yMin <- 10^(allBetas$betas - 1.96*allBetas$stdErr)
    allBetas$yMax <- 10^(allBetas$betas + 1.96*allBetas$stdErr)
    allBetas$betas <- 10^allBetas$betas
    
    motifNames=paste0(motif, " (",glue::collapse(numberOfMotifsPerFrame, sep=" "), " )")
    
    p<- ggplot(allBetas, aes(x = tissueNames, y = betas, col=Frame)) +
      geom_errorbar(aes(ymin = yMin, ymax = yMax)) + 
      geom_point() +
      geom_hline(yintercept = 1, lty = "dashed", col="red")+
      xlab("")+
      theme( axis.text.x = element_text(angle = 90))+
    scale_y_continuous(trans='log2', 
                       breaks=seq(0.4,3,0.6), 
                       labels=seq(0.4,3,0.6))+
      ylab("Effect on PTR ratio")+
      ggtitle(motifNames)
      
    return(p)
  }
  
  
  #################################################################
  ############## This function fits the SNP effect 
  
  
  plot1bpEffect <- function(dmI, responseStr, inputDF, mainStr, xlabStr, ylabStr, effectLog, statCor=F, breaksSeq, limitSeq){
    allCoefs <- list()
    
    responseCols = c()
    for(elem in responseStr){
      responseCols <- c(responseCols, getColumnsOfSubstr(X=inputDF, inputStr=elem))
    }
    responseNames <- colnames(inputDF[,responseCols])
    
    for(correlateType in responseNames){
      correlate = inputDF[[correlateType]]
      lmfI = fitLmModel(dmI,correlate)
      k <- as.data.frame(lmfI[[1]])
      betas <- k[-c(which(rownames(k) %in% c("(Intercept)", "len"))),"estimate"]
      
      if(statCor){
        pvals <- k[-c(which(rownames(k) %in% c("(Intercept)", "len"))),"pvals"]
        betas[pvals > 0.05] <- NA
      }
      
      allCoefs <- lappend(allCoefs, betas)
    }
    
    allCoefs <- do.call(rbind,allCoefs)
    d <- plot_lm(lm_out=lmfI[[1]],lmfI$lm_pos,dmI$num_pos,dmI$posnames, ylabStr, mainStr=mainStr,effectLog, allCoefs, xlabStr, breaksSeq, limitSeq)
    return(d)
  }
  
  ###########################################################
  plot_lm = function(lm_out=NULL, posVec, numPos, PosNames, ylabStr, mainStr, effectLog=T, allCoefs=NULL, xlabStr, breaksSeq, limitSeq)
  {
    
    cm = matrix(NA,nrow = 4,ncol=numPos)
    colnames(cm) = PosNames
    rownames(cm) = c("A","C","G","U")
    
    if(!is.null(allCoefs)){
      allCoefs <- as.data.frame(allCoefs)
      rCoef <- apply(allCoefs, 2,function(x){median(x, na.rm = T)})
      rMin <- apply(allCoefs, 2,function(x){min(x, na.rm = T)})
      rMax <- apply(allCoefs, 2,function(x){max(x, na.rm = T)})
      
      rMin[is.infinite(rMin)] <- NA
      rMax[is.infinite(rMax)] <- NA
      
      if(effectLog){
        f_mat = 10^rCoef
        rMin <- 10^rMin
        rMax <- 10^rMax
      }
    }else{
      rCoef <- coef(lm_out)[c(-1,-2)]
      if(effectLog){
        f_mat = 10^rCoef
      }
    }
    
    ind = as.data.frame(do.call(rbind,strsplit(posVec,'.',fixed=TRUE)))
    ind$V1 <- unlist(lapply(str_split(ind$V1,"P"), function(x){return(x[[2]])}))
    ind$V1 <- gsub("_", "-", ind$V1)
    ind$V1 <- gsub("J", "+", ind$V1)
    
    if(is.null(allCoefs)){
      plotdat = data.frame(coef=f_mat,
                           base=as.factor(ind$V2),
                           positions=as.factor(as.numeric(ind$V1)),
                           stringsAsFactors=F)
      d <- ggplot(data=plotdat,aes(y=coef,x=positions,group=base,colour=base))  +
        geom_point(size=0.5) +  
        geom_line(size=0.2) +
        labs(y=ylabStr, x=xlabStr) + 
        geom_hline(yintercept = 1) + 
        ylim(min(plotdat$coef, na.rm = T)-0.1,max(plotdat$coef, na.rm = T)+0.1)+ 
        scale_colour_manual(values = c("A" = "LawnGreen", "C" = "Blue", "G"="Orange", "U"="Red"))
    }else{
      plotdat = data.frame(coef=f_mat, 
                           rMin = rMin,
                           rMax=rMax,
                           base=as.factor(ind$V2),
                           positions=factor(ind$V1,
                                            levels = unique(ind$V1)),stringsAsFactors=F)
     
       d <- ggplot(data=plotdat,aes(y=coef,x=positions,group=base,colour=base))  +
        geom_point(size=0.5) +  
        geom_errorbar(aes(ymin = rMin, ymax = rMax), size=0.3)+
        labs(y=ylabStr, x="") + 
        theme(legend.position="top") + 
        geom_hline(yintercept = 1, color="darkgrey") + 
        scale_colour_manual(values = c("A" = "LawnGreen", "C" = "Blue", "G"="Orange", "U"="Red"))+
        scale_y_continuous(limits=limitSeq, 
                           trans='log2', 
                           breaks=breaksSeq, 
                           labels=breaksSeq)
    }
    
    
    return(d)
  }
  
  ########################################################
  
  plotGeneCodonScatter <- function(tempDF, xStr, yStr, colStr, xlabStr, ylabStr, labelStr){
    
    p<- ggplot(tempDF, aes_string(x=xStr, y=yStr, col=colStr)) + 
      geom_point(size=0.7, alpha=0.7)+
      xlab(xlabStr) +
      ylab(ylabStr)+
      scale_colour_manual(labels=labelStr, values = c("orange", "blue", "darkgrey"))+
      scale_y_continuous(trans='log10',  
                         breaks = trans_breaks("log10", function(x) 10^x), 
                         labels = trans_format("log10", math_format(10^.x)))+
      scale_x_continuous(trans='log10',  
                         breaks = trans_breaks("log10", function(x) 10^x), 
                         labels = trans_format("log10", math_format(10^.x)))+
      labs(col="Gene type")+
      theme(legend.text = element_text(size=rel(1.2)),
            legend.title = element_text(size=rel(1.5)))
    return(p)
  }
  
  
  plotBoxPlot <- function(inputDF,
                          responseStr,
                          featureStr,
                          xlabStr,
                          ylabStr,
                          colM="lightyellow",
                          labelStr,
                          mainStr="",
                          transLog=T){
    
    k <- data.frame(response=inputDF[[responseStr]],
                    group=inputDF[[featureStr]])
    
    baseMedian <- median(k[ k$group ==0 ,"response"], na.rm = T)
    
    if(transLog){
      k$response <- 10^k$response
    }
    
    p <- ggplot(k, aes(x=factor(group),y=response))+
      geom_boxplot(fill=c("grey", rep(colM,length(unique(factor(k$group)))-1)),
                   coef = 0, outlier.size = 0.1) + 
      labs(x=xlabStr,y=ylabStr)+
      scale_x_discrete(labels=labelStr) +
      stat_summary(fun.data = give.n,
                  geom = "text",
                  fun.y = "median",
                  position = position_dodge(width = 0.5),
                  size=3)+
      theme(legend.position = "none")+
      geom_hline(color="grey", yintercept=baseMedian)+
      ggtitle(mainStr)
    
    if(transLog){
      p <- p+ scale_y_continuous(trans='log10', 
                                 breaks = trans_breaks("log10", function(x) 10^x), 
                                 labels = trans_format("log10", math_format(10^.x)))
    }
    
    return(p)
  }
  
  plotHeatMaps<- function(betasTV, pValsTV, featureSubstr, upThreshold, downThreshold, fileName, sortedNames){
    
    betasTV[pValsTV>0.05] <- NA
    betasTV <- 10^betasTV
    tmp<- unlist(strsplit(colnames(betasTV),featureSubstr))
    tmp<- tmp[tmp!=""]
    colnames(betasTV) <- tmp
    betasTV <- betasTV[,sortedNames]
    rownames(betasTV) <- TissueNamesPrint
    betasTV[betasTV<downThreshold] <- downThreshold
    betasTV[betasTV>upThreshold] <- upThreshold
    
    
    pheatmap(as.matrix(betasTV), cluster_rows = F, cluster_cols = F, treeheight_col = 0,fontsize_col = 15, fontsize_row = 10, na_col = "black", color=colorRampPalette(rev(brewer.pal(n = 7, name = ColDiverging)))(100))
    
  }
  
  
  plotAAAndCodonEffect <- function(inputDF,
                                   codonCols,
                                   aaCols,
                                   y,
                                   responseLabel,
                                   logBase=10){
    
    X <- inputDF[,codonCols]
    results <- getOLSR2(X, 
                        y,
                        featureNames = colnames(X))
    
    tp <- data.frame(codonEffects = results$lmsumDT$Estimate[-1],
                     codon = results$lmsumDT$covariates[-1],
                     codonStdEr = results$lmsumDT$std.error[-1],
                     stringsAsFactors = F)
    
    tp$codon <- sapply(tp$codon, 
                       function(x){strsplit(x,"Codon_")[[1]][2]})
    
    tp$aa <- genCode[tp$codon]
    tp <- as.data.table(tp)
    
    X <- inputDF[,aaCols]
    results <- getOLSR2(X, 
                        y,
                        featureNames = colnames(X))
    
    tpAA <- data.table(AAeffects = results$lmsumDT$Estimate[-1],
                       aa = results$lmsumDT$covariates[-1],
                       AAstEr = results$lmsumDT$std.error[-1],
                       stringsAsFactors = F)
    tpAA$aa <- sapply(tpAA$aa, 
                      function(x)
                        {strsplit(x, "AminoAcid_")[[1]][2]})
    tpAA <- tpAA[order(logBase^tpAA$AAeffects),]
    aaOrder <- tpAA$aa
    tpAA$aa <- factor(tpAA$aa, levels = tpAA$aa)
    
    tpAA[,yMin:=logBase^(AAeffects-1.96*AAstEr),][,yMax:=logBase^(AAeffects+1.96*AAstEr)]
    tpAA[,AAeffects:=logBase^AAeffects,]
    
    p1 <- ggplot(tpAA, aes(x=aa, y=AAeffects, col=aa))+
      geom_point()+
      geom_errorbar(aes(ymin = yMin, ymax = yMax), size=0.3)+
      scale_color_manual(values=rep(c("#999999","#E69F00"),10))+
      theme(panel.spacing = unit(0, "lines"), 
            axis.text.x = element_text(size=15, hjust = 1, vjust=0.5), 
            plot.title = element_text(size=15, hjust = 0.5), 
            legend.position = "None", strip.text = element_text(size=10))+
      geom_hline(yintercept=1, color="darkgrey")+
      scale_y_continuous(trans='log2', breaks = seq(0.8,1.4,0.1) ,
                         labels = seq(0.8,1.4,0.1))+
      ylab(paste0("2 fold amino acid frequency increase \n effect on ", responseLabel))+
      xlab("")
    
    tpAA$aa <- as.character(tpAA$aa)
    tp <- merge(tp, tpAA, by="aa")
    tp <- tp[order(tp$codonEffects),]
    
    tp$codon <- factor(tp$codon, levels=tp$codon)
    tp$aa <- factor(tp$aa, levels=aaOrder)
    
    tp[,yMin:=logBase^(codonEffects-1.96*codonStdEr),][,yMax:=logBase^(codonEffects+1.96*codonStdEr)]
    tp[,codonEffects:=logBase^codonEffects,]
    
    p2 <- ggplot(tp, aes(x=codon, y=codonEffects, col=aa))+
      geom_point()+
      geom_errorbar(aes(ymin = yMin, ymax = yMax), size=0.3)+
      facet_wrap(~aa,  scales = "free_x",  nrow = 1)+
      scale_color_manual(values=rep(c("#999999","#E69F00"),10))+
      theme(panel.spacing = unit(0, "lines"), 
            axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust=0.5), 
            plot.title = element_text(size=12, hjust = 0.5), 
            legend.position = "None", 
            strip.text = element_text(size=15))+
      geom_hline(yintercept=1, color="darkgrey")+
      scale_y_continuous(trans='log2', breaks = seq(0.8,1.4,0.1) ,
                         labels = seq(0.8,1.4,0.1))+
      ylab(paste0("2 fold codon frequency increase \n effect on ", responseLabel))+
      xlab("")
    
    tp$codon <- as.character(tp$codon)
    tp$aa <- as.character(tp$aa)
    
    return(list(p1=p1, p2=p2, codonRes=tp, aaRes=tpAA))
  }
  
  plotAAAndCodonEffectForMultipleTissues <- function(codonBetas,
                                                     quantileNo=NULL,
                                                     aaBetas=NULL,
                                                     responseLabel,
                                                     mainTitle="",
                                                     logBase=10){
   
    if(!is.null(aaBetas)){
      aaBetas <- logBase^aaBetas
      colnames(aaBetas) <-sapply(colnames(aaBetas), 
                                 function(x){strsplit(x,"_")[[1]][2]})
      aaMedianBetas <- apply(aaBetas,2,median)
      aaMedianBetas <-aaMedianBetas[order(aaMedianBetas)]
      aaOrder <- names(aaMedianBetas)
    }else{
      aaMedianBetas <- NULL
    }
    
    codonBetas <- logBase^codonBetas
    colnames(codonBetas) <- sapply(colnames(codonBetas), 
                                   function(x){strsplit(x,"_")[[1]][2]})
    codonMedianBetas <- apply(codonBetas,2,median)
    codonMedianBetas <-codonMedianBetas[order(codonMedianBetas)]
    codonBetasOrder <- names(codonMedianBetas)
    
    codonBetasMelted <- melt(codonBetas)
    codonBetasMelted$aa <- genCode[as.character(codonBetasMelted$variable)]
    
    codonBetasMelted <- as.data.table(codonBetasMelted)
    
    if(!is.null(aaBetas)){
      codonBetasMelted$aa <- factor(codonBetasMelted$aa, levels=aaOrder)
    }else{
      codonBetasMelted[,medianCodonEffect := median(value),by=aa]
      codonBetasMelted <- codonBetasMelted[order(medianCodonEffect),]
      codonBetasMelted$aa <- factor(codonBetasMelted$aa, levels=unique(codonBetasMelted$aa))
    }
    
    codonBetasMelted$variable <- factor(codonBetasMelted$variable,
                                        levels=codonBetasOrder)
    
    
    codonBetasMelted <- codonBetasMelted[order(codonBetasMelted$aa, 
                                               codonBetasMelted$variable),]
    p1=NULL
    
    if(!is.null(aaBetas)){
      aaBetasMelted <- melt(aaBetas)
      aaBetasMelted$variable <- factor(aaBetasMelted$variable, levels = aaOrder)
      
      p1 <- ggplot(aaBetasMelted, aes(x=variable,y=value, fill=variable))+
              geom_boxplot(outlier.shape = NA, width=0.5, lwd=0.3) + 
              scale_fill_manual(values=rep(c("#999999","#E69F00"),10))+
              ylab(paste0("2 fold amino acid frequency increase \n effect on ", responseLabel))+
              xlab("")+
              theme(panel.spacing = unit(0, "lines"), 
                    axis.text.x = element_text(size=7, hjust = 1, vjust=0.5), 
                    plot.title = element_text(size=8, hjust = 0.5), 
                    legend.position = "None", 
                    strip.text = element_text(size=8))+
              geom_hline(yintercept=1, color="darkgrey")+
              scale_y_continuous(trans='log2', breaks = seq(0.2,2.0,0.1) , labels = seq(0.2,2,0.1))
    }
    
    p2 <- ggplot(codonBetasMelted, aes(x=variable,y=value, fill=aa))+
            geom_boxplot(outlier.shape = NA, width=0.5, lwd=0.3) + 
            facet_wrap(~aa,  scales = "free_x",  nrow = 1) +
            scale_fill_manual(values=rep(c("#999999","#E69F00"),10))+
            ylab(paste0("2 fold codon frequency increase \n effect on ", responseLabel))+
            xlab("")+
            theme(panel.spacing = unit(0, "lines"), 
                  axis.text.x = element_text(size=9, angle = 90, hjust = 1, vjust=0.5), 
                  plot.title = element_text(size=7, hjust = 0.5), 
                  legend.position = "None", 
                  strip.text = element_text(size=7))+
            geom_hline(yintercept=1, color="darkgrey")+
            scale_y_continuous(trans='log2', breaks = seq(0.8,1.4,0.1) , labels = seq(0.8,1.4,0.1))
    
      if(mainTitle != "")
      p2 <- p2 + ggtitle(mainTitle)
    
    return(list(p1=p1,
                p2=p2,
                aaMedianBetas=aaMedianBetas,
                codonMedianBetas=codonMedianBetas,
                aaBetas=aaBetas,
                codonBetas=codonBetas))
    
  }
  
  
  ### This function calculates the feature specific R^2 values per tissue
  ### and plots them as a box-plot .
  
  plotFeaturesExplainedVariance <- function(RNA_Protein_Combined, 
                                            isAAIncluded=F, 
                                            responseSubstring, 
                                            ylabStr,
                                            allFeatures,
                                            isTissues=T){
    print(responseSubstring)
    
    allFeatures$ALL <- unlist(allFeatures)
    
    if(isAAIncluded){
      allFeatures$aminoAcidUsage <- grep("AminoAcid_",
                                         colnames(RNA_Protein_Combined))
    }
    
    if(isTissues){
      expVarsOfFeatures <- getExpVarsOfFeatures(allFeatures, 
                                                RNA_Protein_Combined,
                                                responseNames=paste0(TissueNames,responseSubstring)) 
      rownames(expVarsOfFeatures) <- TissueNames
    }else{
      expVarsOfFeatures <- getExpVarsOfFeatures(allFeatures, 
                                                RNA_Protein_Combined,
                                                responseNames=responseSubstring)
    }
    
    mns <- colMedians(expVarsOfFeatures, na.rm=TRUE)
    expVarsOfFeatures <- expVarsOfFeatures[,order(mns)]
    
    
   return(list(expVarsOfFeatures=expVarsOfFeatures))
    
  }
  
  plotCors <- function(inDF, x, y, Xerr=NULL, Yerr=NULL, xLab, yLab, labelStr){
    
    b <- data.frame(x=inDF[,x],
                    y=inDF[,y],
                    aa=inDF[,labelStr])
    
    if(!is.null(Xerr)){
      t <- inDF[,Xerr]
      colnames(t) <- c("xmin", "xmax")
      b <- cbind(b, t)
    }
    
    if(!is.null(Yerr)){  
      t <- inDF[,Yerr]
      colnames(t) <- c("ymin", "ymax")
      b <- cbind(b, t)
    }
  
    k <- cor.test(b$x,
                  b$y,
                  method="spearman")
    print(k)
    
    p <- ggplot(b, aes(x=x, y=y))+
            geom_point(size=0.5, col="red")+
            geom_text(aes(label=aa),
                      hjust=+0.6, 
                      vjust=0,
                      size=2,
                      check_overlap = F)+
            #geom_hline(yintercept=1, color="darkgrey")+
            #geom_vline(xintercept=1, color="darkgrey")+
            #geom_abline(intercept = 0, slope = 1)+
            xlab(xLab)+
            ylab(yLab)+
            ggtitle(paste0(" rho = ",round(k$estimate, digits = 2), " pVal = ",round(k$p.val, digits = 5)))
  
    if(!is.null(Xerr)){
      p <- p+geom_errorbarh(data = b, aes(xmin = xmin,xmax = xmax))
    }
    
    if(!is.null(Yerr)){
      p <- p+geom_errorbar(data = b, aes(ymin = ymin,ymax = ymax))
    }
    
    return(p)
    
  }
  
  
  plotGeneLevels <- function(geneIdentifier,RNA_Protein_Combined){
    mRNALev<-unlist(RNA_Protein_Combined[geneIdentifier, paste0(TissueNames,"_medianExon")])
    proteinLev <- unlist(RNA_Protein_Combined[geneIdentifier, paste0(TissueNames,"_iBAQ")])
    
    tp <- data.frame(Levels= c(mRNALev-median(mRNALev, na.rm = T),
                               proteinLev-median(proteinLev, na.rm = T)),
                     Type= c(rep("mRNA", length(TissueNames)), rep("Protein", length(TissueNames))),
                     Tissues=factor(c(TissueNames, TissueNames ), levels = TissueNames))
    
    tp$Levels <- 10^tp$Levels
    
    print(ggplot(tp, aes(x=Tissues, y=Levels, group= Type, color=Type) )+
      geom_point()+
      geom_line()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      ggtitle(paste0(geneIdentifier,
                     " median centered mRNA and protein levels \n, mRNA S: ",
                     round(RNA_Protein_Combined[geneIdentifier,"mRNATisSScore"], digits = 3),
                     ", prot. S: ",round(RNA_Protein_Combined[geneIdentifier,"proteinTisSScore"], digits = 0)))+
      ylab("Fold change around median level across tissues")+
      scale_y_continuous(trans='log10',  
                         breaks = c(seq(0,3,0.25), seq(3,8,1)), 
                         labels = c(seq(0,3,0.25), seq(3,8,1))))
  
  }
  
  plotDAVIDRes <- function(k){
    if(!is.null(k$plotList$GOTERM_BP_ALL))
      print(k$plotList$GOTERM_BP_ALL)
    
    if(!is.null(k$plotList$GOTERM_CC_ALL))
      print(k$plotList$GOTERM_CC_ALL)
    
    if(!is.null(k$plotList$GOTERM_MF_ALL))
      print(k$plotList$GOTERM_MF_ALL)
    
  }
  
  plotResultFigures <- function(expVarsOfFeatures, mainStr){
    
    ev <- 100*expVarsOfFeatures
    colnames(ev) <- featureCols[colnames(ev), "Feature"]
    
    print(pheatmap(ev, treeheight_row = 0, treeheight_col = 0, fontsize = 7,legend = F))

    ev2 <- as.matrix(apply(ev,2, function(x){x-median(x, na.rm = T)}))
    
    pheatmap(t(ev2), treeheight_row = 0, treeheight_col = 0,
                   main = paste0("Row median centered explained variance in PTR \n by feature groups"), fontsize = 7)
    

    medDP <- apply(expVarsOfFeatures,2, function(x){median(x, na.rm = T)})
    expVarsOfFeatures <- expVarsOfFeatures[,order(medDP)]
    
    colnames(expVarsOfFeatures) <- featureCols[colnames(expVarsOfFeatures), "Feature"]
    fOrder <- colnames(expVarsOfFeatures)
    expVarsOfFeaturesMelted <- melt(expVarsOfFeatures)
    expVarsOfFeaturesMelted$Var2 <- as.character(expVarsOfFeaturesMelted$Var2)
    expVarsOfFeaturesMelted$Var2 <- factor(expVarsOfFeaturesMelted$Var2, levels = fOrder)
    
    p <- ggplot(expVarsOfFeaturesMelted, aes(x=Var2, y=100*value)) +
            ylab(mainStr)+ 
            xlab("")+
            geom_boxplot(lwd=0.3, outlier.size=0.3)+
            coord_flip()+geom_hline(yintercept = 0, col="red")
    
    f <- data.frame(medDP, featureNames=names(medDP))
    
    return(list(f=f,p=p))
  }
  
  
  plotCorrelations <- function(inputDF, xStr, yStr, xLab, yLab, includeLabels=T){
    
    v=100*c(0.001, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35)
    
    b <- data.frame(x=inputDF[,xStr], 
                    y=inputDF[,yStr],
                    featureNames=inputDF$featureNames)
    b[b<0] <- 0
    b[,c("x","y")] <- sqrt(b[,c("x","y")])
    
    
    p <- ggplot(b, aes(x=x, y=y))+
            geom_point(size=0.5, col="red")+
            geom_abline(intercept = 0,slope = 1)+
            xlab(xLab)+
            ylab(yLab)+
            scale_y_continuous(breaks = NULL)+ 
            scale_y_continuous(breaks = sqrt(v), labels = v)+ 
            scale_x_continuous(breaks = NULL)+ 
            scale_x_continuous(breaks = sqrt(v), labels = v) 
    
    if(includeLabels){
      p <- p+ geom_text(aes(label=featureNames),
                        hjust=+0.1,
                        vjust=+0.1,
                        size=3,
                        angle=45,
                        check_overlap = F)
    }
    
    print(p)
    return(p)
  }
  
  
