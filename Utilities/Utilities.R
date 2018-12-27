  require(caret)
  "%ni%" = Negate( "%in%" )
  
  genCode <- GENETIC_CODE
  names(genCode) <- as.vector(unlist(lapply(names(genCode),
                                            function(x){gsub("T","U",x)})))
  
  getPrecisionRecallAUC <- function(perf1){
    auc <- c()
    for(elem in 1:length(perf1@x.values)){
      f <- approxfun(data.frame(unlist(perf1@x.values[[elem]]) ,
                                unlist(perf1@y.values[[elem]]) )) 
      auc <- c(auc, integrate(f, 0, 1, subdivisions = 2000)$value)
    }
    
    return(auc)
  }
  ### Coefficient of variation for log-transformed x values.
  logCV <- function(x){
    return(sqrt(exp(sd(x, na.rm = T)^2)-1))
  }
  
  myPrcomp <- function(m,n=2){
    m <- scale(m, center = T, scale = T)
    ppk<-propack.svd(as.matrix(m),neig=n)
    pca<-t(ppk$d*t(ppk$u))
    return(pca)
  }
  
  lappend <- function (lst, ...){
    lst <- c(lst, list(...))
    return(lst)
  } 
  
  cs1 = ggseqlogo::make_col_scheme(chars=c('A', 'C', 'G','U'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), 
                        cols=c('LawnGreen', 'Blue', 'Orange', 'Red'))
  
  
  myGGSeqLogo <- function(sequences, xlabStr, motifLength, lrOffset, plotYAxis=T, methodS="bits"){
    d <- ggseqlogo( sequences , method = methodS, seq_type='rna', col_scheme=cs1) + 
      theme(legend.position = "None", 
            #axis.text.x = element_text(face = "plain", family = fontFamily, size=7), 
            axis.text.x = element_blank(),
            axis.text.y = element_text( face = "plain", family = fontFamily, size=7),
            plot.margin = unit(c(0.2, -0.1, 0, -0.3), "cm"))
    
    # if(motifLength == 3){
    #   d <- d+ scale_x_continuous(breaks=1:(motifLength+ 2*lrOffset), 
    #                              labels=c(paste0("-",6:1),1:motifLength, paste0("+",1:6)))
    # }else{
    #   d <- d+ scale_x_continuous(breaks=1:(motifLength+ 2*lrOffset), 
    #                      labels=c("-5","","-3","","-1",1:motifLength, "+1","","+3","","+5"))
    # }
    #   
    
    if(plotYAxis){
      return(d+labs(x = xlabStr, y=""))
    }else{
      return(d+ theme(axis.text.y = element_blank())+labs(x = xlabStr, y="")) 
    }
  }
  
  
  computeR2 <- function(response, prediction){
    testMSE = mean((response - prediction)^2, na.rm = T)
    R2 <- (1- testMSE/var(response, na.rm = T))
    return(R2)
  }
  
  
  getExpVarsOfFeatures <- function(allFeatures, RNA_Protein_Combined, responseNames){
    
    expVarsOfFeatures <- list()
    
    for(elem in names(allFeatures)){
      fixedEffectCols <- colnames(RNA_Protein_Combined)[allFeatures[[elem]]]
      X <- as.data.frame(RNA_Protein_Combined[,fixedEffectCols])
      colnames(X) <- fixedEffectCols
      
      expVF <- c()
      for(tissue in responseNames){
        
        if(elem == "codonUsage"){
           aaCols <- grep("AminoAcid_", colnames(RNA_Protein_Combined))
           
           Y <- getYPrime(inputDF=RNA_Protein_Combined, 
                              responseCol=tissue,
                              allFeatureNames=colnames(RNA_Protein_Combined[,aaCols]))
          
        }else if(elem == "codonInitial"){
          
          aaInitialCols <- grep("AminoAcidInitial_", colnames(RNA_Protein_Combined))
          
          Y <- getYPrime(inputDF=RNA_Protein_Combined, 
                         responseCol=tissue,
                         allFeatureNames=colnames(RNA_Protein_Combined[,aaInitialCols]))
        }else{
          Y <- RNA_Protein_Combined[,tissue]
        }
        
        
        res_lm <- compute10FCV_mixedEffectModel(X, 
                                                y = Y, 
                                                rnaCol=NULL, 
                                                fixedEffectCols, 
                                                randomEffectCol=NULL, 
                                                isMixedEffect=F, 
                                                flds = NULL)
        
        expVF <- c(expVF, res_lm$cvR2)
      }
      
      expVarsOfFeatures <- lappend(expVarsOfFeatures, expVF)
    }
    
    names(expVarsOfFeatures) <- names(allFeatures)
    expVarsOfFeatures <- do.call(cbind, expVarsOfFeatures)
    rownames(expVarsOfFeatures) <- responseNames
    return(expVarsOfFeatures)
  }
  
  compute10FCV_mixedEffectModel <- function(X, y, 
                                            rnaCol=NULL,
                                            fixedEffectCols,
                                            randomEffectCol=NULL,
                                            isMixedEffect=F,
                                            flds=NULL,
                                            foldNo=10){
    
    
    X <- as.data.frame(X[!is.na(y),])
    
    if(!is.null(rnaCol)){
      rnaCol <- rnaCol[!is.na(y)]
    }
    
    y <- y[!is.na(y)]
    
    if(length(X) == 1){
      y <- y[!is.na(X)]
      
      if(!is.null(rnaCol)){
        rnaCol <- rnaCol[!is.na(X)]
      }
      
      X <- as.data.frame(X[!is.na(X),])
    }
    
    
    if(is.null(flds)){
      flds <- createFolds(1:nrow(X),
                          k = foldNo,
                          list = TRUE,
                          returnTrain = FALSE)
    }
    
    yPred <- c()
    yVal <- c()
    rnaVal <- c()
    
    for(elem in flds){
      Test <- as.data.frame(X[elem,])
      Train <- as.data.frame(X[-elem,])
      
      if(length(X)== 1){
        colnames(Test) <- fixedEffectCols
        colnames(Train) <- fixedEffectCols
      }
      
      y.test <- y[elem]
      y.train <- y[-elem]
      
      Train$y <- y.train
      Test$y <- y.test
      
      Train <- Train[!is.na(Train$y),]
      Test <- Test[!is.na(Test$y),]
      
      if(isMixedEffect){
        fit <- lme4::lmer(paste0(" y ~ ", paste(fixedEffectCols, collapse = "+"), "+ (1 | ",randomEffectCol,")"),
                          data = Train)
      }else{
        fit <- lm(paste0(" y ~ ", paste(fixedEffectCols, collapse = "+")),
                  data = Train,
                  na.action = na.exclude)
      }
      
      predVal <- stats::predict(fit, newdata=Test)
      
      yVal <- c(yVal, Test$y)
      yPred <- c(yPred, predVal)
      
      if(!is.null(rnaCol)){
        rnaVal <- c(rnaVal,rnaCol[elem])
      }
      
    }
    
    cvR2 <- computeR2(yVal, yPred)
    MAE <- median(abs(yVal-yPred), na.rm = T)
    Cor <- cor(yVal, yPred, method = "spearman")
    return(list(yVal=yVal, yPred=yPred, flds=flds, rnaVal=rnaVal, cvR2=cvR2, MAE=MAE, Cor=Cor))
  }
  
  getNumberOfGenesHavingTheMotifs <- function(motifs, inputDF){
    inputDFTemp <- inputDF[,motifs]
    inputDFTemp[inputDFTemp>1] <-1
    inputDFTemp <- melt(colSums(inputDFTemp))
    inputDFTemp$Motif <- rownames(inputDFTemp)
    return(inputDFTemp)
  }
  ###########################
  ######### This function returns the indices of the columns which has
  ######### the elements of inputStr in its column name.
  getColumnsOfSubstr <- function(X, inputStr , exactMatch=F){
    gg <- c()
    
    for(elem in inputStr){
      if(exactMatch){
        fIndex <- grep(paste0("\\b",elem,"\\b"), colnames(X))
        gg <- c(gg, fIndex)
      }else{
        fIndex <- grep(elem, colnames(X))
        gg <- c(gg, fIndex)
      }
      
    }
    return(gg)
  }
  
  ###############################
  ###### This function returns the column indices in the given data frame
  ###### for the features that are used in our model.
  getAllFeatureColumns <- function(inputDF){
    
    elementCols = list(
      codonUsage=getColumnsOfSubstr(inputDF, "Codon_"),
      kozakSequence=getColumnsOfSubstr(inputDF, "Kozak_"), 
      stopCodonContext=getColumnsOfSubstr(inputDF, "stop_"),
      CDSLength=getColumnsOfSubstr(inputDF, "CDS_length"),
      UTR3Length=getColumnsOfSubstr(inputDF, "UTR3_length"),
      UTR5Length=getColumnsOfSubstr(inputDF, "UTR5_length"),
      CDSGCPercent= getColumnsOfSubstr(inputDF, "CDSGCPercent"), 
      UTR3GCPercent=getColumnsOfSubstr(inputDF, "UTR3GCPercent"), 
      UTR5GCPercent=getColumnsOfSubstr(inputDF, "UTR5GCPercent"),
      UTR3_Motifs=getColumnsOfSubstr(inputDF, "Motif_UTR3"),
      UTR5_ATG= getColumnsOfSubstr(inputDF, "Motif_UTR5_AUG"), 
      UTR5_Motifs=setdiff(getColumnsOfSubstr(inputDF, "Motif_UTR5"), getColumnsOfSubstr(inputDF, "Motif_UTR5_AUG")),
      CDS_AA_Motifs=getColumnsOfSubstr(inputDF, "Motif_AA_"),
      linearProteinMotifs=getColumnsOfSubstr(inputDF, "Motif_LinAA_"), 
      UTR5FE=getColumnsOfSubstr(inputDF, "UTR5FE_"),
      codonPairBias=getColumnsOfSubstr(inputDF, "CPair_"),
      pestMotifs=getColumnsOfSubstr(inputDF, "PESTNo"),
      proteinIP=getColumnsOfSubstr(inputDF, "averagePISO"),
      ntermHydro=getColumnsOfSubstr(inputDF, "ntermHydro")
    )
    
    return(elementCols)
  }
  
  getAllFeatureColumnsWithProteinPTM <- function(inputDF){
    tp <- getAllFeatureColumns(inputDF)
    
    tp$miRNAs=getColumnsOfSubstr(inputDF, "sgn_MiRNA_")
    tp$RBPs=getColumnsOfSubstr(inputDF, "RBP_")
    tp$proteinComplexFormation=getColumnsOfSubstr(inputDF, "ProteinComplex")
    tp$m6ARNAModification=getColumnsOfSubstr(inputDF, "RNAm6A")
    tp$phospSites=getColumnsOfSubstr(inputDF, "numOfPhospSites")
    tp$ubiqSites=getColumnsOfSubstr(inputDF, "numOfUbiSites")
    tp$aceSites=getColumnsOfSubstr(inputDF, "numOfAceSites")
    tp$sumoSites=getColumnsOfSubstr(inputDF, "numOfSumoSites")
    tp$methSites=getColumnsOfSubstr(inputDF, "numOfMethySites")
    
    
    return(tp)
  } 
  
  ################################
  ####### This features returns the residuals for linear fit: response ~ [allFeatureNames- feature]
  getYPrime <- function(inputDF, feature="None", responseCol, allFeatureNames){
    
    otherFeatures <- allFeatureNames[allFeatureNames %ni% feature]
    
    inputDFTemp <- as.data.frame(inputDF[,otherFeatures])
    inputDFTemp$response <- inputDF[,responseCol]
    
    lmout <- lm(response~., data = inputDFTemp, na.action=na.exclude)
    yPrime = residuals(lmout)
    
    return(yPrime=yPrime)
  }
  
  
  correctForGivenFeatures <- function(inputDF, feature="None", responseStr, allFeatureNames){
    responseCols = c()
    for(elem in responseStr){
      responseCols <- c(responseCols, getColumnsOfSubstr(X=inputDF, inputStr=elem))
    }
    
    for(elem in responseCols){
      yP <- getYPrime(inputDF, feature, elem, allFeatureNames)
      inputDF[,elem] <- yP
    }
    
    return(inputDF)
  }
  
  #########################################################
  getOLSR2 <- function(X, y, featureNames){
  
    ###############################################
    #### Estimate the parameters with the fit to the whole data
    tDF3 <- data.frame(cbind(y,X))
    colnames(tDF3) <- c("y", featureNames)
    fitall = lm(y ~ ., data = tDF3, na.action = na.exclude)
    lmsum <- summary(fitall)
    lmsumDT <- data.frame(Estimate=coef(fitall), std.error=NA, tval=NA, pval=NA)
    lmsumDT[rownames(coef(lmsum)),c("std.error","tval","pval")]<- coef(lmsum)[,2:4]
    lmsumDT[,"covariates"] <- rownames(lmsumDT)
    lmsumDT <- as.data.table(lmsumDT)
    
    return(list(lmsumDT=lmsumDT))
  }
  
  ################################################################
  #### This function returns the outputs of a number of OLS fits where
  #### @ inputDF: data frame which contains the feature set and the response variables
  #### @ allFeatureNames: column names of the covariates
  #### @ responseStr: substring of the column names of the response variable(s)
  
  selectBetasAccrossTissues <- function( inputDF, allFeatureNames, responseStr){
    
    ######### Get the responses
    responseCols = c()
    for(elem in responseStr){
      responseCols <- c(responseCols,
                        getColumnsOfSubstr(X=inputDF,
                                           inputStr=elem))
    }
    responseNames <- colnames(inputDF[,responseCols])
    ######## Get the features we are interested in
   
    betas= list()
    pvalues=list()
    covariates = c()
    stdError <- list()
    intercepts <- c()
    
    for(elem in responseCols){
      
      results <- getOLSR2(X=inputDF[,allFeatureNames],
                          y=inputDF[,elem],
                          featureNames = allFeatureNames)
      
      ## Select only the covariates given as arguments
      lmsumDT <- results$lmsumDT[covariates %in% allFeatureNames,]
      intercepts <- c(intercepts,results$lmsumDT[1,Estimate])
      betas <- lappend(betas, lmsumDT[,Estimate,])
      pvalues <- lappend(pvalues, lmsumDT[,pval,])
      
      stdError <- lappend(stdError, lmsumDT[,std.error,])
      covariates <- lmsumDT[,covariates,]
    }
    
    betas <- as.data.frame(do.call(rbind, betas),
                           stringsAsFactor=F)
    
    rownames(betas) <- responseNames
    colnames(betas) <- covariates
    
    pvalues <- as.data.frame(do.call(rbind, pvalues),
                             stringsAsFactor=F)
    pAdjusted <- data.frame(matrix(p.adjust(as.matrix(pvalues),
                                            method = "BH"),
                                   nrow = nrow(pvalues)),
                            stringsAsFactors = F)
    rownames(pAdjusted) <- responseNames
    colnames(pAdjusted) <- covariates
    
    stdError <- as.data.frame(do.call(rbind, stdError),
                              stringsAsFactor=F)
    rownames(stdError) <- responseNames
    colnames(stdError) <- covariates
    
    names(intercepts) <- responseNames
    return(list(betas=betas,
                pvalues=pAdjusted,
                stdError=stdError,
                intercepts=intercepts))
  }
  
  
  getBetasOfFeature <- function(allFeaturesFit, featureSubStr){
    betas <- dplyr::select(allFeaturesFit$betas, matches(featureSubStr))
    pvalues <- dplyr::select(allFeaturesFit$pvalues, matches(featureSubStr))
    stdError <- dplyr::select(allFeaturesFit$stdError, matches(featureSubStr))
    return(list(betas=betas, pvalues=pvalues, stdError=stdError))
  }
  ############################################################################
  #################### This function creates the design matrix our of 
  #################### RNA string elements stored in RNAStringSetInput.
  createDesignMatrixWithGivenSeqs = function(RNAStringSetInput, lroffset=6, msize=3)
  {
    nulllength <- 2*lroffset + msize
    nullgen = rawToChar(as.raw(rep(78,nulllength)))
    res2 <- lapply(RNAStringSetInput , function(x){c(x, nullgen)})
    
    count_base <-  function(x){
      return(data.frame(len=length(x)-1,t(as.vector(consensusMatrix(RNAStringSet(x))[1:4,]))))
    }
    
    temp = foreach(i = res2) %dopar% count_base(i)
    temp = rbindlist(temp)
    
    if(lroffset > 0){
      posnames=c(sapply((lroffset):1,function(x){paste("P_",x,sep='')}),sapply(1:(msize),
                                                                               function(x){paste("P",x,sep='')}), sapply(1:(lroffset),function(x){paste("PJ",x,sep='')}))
    }else{
      posnames=sapply(1:msize,function(x){paste("P",x,sep='')})
    }
    
    
    all_pos = strcat(posnames,c(".A",".C",".G",".U"))
    
    colnames(temp) = c("len",all_pos)
    
    dmI = list(dmatrix.c = temp, seq.c = res2, posnames = posnames, all_pos = all_pos,num_pos = nulllength,nullgen = nullgen)
    return(dmI)
  }
  
  createMatchSeqsForMotif <- function(RNAStringSetInput, lroffset, mismatchNo, motif)
  {
    shortStrOfset=lroffset+2
    RNAStringSetInput_short <- RNAStringSetInput
    RNAStringSetInput_long <- RNAStringSetInput
    msize = nchar(motif)
    nulllength = lroffset*2  + msize
    
    RNAStringSetInput_short <- RNAStringSet(unlist(lapply(RNAStringSetInput_short, 
                                                          function(x){
                                                            return(RNAString(substr(as.character(x),
                                                                                    shortStrOfset,
                                                                                    (nchar(x)-(lroffset+1) ) ) ) )})))
    
    res = vmatchPattern(motif, RNAStringSetInput_short, max.mismatch = mismatchNo)
    
    nullgen = rawToChar(as.raw(rep(78,nulllength)))
    
    ass_func <- function(i) 
    {
      temp = Views(RNAStringSetInput_long[[i]], res[[i]])
      end(temp) = end(temp) +  shortStrOfset-1 + lroffset
      start(temp) = start(temp) + shortStrOfset-1 - lroffset
      return(c(as.character(temp),nullgen))	
    }
    
    ## Parellel execution of ass_func for each sequence
    res2 = foreach(i = seq_along(RNAStringSetInput_short)) %dopar% ass_func(i)
    return(res2)
    
  }  
  
  createDesignMatrixForMotif = function(RNAStringSetInput, lroffset, mismatchNo, motif)
  {
    shortStrOfset=lroffset+2
    RNAStringSetInput_short <- RNAStringSetInput
    RNAStringSetInput_long <- RNAStringSetInput
    msize = nchar(motif)
    nulllength = lroffset*2  + msize
    
    RNAStringSetInput_short <- RNAStringSet(unlist(lapply(RNAStringSetInput_short, function(x){return(RNAString(substr(as.character(x),shortStrOfset,(nchar(x)-(lroffset+1) ) ) ) )})))
    
    res = vmatchPattern(motif, RNAStringSetInput_short, max.mismatch = mismatchNo)
    
    nullgen = rawToChar(as.raw(rep(78,nulllength)))
    
    ass_func <- function(i) 
    {
      temp = Views(RNAStringSetInput_long[[i]], res[[i]])
      end(temp) = end(temp) +  shortStrOfset-1 + lroffset
      start(temp) = start(temp) + shortStrOfset-1 - lroffset
      return(c(as.character(temp),nullgen))	
    }
    
    ## Parellel execution of ass_func for each sequence
    res2 = foreach(i = seq_along(RNAStringSetInput_short)) %dopar% ass_func(i)
    
    
    #make design matrix
    #count bases
    count_base <-  function(x){
      return(data.frame(len=length(x)-1,t(as.vector(consensusMatrix(RNAStringSet(x))[1:4,]))))
    }
    
    temp = foreach(i = res2) %dopar% count_base(i)
    temp = rbindlist(temp)
    
    if(lroffset > 0){
      posnames=c(sapply((lroffset):1,function(x){paste("P_",x,sep='')}),sapply(1:(msize),function(x){paste("P",x,sep='')}), sapply(1:(lroffset),function(x){paste("PJ",x,sep='')}))
    }else{
      posnames=sapply(1:msize,function(x){paste("P",x,sep='')})
    }
    
    all_pos = strcat(posnames,c(".A",".C",".G",".U"))
    
    colnames(temp) = c("len",all_pos)
    
    dmI = list(dmatrix.c = temp, seq.c = res2, posnames = posnames, all_pos = all_pos,num_pos = nulllength,nullgen = nullgen)
    return(dmI)
  }
  #############################################################################
  ################ This function fits the linear model on motif base pair variances. 
  
  fitLmModel = function(dmI, correlate)
  {
    motif = RNAStringSet(unlist(dmI$seq.c))
    background = apply(consensusMatrix(motif,as.prob=TRUE)[1:4,],2,which.max)
    ind_rm = sapply(1:length(background),function(x){(x-1)*4+background[x]})
    lm_pos = dmI$all_pos[-ind_rm]
    
    formula.txt = paste("correlate","~len+",paste(lm_pos,collapse="+"),sep='')
    lmFrame <- as.data.frame(dmI$dmatrix.c)
    lmFrame[["correlate"]] <- correlate
    lmFrame <- lmFrame[!is.na(lmFrame$correlate),]
    estimate=coef((lm(as.formula(formula.txt),data=lmFrame)))
    pvals=estimate
    pvals[rownames(coef(summary(lm(as.formula(formula.txt),data=lmFrame))))] <- coef(summary(lm(as.formula(formula.txt),data=lmFrame)))[,4]
    
    
    lmfI = list(lm_out_c = data.frame(estimate=estimate, pvals=pvals) ,lm_pos=lm_pos,ind_rm=ind_rm)
    return(lmfI)
    
  }
  ############################################################
  #####
  density_kmer_ggplot <- function(align.right, kmers,colors,
                                  dat,seqcol,responsecol,window,max.length,
                                  seqlengthcol, xlab, ylab="Motif density"){
    require(data.table)
    maxL = max.length-50
    k_n_list <- density_kmer(kmers,
                             colors,
                             dat,
                             seqcol,
                             responsecol,
                             genenamecol='transID',
                             seqlengthcol ,max.length ,
                             align.right, window)
    dat <- rbindlist(k_n_list, idcol='motif')
    dat[ ,density := window*(k/n/theo.w)]
    dat[ ,sd := sqrt(density*(1-density)/n)]
    dat <- dat[,.SD[1:maxL],by=motif]
    #dat[ ,position := rev(position)]
    density_plot <- ggplot(dat, aes(plot_position, density, ymax=density+sd, ymin=density-sd, color=color, fill=color)) + 
      geom_smooth(stat='identity', alpha=0.2)
    if (!isTRUE(all.equal(dat$motif, dat$color))){
      density_plot <- density_plot + scale_fill_identity() + scale_color_identity()
    }
    breaks <- seq(0,max(dat$position),50)
    
    if (align.right){
      density_plot <- density_plot + scale_x_continuous(breaks=breaks,labels=-rev(breaks))
    }
    
    density_plot <- density_plot + 
      labs(x=xlab, y=ylab) +
      theme_minimal()+
      theme(axis.title=element_text(size=1.2, face='bold'))
    return(density_plot)
  }
  ############################################################
  #######
  density_kmer <- function(kmers,colors=NULL,dat, seqcol, responsecol, genenamecol='transID',
                           seqlengthcol ,max.length ,align.right, window){
    # Check colors
    if (is.null(colors)){
      colors <- kmers
    }else{
      if (length(kmers)!=length(colors)){
        stop('Length of kmers not equal to length of color vector.')
      }
    }
    # kill factor
    dat <- data.frame(lapply(dat, as.character), stringsAsFactors=FALSE)
    rownames(dat) = dat[,genenamecol]
    # Compute statistics
    k_n_list = list()
    for (k in seq(length(kmers))){
      motif <- kmers[k]
      ## Create data frame with motif, gene, rates, UTR.length, dist
      # Find motif positions at each gene
      positions = start(vmatchPattern(motif, RNAStringSet(dat[,seqcol]), fixed=FALSE))
      names(positions) = rownames(dat)
      # Find the index of genes which have target motif
      indx =c()
      for(i in 1:length(positions)){
        if (length(positions[[i]]) != 0){
          indx = c(indx,i)
        }
      }
      
      tab = t(as.data.frame(sapply(rownames(dat)[indx],
                                   function(i) sapply(positions[[i]],
                                                      function(j) c(i,motif, dat[i,seqlengthcol],
                                                                    as.numeric(j)+nchar(motif)/2)
                                   )
      )
      )
      )
      
      
      
      colnames(tab) = c('geneID','motif','length','dist')
      
      tab = as.data.frame(tab,stringsAsFactors=F)
      
      if(align.right){
        tab$dist = as.numeric(tab$length)-as.numeric(tab$dist)
      }
      
      ## Compute number of motifs at each positon and number of sequence coverage
      # normalization vector
      n=round(nchar(motif)/2,0)
      
      theo = c(rep(0,n),rep(1,max.length-2*n),rep(0,n)) # Theoretical possible kmer numbers at each position
      
      theo.w = c()  # Theoretical kmer numbers at each window
      
      for(i in 1:(max.length)){theo.w = c(theo.w, sum( theo[ max((i-(window/2)+1),1):(min((i+window/2),max.length)) ] ))}
      
      # k: number of motifs, n: sequence coverage
      k_n = t(sapply( # number of motifs in the window
        1:max.length,
        function(i){
          c(sum(as.numeric(tab$dist) >= (i-window/2) & as.numeric(tab$dist) <= i+(window/2)), sum(as.numeric(tab$length) >= i+(window/2))) 
        }
      ))
      
      k_n = data.table(k_n,theo.w,1:length(theo.w),colors[k]) # the columns are: oberseved # of motifs in the window, number of genes that UTR/promter longer than the extend of the window, theoretical number of possible motifs
      colnames(k_n) <- c('k','n','theo.w','position','color')
      if (align.right){
        k_n$plot_position <- rev(k_n$position)
      }else{
        k_n$plot_position <- k_n$position
      }
      k_n_list[[motif]] = k_n
    }
    return(k_n_list)
  }
  
  
  addMotifsToDF <- function(DF, motifDB, consOrUpdated="updatedMotif", allCount=F){
    
    for(i in 1:nrow(motifDB)){
      
      counts <- vcountPattern(as.character(motifDB[i,consOrUpdated]),
                                (RNAStringSet(DF[,paste0(motifDB[i,"motifType"],"Seq")])), 
                                fixed=F,
                                max.mismatch=0)
      
      # counts <- sapply( DF[,paste0(motifDB[i,"motifType"],"Seq")],
      #                     function(x){countPWM(motifPWMs[[motifDB[i,consOrUpdated]]],
      #                                          subject=DNAString(x),
      #                                          min.score="90%")})
      
      if(!allCount & motifDB[i,"binary"]){
        counts[counts>1]<-1
      }
      
      DF[paste0("Motif_",motifDB[i,"motifType"],"_",motifDB[i,consOrUpdated])] <- counts
      
    }
    
    return(DF)
  }
  
  addMotifCounts <- function(RNA_Protein_Combined, 
                             motifList,
                             sequenceVector,
                             columnPrefix,
                             binary,
                             fixedP=T,
                             isAA=F,
                             maxMismatch=0){
    
    for(motif in motifList){
      
      if(isAA){
        counts = vcountPattern(motif,(AAStringSet(sequenceVector)),
                               fixed=fixedP,
                               max.mismatch=maxMismatch)
      }else{
        counts = vcountPattern(motif,(RNAStringSet(sequenceVector)),
                               fixed=fixedP,
                               max.mismatch=maxMismatch)
      }
      if(binary){
        counts[counts > 1] <- 1
      }
      RNA_Protein_Combined[[paste0(columnPrefix,motif)]] <- counts
    }
    
    return(RNA_Protein_Combined)
  }
  
  addCDSMotifsFrame <- function( dF, motif, motifColStr, regionStr){
    
    dF[[paste0(motif,"Pos")]] <- apply(dF, 1, function(x){
      k<- as.numeric(unlist(start(vmatchPattern(motif,RNAStringSet(as.character(x[regionStr])),
                                                max.mismatch = 0,
                                                fixed = T))))
      return(unlist(tail(unlist(k), x[[motifColStr]])))
    })
    
    dF[dF[[motifColStr]]==0,paste0(motif,"Pos")]<- -1
    dF[[paste0(motif,"Frame")]] <- lapply(dF[[paste0(motif,"Pos")]], function(x){ k<- as.numeric(x) %% 3; k[k==0]  <- 3; return(unlist(k))})
    dF[dF[[motifColStr]]==0,paste0(motif,"Frame")]<- -1
    return(data.frame(dF, stringsAsFactors = F))
  }
  
  getTopForallBranches <- function(subsetGeneNames, 
                                   genePopulation, 
                                   explanationStr="",
                                   geneIdentifiesStr="geneSymbol", 
                                   motif=""){
    nTerms=15
    
    print("###############################################################")
    print( paste0(" Biological Processes "))
    mgsaFit1 <- compute_mgsa_and_pv(study_set=as.vector(subsetGeneNames), 
                                    population=genePopulation,
                                    go_branch='bp',
                                    organism = "human",
                                    folderPath=GoDir)
    mgsaFit1 <- mgsaFit1[p.adjust< 0.01,]
    mgsaFit1 <- mgsaFit1[order(estimate, decreasing=T)]
    if(nrow(mgsaFit1) < 15){
      nTerms = nrow(mgsaFit1)
    }else{
      nTerms=15
    }
    if(nTerms != 0){
        q <- qplot(reorder(mgsaFit1[1:nTerms,]$term, 
                          mgsaFit1[1:nTerms,]$estimate),
                  mgsaFit1[1:nTerms,]$estimate,
                  ylab="Pathway Estimate",
                  xlab="",
                  main = "Biological process") + 
              coord_flip()+
              theme(axis.title = element_text(face="bold",size=rel(1)), 
                    axis.text = element_text(size=rel(1)), 
                    plot.title = element_text(face="bold",size=rel(1), hjust = 0.5))
        if(motif != ""){
          q <- q + ggtitle(paste0(" Biological processes GO terms for motif: ", motif))
        }
          
        print(q)   
        
    }
    
    print( paste0(" Molecular Function "))
    mgsaFit2 <- compute_mgsa_and_pv(study_set=as.vector(subsetGeneNames),
                                    population=as.vector(genePopulation),
                                    go_branch='mf',
                                    organism = "human",
                                    folderPath=GoDir)
    mgsaFit2 <- mgsaFit2[p.adjust< 0.01,]
    mgsaFit2 <- mgsaFit2[order(estimate, decreasing=T)]
    
    if(nrow(mgsaFit2) < 15){
      nTerms = nrow(mgsaFit2)
    }else{
      nTerms=15
    }
    if(nTerms !=0){
      q <- qplot(reorder(mgsaFit2[1:nTerms,]$term,
                          mgsaFit2[1:nTerms,]$estimate),
                  mgsaFit2[1:nTerms,]$estimate,
                  ylab="Pathway Estimate",
                  xlab="",
                  main="Molecular function") + 
              coord_flip()+
              theme(axis.title = element_text(face="bold",size=rel(1)),
                    axis.text = element_text(size=rel(1)),
                    plot.title = element_text(face="bold",size=rel(1), hjust = 0.5))
      if(motif != ""){
        q <- q + ggtitle(paste0("Molecular function GO terms for motif: ", motif))
      }
              
      print(q)
    }
    
    
    print( paste0(" Cellular Component "))
    mgsaFit3 <- compute_mgsa_and_pv(study_set=as.vector(subsetGeneNames),
                                    population=as.vector(genePopulation),
                                    go_branch='cc', organism = "human",
                                    folderPath=GoDir)
    mgsaFit3 <- mgsaFit3[p.adjust< 0.01,]
    
    mgsaFit3 <- mgsaFit3[order(estimate, decreasing=T)]
    
    if(nrow(mgsaFit3) < 15){
      nTerms = nrow(mgsaFit3)
    }else{
      nTerms=15
    }
    
    if(nTerms !=0){
      q <- qplot(reorder(mgsaFit3[1:nTerms,]$term,  
                        mgsaFit3[1:nTerms,]$estimate),
                mgsaFit3[1:nTerms,]$estimate, 
                ylab="Pathway Estimate", 
                xlab="", 
                main="Cellular component") + 
      coord_flip()+
      theme(axis.title = element_text(face="bold",size=rel(1)), 
            axis.text = element_text(size=rel(1)), 
            plot.title = element_text(face="bold",size=rel(1), 
                                      hjust = 0.5))
      if(motif != ""){
        q <-q + ggtitle(paste0("Cellular component GO terms for motif: ", motif))
      }
        
      print(q)
    }
    
    print("############################################################")
    return(list(mgsaFit1=mgsaFit1, mgsaFit2=mgsaFit2,mgsaFit3=mgsaFit3))
  }
  
  sumStat <- function(yValues, yPredVlues, includeCV=T){
    scor <- round(cor(yValues, yPredVlues, method = "spearman"), digits = 2)
    pcor <- round(cor(yValues, yPredVlues, method = "pearson"), digits = 2)
    mae <- round(median(abs(yValues - yPredVlues)), digits = 2)
    CVR2 <- round(computeR2(yValues, yPredVlues), digits = 2)
    
    if(includeCV){
      summary=paste0("MAE = ", mae, " ,Rs = ", scor, " ,Rp = ", pcor, " ,R^2 = ", CVR2)
    }else{
      summary=paste0( "Rs = ", scor, " ,Rp = ", pcor)
    }
    
    return(list(summary=summary, scor=scor, pcor=pcor, mae=mae, CVR2=CVR2))
  }
  
  plotHydroBoxplot <- function(Indv, positions, ylabStr){
    Indv <- melt(Indv)
    
    p <- ggplot(Indv, aes(x=variable, y=value))+
            geom_boxplot()+
            scale_x_discrete(breaks=positions, labels=positions+1)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            ylab(ylabStr)
    return(p)
  }
  
  getMovingAverage <- function(Indv, rloffSet){
    
    slidingMean15 <- apply(Indv,1, function(x){
      x <- as.numeric(x)
      
      k=c()
      for(i in (rloffSet+1):(ncol(Indv)-rloffSet)){
        k = c(k,mean(as.numeric(x)[(i-rloffSet):(i+rloffSet)], na.rm=T))
      }
      
      return(as.list(k))
    })
    
    slidingMean15 <- as.data.frame(do.call(rbind,slidingMean15))
    
    for(i in 1:ncol(slidingMean15)){
      slidingMean15[,i] <- unlist(slidingMean15[,i])
    }
    
    return(slidingMean15)
  }
  
  
  # getMotifSignificance<- function(dfTemp, motif, side){
  #   print(motif)
  #   
  #   dfTempK <- dfTemp[dfTemp$motif %in% c(paste0("SCRAMBLED_",motif), motif),]
  #   dfTempK$motif <- drop.levels(dfTempK$motif)
  #   dfTempK$motif <- factor(dfTempK$motif, levels = c(motif, paste0("SCRAMBLED_",motif)))
  #   dfTempK$TB <- paste0(dfTempK$motif, "_",dfTempK$TB)
  #   
  #   wTestResPValues <- c()
  #   fc <- c()
  #   for(elem in levels(dfTempK$Experiment)){
  #     
  #     dfTempK_expTemp <- dfTempK[dfTempK$Experiment == elem,]
  #     
  #     fit <- lm(dfTempK_expTemp$intensity ~ dfTempK_expTemp$time.point)
  #     dfTempK_expTemp$intensityRes <- (residuals(fit) + coef(fit)["(Intercept)"])
  #     dfTempK_expTemp <- as.data.table(dfTempK_expTemp)
  #     dfTempK_expTemp[,medianTBIntRes := median(intensityRes), by = TB]
  #     
  #     dfTempK_expTemp[ ,c("time.point","intensity","intensityRes", "GLUC_measurement", "SEAP_measurement")] <- NULL
  #     groupMedians <- as.data.frame(dfTempK_expTemp[,median(medianTBIntRes),by=motif])
  #     dfTempK_expTemp <- as.data.frame(unique(dfTempK_expTemp))
  #     
  #     levels(dfTempK_expTemp$motif)
  #     wTestRes <- wilcox.test(dfTempK_expTemp$medianTBIntRes ~ dfTempK_expTemp$motif, alternative=side)
  #     wTestResPValues <- c(wTestResPValues, wTestRes$p.value)
  #     
  #     fc <- c(fc, 10^(groupMedians[groupMedians$motif == motif, "V1"] - groupMedians[groupMedians$motif == paste0("SCRAMBLED_",motif), "V1"]))
  #   }
  #   
  #    print(wTestResPValues)
  #    
  #   return(list(pVal = wTestResPValues, fc=fc))
  # }
  
  generateCodonFrequency <- function(seqs){
    names(seqs) <- NULL
    
    codonOccList <- lapply(seqs, function(x)
    {trinucleotideFrequency(RNAString(x), step=3)
    })
    
    names(codonOccList) <- 1:length(codonOccList)
    codonFreqDFNoBias <- as.data.frame(do.call(rbind, codonOccList),
                                       stringsAsFactors=F)
    codonFreqDFNoBias[, c("UAA","UGA","UAG")] <- NULL
    
    codonFreqDFNoBias <- as.data.frame(as.matrix(codonFreqDFNoBias)+1)
    totalNo <- rowSums(codonFreqDFNoBias)
    
    codonFreqDFNoBias <- sweep(codonFreqDFNoBias, 1, totalNo, `/`)
    codonFreqDFNoBias <- log2(codonFreqDFNoBias)
    
    return(codonFreqDFNoBias)
  }
  
  generateAminoAcidFrequency <- function(seqs){
    names(seqs) <- NULL
    
    aminoacids <- as.vector(unique(GENETIC_CODE))
    aminoacids <- aminoacids[-5]
    
    aaOccu = function(proteinSeq){
      aaNo = c()
      m = sapply(proteinSeq, 
                 function(i) 
                   substring(i,seq(1,nchar(i),1),seq(1,nchar(i),1)))
      for(elem in aminoacids){
        aaNo= c(aaNo, sum(str_count(m, elem)))
      }
      return(aaNo)
    }
    
    aaOccList <- lapply(seqs, 
                        function(x){aaOccu(x)})
    
    aaFreqDF <- do.call(rbind, aaOccList)
    aaFreqDF <- (aaFreqDF+1)
    colnames(aaFreqDF) <- aminoacids
    totalNo <- rowSums(aaFreqDF)
    aaFreqDFFraction <- as.data.frame(sweep(aaFreqDF, 1, totalNo, `/`))
    aaFreqDFFraction <- log2(aaFreqDFFraction)
    
    return(aaFreqDFFraction)
    
  }
  
  getCodonEffectGCContentCor <- function(codonMedianBetas, yLab){
    
    Codon_GCContent <- unlist(round(sapply(names(codonMedianBetas),
                                     function(x){return(letterFrequency(RNAString(x), letters = "GC")) }) / 3, digits = 2))
    
    
    tp <- data.frame(codonEfect=codonMedianBetas,
                     codon=names(codonMedianBetas),
                     GCContent = Codon_GCContent)
    tp <- tp[order(tp$codonEfect),]
    tp$codon <- factor(tp$codon, levels = tp$codon)
    p2 <- ggplot(tp, aes(x=codon, y=codonEfect))+
      geom_bar(stat="identity")+coord_flip()+theme_bw()
    
    
    tp$GCContent <- as.factor(tp$GCContent)
    
    p <- ggplot(tp, aes(x=GCContent, y=codonEfect))+
      geom_boxplot()+
      theme_bw()+
      ylab(paste0("Codon effect ", yLab))+
      xlab(" Codon GC Content ")
    
    p <- add_pval(p, pairs = list(c(1,3)), log = F, fold_change = T)
    
    tp$aa <- genCode[as.character(tp$codon)]
    
    return(list(p=p, tp=tp, p2=p2))
    
  }
  
  getCodonAAVals <- function(RNA_Protein_Combined,
                             responseStr,
                             rLab,
                             includeAA=T,
                             allfeatures){
    
    codonFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                           allFeatureNames=colnames(RNA_Protein_Combined[,unlist(allfeatures)]),
                                           responseStr=responseStr)
    codon_TV <- getBetasOfFeature(allFeaturesFit = codonFit,
                                  featureSubStr="Codon_")
    
    if(includeAA){
      allfeaturesAA <- allfeatures
      allfeaturesAA$codonUsage <- NULL
      allfeaturesAA$aaUsage <- grep("AminoAcid_",colnames(RNA_Protein_Combined))
      allFeatureNamesAA <- colnames(RNA_Protein_Combined[,unlist(allfeaturesAA)])
      
      aaFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                          allFeatureNamesAA,
                                          responseStr=responseStr)
      
      aa_TV <- getBetasOfFeature(aaFit, featureSubStr="AminoAcid_")
      
      ptrEf <- plotAAAndCodonEffectForMultipleTissues(codonBetas=codon_TV$betas, 
                                                      aaBetas=aa_TV$betas,
                                                      responseLabel= rLab)
      
    }else{
      ptrEf <- plotAAAndCodonEffectForMultipleTissues(codonBetas=codon_TV$betas, 
                                                      aaBetas=NULL,
                                                      responseLabel= rLab)
      
    }
    
    return(ptrEf)
  }
  
  getMergedTable <- function(firstVals, secondVals, firstName, secondName, byName){
    
    x <- data.table(a = firstVals,
                    c=names(firstVals))
    
    y <- data.table(b = secondVals,
                    c=names(secondVals))
    
    t <- as.data.frame(merge(x,y, by="c"))
    colnames(t) <- c(byName,firstName, secondName)
    return(t)
  }
  
  
  calculateGerstbergerScore <- function(inputDF){
    N <- ncol(inputDF)
    inputDF <- 10^inputDF
    tpX <- t(apply(as.matrix(inputDF),1,function(x){ return(x/sum(x, na.rm = T))}))
    k <- apply(tpX,1,function(x){return(-sum(x*log2(x), na.rm = T))})
    S=log2(N)-k
    return(S)
  }
  
  correlateRBPEffectWithItsExpression <- function(responseStr, responseType, RNA_Protein_Combined){
    
    allFeatures <- getAllFeatureColumnsWithProteinPTM(RNA_Protein_Combined)
    rownames(RNA_Protein_Combined) <- RNA_Protein_Combined$GeneName
    rbpFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                         allFeatureNames=colnames(RNA_Protein_Combined[,unlist(allFeatures)]),
                                         responseStr=responseStr)
    rbp_tv <- getBetasOfFeature(rbpFit, featureSubStr="RBP_")
    
    betas <- 10^rbp_tv$betas
    colnames(betas) <- RBPNames
    
    contRBPs <- RBPNames[which(RBPNames %in% RNA_Protein_Combined$GeneName)]
    
    ExpressionLevels <- RNA_Protein_Combined[contRBPs,paste0(TissueNames,"_iBAQ")]
    betas <- betas[,contRBPs]
    
    cors <- c()
    pVals <- c()
    
    for(elem in contRBPs){
      cT <- cor.test(unlist(betas[,elem]),
                     unlist(ExpressionLevels[elem,]),
                     method = "spearman", na.action=na.exclude)
      cors <- c(cors, cT$estimate)
      pVals <- c(pVals, cT$p.value)
    }
    
    pVals <- p.adjust(pVals, method = "BH")
    names(cors) <- contRBPs
    print("Spearman's correlation coefficients: ")
    print(cors[which(pVals < 0.05)])
    print("P-values: ")
    print(pVals[which(pVals < 0.05)])
    
    tp <- data.frame(RBP=contRBPs, 
                     cors=cors,       
                     FDR=round(pVals, digits = 2))
    
    #correlatingRBPs <- names(cors[which(pVals < 0.05)])
    correlatingRBPs <- names(cors)
    
    betas <- as.matrix(betas[,correlatingRBPs])
    rownames(betas) <- paste0(TissueNames, "_beta")
    betas <- as.data.table(melt(betas))
    betas[,c("Tissue","Var1") := tstrsplit(Var1,"_"),][,Var1:=NULL,]
    colnames(betas) <- c("RBP", "Beta", "Tissue")
    betas <- merge(betas, tp, by="RBP")
    
    ExpressionLevels <- t(ExpressionLevels[correlatingRBPs,])
    ExpressionLevels <- as.data.table(melt(ExpressionLevels))
    ExpressionLevels[,c("Tissue","Var1") := tstrsplit(Var1,"_"),][,Var1:=NULL,]
    colnames(ExpressionLevels) <- c("RBP", "Expression", "Tissue")
    
    dfTemp <- merge(betas, ExpressionLevels, by=c("RBP", "Tissue"))
    
    dfTemp$Expression <- 10^dfTemp$Expression
    dfTemp$RBP <- paste0(dfTemp$RBP, ", \n FDR =", dfTemp$FDR)
    p <-ggplot(dfTemp, aes(x=Beta, y=Expression))+
            geom_point(size=0.1, col="red")+
            facet_wrap(~ RBP, ncol = 10)+
            geom_hline(yintercept=1, color="darkgrey")+
            geom_vline(xintercept=1, color="darkgrey")+
            xlab(paste0("Tissue-specific effect size on ", responseType))+
            ylab("Tissue-specific protein expression level of the RBP")+
            scale_y_continuous(trans='log10', limits = c(10^6,10^10))
    
    ggsave("./PaperFigures/Figures/Appendix/13_2.pdf", 
           plot = p, 
           width = 8.0,
           height = 10)
    
    return(list(rbpFit=rbpFit, signRBPs=correlatingRBPs))
    
  }
  
  
  getExpVarOneMatrixByAnother <- function(M1, M2){
    
    M1[is.na(M1)] <- 0
    M2NA <- is.na(t(as.matrix(M2)))
    
    M2[is.na(M2)] <- 0
    M2Means <- colMeans(t(as.matrix(M2)), na.rm = T)
    
    M1 <- scale(t(M1), center = T, scale = F)
    M2 <- scale(t(M2), center = T, scale = F)
    
    m1.pca <- prcomp(M1)
    
    L <- m1.pca$rotation
    
    ### Predicted coordinates of protein variables (tissues)
    ## s <- protein %*% L
    m2Pre <- stats::predict(m1.pca,
                            newdata=M2)
    m2Recons <- m2Pre %*% t(L)
    
    M2[M2NA] <- NA
    m2Recons[M2NA] <- NA
    
    M2Var <- sum(M2^2, na.rm = TRUE)
    r2 <- round(1 - (sum((M2 - m2Recons)^2, na.rm = TRUE)/M2Var), digits = 2)
    return(r2)
    
  }
  
  
  readAttractDBPWM = function(id)
  {
    pwm = read.table(paste0(tempFileDir, "/UnsharedData/attract/processed/pwm/",id,".txt"))
    colnames(pwm) = c('A','C','G','T')
    pwm = log2(pwm)-log2(0.25)
    pwm = t(pwm)
    return(list(id, as.matrix(pwm)))
  }
  
  readRBPDBPWM = function(id, geneName)
  {
    pwm = read.table(paste0(tempFileDir,"/UnsharedData/RBPDB/PWMDir/",id))
    rownames(pwm) = c('A','C','G','T')
    return(list(id, as.matrix(pwm)))
  }
  
  
  pwmScore = function(pwmL,seqM)
  {
    id = pwmL[[1]][1]
    RBP_ensemblID = pwmL[[1]][2]
    pwm = pwmL[[2]]
    k = nchar(seqM[1])
    org = seqM
    if(k - ncol(pwm) + 1 < 1)
    {
      r = ncol(pwm) - k 
      seqM = paste0(paste0(rep('N',r),collapse=''),seqM,paste0(rep('N',r),collapse=''))
      k = nchar(seqM[1])
    }
    index = 1:(k - ncol(pwm) + 1)
    eax = data.table(kmer=org,pat=seqM)
    eax[,pwmScore:=max(PWMscoreStartingAt(pwm,pat, starting.at=index)),by=pat]
    eax[,pwmId:=id]
    eax[,RBP_ensemblID := RBP_ensemblID,]
    
    #maxScore = sum(apply(pwm,2,max)) #- log2(0.25)*ncol(pwm)
    
    k = nchar(org)
    
    if(k>=ncol(pwm))
    {
      maxScore = sum(apply(pwm,2,max))
    }else{
      maxScore = apply(pwm,2,max)
      r = ncol(pwm)-k
      
      maxScore = min(sapply(1:(r+1),function(x){sum(maxScore[x:(x+k-1)])}))
    }
    
    eax[,rScore:=pwmScore/maxScore]
    
    
    cseq = c('A','C','G','T')[apply(pwm,2,which.max)]%>%paste0(collapse='')
    
    eax[,consensus:=cseq]
    eax[,tLen:=nchar(kmer),]
    eax[,mLen:=nchar(consensus),]
    
    return(eax)
  }
  
