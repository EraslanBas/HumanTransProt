#### Query DAVID database for enrichment of GO terms associated with
#### foreground genes in a population of background genes.
####
#### @foregroundGenes 
#### @backgroundGenes - should include foreground genes
#### @idType
#### @fileNamePrefix - should be unique per query

queryDavidGO <- function(foregroundGenes, 
                       backgroundGenes,
                       idType="ENSEMBL_GENE_ID",
                       filePath="/Users/basak/WORKSPACE/humanTransProt/DAVID_QUERIES/",
                       fileNamePrefix)
{
  # Load RDAVIDWebService.
  library("RDAVIDWebService")
  
  # Create a DAVIDWebService object connected to David, using your registration email.
  # To register, go to: http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html.
  david <- DAVIDWebService$new(email='eraslan@in.tum.de', 
                               url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  
  # Define foreground and background gene lists.
  # The foreground list should be contained within the background list.
  FG <- addList(david, foregroundGenes, idType="ENSEMBL_GENE_ID", listName="isClass", listType="Gene")
  BG <- addList(david, backgroundGenes, idType="ENSEMBL_GENE_ID", listName="all", listType="Background")
  

  # Specifiy annotation categories.
  setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
  
  # Get functional annotation chart as R object.
  FuncAnnotChart <- getFunctionalAnnotationChart(david)
  
  # Print functional annotation chart to file.
  getFunctionalAnnotationChartFile(david, file.path(filePath, paste0(fileNamePrefix, "_chart.tsv")))
  
  # Get functional annotation clustering (limited to 3000 genes).
  FuncAnnotClust <- getClusterReport(david)
  
  # Print functional annotation clustering to file (limited to 3000 genes).
  getClusterReportFile(david, file.path(filePath, paste0(fileNamePrefix, "_clust.tsv")))
  
  if(file.exists(file.path( filePath, paste0(fileNamePrefix, "_chart.tsv")))){
    return(T)
  }else
    {return(F)}
  
}


getDAVIDGO <- function(foregroundGenes, 
                        backgroundGenes,
                        idType="ENSEMBL_GENE_ID",
                        filePath="/Users/basak/WORKSPACE/humanTransProt/DAVID_QUERIES/",
                        fileNamePrefix,
                        mainStr="")
{
  
  if(!file.exists(file.path(filePath, paste0(fileNamePrefix, "_chart.tsv")))){
      res <- queryDavidGO(foregroundGenes, 
               backgroundGenes,
               idType,
               filePath,
               fileNamePrefix)
      if(!res){
        print(" Query to DAVID database could not be accomplished.")
        return(NULL)
      }
  }
  
  chartTable <- as.data.table(read.csv(file.path(filePath, paste0(fileNamePrefix, "_chart.tsv")), 
                                       sep = "\t",
                                       stringsAsFactors = F) )   
  setorder(chartTable, Category, FDR)
  
  chartTable <- chartTable[FDR<0.1,.(Term, Count, Fold.Enrichment, FDR),by=Category]
  
  if(nrow(chartTable) > 0){
    chartTable[,c("Termx","Term") := tstrsplit(Term,"~"),]
    plotList <- plotDAVIDGO(chartTable, mainStr)
  }else{
    chartTable <- NULL
    plotList <- NULL
  }
  
  return(list(chartTable=chartTable, plotList=plotList))
  
}

plotDAVIDGO <- function(inputTB, mainStr=""){
  x <- data.frame(category=c("GOTERM_BP_ALL", "GOTERM_CC_ALL","GOTERM_MF_ALL"),
                  names= c("Biological Process",
                           "Cellular Component",
                           "Metabolic Function"),
                  stringsAsFactors = F)
  rownames(x) <- x$category
  
  plotList <- list()
  for(elem in x$category){
    tp <- as.data.frame(inputTB[Category == elem,])
    tp <- tp[order(-tp$FDR),]
    
    if(nrow(tp)==0){
      p <- NULL
    }else{
      
      tp <- tp[1:15,]
      tp$Term <- factor(tp$Term, levels=tp$Term)
      p <- ggplot(tp, aes(x=Term, y=-log10(FDR)))+
        geom_bar(stat="identity", fill="lightgrey")+
        coord_flip()+
        geom_text(aes(x=Term, y=0.05,label=Term),
                  #position = position_dodge(width = 1),
                  hjust = 0,
                  size=2.5) +
        theme(axis.text.y = element_blank(),
              axis.title.y=element_blank())+
        ggtitle(paste0( x[elem,"names"] ))
      
    }
    plotList <- lappend(plotList, p)
  }
  names(plotList) <- x$category
  
  return(plotList)
}


plotForegroundGeneLevels <- function(inputDF, foregroundGenes, expStr){
  
  tp <- inputDF[foregroundGenes,paste0(TissueNames, expStr)]
  colnames(tp) <- TissueNames
  tp <- t(as.matrix(apply(tp, 1, function(x){x-median(x, na.rm = T)})))
  tp[tp< -1] <- -1
  tp[tp> 1] <- 1
  print(pheatmap(tp[,], show_rownames = F, treeheight_col = 0, treeheight_row = 0,border_color = NA))
}
