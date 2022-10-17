#' Non-additive expression analysis
#'
#' @description About the classification method based on |d/a|, the additive (a) and dominant (d) values were calculated by the expression level of each miRNA. Edwards et al. proposed that the "|d/a|" can be used as the criterion to estimate the expression patterns of miRNAs. Specific classification criteria are as follows, |d/a| <= 0.2, additivity; |d/a| > 0.2 and |d/a| <= 0.8, partial dominance; |d/a| > 0.8 and |d/a| <= 1.2, dominance; |d/a| > 1.2, overdominance.
#'
#' @param P1_RPM A data frame. The RPM table of genes in P1 species. For the RPM table, the first column is the gene identifier, and other columns are the RPM values of the genes in each biological replicate.
#' @param P2_RPM A data frame. The RPM table of genes in P2 species.
#' @param F1_RPM A data frame. The RPM table of genes in F1 species.
#' @param type A character. "sRNA" or "mRNA".
#' @param homoeologs A data frame. Orthologous relationships of genes in the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param rpm_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 1 (the average RPM of all replicates).
#' @details The 'homoeologs' table contains the orthologs pairs. In detail, the first column is the group name (unique) of homoeologs among three species (Parents: P1; P2, Progeny: F1), the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column and fifth columns are the identifier of F1 orthologs derived from P1 and P2 ancestors, respectively (e.g. "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1").
#' @references Edwards MD, Stuber CW, Wendel JF. Molecular-marker-facilitated investigations of quantitative-trait loci in maize. I. Numbers, genomic distribution and types of gene action. Genetics. 1987 May;116(1):113-25.
#' @return A data frame. Classification results of non-additive expression analysis based on |d/a|.
#' @export
#' @examples
#' DAresult <- GetDAtable(P1_RPM = P1_miRNA_rpm,
#'                        P2_RPM = P2_miRNA_rpm,
#'                        F1_RPM = F1_miRNA_rpm,
#'                        type = "sRNA", rpm_threshold = 1)
GetDAtable <- function(P1_RPM,P2_RPM,F1_RPM,type,homoeologs,rpm_threshold = 1){
  if(ncol(P1_RPM) == ncol(P2_RPM) & ncol(P2_RPM) == ncol(F1_RPM)){
    colnum = ncol(P1_RPM)
    sample_number = colnum-1
  }else{
    message("Error!!! Inconsistent biological replicates between different samples.")
  }
  func_mirnafiliter <- function(data1,threshold,group){
    Average_rpm <- as.data.frame(apply(data1[,c(-1)],1,mean))
    data1_value <- Average_rpm >= threshold
    data2 <- data1[data1_value,]
    result <- cbind(data2[,1],as.data.frame(apply(data2[,c(-1)],1,mean)))
    colnames(result) <- c("sequence",group)
    return(result)
  }
  if(type=="mRNA"){
    P1_colname <- "P1_id"
    P2_colname <- "P2_id"
    F1_colname <- "F1_id"
    for (i in (1:sample_number)){
      P1_id <- paste("P1_",i,sep = "")
      P2_id <- paste("P2_",i,sep = "")
      F1_id <- paste("F1_",i,sep = "")
      P1_colname <- c(P1_colname,P1_id)
      P2_colname <- c(P2_colname,P2_id)
      F1_colname <- c(F1_colname,F1_id)
    }
    names(P1_RPM) <- P1_colname
    names(P2_RPM) <- P2_colname
    names(F1_RPM) <- F1_colname
    ###
    P1_homo <- homoeologs[,c(1,2)]
    P2_homo <- homoeologs[,c(1,3)]
    colnames(P1_homo) <- c("orthogroup","P1_id")
    colnames(P2_homo) <- c("orthogroup","P2_id")
    F1P1_homo <- homoeologs[,c(1,4)]
    F1P2_homo <- homoeologs[,c(1,5)]
    colnames(F1P1_homo) <- c("orthogroup","F1_id")
    colnames(F1P2_homo) <- c("orthogroup","F1_id")
    ###
    P1_gene <- merge.data.frame(P1_homo,P1_RPM,by="P1_id")
    P2_gene <- merge.data.frame(P2_homo,P2_RPM,by="P2_id")
    F1P1_gene <- merge.data.frame(F1P1_homo,F1_RPM,by="F1_id")
    F1P2_gene <- merge.data.frame(F1P2_homo,F1_RPM,by="F1_id")
    ###
    F1_res <- merge.data.frame(F1P1_gene,F1P2_gene,by="orthogroup")
    F1_res$F1_1 <- F1_res$F1_1.x+F1_res$F1_1.y
    F1_res$F1_2 <- F1_res$F1_2.x+F1_res$F1_2.y
    F1_res$F1_3 <- F1_res$F1_3.x+F1_res$F1_3.y
    F1_gene <- F1_res[,c("orthogroup","F1_1","F1_2","F1_3")]
    P1_gene <- subset(P1_gene,select=-P1_id)
    P2_gene <- subset(P2_gene,select=-P2_id)
    P1_mirna <- func_mirnafiliter(data1 = P1_gene,threshold = rpm_threshold,group = "P1_average_rpm")
    P2_mirna <- func_mirnafiliter(data1 = P2_gene,threshold = rpm_threshold,group = "P2_average_rpm")
    F1_mirna <- func_mirnafiliter(data1 = F1_gene,threshold = rpm_threshold,group = "F1_average_rpm")
  }else if(type=="sRNA"){
    P1_colname <- "sequence"
    P2_colname <- "sequence"
    F1_colname <- "sequence"
    for (i in (1:sample_number)){
      P1_id <- paste("P1_",i,sep = "")
      P2_id <- paste("P2_",i,sep = "")
      F1_id <- paste("F1_",i,sep = "")
      P1_colname <- c(P1_colname,P1_id)
      P2_colname <- c(P2_colname,P2_id)
      F1_colname <- c(F1_colname,F1_id)
    }
    names(P1_RPM) <- P1_colname
    names(P2_RPM) <- P2_colname
    names(F1_RPM) <- F1_colname
    P1_mirna <- func_mirnafiliter(data1 = P1_RPM,threshold = rpm_threshold,group = "P1_average_rpm")
    P2_mirna <- func_mirnafiliter(data1 = P2_RPM,threshold = rpm_threshold,group = "P2_average_rpm")
    F1_mirna <- func_mirnafiliter(data1 = F1_RPM,threshold = rpm_threshold,group = "F1_average_rpm")
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
  ####################################
  Parent <- merge.data.frame(P1_mirna,P2_mirna,all=TRUE)
  filter_result <- merge.data.frame(Parent,F1_mirna,all=TRUE)
  fina_result <- na.omit(filter_result)
  #####
  fina_result$d <- fina_result$F1_average-((fina_result$P1_average+fina_result$P2_average)/2)
  fina_result$a <- (fina_result$P1_average-fina_result$P2_average)/2
  fina_result$`abs(d/a)` <- abs(fina_result$d/fina_result$a)
  da <- fina_result$`abs(d/a)`
  if (length(da[da<=0.2]>0)){
    additivity <- fina_result[fina_result$`abs(d/a)`<=0.2,]
    additivity$condition <- "Additivity"
  }else{
    additivity <- data.frame()
  }
  if (length(da[(da>0.2)&(da<=0.8)]>0)){
    partialDom <- fina_result[(fina_result$`abs(d/a)`>0.2)&(fina_result$`abs(d/a)`<=0.8),]
    partialDom$condition <- "Partial Dominance"
  }else{
    partialDom <- data.frame()
  }
  if (length(da[(da>0.8)&(da<=1.2)]>0)){
    Dominance <- fina_result[(fina_result$`abs(d/a)`>0.8)&(fina_result$`abs(d/a)`<=1.2),]
    Dominance$condition <- "Dominance"
  }else{
    Dominance <- data.frame()
  }
  if (length(da[da>1.2]>0)){
    Overdominance <- fina_result[fina_result$`abs(d/a)`>1.2,]
    Overdominance$condition <- "Overdominance"
  }else{
    Overdominance <- data.frame()
  }
  result <- rbind(additivity,partialDom,Dominance,Overdominance)
  row.names(result)<- 1:nrow(result)
  if(type=="mRNA"){
    colnames(result)[1] <- "orthogroup"
    return(result)
  }else if(type=="sRNA"){
    return(result)
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
}
