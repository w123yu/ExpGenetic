#' @title Filtering out lowly expressed genes based on count
#' @description Regarding the criteria for filtering out lowly expressed genes, no less than the count threshold in all replicates.
#' @param P1_count A data frame. The count table of genes in P1 species. For the count table, the first column is the gene identifier, and other columns are read counts of the gene in each biological replicate.
#' @param P2_count A data frame. The count table of genes in P2 species.
#' @param F1_count A data frame. The count table of genes in F1 species.
#' @param type A character. "sRNA" or "mRNA".
#' @param homoeologs A data frame. Orthologous relationships of genes within the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param count_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 5 (the count values in all replicates).
#' @return A data frame.
#' @export
#' @details The 'homoeologs' table contains the orthologs pairs. In detail, the first column is the group name (unique) of homoeologs among three species (Parents: P1; P2, Progeny: F1), the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column and fifth columns are the identifier of F1 orthologs derived from P1 and P2 ancestors, respectively (e.g. "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1").
#' @examples
#' Count5result <- Countfilter(P1_count = P1_miRNA_count,
#'                             P2_count = P2_miRNA_count,
#'                             F1_count = F1_miRNA_count,
#'                             type = "sRNA", count_threshold = 5)
Countfilter <- function(P1_count,P2_count,F1_count,type,homoeologs,count_threshold = 5){
  if(ncol(P1_count) == ncol(P2_count) & ncol(P1_count) == ncol(F1_count)){
    colnum = ncol(P2_count)
    sample_number = colnum-1
  }else{
    message("Error!!! Inconsistent biological replicates between different samples.")
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
    names(P1_count) <- P1_colname
    names(P2_count) <- P2_colname
    names(F1_count) <- F1_colname
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
    P1_gene <- merge.data.frame(P1_homo,P1_count,by="P1_id")
    P2_gene <- merge.data.frame(P2_homo,P2_count,by="P2_id")
    F1P1_gene <- merge.data.frame(F1P1_homo,F1_count,by="F1_id")
    F1P2_gene <- merge.data.frame(F1P2_homo,F1_count,by="F1_id")
    ###
    F1_res <- merge.data.frame(F1P1_gene,F1P2_gene,by="orthogroup")
    F1_res$F1_1 <- F1_res$F1_1.x+F1_res$F1_1.y
    F1_res$F1_2 <- F1_res$F1_2.x+F1_res$F1_2.y
    F1_res$F1_3 <- F1_res$F1_3.x+F1_res$F1_3.y
    F1_gene <- F1_res[,c("orthogroup","F1_1","F1_2","F1_3")]
    P1_gene <- subset(P1_gene,select=-P1_id)
    P2_gene <- subset(P2_gene,select=-P2_id)
    input_data <- merge.data.frame(merge.data.frame(P1_gene,P2_gene,by="orthogroup"),F1_gene,by="orthogroup")
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
    names(P1_count) <- P1_colname
    names(P2_count) <- P2_colname
    names(F1_count) <- F1_colname
    input_data <- merge.data.frame(merge.data.frame(P1_count,P2_count,all=TRUE),F1_count,all=TRUE)
    input_data[is.na(input_data)] <- 0
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
  ############
  sum <- 3*sample_number
  row.names(input_data) <- input_data[,1]
  input_data <- input_data[,-1]
  keep <- as.data.frame(input_data >= count_threshold)
  value <- rowSums(keep[,c(1:sample_number)])==sample_number | rowSums(keep[,c((sample_number+1):(2*sample_number))])==sample_number | rowSums(keep[,c((2*sample_number+1):(sample_number*3))])==sample_number
  filter_result<-input_data[value,]
  filter_result<-cbind(row.names(filter_result),data.frame(filter_result,row.names = NULL))
  if(type=="mRNA"){
    colnames(filter_result)[1] <- "orthogroup"
  }else if(type=="sRNA"){
    colnames(filter_result)[1] <- "sequence"
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
  return(filter_result)
}

#' @title Filtering out lowly expressed genes based on RPM
#' @description Regarding the criteria for filtering out lowly expressed genes, no less than the RPM threshold in all replicates.
#' @param P1_RPM A data frame. The RPM table of genes in P1 species. For the RPM table, the first column is the gene identifier (e.g. sequences of sRNA, Gene ID), and other columns are the RPM values of the gene in each biological replicate.
#' @param P2_RPM A data frame. The RPM table of genes in P2 species.
#' @param F1_RPM A data frame. The RPM table of genes in F1 species.
#' @param type A character. "sRNA" or "mRNA".
#' @param homoeologs A data frame. Orthologous relationships of genes within the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param rpm_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 1 (the average RPM of all replicates).
#' @details The 'homoeologs' table contains the orthologs pairs. In detail, the first column is the group name (unique) of homoeologs among three species (Parents: P1; P2, Progeny: F1), the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column and fifth columns are the identifier of F1 orthologs derived from P1 and P2 ancestors, respectively (e.g. "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1").
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' Rpm1result <- Rpmfilter(P1_RPM = P1_miRNA_rpm,
#'                          P2_RPM = P2_miRNA_rpm,
#'                          F1_RPM = F1_miRNA_rpm,
#'                          type = "sRNA", rpm_threshold = 1)
Rpmfilter <-function(P1_RPM,P2_RPM,F1_RPM,type,homoeologs,rpm_threshold = 1){
  if(ncol(P1_RPM) == ncol(P2_RPM) & ncol(P1_RPM) == ncol(F1_RPM)){
    colnum = ncol(P2_RPM)
    sample_number = colnum-1
  }else{
    message("Error!!! Inconsistent biological replicates between different samples.")
  }
  func_mirnafiliter <- function(data1,threshold){
    Average_rpm <- as.data.frame(apply(data1[,c(-1)],1,mean))
    data1_value <- Average_rpm >= threshold
    data2 <- data1[data1_value,]
    return(data2)
  }
  if(type=="sRNA"){
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
    ###########################
    P1_mirna <- func_mirnafiliter(data1 = P1_RPM, threshold = rpm_threshold)
    P2_mirna <- func_mirnafiliter(data1 = P2_RPM, threshold = rpm_threshold)
    F1_mirna <- func_mirnafiliter(data1 = F1_RPM, threshold = rpm_threshold)
    ####################################
    Parent <- merge.data.frame(P1_mirna,P2_mirna,all=TRUE)
    filter_result <- merge.data.frame(Parent,F1_mirna,all=TRUE)
    fina_result <- na.omit(filter_result)
    ####################################
    rownames(fina_result) <- 1:nrow(fina_result)
    return(fina_result)
  }else if(type=="mRNA"){
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
    ###########################
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
    P1_mRNA <- func_mirnafiliter(data1 = P1_gene, threshold = rpm_threshold)
    P2_mRNA <- func_mirnafiliter(data1 = P2_gene, threshold = rpm_threshold)
    F1_mRNA <- func_mirnafiliter(data1 = F1_gene, threshold = rpm_threshold)
    ####################################
    Parent <- merge.data.frame(P1_mRNA,P2_mRNA,all = TRUE)
    filter_result <- merge.data.frame(Parent,F1_mRNA,all = TRUE)
    fina_result <- na.omit(filter_result)
    ####################################
    rownames(fina_result) <- 1:nrow(fina_result)
    return(fina_result)
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
}

