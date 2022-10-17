#' @title Get the results of differential expression analysis.
#' @description Extract the results of differential expression analysis.
#' @param P1_count A data frame. The count table of genes in P1 species. For the count table, the first column is the gene identifier, and other columns are the corresponding expression levels of the genes in each biological replicate.
#' @param P2_count A data frame. The count table of genes in P2 species.
#' @param F1_count A data frame. The count table of genes in F1 species.
#' @param type A character. "sRNA" or "mRNA".
#' @param output_type A character. "F1_vs_P1", "F1_vs_P2" or "P2_vs_P1".
#' @param homoeologs A data frame. Orthologous relationships of genes in the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param count_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 5 (the count values in all replicates).
#' @details F1_vs_P1: Results of differential expression analysis using DESeq2. Parental P1 was used as the control group and F1 was used as the treatment group. If the log2FoldChange of a gene is positive, it means that the expression level of the gene in F1 is higher than that in P1. F1_vs_P2: Results of differential expression analysis using DESeq2. Parental P2 was used as the control group and F1 was used as the treatment group. P2_vs_P1: Results of differential expression analysis using DESeq2. Parental P1 was used as the control group and P2 was used as the treatment group.
#' @return A data frame. Differential expression analysis results.
#' @export
#'
#' @examples
#' \donttest{
#' P2_vs_P1 <- GetDEdata(P1_count = P1_miRNA_count,
#'                       P2_count = P2_miRNA_count,
#'                       F1_count = F1_miRNA_count,
#'                       output_type = "P2_vs_P1", type="sRNA")
#' }
GetDEdata <- function(P1_count,P2_count,F1_count,output_type,type,homoeologs,count_threshold = 5){
  if(ncol(P1_count) == ncol(P2_count) & ncol(P1_count) == ncol(F1_count)){
    colnum = ncol(P2_count)
    sample_number = colnum-1
  }else{
    message("Error!!! Inconsistent biological replicates between different samples.")
  }
  if(output_type == "F1_vs_P1" | output_type == "F1_vs_P2" | output_type == "P2_vs_P1"){
    message("Run start:","\n")
  }else{
    message("Error!!! Please select the output option of 'F1_vs_P1', 'F1_vs_P2' or 'P2_vs_P1'.")
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
  ####################################
  P1.tab <- filter_result[,c(1:sample_number)]
  P2.tab <- filter_result[,c((sample_number+1):(2*sample_number))]
  F1.tab <- filter_result[,c((2*sample_number+1):(sample_number*3))]
  #############################
  P1_condition <- factor(c(rep("P1",sample_number),rep("F1",sample_number)),levels=c("P1","F1"))
  P2_condition <- factor(c(rep("P2",sample_number),rep("F1",sample_number)),levels=c("P2","F1"))
  P12_condition <- factor(c(rep("P1",sample_number),rep("P2",sample_number)),levels=c("P1","P2"))
  #############################
  P1F1.tab <- cbind(P1.tab,F1.tab)
  P2F1.tab <- cbind(P2.tab,F1.tab)
  P12.tab <- cbind(P1.tab,P2.tab)
  ##################################
  P1_colData <- data.frame(colnames(P1F1.tab),P1_condition)
  P2_colData <- data.frame(colnames(P2F1.tab),P2_condition)
  P12_colData <- data.frame(colnames(P12.tab),P12_condition)
  #############################
  P1_dds <- DESeq2::DESeqDataSetFromMatrix(countData = P1F1.tab,colData = P1_colData,design = ~P1_condition)
  P2_dds <- DESeq2::DESeqDataSetFromMatrix(countData = P2F1.tab,colData = P2_colData,design = ~P2_condition)
  P12_dds <- DESeq2::DESeqDataSetFromMatrix(countData = P12.tab,colData = P12_colData,design = ~P12_condition)
  ###################################
  ###################################
  ###
  if(type=="sRNA"){
    if(output_type == "F1_vs_P1"){
      P1_dds <- DESeq2::DESeq(P1_dds)
      F1_vs_P1 <- as.data.frame(DESeq2::results(P1_dds))
      #
      F1_vs_P1<-cbind(row.names(F1_vs_P1),data.frame(F1_vs_P1,row.names = NULL))
      colnames(F1_vs_P1)[1] <- "sequence"
      F1_vs_P1 <- transform(F1_vs_P1,Treat_vs_Control <- "F1_vs_P1")
      return(F1_vs_P1)
    }else if(output_type == "F1_vs_P2"){
      P2_dds <- DESeq2::DESeq(P2_dds)
      F1_vs_P2 <- as.data.frame(DESeq2::results(P2_dds))
      #
      F1_vs_P2<-cbind(row.names(F1_vs_P2),data.frame(F1_vs_P2,row.names = NULL))
      colnames(F1_vs_P2)[1] <- "sequence"
      F1_vs_P2 <- transform(F1_vs_P2,Treat_vs_Control <- "F1_vs_P2")
      return(F1_vs_P2)
    }else if(output_type == "P2_vs_P1"){
      P12_dds <- DESeq2::DESeq(P12_dds)
      P2_vs_P1 <- as.data.frame(DESeq2::results(P12_dds))
      #
      P2_vs_P1<-cbind(row.names(P2_vs_P1),data.frame(P2_vs_P1,row.names = NULL))
      colnames(P2_vs_P1)[1] <- "sequence"
      P2_vs_P1 <- transform(P2_vs_P1,Treat_vs_Control <- "P2_vs_P1")
      return(P2_vs_P1)
    }else{
      message("Error!!! Please select the output option of 'F1_vs_P1', 'F1_vs_P2' or 'P2_vs_P1'.")
    }
  }else if(type=="mRNA"){
    if(output_type == "F1_vs_P1"){
      P1_dds <- DESeq2::DESeq(P1_dds)
      F1_vs_P1 <- as.data.frame(DESeq2::results(P1_dds))
      ##
      F1_vs_P1_re<-cbind(row.names(F1_vs_P1),data.frame(F1_vs_P1,row.names = NULL))
      colnames(F1_vs_P1_re)[1] <- "orthogroup"
      F1_vs_P1 <- transform(F1_vs_P1_re,Treat_vs_Control <- "F1_vs_P1")
      return(F1_vs_P1)
    }else if(output_type == "F1_vs_P2"){
      P2_dds <- DESeq2::DESeq(P2_dds)
      F1_vs_P2 <- as.data.frame(DESeq2::results(P2_dds))
      ##
      F1_vs_P2_re<-cbind(row.names(F1_vs_P2),data.frame(F1_vs_P2,row.names = NULL))
      colnames(F1_vs_P2_re)[1] <- "orthogroup"
      F1_vs_P2 <- transform(F1_vs_P2_re,Treat_vs_Control <- "F1_vs_P2")
      return(F1_vs_P2)
    }else if(output_type == "P2_vs_P1"){
      P12_dds <- DESeq2::DESeq(P12_dds)
      P2_vs_P1 <- as.data.frame(DESeq2::results(P12_dds))
      ##
      P2_vs_P1_re<-cbind(row.names(P2_vs_P1),data.frame(P2_vs_P1,row.names = NULL))
      colnames(P2_vs_P1_re)[1] <- "orthogroup"
      P2_vs_P1 <- transform(P2_vs_P1_re,Treat_vs_Control <- "P2_vs_P1")
      return(P2_vs_P1)
    }else{
      message("Error!!! Please select the output option of 'F1_vs_P1', 'F1_vs_P2' or 'P2_vs_P1'.")
    }
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
}








