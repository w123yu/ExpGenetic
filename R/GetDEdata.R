#' GetDEdata: Get specific results of differential expressed genes.
#'
#' @param P1_count A dataframe. The count data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the count of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_count A dataframe. Similar with P1_count, the count data of miRNA from the P2 species.
#' @param F1_count A dataframe. Similar with P1_count, the count data of miRNA from the F1 species.
#' @param P1_name A character. The text in the figure represents the name of one of the parents (P1).
#' @param P2_name A character. The text in the figure represents the name of one of the parents (P2).
#' @param F1_name A character. The text in the figure represents the name of hybrid offspring (F1).
#' @param output_type A character. 'F1_vs_P1', 'F1_vs_P2' or 'P2_vs_P1'.
#' @param type A character. "sRNA" or "mRNA".
#' @param homoeologs A dataframe. if the "type" is "mRNA", the option is required. The first column is the group name (unique) of homoeologs among three species, the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column is the ID in the F1 generation that is consistent with the P1 genome type, while the fifth column is the ID in the F1 generation that is consistent with the P2 genome type. Such as "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1".
#' @param count_threshold A numeric. Among the three species(P1,P2,F1), all replicates in at least one species with count values greater than or equal to 5 were retained. By default, the count value more than or equal to 5 is retained.
#' @param Pvalue A numeric. The threshold of significance test among different groups. Default is 0.05.
#'
#' @return A dataframe. Differential expression analysis results of miRNA/siRNA/gene expressed in each two species (count >= count_threshold).
#' @export
#'
#' @examples
#' \dontrun{
#' P2_vs_P1 <- GetDEdata(P1_count = P1_miRNA_count,
#'                       P2_count = P2_miRNA_count,
#'                       F1_count = F1_miRNA_count,
#'                       P1_name = "B.napus(AACC)",
#'                       P2_name = "B.rapa(AA)",
#'                       F1_name = "B.napus x B.rapa (AAAACC)",
#'                       output_type = "P2_vs_P1", type="sRNA")
#' }
GetDEdata <- function(P1_count,P2_count,F1_count,P1_name,P2_name,F1_name,output_type,type,homoeologs,count_threshold = 5,Pvalue = 0.05){
  if(ncol(P1_count) == ncol(P2_count) & ncol(P1_count) == ncol(F1_count)){
    colnum = ncol(P2_count)
    sample_number = colnum-1
  }else{
    cat("Error!!! Inconsistent biological replicates between different samples.")
  }
  if(output_type == "F1_vs_P1" | output_type == "F1_vs_P2" | output_type == "P2_vs_P1"){
    cat("Run start:","\n")
  }else{
    cat("Error!!! Please select the output option of 'F1_vs_P1', 'F1_vs_P2' or 'P2_vs_P1'.")
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
    colnames(P1_homo) <- c("orthgroup","P1_id")
    colnames(P2_homo) <- c("orthgroup","P2_id")
    F1P1_homo <- homoeologs[,c(1,4)]
    F1P2_homo <- homoeologs[,c(1,5)]
    colnames(F1P1_homo) <- c("orthgroup","F1_id")
    colnames(F1P2_homo) <- c("orthgroup","F1_id")
    ###
    P1_gene <- merge.data.frame(P1_homo,P1_count,by="P1_id")
    P2_gene <- merge.data.frame(P2_homo,P2_count,by="P2_id")
    F1P1_gene <- merge.data.frame(F1P1_homo,F1_count,by="F1_id")
    F1P2_gene <- merge.data.frame(F1P2_homo,F1_count,by="F1_id")
    ###
    F1_res <- merge.data.frame(F1P1_gene,F1P2_gene,by="orthgroup")
    F1_res$F1_1 <- F1_res$F1_1.x+F1_res$F1_1.y
    F1_res$F1_2 <- F1_res$F1_2.x+F1_res$F1_2.y
    F1_res$F1_3 <- F1_res$F1_3.x+F1_res$F1_3.y
    F1_gene <- F1_res[,c("orthgroup","F1_1","F1_2","F1_3")]
    P1_gene <- subset(P1_gene,select=-P1_id)
    P2_gene <- subset(P2_gene,select=-P2_id)
    input_data <- merge.data.frame(merge.data.frame(P1_gene,P2_gene,by="orthgroup"),F1_gene,by="orthgroup")
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
    cat("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
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
      cat("Error!!! Please select the output option of 'F1_vs_P1', 'F1_vs_P2' or 'P2_vs_P1'.")
    }
  }else if(type=="mRNA"){
    if(output_type == "F1_vs_P1"){
      P1_dds <- DESeq2::DESeq(P1_dds)
      F1_vs_P1 <- as.data.frame(DESeq2::results(P1_dds))
      ##
      F1_vs_P1_re<-cbind(row.names(F1_vs_P1),data.frame(F1_vs_P1,row.names = NULL))
      colnames(F1_vs_P1_re)[1] <- "orthgroup"
      F1_vs_P1 <- transform(F1_vs_P1_re,Treat_vs_Control <- "F1_vs_P1")
      return(F1_vs_P1)
    }else if(output_type == "F1_vs_P2"){
      P2_dds <- DESeq2::DESeq(P2_dds)
      F1_vs_P2 <- as.data.frame(DESeq2::results(P2_dds))
      ##
      F1_vs_P2_re<-cbind(row.names(F1_vs_P2),data.frame(F1_vs_P2,row.names = NULL))
      colnames(F1_vs_P2_re)[1] <- "orthgroup"
      F1_vs_P2 <- transform(F1_vs_P2_re,Treat_vs_Control <- "F1_vs_P2")
      return(F1_vs_P2)
    }else if(output_type == "P2_vs_P1"){
      P12_dds <- DESeq2::DESeq(P12_dds)
      P2_vs_P1 <- as.data.frame(DESeq2::results(P12_dds))
      ##
      P2_vs_P1_re<-cbind(row.names(P2_vs_P1),data.frame(P2_vs_P1,row.names = NULL))
      colnames(P2_vs_P1_re)[1] <- "orthgroup"
      P2_vs_P1 <- transform(P2_vs_P1_re,Treat_vs_Control <- "P2_vs_P1")
      return(P2_vs_P1)
    }else{
      cat("Error!!! Please select the output option of 'F1_vs_P1', 'F1_vs_P2' or 'P2_vs_P1'.")
    }
  }else{
    cat("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
}








