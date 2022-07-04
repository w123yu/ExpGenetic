#' @title Differential expression analysis
#'
#' @param P1_count A dataframe. The count data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the count of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_count A dataframe. Similar with P1_count, the count data of miRNA from the P2 species.
#' @param F1_count A dataframe. Similar with P1_count, the count data of miRNA from the F1 species.
#' @param P1_name A character. The text in the figure represents the name of one of the parents (P1).
#' @param P2_name A character. The text in the figure represents the name of one of the parents (P2).
#' @param F1_name A character. The text in the figure represents the name of hybrid offspring (F1).
#' @param type A character. "sRNA" or "mRNA".
#' @param homoeologs A dataframe. if the "type" is "mRNA", the option is required. The first column is the group name (unique) of homoeologs among three species, the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column is the ID in the F1 generation that is consistent with the P1 genome type, while the fifth column is the ID in the F1 generation that is consistent with the P2 genome type. Such as "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1".
#' @param count_threshold A numeric. Among the three species(P1,P2,F1), all replicates in at least one species with count values greater than or equal to 5 were retained. By default, the count value more than or equal to 5 is retained.
#' @param Pvalue A numeric. The threshold of significance test among different groups. Default is 0.05.
#'
#' @return A dataframe. Differential expression analysis results of miRNAs/siRNAs/genes expressed in each two species (count >= count_threshold).
#' @export
#' @examples
#' \dontrun{
#' polyDESeq(P1_count = P1_miRNA_count,
#'           P2_count = P2_miRNA_count,
#'           F1_count = F1_miRNA_count,
#'           P1_name = "B.napus(AACC)",
#'           P2_name = "B.rapa(AA)",
#'           F1_name = "B.napus x B.rapa (AAAACC)",type="sRNA")
#' }
polyDESeq <- function(P1_count,P2_count,F1_count,P1_name,P2_name,F1_name,type,homoeologs,count_threshold = 5,Pvalue = 0.05){
  if(ncol(P1_count) == ncol(P2_count) & ncol(P1_count) == ncol(F1_count)){
    colnum = ncol(P2_count)
    sample_number = colnum-1
  }else{
    cat("Error!!! Inconsistent biological replicates between different samples.")
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
  ##############################
  cat("F1_vs_P1","\n")
  P1_dds <- DESeq2::DESeq(P1_dds)
  cat("F1_vs_P2","\n")
  P2_dds <- DESeq2::DESeq(P2_dds)
  cat("P2_vs_P1","\n")
  P12_dds <- DESeq2::DESeq(P12_dds)
  ########
  F1_vs_P1 <- as.data.frame(DESeq2::results(P1_dds))
  F1_vs_P2 <- as.data.frame(DESeq2::results(P2_dds))
  P2_vs_P1 <- as.data.frame(DESeq2::results(P12_dds))
  ###########################
  F1_vs_P1 <- na.omit(F1_vs_P1)
  F1_vs_P2 <- na.omit(F1_vs_P2)
  P2_vs_P1 <- na.omit(P2_vs_P1)
  ####F1_vs_P1
  F1_vs_P1_up <-  F1_vs_P1[(F1_vs_P1$log2FoldChange > 1 & F1_vs_P1$pvalue < Pvalue), ]
  P1_vs_F1_up <-  F1_vs_P1[(F1_vs_P1$log2FoldChange < -1 & F1_vs_P1$pvalue < Pvalue), ]
  F1_vs_P1_num <- nrow(F1_vs_P1_up)
  P1_vs_F1_num <- nrow(P1_vs_F1_up)
  ####F1_vs_P2
  F1_vs_P2_up <-  F1_vs_P2[(F1_vs_P2$log2FoldChange > 1 & F1_vs_P2$pvalue < Pvalue), ]
  P2_vs_F1_up <-  F1_vs_P2[(F1_vs_P2$log2FoldChange < -1 & F1_vs_P2$pvalue < Pvalue), ]
  F1_vs_P2_num <- nrow(F1_vs_P2_up)
  P2_vs_F1_num <- nrow(P2_vs_F1_up)
  ####P1_vs_P2
  P2_vs_P1_up <-  P2_vs_P1[(P2_vs_P1$log2FoldChange > 1 & P2_vs_P1$pvalue < Pvalue), ]
  P1_vs_P2_up <-  P2_vs_P1[(P2_vs_P1$log2FoldChange < -1 & P2_vs_P1$pvalue < Pvalue), ]
  P2_vs_P1_num <- nrow(P2_vs_P1_up)
  P1_vs_P2_num <- nrow(P1_vs_P2_up)
  ###
  grid::grid.newpage()
  grid::grid.lines(x = grid::unit(c(1.5/5, 2.5/5), "npc"), y = grid::unit(c(1.5/5, 3.5/5), "npc"),
                   default.units = "npc",arrow = NULL, name = NULL,gp=grid::gpar(), draw = TRUE, vp = NULL)
  grid::grid.lines(x = grid::unit(c(1.5/5, 3.5/5), "npc"),y = grid::unit(c(1.5/5, 1.5/5), "npc"),
                   default.units = "npc",arrow = NULL, name = NULL,gp=grid::gpar(), draw = TRUE, vp = NULL)
  grid::grid.lines(x = grid::unit(c(3.5/5, 2.5/5), "npc"),y = grid::unit(c(1.5/5, 3.5/5), "npc"),
                   default.units = "npc",arrow = NULL, name = NULL, gp=grid::gpar(), draw = TRUE, vp = NULL)
  grid::grid.circle(x=1.5/5, y=1.5/5, r=0.1, default.units="npc", name=NULL, gp=grid::gpar(fill="#0099FF"), draw=TRUE, vp=NULL)
  grid::grid.circle(x=3.5/5, y=1.5/5, r=0.1, default.units="npc", name=NULL, gp=grid::gpar(fill="#FF9900"), draw=TRUE, vp=NULL)
  grid::grid.circle(x=2.5/5, y=3.5/5, r=0.1, default.units="npc", name=NULL, gp=grid::gpar(fill="#00CC00"), draw=TRUE, vp=NULL)
  grid::grid.text(P1_name,x=1.5/5,y=1.5/5)
  grid::grid.text(P2_name,x=3.5/5,y=1.5/5)
  grid::grid.text(F1_name,x=2.5/5,y=3.5/5)
  grid::grid.text(P1_vs_F1_num,x=1.65/5,y=2.15/5)
  grid::grid.text(P1_vs_P2_num,x=2/5,y=1.3/5)
  grid::grid.text(P2_vs_P1_num,x=3/5,y=1.3/5)
  grid::grid.text(P2_vs_F1_num,x=3.35/5,y=2.15/5)
  grid::grid.text(F1_vs_P1_num,x=2/5,y=3/5)
  grid::grid.text(F1_vs_P2_num,x=3/5,y=3/5)
}
