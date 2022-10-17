#' @title Make a Triangle Diagram
#' @description The count matrix of different species as the input data to perform differential expression analysis using DESeq2. And the number of differentially expressed genes between any two species is marked on the triangle diagram.
#' @param P1_count A data frame. The count table of genes in P1 species. For the count table, the first column is the gene identifier, and other columns are the read counts of the genes in each biological replicate.
#' @param P2_count A data frame. The count table of genes in P2 species.
#' @param F1_count A data frame. The count table of genes in F1 species.
#' @param P1_name A character. Category names of P1 species.
#' @param P2_name A character. Category names of P2 species.
#' @param F1_name A character. Category names of F1 species.
#' @param type A character. "sRNA" or "mRNA".
#' @param homoeologs A data frame. Orthologous relationships of genes in the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param count_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 5 (the count values in all replicates).
#' @param Pvalue A numeric. Threshold for significance test in differential expression analysis. Default is 0.05.
#' @details The 'homoeologs' table contains the orthologs pairs. In detail, the first column is the group name (unique) of homoeologs among three species (Parents: P1;P2, Progeny: F1), the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column and fifth columns are the identifier of F1 orthologs derived from P1 and P2 ancestors, respectively (e.g. "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1").
#' @return Triangle Diagram
#' @export
#' @examples
#' \donttest{
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
  message("F1_vs_P1","\n")
  P1_dds <- DESeq2::DESeq(P1_dds)
  message("F1_vs_P2","\n")
  P2_dds <- DESeq2::DESeq(P2_dds)
  message("P2_vs_P1","\n")
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
