#' Make a three-set Venn Diagram
#'
#' @description This function creates a Venn Diagram to display the overlap of expressed genes between three sets (parents and progeny).
#'
#' @param P1_RPM A data frame. The RPM table of genes in P1 species. For the RPM table, the first column is the gene identifier, and other columns are the RPM values of the genes in each biological replicate.
#' @param P2_RPM A data frame. The RPM table of genes in P2 species.
#' @param F1_RPM A data frame. The RPM table of genes in F1 species.
#' @param P1_name Character. Category names of P1 species.
#' @param P2_name Character. Category names of P2 species.
#' @param F1_name Character. Category names of F1 species.
#' @param type Character. "sRNA" or "mRNA".
#' @param homoeologs A data frame. Orthologous relationships of genes in the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param rpm_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 1 (the average RPM of all replicates).
#' @details The 'homoeologs' table contains the orthologs pairs. In detail, the first column is the group name (unique) of homoeologs among three species (Parents: P1; P2, Progeny: F1), the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column and fifth columns are the identifier of F1 orthologs derived from P1 and P2 ancestors, respectively (e.g. "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1").
#' @return Venn diagram.
#' @export
#' @examples
#' #miRNA
#' VennPlot(P1_RPM = P1_miRNA_rpm,
#'          P2_RPM = P2_miRNA_rpm,
#'          F1_RPM = F1_miRNA_rpm,
#'          P1_name = "B.napus(AACC)",
#'          P2_name = "B.rapa(AA)",
#'          F1_name = "B.napus x B.rapa(AAAACC)",type="sRNA")
VennPlot <- function(P1_RPM,P2_RPM,F1_RPM,P1_name,P2_name,F1_name,type,homoeologs,rpm_threshold = 1){
  if(ncol(P1_RPM) == ncol(P2_RPM) & ncol(P1_RPM) == ncol(F1_RPM)){
    colnum = ncol(P2_RPM)
    sample_number = colnum-1
  }else{
    message("Error!!! Inconsistent biological replicates between different samples.")
  }
  ############################
  func_mirnafiliter <- function(data1,threshold){
    Average_rpm <- as.data.frame(apply(data1[,c(-1)],1,mean))
    data1_value <- Average_rpm >= threshold
    data2 <- data1[data1_value,]
    result <- as.vector(as.matrix(data2[,1]))
    return(result)
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
    P1_mirna <- func_mirnafiliter(data1 = P1_RPM, threshold = rpm_threshold)
    P2_mirna <- func_mirnafiliter(data1 = P2_RPM, threshold = rpm_threshold)
    F1_mirna <- func_mirnafiliter(data1 = F1_RPM, threshold = rpm_threshold)
    x = list(P1=P1_mirna,P2=P2_mirna,F1=F1_mirna)
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
    P1_mRNA <- func_mirnafiliter(data1 = P1_gene, threshold = rpm_threshold)
    P2_mRNA <- func_mirnafiliter(data1 = P2_gene, threshold = rpm_threshold)
    F1_mRNA <- func_mirnafiliter(data1 = F1_gene, threshold = rpm_threshold)
    x = list(P1=P1_mRNA,P2=P2_mRNA,F1=F1_mRNA)
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
  grid::grid.newpage()
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.plot <- VennDiagram::venn.diagram(
    x,category.names = c(P1_name,P2_name,F1_name),
    euler.d = TRUE,filename = NULL,fontfamily = "serif",col="white",
    fill=c(colors()[616], colors()[38], colors()[468]),
    alpha=c(0.4, 0.4, 0.4),lwd=c(0.1, 0.1, 0.1),
    cat.dist = c(0.03, 0.03, 0.03), #调整分类条目位置
    cex = 2,cat.cex = 2,reverse = TRUE)
  grid::grid.draw(venn.plot)
}

#' Get the details of the Venn Diagram
#'
#' @description Get the information for each region of the venn diagram.
#'
#' @param P1_RPM A data frame. The RPM table of genes in P1 species. For the RPM table, the first column is the gene identifier, and other columns are the RPM values of the genes in each biological replicate.
#' @param P2_RPM A data frame. The RPM table of genes in P2 species.
#' @param F1_RPM A data frame. The RPM table of genes in P2 species.
#' @param type Character. "sRNA" or "mRNA".
#' @param homoeologs A data frame. Orthologous relationships of genes in the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param rpm_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 1 (the average RPM of all replicates).
#' @param output_file "venn_list", "P1_specific", "P2_specific", "F1_specific", or "all_common".
#' @details The 'homoeologs' table contains the orthologs pairs. In detail, the first column is the group name (unique) of homoeologs among three species (Parents: P1; P2, Progeny: F1), the second column is the Gene ID of P1, the third column is the Gene ID of P2. And the fourth column and fifth columns are the identifier of F1 orthologs derived from P1 and P2 ancestors, respectively (e.g. "Homoeolog1	BraA01t00004Z	BolC01g000040.2J BnA01g0000030.1 BnC01g0424620.1").
#' @return A data frame.
#' @export
#'
#' @examples
#' #output_file = "venn_list"
#' venn_list <- VennData(P1_RPM = P1_miRNA_rpm,
#'                       P2_RPM = P2_miRNA_rpm,
#'                       F1_RPM = F1_miRNA_rpm,
#'                       type="sRNA",rpm_threshold = 1,
#'                       output_file = "venn_list")
#' ##output_file = "P1_specific"
#' P1_specific <- VennData(P1_RPM = P1_miRNA_rpm,
#'                           P2_RPM = P2_miRNA_rpm,
#'                           F1_RPM = F1_miRNA_rpm,
#'                           type="sRNA",rpm_threshold = 1,
#'                           output_file = "P1_specific")
#' ##output_file = "P2_specific"
#' P2_specific <- VennData(P1_RPM = P1_miRNA_rpm,
#'                           P2_RPM = P2_miRNA_rpm,
#'                           F1_RPM = F1_miRNA_rpm,
#'                           type="sRNA",rpm_threshold = 1,
#'                           output_file = "P2_specific")
#' ##output_file = "F1_specific"
#' F1_specific <- VennData(P1_RPM = P1_miRNA_rpm,
#'                           P2_RPM = P2_miRNA_rpm,
#'                           F1_RPM = F1_miRNA_rpm,
#'                           type="sRNA",rpm_threshold = 1,
#'                           output_file = "F1_specific")
#' ##output_file = "all_common"
#' all_common <- VennData(P1_RPM = P1_miRNA_rpm,
#'                          P2_RPM = P2_miRNA_rpm,
#'                          F1_RPM = F1_miRNA_rpm,
#'                          type="sRNA",rpm_threshold = 1,
#'                          output_file = "all_common")
VennData <- function(P1_RPM,P2_RPM,F1_RPM,type,homoeologs,rpm_threshold = 1,output_file="venn_list"){
  if(ncol(P1_RPM) == ncol(P2_RPM) & ncol(P1_RPM) == ncol(F1_RPM)){
    colnum = ncol(P2_RPM)
    sample_number = colnum-1
  }else{
    message("Error!!! Inconsistent biological replicates between different samples.")
  }
  #######################
  func_mirnafiliter <- function(data1,threshold){
    Average_rpm <- as.data.frame(apply(data1[,c(-1)],1,mean))
    data1_value <- Average_rpm >= threshold
    data2 <- data1[data1_value,]
    result <- as.vector(as.matrix(data2[,1]))
    return(result)
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
    P1_mirna <- func_mirnafiliter(data1 = P1_RPM, threshold = rpm_threshold)
    P2_mirna <- func_mirnafiliter(data1 = P2_RPM, threshold = rpm_threshold)
    F1_mirna <- func_mirnafiliter(data1 = F1_RPM, threshold = rpm_threshold)
    x = list(P1=P1_mirna,P2=P2_mirna,F1=F1_mirna)
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
    P1_mRNA <- func_mirnafiliter(data1 = P1_gene, threshold = rpm_threshold)
    P2_mRNA <- func_mirnafiliter(data1 = P2_gene, threshold = rpm_threshold)
    F1_mRNA <- func_mirnafiliter(data1 = F1_gene, threshold = rpm_threshold)
    x = list(P1=P1_mRNA,P2=P2_mRNA,F1=F1_mRNA)
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
  inter <- VennDiagram::get.venn.partitions(x)
  for (i in 1:nrow(inter)){
    inter[i,'..values..'] <- paste(inter[[i,'..values..']], collapse = ',')
  }
  for (i in 1:nrow(inter)){
    if (inter[i,1]=="TRUE" & inter[i,2]=="TRUE" & inter[i,3]=="TRUE"){
      P1_P2_F1 <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(P1_P2_F1) <- "P1_P2_F1"
    }else if (inter[i,1]=="TRUE" & inter[i,2]=="FALSE" & inter[i,3]=="FALSE"){
      P1_only <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(P1_only) <- "P1_only"
    }else if (inter[i,2]=="TRUE" & inter[i,1]=="FALSE" & inter[i,3]=="FALSE"){
      P2_only <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(P2_only) <- "P2_only"
    }else if (inter[i,3]=="TRUE" & inter[i,1]=="FALSE" & inter[i,2]=="FALSE"){
      F1_only <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(F1_only) <- "F1_only"
    }else{
      next
    }
  }
  if (output_file == "all_common"){
    return(P1_P2_F1)
  }else if (output_file == "P1_specific"){
    return(P1_only)
  }else if (output_file == "P2_specific"){
    return(P2_only)
  }else if (output_file == "F1_specific"){
    return(F1_only)
  }else{
    inter$..values.. <- unlist(inter$..values..)
    colnames(inter) <- c("P1","P2","F1","set","values","count")
    return(inter)
  }
}
