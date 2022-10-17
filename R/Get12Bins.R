#' Non-additive expression analysis
#' @description Rapp et al. proposed the classification of 12 expression patterns in allopolyploids, including additivity (I, XII), ELD (II, XI, IV, IX), transgressive down-regulation (III, VII, X) and transgressive up-regulation (V, VI, VIII).
#' @param P1_count A data frame. The count table of genes in P1 species. For the count table, the first column is the gene identifier, and other columns are the corresponding expression levels of the genes in each biological replicate.
#' @param P2_count A data frame. The count table of genes in P2 species.
#' @param F1_count A data frame. The count table of genes in F1 species.
#' @param type A character. "sRNA" or "mRNA".
#' @param homoeologs A data frame. Orthologous relationships of genes in the parental species and their progeny. Only required when the 'type' is 'mRNA'.
#' @param count_threshold A numeric. Threshold for filtering out the lowly expressed genes. The default is 5 (the count values in all replicates).
#' @param Pvalue A numeric. The P value of differential expression analysis using DESeq2. Default is 0.05.
#' @param log2FC A numeric. The log2-transformed expression fold of differential expression analysis using DESeq2. Default is 1.
#' @references Rapp RA, Udall JA, Wendel JF. Genomic expression dominance in allopolyploids. BMC Biol. 2009 May 1;7:18.
#' @return A data frame. Classification results of non-additive analysis based on the ELD method.
#' @export
#' @details pv11: P value of differential expression analysis using DESeq2. Parental P1 was used as the control group and F1 was used as the treatment group. pv12: P value of differential expression analysis using DESeq2. Parental P2 was used as the control group and F1 was used as the treatment group. pv21: P value of differential expression analysis using DESeq2. Parental P1 was used as the control group and P2 was used as the treatment group. Besides, "fc" represents the log2FoldChange of differential expression analysis.
#' @examples
#' \donttest{
#' miRNA_12bin <- Get12Bins(P1_count = P1_miRNA_count,
#'                          P2_count = P2_miRNA_count,
#'                          F1_count = F1_miRNA_count,type = "sRNA")
#'}
Get12Bins <- function(P1_count,P2_count,F1_count,type,homoeologs,count_threshold = 5,Pvalue = 0.05,log2FC = 1){
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
  new_F1_vs_P1 <- cbind(as.data.frame(row.names(F1_vs_P1)),F1_vs_P1$pvalue,F1_vs_P1$log2FoldChange)
  new_F1_vs_P2 <- cbind(as.data.frame(row.names(F1_vs_P2)),F1_vs_P2$pvalue,F1_vs_P2$log2FoldChange)
  new_P2_vs_P1 <- cbind(as.data.frame(row.names(P2_vs_P1)),P2_vs_P1$pvalue,P2_vs_P1$log2FoldChange)
  if(type=="sRNA"){
    names(new_F1_vs_P1) <- c("sequence","pv11","fc11")
    names(new_F1_vs_P2) <- c("sequence","pv12","fc12")
    names(new_P2_vs_P1) <- c("sequence","pv21","fc21")
  }else if(type=="mRNA"){
    names(new_F1_vs_P1) <- c("orthogroup","pv11","fc11")
    names(new_F1_vs_P2) <- c("orthogroup","pv12","fc12")
    names(new_P2_vs_P1) <- c("orthogroup","pv21","fc21")
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
  #####################
  SeGedata <- merge.data.frame(merge.data.frame(new_F1_vs_P1,new_F1_vs_P2,all=TRUE),new_P2_vs_P1,all=TRUE)
  SeGedata <- na.omit(SeGedata)
  flg2FC <- -log2FC
  c_i <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11>=log2FC & SeGedata$fc12<=flg2FC & SeGedata$fc21>=log2FC),]
  c_xii <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11<=flg2FC & SeGedata$fc12>=log2FC & SeGedata$fc21<=flg2FC),]
  c_ii <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12>=Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11>=log2FC & abs(SeGedata$fc12)<log2FC & SeGedata$fc21>=log2FC),]
  c_xi <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12>=Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11<=flg2FC & abs(SeGedata$fc12)<log2FC & SeGedata$fc21<=flg2FC),]
  c_iv <- SeGedata[which(SeGedata$pv11>=Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & abs(SeGedata$fc11)<log2FC & SeGedata$fc12>=log2FC & SeGedata$fc21<=flg2FC),]
  c_ix <- SeGedata[which(SeGedata$pv11>=Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & abs(SeGedata$fc11)<log2FC & SeGedata$fc12<=flg2FC & SeGedata$fc21>=log2FC),]
  c_iii <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11<=flg2FC & SeGedata$fc12<=flg2FC & SeGedata$fc21>=log2FC),]
  c_vii <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21>=Pvalue & SeGedata$fc11<=flg2FC & SeGedata$fc12<=flg2FC & abs(SeGedata$fc21)<log2FC),]
  c_x <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11<=flg2FC & SeGedata$fc12<=flg2FC & SeGedata$fc21<=flg2FC),]
  c_v <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11>=log2FC & SeGedata$fc12>=log2FC & SeGedata$fc21>=log2FC),]
  c_vi <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21<Pvalue & SeGedata$fc11>=log2FC & SeGedata$fc12>=log2FC & SeGedata$fc21<=flg2FC),]
  c_viii <- SeGedata[which(SeGedata$pv11<Pvalue & SeGedata$pv12<Pvalue & SeGedata$pv21>=Pvalue & SeGedata$fc11>=log2FC & SeGedata$fc12>=log2FC & abs(SeGedata$fc21)<log2FC),]
  #Additivity
  if(nrow(c_i)==0){
    c_i <- data.frame()
  }else{
    c_i$Categories <- "I"
    c_i$Group <- "Additivity"
  }
  if(nrow(c_xii)==0){
    c_xii <- data.frame()
  }else{
    c_xii$Categories <- "XII"
    c_xii$Group <- "Additivity"
  }
  #P2-expression Level dominance
  if(nrow(c_ii)==0){
    c_ii <- data.frame()
  }else{
    c_ii$Categories <- "II"
    c_ii$Group <- "P2-expression level dominance"
  }
  if(nrow(c_xi)==0){
    c_xi <- data.frame()
  }else{
    c_xi$Categories <- "XI"
    c_xi$Group <- "P2-expression level dominance"
  }
  #P1-expression level dominance
  if(nrow(c_iv)==0){
    c_iv <- data.frame()
  }else{
    c_iv$Categories <- "IV"
    c_iv$Group <- "P1-expression level dominance"
  }
  if(nrow(c_ix)==0){
    c_ix <- data.frame()
  }else{
    c_ix$Categories <- "IX"
    c_ix$Group <- "P1-expression level dominance"
  }
  #Transgressive down-regulation
  if(nrow(c_iii)==0){
    c_iii <- data.frame()
  }else{
    c_iii$Categories <- "III"
    c_iii$Group <- "Transgressive down-regulation"
  }
  if(nrow(c_vii)==0){
    c_vii <- data.frame()
  }else{
    c_vii$Categories <- "VII"
    c_vii$Group <- "Transgressive down-regulation"
  }
  if(nrow(c_x)==0){
    c_x <- data.frame()
  }else{
    c_x$Categories <- "X"
    c_x$Group <- "Transgressive down-regulation"
  }
  #Transgressive up-regulation
  if(nrow(c_v)==0){
    c_v <- data.frame()
  }else{
    c_v$Categories <- "V"
    c_v$Group <- "Transgressive up-regulation"
  }
  if(nrow(c_vi)==0){
    c_vi <- data.frame()
  }else{
    c_vi$Categories <- "VI"
    c_vi$Group <- "Transgressive up-regulation"
  }
  if(nrow(c_viii)==0){
    c_viii <- data.frame()
  }else{
    c_viii$Categories <- "VIII"
    c_viii$Group <- "Transgressive up-regulation"
  }
  ##################
  result <- rbind(c_i,c_xii,c_ii,c_xi,c_iv,c_ix,c_iii,c_vii,c_x,c_v,c_vi,c_viii)
  if(type=="sRNA"){
    eld <- result$sequence
    nc <- SeGedata[!(SeGedata$sequence %in% eld),]
    if(nrow(nc)==0){
      nc <- data.frame()
    }else{
      nc$Categories <- "No Change"
      nc$Group <- "No change"
    }
    fina_result <- rbind(result,nc)
    row.names(fina_result) <- 1:nrow(fina_result)
    return(fina_result)
  }else if(type=="mRNA"){
    eld <- result$orthgroup
    nc <- SeGedata[!(SeGedata$orthgroup %in% eld),]
    if(nrow(nc)==0){
      nc <- data.frame()
    }else{
      nc$Categories <- "No Change"
      nc$Group <- "No change"
    }
    fina_result <- rbind(result,nc)
    row.names(fina_result) <- 1:nrow(fina_result)
    return(fina_result)
  }else{
    message("Error!!! Please select the type option of 'mRNA' or 'sRNA'.")
  }
}
