#' Generate the data of miRNA base frequency in each position
#'
#' Generally, the "T" base account for the highest percentage of miRNA in the first position.The function of "mirnapredata" can provide the input data for the next drawing of miRNA base distribution in each position.
#'
#' @param miRNAseq A Character. All miRNA sequences.
#'
#' @return A dataframe. About the output results, the first column is the base, the second column is the base frequency, the third column is the position.
#' @export
#' @examples
#' #F1
#' F1_miRNA <- F1_miRNA_count[,1]
#' F1_bf <- mirnapredata(miRNAseq = F1_miRNA)
mirnapredata <- function(miRNAseq){
  max_len <- max(nchar(as.vector(as.matrix(miRNAseq))))
  mirna_list <- strsplit(as.vector(as.matrix(miRNAseq)),"")
  new_list <- lapply(mirna_list, function(x) {c(x, rep(NA, max_len - length(x)))})
  mirna_matrix <- do.call(rbind,new_list)
  res <- data.frame()
  for(i in (1:max_len)){
    data <- as.data.frame(table(mirna_matrix[,i]))
    data$Position <- i
    res <- rbind(res,data)
  }
  names(res) <- c("Base","Frequency","Position")
  return(res)
}
#' Generate the base frequency plot of miRNA
#'
#' @param miRNAdata A dataframe. The output result after running mirnapredata.
#' @param width A numeric. The width of the output bar plot, and default is 0.6.
#' @param text_size A numeric. The size of text in the output plot, and default is 10.
#' @param title_size A numeric. The size of title text in the output plot, and default is 12.
#'
#' @return The miRNA base frequency plot.
#' @export
#' @examples
#' #F1
#' F1_miRNA <- F1_miRNA_count[,1]
#' F1_bf <- mirnapredata(miRNAseq = F1_miRNA)
#' basepreplot(miRNAdata = F1_bf)
basepreplot <- function(miRNAdata,width = 0.6,text_size=10,title_size=12){
  colnames(miRNAdata) <- c("Base","Frequency","Position")
  mirna <- plyr::ddply(miRNAdata,"Position",transform,Percent = Frequency/sum(Frequency)*100)
  ggplot2::ggplot(mirna,ggplot2::aes(Position,Percent,fill=Base))+
    ggplot2::scale_x_continuous(breaks = seq(min(mirna$Position),max(mirna$Position),1))+
    ggplot2::geom_bar(stat = "identity",position = "stack",width = width)+
    ggsci::scale_fill_npg()+ggplot2::theme_classic()+ggplot2::xlab("Position")+ggplot2::ylab("Percent")+
    ggplot2::theme(axis.text = ggplot2::element_text(size = text_size,family="serif"),
                   axis.title = ggplot2::element_text(size = title_size,family="serif"),
                   legend.text = ggplot2::element_text(size = text_size,family = "serif"),
                   legend.title = ggplot2::element_text(size = title_size,family = "serif"),
                   panel.border = ggplot2::element_rect(fill=NA,color="black", size=1, linetype="solid")
    )
}
