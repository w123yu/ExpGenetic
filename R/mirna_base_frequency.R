#' @title Base frequency distribution of small RNA (sRNA)
#' @description Get the base frequency distribution table.
#' @param sRNAseq Character. All sRNA sequences in vector format.
#' @return A data frame. The output consists of three columns, i.e., base, base frequency and position.
#' @export
#' @examples
#' #F1
#' F1_miRNA <- F1_miRNA_count[,1]
#' F1_bf <- mirnapredata(sRNAseq = F1_miRNA)
#' #output result
#' head(F1_bf)
#' #  Base Frequency Position
#' #1    A        32        1
#' #2    C        27        1
#' #3    G        31        1
#' #4    T       115        1
#' #5    A        27        2
#' #6    C        50        2

mirnapredata <- function(sRNAseq){
  max_len <- max(nchar(as.vector(as.matrix(sRNAseq))))
  mirna_list <- strsplit(as.vector(as.matrix(sRNAseq)),"")
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

#' @title Plot the base frequency distribution diagram for small RNA (sRNA)
#' @param sRNAdata A data frame. Base frequency distribution of sRNAs.
#' @param width A numeric. Bar width, and default is 0.6.
#' @param font_size A numeric. Size of axis ticks and legend item labels, and default is 10.
#' @param title_size A numeric. Size of axis titles and legend titles, and default is 12.
#'
#' @return Base frequency distribution plot of sRNAs.
#' @export
#' @examples
#' #F1
#' F1_miRNA <- F1_miRNA_count[,1]
#' F1_bf <- mirnapredata(sRNAseq = F1_miRNA)
#' basepreplot(sRNAdata = F1_bf)
basepreplot <- function(sRNAdata,width = 0.6,font_size=10,title_size=12){
  colnames(sRNAdata) <- c("Base","Frequency","Position")
  mirna <- plyr::ddply(sRNAdata,"Position",transform,Percent = Frequency/sum(Frequency)*100)
  ggplot2::ggplot(mirna,ggplot2::aes(Position,Percent,fill=Base))+
    ggplot2::scale_x_continuous(breaks = seq(min(mirna$Position),max(mirna$Position),1))+
    ggplot2::geom_bar(stat = "identity",position = "stack",width = width)+
    ggsci::scale_fill_npg()+ggplot2::theme_classic()+ggplot2::xlab("Position")+ggplot2::ylab("Percent")+
    ggplot2::theme(axis.text = ggplot2::element_text(size = font_size,family="serif"),
                   axis.title = ggplot2::element_text(size = title_size,family="serif"),
                   legend.text = ggplot2::element_text(size = font_size,family = "serif"),
                   legend.title = ggplot2::element_text(size = title_size,family = "serif"),
                   panel.border = ggplot2::element_rect(fill=NA,color="black", size=1, linetype="solid")
    )
}
