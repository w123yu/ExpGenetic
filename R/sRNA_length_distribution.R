#' Generate the data of sRNA length distribution
#'
#' Generally, the length interval of sRNA is 21-24. The function of "srnapredata" can provide the input data for the next drawing of sRNA length distribution among different species.
#'
#' @param sRNAseq Character. All sRNA sequences.
#' @param Group A character. You an select a representative group name for next drawing.
#'
#' @return A dataframe. The output results are consist of three columns, the first column is the length of sRNA, the second column id the frequency, and the third column is the group name.
#' @export
#' @examples
#' #P1(B.napus)
#' B.napu_sRNA <- srnapredata(sRNAseq = P1_sRNA_seq,Group = "B.napus(AACC)")
#' #P2(B.rapa)
#' B.rapa_sRNA <- srnapredata(sRNAseq = P2_sRNA_seq,Group = "B.rapa(AA)")
#' #F1(B.napus X B.rapa)
#' B.nr_sRNA <- srnapredata(sRNAseq = F1_sRNA_seq,Group = "B.napus x B.rapa(AAAACC)")
#' #intergrate these data for length distribution plot
#' sRNA_data <- rbind(B.napu_sRNA,B.rapa_sRNA,B.nr_sRNA)
srnapredata <- function(sRNAseq,Group){
  length <- nchar(as.vector(as.matrix(sRNAseq)))
  dataset1 <- as.data.frame(table(length))
  dataset1$group <- Group
  colnames(dataset1) <- c("Length","Frequency","Group")
  return(dataset1)
}
#' Generate the sRNA length distribution plot
#'
#' @param sRNAdata A dataframe. The output result after running "srnapredata".
#' @param type A character. "bar" or "line", and specific option will generate specific image.
#' @param width A numeric. The width of the output bar plot, and default is 0.6.
#' @param text_size A numeric. The size of text in the output plot, and default is 10.
#' @param title_size A numeric. The size of title text in the output plot, and default is 12.
#'
#' @return The sRNA length distribution plot.
#' @export
#' @examples
#' #P1(B.napus)
#' B.napu_sRNA <- srnapredata(sRNAseq = P1_sRNA_seq,Group = "B.napus(AACC)")
#' #P2(B.rapa)
#' B.rapa_sRNA <- srnapredata(sRNAseq = P2_sRNA_seq,Group = "B.rapa(AA)")
#' #F1(B.napus X B.rapa)
#' B.nr_sRNA <- srnapredata(sRNAseq = F1_sRNA_seq,Group = "B.napus x B.rapa(AAAACC)")
#' #intergrate these data for length distribution plot
#' sRNA_data <- rbind(B.napu_sRNA,B.rapa_sRNA,B.nr_sRNA)
#' #plot
#' lenplot(sRNAdata = sRNA_data,type = "line")
#' lenplot(sRNAdata = sRNA_data,type = "bar")
lenplot <- function(sRNAdata,type,width=0.6,text_size=10,title_size=12){
  colnames(sRNAdata) <- c("Length","Frequency","Group")
  if (type == "bar"){
    srna <- plyr::ddply(sRNAdata,"Group",transform,Percent = Frequency/sum(Frequency)*100)
    ggplot2::ggplot(srna,ggplot2::aes(Length,Percent,fill=Group))+
      ggplot2::geom_bar(stat = "identity",position = "dodge",width = width)+
      ggsci::scale_fill_npg()+ggplot2::theme_classic()+ggplot2::xlab("Position")+ggplot2::ylab("Percent")+
      ggplot2::theme(axis.text = ggplot2::element_text(size = text_size,family="serif"),
                     axis.title = ggplot2::element_text(size = title_size,family="serif"),
                     legend.text = ggplot2::element_text(size = text_size,family = "serif"),
                     legend.title = ggplot2::element_text(size = title_size,family = "serif"),
                     panel.border = ggplot2::element_rect(fill=NA,color="black", size=1, linetype="solid")
      )
  }else if(type == "line"){
    sRNAdata$Length <- as.numeric(as.character(sRNAdata$Length))
    options(scipen = 200)
    ggplot2::ggplot(sRNAdata,ggplot2::aes(x=Length,y=Frequency,colour = Group))+
      ggplot2::geom_line(size=1)+ggplot2::geom_point(size=2)+
      ggplot2::xlab("Length")+ggplot2::ylab("Frequency")+
      ggplot2::theme_classic()+ggplot2::scale_x_continuous(breaks = sRNAdata$Length)+
      ggplot2::theme(axis.text = ggplot2::element_text(size = text_size,family="serif"),
                     axis.title = ggplot2::element_text(size = title_size,family="serif"),
                     legend.text = ggplot2::element_text(size = text_size,family = "serif"),
                     legend.title = ggplot2::element_text(size = title_size,family = "serif"),
                     axis.line = ggplot2::element_line(colour = "black",size=1))
  }else{
    cat("Eorr! Please select the pattern of output image using the option of type (type = bar or type = line).")
  }
}

utils::globalVariables(
  c("Length","Frequency","Group","Percent","Base","Position","colors","na.omit","write.csv")
)
