#' Length distribution of small RNAs (sRNAs)
#'
#' @description Get the length distribution of sRNAs.
#' @param sRNAseq Character. All sRNA sequences in vector format.
#' @param Group Character. Group name.
#'
#' @return A data frame. The output consists of three columns, i.e., length, frequency and group name.
#' @export
#' @examples
#' #P1(B.napus)
#' B.napu_sRNA <- srnapredata(sRNAseq = P1_sRNA_seq, Group = "B.napus(AACC)")
#' #P2(B.rapa)
#' B.rapa_sRNA <- srnapredata(sRNAseq = P2_sRNA_seq, Group = "B.rapa(AA)")
#' #F1(B.napus X B.rapa)
#' B.nr_sRNA <- srnapredata(sRNAseq = F1_sRNA_seq, Group = "B.napus x B.rapa(AAAACC)")
#' #intergrate these data for length distribution plot
#' sRNA_data <- rbind(B.napu_sRNA, B.rapa_sRNA, B.nr_sRNA)
#' #output result
#' head(sRNA_data)
#' #  Length Frequency         Group
#' #1     15         8 B.napus(AACC)
#' #2     16         7 B.napus(AACC)
#' #3     17        13 B.napus(AACC)
#' #4     18        16 B.napus(AACC)
#' #5     19        25 B.napus(AACC)
#' #6     20        33 B.napus(AACC)

srnapredata <- function(sRNAseq,Group){
  length <- nchar(as.vector(as.matrix(sRNAseq)))
  dataset1 <- as.data.frame(table(length))
  dataset1$group <- Group
  colnames(dataset1) <- c("Length","Frequency","Group")
  return(dataset1)
}

#' @title Plot the length distribution diagram for small RNAs (sRNAs)
#' @description There are two types of pictures: bar plot (type = "bar") and line plot (type = "line"). For the bar plot, the Y-axis displays the proportion of sRNAs in a certain length, the X-axis represents sRNAs in different length. And for line plot, the Y-axis displays the abundance of sRNAs in a certain length, the X-axis represents sRNAs in different length.
#'
#' @param sRNAdata A data frame. Frequency distribution of sRNAs in different length.
#' @param type A character. "bar" or "line".
#' @param width A numeric. Bar width, and default is 0.6. if the type is "line", the parameter does not need to be given.
#' @param font_size A numeric. Size of axis ticks and legend item labels, and default is 10.
#' @param title_size A numeric. Size of axis titles and legend titles, and default is 12.
#'
#' @return Length distribution plot of sRNAs.
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
lenplot <- function(sRNAdata,type,width=0.6,font_size=10,title_size=12){
  colnames(sRNAdata) <- c("Length","Frequency","Group")
  if (type == "bar"){
    srna <- plyr::ddply(sRNAdata,"Group",transform,Percent = Frequency/sum(Frequency)*100)
    ggplot2::ggplot(srna,ggplot2::aes(Length,Percent,fill=Group))+
      ggplot2::geom_bar(stat = "identity",position = "dodge",width = width)+
      ggsci::scale_fill_npg()+ggplot2::theme_classic()+ggplot2::xlab("Position")+ggplot2::ylab("Percent")+
      ggplot2::theme(axis.text = ggplot2::element_text(size = font_size,family="serif"),
                     axis.title = ggplot2::element_text(size = title_size,family="serif"),
                     legend.text = ggplot2::element_text(size = font_size,family = "serif"),
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
      ggplot2::theme(axis.text = ggplot2::element_text(size = font_size,family="serif"),
                     axis.title = ggplot2::element_text(size = title_size,family="serif"),
                     legend.text = ggplot2::element_text(size = font_size,family = "serif"),
                     legend.title = ggplot2::element_text(size = title_size,family = "serif"),
                     axis.line = ggplot2::element_line(colour = "black",size=1))
  }else{
    message("Eorr! Please select the pattern of output image using the option of type (type = bar or type = line).")
  }
}

utils::globalVariables(
  c("Length","Frequency","Group","Percent","Base","Position","colors","na.omit","write.csv")
)
