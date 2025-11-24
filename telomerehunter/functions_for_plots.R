# Description: contains functions for plotting scripts of TelomereHunter

# Copyright 2015 Lina Sieverling, Philip Ginsbach, Lars Feuerbach

# This file is part of TelomereHunter.

# TelomereHunter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TelomereHunter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TelomereHunter.  If not, see <http://www.gnu.org/licenses/>.



library(ggplot2, quietly=TRUE, warn.conflicts=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE)
library(grid, quietly=TRUE, warn.conflicts=FALSE)
library(gridExtra, quietly=TRUE, warn.conflicts=FALSE)
library(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE)


# zipFastener for TWO dataframes of unequal length
zipFastener <- function(df1, df2, along=2)
{
  # parameter checking
  if(!is.element(along, c(1,2))){
    stop("along must be 1 or 2 for rows and columns
         respectively")
  }
  # if merged by using zip feeding along the columns, the
  # same no. of rows is required and vice versa
  if(along==1 & (ncol(df1)!= ncol(df2))) {
    stop ("the no. of columns has to be equal to merge
          them by zip feeding")
    }
    if(along==2 & (nrow(df1)!= nrow(df2))) {
        stop ("the no. of rows has to be equal to merge them by
              zip feeding")
    }

    # zip fastener preperations
    d1 <- dim(df1)[along]
    d2 <- dim(df2)[along]
    i1 <- 1:d1           # index vector 1
    i2 <- 1:d2 + d1      # index vector 2

    # set biggest dimension dMax
    if(d1==d2) {
        dMax <- d1
    } else if (d1 > d2) {
        length(i2) <- length(i1)    # make vectors same length, 
        dMax <- d1                  # fill blanks with NAs   
    } else  if(d1 < d2){
        length(i1) <- length(i2)    # make vectors same length,
        dMax <- d2                  # fill blanks with NAs   
    }
    
    # zip fastener operations
    index <- as.vector(matrix(c(i1, i2), ncol=dMax, byrow=T))
    index <- index[!is.na(index)]         # remove NAs
    
    if(along==1){
        colnames(df2) <- colnames(df1)   # keep 1st colnames                  
        res <- rbind(df1,df2)[ index, ]  # reorder data frame
    }
    if(along==2) res <- cbind(df1,df2)[ , index]           

    return(res)
}


# function for making bar plots
barplot_repeattype <- function (height,
                                plot_file_prefix, plot_file_format=c("pdf", "png", "svg"), width=550, mar=c(3.1, 5.1, 5.1, 9.3),
                                main, sub_title="", ylab="Telomere Reads (per Million Reads)",
                                plot_reverse_complement=FALSE,
                                repeat_threshold, count_type, mapq_threshold,
                                inset_legend=c(-0.34,0),    #legend.text=row.names(height)
                                sample_label=FALSE, sample_label_text=NA,
                                axis=TRUE, axis_simple=FALSE, labels=NA, xlas=2, tick=TRUE, cex.axis=1, cex.lab=1, cex.main=1.2, cex.names=1) {

  if (plot_reverse_complement==FALSE){
    height_patterns = as.matrix(height[1:(dim(height)[1]-1), ])
    height_summed = rowsum(height_patterns, group=rep(row.names(height_patterns)[seq(1, dim(height_patterns)[1], by=2)], each=2), reorder=FALSE)   
    height_summed = rbind(height_summed, height[dim(height)[1],])
    row.names(height_summed)[dim(height_summed)[1]] = "other"
    height = height_summed                         
  }

  searched_patterns = row.names(height)
  
  #get colors
  col_all = c("#CB181D", "#FF656A", "#74C476", "#C1FFC3", "#2171B5", "#6EBEFF", "#FFA500", "#FFF24D", "#9370DB", "#E0BDFF", "#000000")
  names(col_all) = c("TTAGGG", "CCCTAA", "TGAGGG", "CCCTCA", "TCAGGG", "CCCTGA", "TTGGGG", "CCCCAA",  "TTCGGG", "CCCGAA", "other")
  col=c(col_all[searched_patterns])  
  palette = colorRampPalette(brewer.pal(9, "Greys"))
  col[is.na(col)] = palette(sum(is.na(col))+2)[-c(1,sum(is.na(col))+2)]
  
  for (plot_type in plot_file_format){
    
    if (plot_type=="png"){png(paste0(plot_file_prefix,".png"), width=width, height=14, units="cm", res=300)}  
    if (plot_type=="pdf"){pdf(paste0(plot_file_prefix,".pdf"), width=width*0.4, height=14*0.4)}
    if (plot_type=="svg"){svg(paste0(plot_file_prefix, ".svg"), width=width/2.56, height=14/2.56)}
        
    par(mar=mar, xpd=TRUE)
    mp = barplot(height,
                 col=col,
                 main=main,
                 sub=sub_title,
                 ylab=ylab,
                 las=xlas,
                 cex.main=cex.main,
                 cex.axis=cex.axis,
                 cex.lab=cex.lab,
                 cex.names=cex.names,
                 border=NA)
    
    legend.title = paste(repeat_threshold, count_type, "Repeats", "\nMapq Threshold =", mapq_threshold)
    legend_pos = legend("topright", inset=inset_legend, title=legend.title, legend=searched_patterns, fill=col, bty="n", cex= 0.75, border="white")  
    
    xleft = legend_pos$rect[["left"]]
    ytop = legend_pos$rect[["top"]]
    ybottom = ytop - legend_pos$rect[["h"]]
    xright = xleft + legend_pos$rect[["w"]]
    rect(xleft, ybottom, xright, ytop+(ytop-ybottom)/(1.25*dim(height)[1]))
    
    if (sample_label==TRUE){
      text(x=mp, y=0, labels=sample_label_text, pos=1, cex=0.7)
    }

    
    if (axis==TRUE){
      mp2 = (mp[seq(1,length(mp),2)] + mp[seq(2,length(mp),2)]) / 2
      axis(1, labels=labels, at=mp2, las=xlas, tick=tick, cex.axis=cex.axis)
    } else if (axis_simple==TRUE){
      axis(1, labels=labels, at=mp, las=xlas, tick=tick, cex.axis=cex.axis)
    }
    
    garbage = dev.off()   # message from dev.off is not printed
  }
}


make_plots <- function(filename,
                       plot,
                       width=18,
                       height=18,
                       plot_types=c("png", "pdf", "svg"),
                       arrange_grob=FALSE){
  
  for (plot_type in plot_types){    
    if (plot_type=="png"){png(paste0(filename, ".png"), width=width, height=height, units="cm", res=300)}  
    if (plot_type=="pdf"){pdf(paste0(filename, ".pdf"), width=width*0.4, height=height*0.4)}     
    if (plot_type=="svg"){svg(paste0(filename, ".svg"), width=width/2.56, height=height/2.56)}
    
    if (arrange_grob==TRUE){
      grid::grid.draw(plot)
    }else{
      print(plot)  
    }
    
    dev.off()
  }
}