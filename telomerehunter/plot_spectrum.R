# Usage: R --no-save --slave --args <FUNCTION_DIR> <PID> <SPECTRUM_SUMMARY_DIR> <REPEAT_THRESHOLD> <CONSECUTIVE_FLAG> <MAPQ> <BANDING_FILE> <PLOT_FILE_FORMAT> < plot_spectrum.R
# Description: Makes a bar plot for each chromosome. The chromosome bands are plotted on the x-axis and the specific telomere reads per million bases 
#              and per billion reads are plotted on the y-axis. Attention: the "junction bands" are not divided by band length (for these bands
#              the y-axis is "specific telomere reads per billion reads"), so the scaling is not correct in comparison to "normal" bands!!! 

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

# get commandline arguments
commandArgs = commandArgs()
pipeline_dir = commandArgs[5]
pid = commandArgs[6]
spectrum_dir = paste0(commandArgs[7], "/", pid, "/")
repeat_threshold = commandArgs[8]
consecutive_flag = commandArgs[9]
mapq_threshold = commandArgs[10]
banding_file = commandArgs[11]
plot_reverse_complement = commandArgs[12]
plot_file_format = commandArgs[13]


if (consecutive_flag == "True"){
  count_type = "Consecutive"
}else{
  count_type = "Non-consecutive"
}

if (plot_reverse_complement == "True"){
  plot_reverse_complement = TRUE
}else{
  plot_reverse_complement = FALSE
}

if (plot_file_format=="all"){
  plot_file_format=c("pdf", "png", "svg")
}


#------------------------------------------------------------------------------------------------------------------------

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




#------------------------------------------------------------------------------------------------------------------------

# get samples
spectrum_tumor_file = paste0(spectrum_dir, "/tumor_TelomerCnt_", pid,"/", pid, "_spectrum.tsv")
spectrum_control_file = paste0(spectrum_dir, "/control_TelomerCnt_", pid,"/", pid, "_spectrum.tsv")

samples = c()
if (file.exists(spectrum_tumor_file)){samples = c(samples, "tumor"); names(samples)[samples=="tumor"]="Tumor"}
if (file.exists(spectrum_control_file)){samples = c(samples, "control"); names(samples)[samples=="control"]="Control"}


# get band lengths
band_info = read.table(banding_file)
colnames(band_info) = c("chr", "start", "end", "band_name", "stain")
band_info[,"chr"] = gsub("chr", "", band_info[,1])
band_info[,"length"] = band_info[,"end"] - band_info[,"start"]

# get plot directory
plot_dir = file.path(spectrum_dir, "plots")
if (!(file.exists(plot_dir))){dir.create(plot_dir, recursive=TRUE)}



spectrum_list = list()

for (sample in samples){
  #get spectrum
  spectrum = read.table(paste0(spectrum_dir, "/",sample, "_TelomerCnt_", pid, "/", pid, "_spectrum.tsv"), header=TRUE, comment.char="")
 
  # get total number of reads 
  read_count_file = paste0(spectrum_dir, "/",sample, "_TelomerCnt_", pid, "/", pid, "_readcount.tsv")
  readcount = read.table(read_count_file, header=TRUE, comment.char="")
  total_reads = sum(as.numeric(readcount$reads))
    
  #normalize
  spectrum[ ,4:ncol(spectrum)] = spectrum[ ,4:ncol(spectrum)] * (spectrum$reads_with_pattern / apply(spectrum[ ,4:ncol(spectrum)], 1,sum))
  spectrum[ ,4:ncol(spectrum)] = spectrum[ ,4:ncol(spectrum)]*1000000000 / total_reads
  
  spectrum_list[[sample]] = spectrum
}
  
  
for (chr in c(1:22, "X", "Y")){
  
  spectrum_chr_list =list()
  
  for (sample in samples){
    spectrum = spectrum_list[[sample]]
    spectrum_chr = spectrum[spectrum$chr==chr,]
        
    #normalize by band length
    band_info_chr = band_info[band_info$chr==chr,]
    spectrum_chr[2:(nrow(spectrum_chr)-1),3:ncol(spectrum_chr)] = spectrum_chr[2:(nrow(spectrum_chr)-1),3:ncol(spectrum_chr)] *(1000000/ band_info_chr$length) 
    
    height = t(spectrum_chr[ ,4:ncol(spectrum_chr)])
    spectrum_chr_list[[sample]] = spectrum_chr
  }
  
  bands=spectrum_chr_list[[1]][ , "band"]
  
  if (length(samples)==2){
    height = zipFastener(t(spectrum_chr_list[["tumor"]][ ,4:ncol(spectrum_chr)]), t(spectrum_chr_list[["control"]][ ,4:ncol(spectrum_chr)]))
    sub_title = "Left: Tumor, Right: Control"
    axis=TRUE
    axis_simple=FALSE
  }else{
    height = t(spectrum_chr_list[[1]][ ,4:ncol(spectrum_chr)])
    sub_title = paste0(names(samples), " Sample")
    axis=FALSE
    axis_simple=TRUE
  }
  
  colnames(height)=NULL
  
  main = paste0(pid,": Telomere Repeat Types in Chr", chr)
  plot_file_prefix = paste0(plot_dir,"/", pid, "_", chr)
  
  barplot_repeattype(height=height,
                     plot_file_prefix=plot_file_prefix, plot_file_format=plot_file_format, width=28, mar=c(5.1, 5.1, 5.1, 9.3),
                     main=main, ylab="Normalized Number of Telomere Reads", #bands: Telomere Reads (per Million Bases and per Billion Reads), junctions: Telomere Reads (per Billion Reads)
                     plot_reverse_complement=plot_reverse_complement,
                     repeat_threshold=repeat_threshold, count_type=count_type, mapq_threshold=mapq_threshold,
                     inset_legend=c(-0.18,0),
                     axis=axis, axis_simple=axis_simple, labels=bands, cex.axis=0.75, cex.lab=0.9, sub_title=sub_title) 
}

  
  
  