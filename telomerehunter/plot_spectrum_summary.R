# Usage: R --no-save --slave --args <PID> <SPECTRUM_SUMMARY_DIR> <REPEAT_THRESHOLD> <CONSECUTIVE_FLAG> <MAPQ> <PLOT_FILE_FORMAT> < plot_spectrum_summary.R
# Description: Makes a bar plot of the number of telomere reads (per million reads) in the different fractions
#               of the tumor and/or control sample.

# Copyright 2018 Lina Sieverling

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
#pipeline_dir = commandArgs[5]
pid = commandArgs[5]
spectrum_dir = paste0(commandArgs[6], "/", pid, "/")
repeat_threshold = commandArgs[7]
consecutive_flag = commandArgs[8]
mapq_threshold = commandArgs[9]
plot_reverse_complement = commandArgs[10]
plot_file_format = commandArgs[11]

library(reshape2, quietly=TRUE, warn.conflicts=FALSE)
library(dplyr, quietly=TRUE, warn.conflicts=FALSE)
library(cowplot, quietly=TRUE, warn.conflicts=FALSE)
library(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE)

# source(file.path(pipeline_dir, "/functions_for_plots.R"))

if (consecutive_flag == "True"){
  count_type = "consecutive"
}else{
  count_type = "non-consecutive"
}

if (plot_reverse_complement == "True"){
  plot_reverse_complement = TRUE
}else{
  plot_reverse_complement = FALSE
}

if (plot_file_format=="all"){
  plot_file_format=c("pdf", "png", "svg")
}


spectrum_tumor_file = paste0(spectrum_dir, "/tumor_TelomerCnt_", pid,"/", pid, "_spectrum.tsv")
spectrum_control_file = paste0(spectrum_dir, "/control_TelomerCnt_", pid,"/", pid, "_spectrum.tsv")

samples = c()
samples_short = c()
if (file.exists(spectrum_tumor_file)){samples = c(samples, "tumor"); samples_short = c(samples_short, "T")}
if (file.exists(spectrum_control_file)){samples = c(samples, "control"); samples_short = c(samples_short, "C")}

col_repeat_types = c("#CB181D", "#FF656A", "#74C476", "#C1FFC3", "#2171B5", "#6EBEFF", "#FFA500", "#FFF24D", "#9370DB", "#E0BDFF", "#000000")
names(col_repeat_types) = c("TTAGGG", "CCCTAA", "TGAGGG", "CCCTCA", "TCAGGG", "CCCTGA", "TTGGGG", "CCCCAA",  "TTCGGG", "CCCGAA", "other")

####################################################################################################################################

spectrum_summary_norm_list = list()

for (sample in samples){
  read_count_file = paste0(spectrum_dir, "/",sample, "_TelomerCnt_", pid, "/", pid, "_readcount.tsv")
  
  #get spectrum
  spectrum = read.table(paste0(spectrum_dir, "/",sample, "_TelomerCnt_", pid, "/", pid, "_spectrum.tsv"), header=TRUE, comment.char="")
  
  # make empty spectrum summary tables
  spectrum_summary = matrix(data=0, nrow=4, ncol=(ncol(spectrum)-2), byrow=TRUE)
  row.names(spectrum_summary)=c("intra_chromosomal", "subtelomeric", "junction_spanning", "intra_telomeric")
  colnames(spectrum_summary) = colnames(spectrum)[3:ncol(spectrum)]
  spectrum_summary = data.frame(spectrum_summary)
  
  
  # get intra-telomeric reads (= unmapped reads)
  spectrum_summary["intra_telomeric",] = spectrum[spectrum$chr=="unmapped",3:ncol(spectrum)]
  
  
  for (chr in c(1:22, "X", "Y")){
    spectrum_chr = spectrum[spectrum$chr==chr,3:ncol(spectrum)]
    
    # get junction spanning reads
    spectrum_summary["junction_spanning", ] = spectrum_summary["junction_spanning", ] + spectrum_chr[1,] + spectrum_chr[nrow(spectrum_chr),]
    
    # get subtelomeric reads (= first and last band)
    spectrum_summary["subtelomeric", ] = spectrum_summary["subtelomeric", ] + spectrum_chr[2,] + spectrum_chr[nrow(spectrum_chr)-1,]
    
    # get intra-chromosomal reads (= all other bands) 
    spectrum_summary["intra_chromosomal", ] = spectrum_summary["intra_chromosomal", ] + apply(spectrum_chr[3:(nrow(spectrum_chr)-2),], 2, sum)   
    
  }
  
  
  # get total number of reads 
  readcount = read.table(read_count_file, header=TRUE, comment.char="")
  total_reads = sum(as.numeric(readcount$reads))
  
  
  #normalize
  spectrum_summary_norm = spectrum_summary[,2:ncol(spectrum_summary)] * (spectrum_summary$reads_with_pattern / apply(spectrum_summary[ ,2:ncol(spectrum_summary)], 1,sum))
  spectrum_summary_norm = spectrum_summary_norm * (1000000 / total_reads)
  
  
  #bring into correct format
  spectrum_summary_norm$region = row.names(spectrum_summary_norm)
  
  spectrum_summary_norm.m = melt(spectrum_summary_norm, id.vars = "region", variable.name="repeat_type")
  
  spectrum_summary_norm.m$sample = sample
  
  spectrum_summary_norm_list[[sample]] = spectrum_summary_norm.m
}


spectrum_summary = do.call("rbind", spectrum_summary_norm_list)


#combine forward and reverse patterns
if(!plot_reverse_complement){
  spectrum_summary = spectrum_summary %>%
    mutate(repeat_type = as.character(repeat_type)) %>%
    mutate(repeat_type_forward = ifelse(grepl("GGG", repeat_type) | repeat_type=="other", 
                                        repeat_type,
                                        sapply(lapply(strsplit(chartr("ATGC","TACG",repeat_type), NULL), rev), paste, collapse=""))) %>%
    group_by(sample, region, repeat_type_forward) %>%
    summarise(sum = sum(value)) %>%
    rename(repeat_type = repeat_type_forward)
}else{
  spectrum_summary = spectrum_summary %>%
    rename(sum = value)
}



#change ordering of samples and region
spectrum_summary$sample = factor(spectrum_summary$sample, 
                                 levels=c("tumor", "control"))
spectrum_summary$region = factor(spectrum_summary$region, 
                                 levels=c("intra_chromosomal", "subtelomeric",
                                          "junction_spanning", "intra_telomeric"))
levels(spectrum_summary$region) = c("intrachromosomal", "subtelomeric", "junction spanning", "intratelomeric")


#change ordering of patterns and add colors
repeats_no_color = unique(as.character(spectrum_summary$repeat_type[!spectrum_summary$repeat_type %in% names(col_repeat_types)]))
n_repeats_no_color = length(repeats_no_color)
spectrum_summary$repeat_type = factor(spectrum_summary$repeat_type,
                            levels = c(names(col_repeat_types)[names(col_repeat_types)!="other"], 
                                       repeats_no_color,
                                       "other"))

palette = colorRampPalette(brewer.pal(9, "Greys"))

col_repeat_types[repeats_no_color] = rev(palette(n_repeats_no_color+2)[-c(1,n_repeats_no_color+2)])



#title
if (nchar(pid)<11){
  main = paste0(pid,": Repeat types in telomeric reads")
}else{
  main = paste0(pid,":\nRepeat types in telomeric reads")
}


#replace NaN with 0
spectrum_summary = spectrum_summary %>%
  mutate(sum = ifelse(is.nan(sum), 0, sum))


#plot
p = ggplot(spectrum_summary, aes(sample, sum, fill=repeat_type)) +
  geom_bar(stat="identity") +
  facet_grid(.~region) +
  scale_fill_manual(name = paste(repeat_threshold, count_type, "repeats\nmapq threshold =", mapq_threshold),
                    values=col_repeat_types) +
  theme(axis.title.x = element_blank()) +
  ylab("Telomeric reads (per million reads)") +
  ggtitle(main)


#save
plot_dir = file.path(spectrum_dir, "plots")

if (!(file.exists(plot_dir))){
  dir.create(plot_dir, recursive=TRUE)
}

plot_file_prefix = paste0(plot_dir,"/", pid, "_sorted_telomere_read_counts")

for (plot_type in plot_file_format){
  ggsave(paste0(plot_file_prefix, ".", plot_type), p, width=10, height=5)
}

remove = file.remove("Rplots.pdf")