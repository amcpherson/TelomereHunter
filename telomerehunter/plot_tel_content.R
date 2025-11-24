# Usage: R --no-save --slave --args <PID> <SPECTRUM_SUMMARY_DIR> <REPEAT_THRESHOLD> <CONSECUTIVE_FLAG> <MAPQ> <PLOT_FILE_FORMAT> <GC_LOWER_LIMIT> <GC_UPPER_LIMIT> < plot_tel_content.R
# Description: makes a bar plot of the number of intratelomeric reads (per million reads) in the tumor and control sample.

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
pid = commandArgs[5]
spectrum_dir = paste0(commandArgs[6], "/", pid, "/")
repeat_threshold = commandArgs[7]
consecutive_flag = commandArgs[8]
mapq_threshold = commandArgs[9]
plot_reverse_complement = commandArgs[10]
plot_file_format = commandArgs[11]
gc_lower_limit = as.numeric(commandArgs[12])
gc_upper_limit = as.numeric(commandArgs[13])

library(cowplot, quietly=TRUE, warn.conflicts=FALSE)
library(dplyr, quietly=TRUE, warn.conflicts=FALSE)
library(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE)

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

col_repeat_types = c("#CB181D", "#FF656A", "#74C476", "#C1FFC3", "#2171B5", "#6EBEFF", "#FFA500", "#FFF24D", "#9370DB", "#E0BDFF", "#000000")
names(col_repeat_types) = c("TTAGGG", "CCCTAA", "TGAGGG", "CCCTCA", "TCAGGG", "CCCTGA", "TTGGGG", "CCCCAA",  "TTCGGG", "CCCGAA", "other")

############################################################################################################

spectrum_tumor_file = paste0(spectrum_dir, "/tumor_TelomerCnt_", pid,"/", pid, "_spectrum.tsv")
spectrum_control_file = paste0(spectrum_dir, "/control_TelomerCnt_", pid,"/", pid, "_spectrum.tsv")

samples = c()
if (file.exists(spectrum_tumor_file)){samples = c(samples, "tumor")}
if (file.exists(spectrum_control_file)){samples = c(samples, "control")}

heights_samples_list = list()
for (sample in samples){
  #get spectrum
  spectrum = read.table(paste0(spectrum_dir, "/", sample, "_TelomerCnt_", pid,"/", pid, "_spectrum.tsv"), header=TRUE, comment.char="")
  spectrum = spectrum[spectrum$chr=="unmapped",]
  
  #get gc corrected telomere content
  summary = read.table(paste0(spectrum_dir, "/", sample, "_TelomerCnt_", pid,"/", pid, "_", sample, "_summary.tsv"), header=TRUE, sep="\t")
  tel_content = summary$tel_content
  
  #get relative telomere repeat type occurrences in telomere content
  spectrum[ ,4:ncol(spectrum)] = spectrum[ ,4:ncol(spectrum)] / sum(spectrum[ ,4:ncol(spectrum)]) * tel_content
  
  height = t(spectrum[ ,4:ncol(spectrum)])
  colnames(height) = "value"
  
  height = as.data.frame(height)
  height$repeat_type = row.names(height)
  height$sample = sample
  
  heights_samples_list[[sample]] = height
}


height = do.call("rbind", heights_samples_list)

#combine forward and reverse patterns
if(!plot_reverse_complement){
  height = height %>%
    mutate(repeat_type_forward = ifelse(grepl("GGG", repeat_type) | repeat_type=="other", 
                                        repeat_type,
                                        sapply(lapply(strsplit(chartr("ATGC","TACG",repeat_type), NULL), rev), paste, collapse=""))) %>%
    group_by(sample, repeat_type_forward) %>%
    summarise(sum = sum(value)) %>%
    rename(repeat_type = repeat_type_forward)
}else{
  height = height %>%
    rename(sum=value)
}


#change ordering of samples
height$sample = factor(height$sample, 
                                 levels = c("tumor", "control"))

#change ordering of patterns and add colors
repeats_no_color = unique(as.character(height$repeat_type[!height$repeat_type %in% names(col_repeat_types)]))
n_repeats_no_color = length(repeats_no_color)
height$repeat_type = factor(height$repeat_type,
                            levels = c(names(col_repeat_types)[names(col_repeat_types)!="other"], 
                                       repeats_no_color,
                                       "other"))

palette = colorRampPalette(brewer.pal(9, "Greys"))

col_repeat_types[repeats_no_color] = rev(palette(n_repeats_no_color+2)[-c(1,n_repeats_no_color+2)])


#plot
p_tel_content = ggplot(height, aes(sample, sum, fill=repeat_type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name = paste(repeat_threshold, count_type, "repeats\nmapq threshold =", mapq_threshold),
                    values=col_repeat_types) +
  theme(axis.title.x = element_blank()) +
  ylab(paste0("Telomere content (intratelomeric reads per\nmillion reads with GC Content of ", gc_lower_limit, "-", gc_upper_limit, "%)")) +
  ggtitle(paste0(pid,": GC corrected telomere content"))

plot_dir = file.path(spectrum_dir, "plots")
if (!(file.exists(plot_dir))){dir.create(plot_dir, recursive=TRUE)}

plot_file_prefix = paste0(plot_dir,"/", pid, "_telomere_content")

for (plot_type in plot_file_format){
  ggsave(paste0(plot_file_prefix, ".", plot_type), p_tel_content, width=6, height=5)
}

temp_dir = paste0(plot_dir, "/temp")
if (!(file.exists(temp_dir))){dir.create(temp_dir, recursive=TRUE)}
save(p_tel_content, col_repeat_types, file = paste0(temp_dir, "/", pid, "_plot_tel_content.rds"))

remove = file.remove("Rplots.pdf")