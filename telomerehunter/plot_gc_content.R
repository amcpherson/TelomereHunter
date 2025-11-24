# Usage: R --no-save --slave --args <FUNCTION_FILE> <PID> <OUTPUT_DIR> <PLOT_FILE_FORMAT> <GC_LOWER_LIMIT> <GC_UPPER_LIMIT> < plot_gc_content.R
# Description: Makes a summary plot of the GC content per read for all reads and for intratelomeric reads in the tumor and control sample

# Copyright 2015 Lina Sieverling

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
out_dir = paste0(commandArgs[7], "/", pid, "/")
plot_file_format = commandArgs[8]
gc_lower_limit = as.numeric(commandArgs[9])
gc_upper_limit = as.numeric(commandArgs[10])

library(cowplot, quietly=TRUE, warn.conflicts=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE)

#source(file.path(pipeline_dir, "functions_for_plots.R"))

if (plot_file_format=="all"){
  plot_file_format=c("pdf", "png", "svg")
}

colors = c("blue", "darkorange")
names(colors) = c("tumor", "control")

samples = c()
df_all = data.frame(gc_content_percent=c(0:100))
df_intratel = data.frame(gc_content_percent=c(0:100))

# get gc content tables
gc_content_file_T = paste0(out_dir, "/tumor_TelomerCnt_", pid, "/", pid, "_tumor_gc_content.tsv")
gc_content_file_C = paste0(out_dir, "/control_TelomerCnt_", pid, "/", pid, "_control_gc_content.tsv")
gc_content_file_intratel_T = paste0(out_dir, "/tumor_TelomerCnt_", pid, "/", pid, "_intratelomeric_tumor_gc_content.tsv")
gc_content_file_intratel_C = paste0(out_dir, "/control_TelomerCnt_", pid, "/", pid, "_intratelomeric_control_gc_content.tsv")

if (file.exists(gc_content_file_T)){
  samples=c(samples, "tumor")
   
  gc_content_T = read.table(gc_content_file_T, header=TRUE)
  gc_content_intratel_T = read.table(gc_content_file_intratel_T, header=TRUE)
  
  # add percentage of reads to gc content tables
  gc_content_T$fraction_of_reads = gc_content_T$read_count/sum(as.numeric(gc_content_T$read_count))
  gc_content_intratel_T$fraction_of_reads = gc_content_intratel_T$read_count/sum(as.numeric(gc_content_intratel_T$read_count))
  df_all$tumor = gc_content_T$fraction_of_reads
  df_intratel$tumor = gc_content_intratel_T$fraction_of_reads
}


if(file.exists(gc_content_file_C)){
  samples=c(samples, "control")
    
  gc_content_C = read.table(gc_content_file_C, header=TRUE)
  gc_content_intratel_C = read.table(gc_content_file_intratel_C, header=TRUE)
  
  # add percentage of reads to gc content tables
  gc_content_C$fraction_of_reads = gc_content_C$read_count/sum(as.numeric(gc_content_C$read_count))
  gc_content_intratel_C$fraction_of_reads = gc_content_intratel_C$read_count/sum(as.numeric(gc_content_intratel_C$read_count)) 
  df_all$control = gc_content_C$fraction_of_reads
  df_intratel$control = gc_content_intratel_C$fraction_of_reads
}


dfm_all = melt(df_all, id.var = c("gc_content_percent"), variable.name="sample", value.name="fraction_of_reads")
dfm_intratel = melt(df_intratel, id.var = c("gc_content_percent"), variable.name="sample", value.name="fraction_of_reads")

gc_bins = c(gc_lower_limit:gc_upper_limit)


#merge
dfm_all$read_type = "All reads"
dfm_intratel$read_type = "Intratelomeric reads"

dfm = rbind(dfm_all, dfm_intratel)


#################
### make plot ###
#################

plot_dir = file.path(out_dir, "plots")
plot_file_prefix = paste0(plot_dir, "/",  pid, "_gc_content")

if (!(file.exists(plot_dir))){dir.create(plot_dir, recursive=TRUE)}


#sort: control should be plotted first
dfm$group = factor(dfm$sample, levels = c("control", "tumor"))


#plot
plot = ggplot(dfm) +
  geom_rect(data=data.frame(read_type=factor(c("All reads"))),
            xmin = min(gc_bins-0.5), xmax = max(gc_bins+0.5),
            ymin = 0, ymax = max(dfm_all$fraction_of_reads),
            alpha = 0.5, fill = "grey") +
  geom_line(aes(gc_content_percent, fraction_of_reads, group=group, color = sample)) +
  theme(legend.position = "top", legend.justification = "center") +
  theme(legend.title=element_blank()) +
  ggtitle(paste0(pid, ": GC content")) +
  xlab("GC content [%]") +
  ylab("Fraction of reads") +
  scale_color_manual(values=colors[samples]) +
  geom_text(data=data.frame(read_type=factor(c("All reads"))),
            x = min(gc_upper_limit, 88), y = max(dfm_all$fraction_of_reads), 
            label = "Bins used for \nGC correction", color="grey40", size=4, hjust = 0, vjust=2) +
  facet_wrap(~read_type,scales="free")


for (plot_type in plot_file_format){
  ggsave(paste0(plot_file_prefix, ".", plot_type), plot, width = 10, height = 5)
}

remove = file.remove("Rplots.pdf")