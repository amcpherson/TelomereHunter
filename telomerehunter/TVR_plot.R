# usage: R --no-save --slave --args <MAIN_PATH> <PID> <PLOT_FILE_FORMAT> < TVR_plot.R
# description: uses the output of normalize_TVR_counts.py to make a plot containing:
#              1) a barplot with mean pattern occurrences (variant repeats) per intratelomeric read in tumor and control sample
#              2) a barplot with the log2 ratio (tumor/control) of pattern occurrences per intratelomeric read


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
main_path = commandArgs[5]
pid = commandArgs[6]
plot_file_format = commandArgs[7]

library(cowplot, quietly=TRUE, warn.conflicts=FALSE)
library(dplyr, quietly=TRUE, warn.conflicts=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE)

if (plot_file_format=="all"){
  plot_file_format=c("pdf", "png", "svg")
}

col_repeat_types = c("#CB181D", "#FF656A", "#74C476", "#C1FFC3", "#2171B5", "#6EBEFF", "#FFA500", "#FFF24D", "#9370DB", "#E0BDFF", "#000000")
names(col_repeat_types) = c("TTAGGG", "CCCTAA", "TGAGGG", "CCCTCA", "TCAGGG", "CCCTGA", "TTGGGG", "CCCCAA", "TTCGGG", "CCCGAA", "other")

########################################################################################################################################################

table_file = paste0(main_path, "/", pid, "_normalized_TVR_counts.tsv")
table_merged = read.table(table_file, header=TRUE, stringsAsFactors = FALSE)

plot_dir = paste0(main_path, "/plots/")
temp_dir = paste0(plot_dir, "/temp")
if (!(file.exists(temp_dir))){dir.create(temp_dir, recursive=TRUE)}

#only keep patterns if the normalized pattern count is bigger than 0.01 in tumor or control
table_merged_filter = table_merged %>%
  filter(Count_norm_by_intratel_reads_T>=0.01 | Count_norm_by_intratel_reads_C>=0.01)

####################################################
### make pattern screen bar plot absolute counts ###
####################################################

table_norm_by_intratel = table_merged_filter %>%
  select(Pattern, Count_norm_by_intratel_reads_T, Count_norm_by_intratel_reads_C) %>%
  rename(tumor = Count_norm_by_intratel_reads_T,
         control = Count_norm_by_intratel_reads_C)

table_norm_by_intratel = melt(table_norm_by_intratel,
                              id.vars = "Pattern",
                              variable.name = "sample",
                              value.name = "count_norm_by_intratel_reads")

#order patterns
table_norm_by_intratel$Pattern = factor(table_norm_by_intratel$Pattern, levels = table_merged_filter$Pattern)

#define colors
col_repeat_types["GTTGGG"] = col_repeat_types["TTGGGG"]


col_repeat_types_plot = col_repeat_types

for(pattern in table_norm_by_intratel$Pattern){
  if(! pattern %in% names(col_repeat_types)){
    col_repeat_types_plot[pattern] = "#454545"
  }
}

#plot with TTAGGG
p_TVR_counts_TTAGGG = ggplot(table_norm_by_intratel, aes(Pattern, count_norm_by_intratel_reads, alpha=sample, fill=Pattern)) +
  geom_bar(stat="identity", position="dodge") +
  scale_alpha_manual(values = c(1, 0.25)) +
  scale_fill_manual(values = col_repeat_types_plot) +
  theme(axis.title=element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  theme(legend.position = "none") +
  ggtitle("Including TTAGGG") +
  theme(plot.title = element_text(size = 10, face="plain"))

#plot without TTAGGG
p_TVR_counts = table_norm_by_intratel %>%
  filter(Pattern != "TTAGGG") %>%
  ggplot(aes(Pattern, count_norm_by_intratel_reads, alpha=sample, fill=Pattern)) +
  geom_bar(stat="identity", position="dodge") +
  scale_alpha_manual(values = c(1, 0.25)) +
  scale_fill_manual(values = col_repeat_types_plot, guide = FALSE) +
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) +  # rotate x-axis labels
  ylab("Mean TVR counts per intratelomeric read")


p_TVR_counts_all = ggdraw() +
  draw_plot(p_TVR_counts) +
  draw_plot(p_TVR_counts_TTAGGG, x = 0.7, y = 0.65, width = 0.3, height = 0.3) 



####################################################
### make pattern screen bar plot log2 ratio      ###
####################################################

if(!is.na(table_merged_filter$log2_ratio_count_norm_by_intratel_reads)[1]){
  
  # if a absolute log2 ratio is larger than 5, add the number to the plot (the plot is cut off at 5)
  table_merged_filter = table_merged_filter %>%
    mutate(label_top = ifelse(log2_ratio_count_norm_by_intratel_reads > 5,
                              formatC(round(log2_ratio_count_norm_by_intratel_reads, 1), format='f', digits=1), NA)) %>%
    mutate(label_bottom = ifelse(log2_ratio_count_norm_by_intratel_reads < -5,
                                 formatC(round(log2_ratio_count_norm_by_intratel_reads, 1), format='f', digits=1), NA))
  
  #sort
  table_merged_filter$Pattern = factor(table_merged_filter$Pattern, levels=table_merged_filter$Pattern)
  
  #colors
  col_repeat_types_log2 = col_repeat_types
  for(pattern in table_norm_by_intratel$Pattern){
    if(! pattern %in% names(col_repeat_types)){
      col_repeat_types_log2[pattern] = "grey50"
    }
  }
  
  p_TVR_log2 = ggplot(table_merged_filter, aes(Pattern, log2_ratio_count_norm_by_intratel_reads, fill=Pattern)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = col_repeat_types_log2, guide=FALSE) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) +  # rotate x-axis labels
    ylab("Normalized TVR counts tumor/control (log2)") +
    coord_cartesian(ylim=c(-5, 5)) +
    geom_text(y=5, aes(x=Pattern, label=label_top), na.rm=TRUE) +
    geom_text(y=-5, aes(x=Pattern, label=label_bottom), na.rm=TRUE)
}

####################################################
### merge bar plots and save                     ###
####################################################

#merge and save R objects
if(exists("p_TVR_log2")){
  p_no_title = plot_grid(p_TVR_counts_all,
                p_TVR_log2 + theme(plot.margin=margin(l=15, r=7, t=7, b=7, unit="pt")))
  width = 10
  
  save(p_TVR_log2, file = paste0(temp_dir, "/", pid, "_TVR_plots.rds"))
  
}else{
  p_no_title = p_TVR_counts_all
  width = 5
  
  save(p_TVR_counts, p_TVR_counts_TTAGGG, file = paste0(temp_dir, "/", pid, "_TVR_plots.rds"))
}

#add title
main_title = paste0(pid,": TVRs found in intratelomeric reads")

p = plot_grid(ggdraw() + draw_label(main_title, fontface='bold'), 
              p_no_title,
              ncol=1, rel_heights=c(0.1, 1))

#save plots
for (plot_type in plot_file_format){
  ggsave(paste0(plot_dir, pid, "_TVR_barplot.", plot_type),
         p, width=width, height=5)
}



#######################################
### make pattern screen scatterplot ###
#######################################

row.names(table_merged) = table_merged$Pattern 

# skip if counts for tumor or control sample don't exist
if(!is.na(table_merged["TTAGGG", "Count_norm_by_intratel_reads_T"]) & !is.na(table_merged["TTAGGG", "Count_norm_by_intratel_reads_C"])){

  # set pattern counts of 0 to 0.000001 for logarithmic scale of plot
  table_merged[table_merged$Count_norm_by_intratel_reads_T==0, "Count_norm_by_intratel_reads_T"] = 0.000001
  table_merged[table_merged$Count_norm_by_intratel_reads_C==0, "Count_norm_by_intratel_reads_C"] = 0.000001

  table_merged$label="other"
  table_merged["TTAGGG", "label"] = "TTAGGG"
  table_merged["TGAGGG", "label"] = "TGAGGG"
  table_merged["TCAGGG", "label"] = "TCAGGG"
  table_merged["GTTGGG", "label"] = "TTGGGG"
  table_merged["TTCGGG", "label"] = "TTCGGG"

  table_merged$label = factor(table_merged$label, levels = c("TTAGGG", "TGAGGG", "TCAGGG", "TTGGGG", "TTCGGG", "other"))

  plot_file_prefix = paste0(plot_dir, pid, "_pattern_screen_scatterplot")
     
  plot = ggplot(table_merged, aes(Count_norm_by_intratel_reads_C, Count_norm_by_intratel_reads_T, colour = label)) +
    geom_point(size = 3) +
    scale_y_log10() +
    scale_x_log10() +
    ylab("Mean TVR counts per intratel. read (tumor)")  +             # set y-axis label  
    xlab("Mean TVR counts per intratel. read (control)")  +             # set x-axis label  
    geom_abline(intercept = 0, slope = 1) + 
    geom_abline(intercept = 0.25, slope = 1, colour = "darkgrey") + 
    geom_abline(intercept = -0.25, slope = 1, colour = "grey") + 
    geom_abline(intercept = 0.5, slope = 1, colour = "grey") + 
    geom_abline(intercept = -0.5, slope = 1, colour = "grey") + 
    scale_colour_manual(values=c(unname(col_repeat_types[c("TTAGGG", "TGAGGG", "TCAGGG", "TTGGGG", "TTCGGG")]) , "black")) +
    theme(legend.title=element_blank()) +    # remove legend title
    theme(legend.position="top") +  
    ggtitle(paste0(pid, ": TVRs per intratelomeric read"))  +
    theme(plot.margin=margin(l=7, r=20, t=7, b=7, unit="pt"))  #add margin so that x axis label is not cut off

    
  for (plot_type in plot_file_format){
    ggsave(paste0(plot_dir, pid, "_TVR_scatterplot.", plot_type), plot, width=5, height=5)
  }
}


remove = file.remove("Rplots.pdf")



