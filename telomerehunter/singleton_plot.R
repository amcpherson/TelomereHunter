# usage: R --no-save --slave --args <PIPELINE_DIR> <MAIN_PATH> <PID> <PLOT_FILE_FORMAT> < singleton_plot.R
# description: uses the output of TVR_context_summary_tables.R to make a plot containing:
#              - the raw singleton counts in the tumor and in the control sample
#              - the normalized singleton counts
#              - the normalized singleton count log2 ratio
#              - the distance to the expected singleton count log2 ratio

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

library("ggplot2", quietly=TRUE, warn.conflicts=FALSE)
library("reshape2", quietly=TRUE, warn.conflicts=FALSE)
library("cowplot", quietly=TRUE, warn.conflicts=FALSE)
library("dplyr", quietly=TRUE, warn.conflicts=FALSE)


# get commandline arguments
commandArgs = commandArgs()
pipeline_dir = commandArgs[5]
main_path = commandArgs[6]
pid = commandArgs[7]
plot_file_format = commandArgs[8]

if (plot_file_format=="all"){
  plot_file_format=c("pdf", "png", "svg")
}


col_repeat_types = c("#CB181D", "#FF656A", "#74C476", "#C1FFC3", "#2171B5", "#6EBEFF", "#FFA500", "#FFF24D", "#9370DB", "#E0BDFF", "#000000")
names(col_repeat_types) = c("TTAGGG", "CCCTAA", "TGAGGG", "CCCTCA", "TCAGGG", "CCCTGA", "TTGGGG", "CCCCAA", "TTCGGG", "CCCGAA", "other")


plot_dir = paste0(main_path, "/plots/")
temp_dir = paste0(plot_dir, "/temp")
if (!(file.exists(temp_dir))){dir.create(temp_dir, recursive=TRUE)}


plot_file_prefix = paste0(plot_dir, pid, "_singletons")

##############################################################################################################

singletons = read.table(paste0(main_path, "/", pid, "_singletons.tsv"),
           sep="\t", header=TRUE, stringsAsFactors = FALSE)

singletons = singletons %>%
  arrange(-singleton_count_tumor) %>%
  mutate(pattern = factor(pattern, levels = pattern))


#-------------------------------------------------------------
# plot counts
#-------------------------------------------------------------

counts = melt(singletons, id.vars = c("PID", "pattern"), measure.vars=c("singleton_count_tumor", "singleton_count_control"),
              variable.name = "sample", value.name = "counts")

counts$sample = gsub("singleton_count_", "", counts$sample)
counts$sample = factor(counts$sample, levels=c("tumor", "control"))

counts = counts[!is.na(counts$counts),]

p_counts = ggplot(counts, aes(pattern, counts)) +
  geom_bar(aes(fill=pattern, alpha=sample), stat="identity", position = "dodge", width=0.5) +
  #theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #scale_fill_manual(values=col_repeat_types[unique(counts$pattern)], na.value="grey27", guide=FALSE) +
  scale_fill_manual(values=col_repeat_types, na.value="grey27", guide=FALSE) +
  scale_alpha_manual(values=c(1, 0.3)) +
  theme(axis.title.x = element_blank()) +
  ylab("Raw singleton counts") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top", legend.justification = "center") +
  theme(legend.margin=margin(t = 0, unit='cm')) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) # rotate x-axis labels

#-------------------------------------------------------------
# plot normalized counts (raw counts/total read count)
#-------------------------------------------------------------
counts_norm = melt(singletons, id.vars = c("PID", "pattern"), measure.vars=c("singleton_count_tumor_norm", "singleton_count_control_norm"),
                   variable.name = "sample", value.name = "counts_norm")

counts_norm$sample = gsub("singleton_count_", "", counts$sample)
counts_norm$sample = factor(counts_norm$sample, levels=c("tumor", "control"))

counts_norm = counts_norm[!is.na(counts_norm$counts_norm),]

p_singletons_norm = ggplot(counts_norm, aes(pattern, counts_norm)) +
  geom_bar(aes(fill=pattern, alpha=sample), stat="identity", position = "dodge", width=0.5) +
  #theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values=col_repeat_types, na.value="grey27", guide=FALSE) +
  #scale_fill_manual(values=col_repeat_types[unique(counts$pattern)], na.value="grey27", guide=FALSE) +
  scale_alpha_manual(values=c(1, 0.3)) +
  theme(axis.title.x = element_blank()) +
  ylab("Normalized singleton counts\n(singleton count / total read count)") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top", legend.justification = "center") +
  theme(legend.margin=margin(t = 0, unit='cm')) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) # rotate x-axis labels



title <- ggdraw() + draw_label(paste0(pid, "\nSingletons (TTAGGG)3-NNNGGG-(TTAGGG)3 in intratelomeric reads"),
                               fontface='bold')

if(all(is.na(singletons$singleton_count_log2_ratio))){
  
  # only plot counts if only one sample exists
  
  p_grid = plot_grid(p_counts, p_singletons_norm, nrow=1, align="vh", axis="tb", scale=0.93)
  
  p_grid_title = plot_grid(title, p_grid, ncol=1, rel_heights=c(0.1, 1)) 
  
  for (plot_type in plot_file_format){
    ggsave(paste0(plot_file_prefix, ".", plot_type), p_grid_title, width=10, height=5)
  }
  
  
  save(p_singletons_norm, file = paste0(temp_dir, "/", pid, "_singletons.rds"))
  
}else{
  
  #-------------------------------------------------------------
  # plot log2 ratio
  #-------------------------------------------------------------
  
  #set NA to 0 to prevent warning
  singletons = singletons %>%
    mutate(singleton_count_log2_ratio_norm = ifelse(is.na(singleton_count_log2_ratio_norm),
                                                    0,
                                                    singleton_count_log2_ratio_norm))
  
  # if a absolute log2 ratio is larger than 5, add the number to the plot (the plot is cut off at 5)
  singletons = singletons %>%
    mutate(label_top_log2 = ifelse(singleton_count_log2_ratio_norm > 5,
                                   formatC(round(singleton_count_log2_ratio_norm, 1), format='f', digits=1), NA)) %>%
    mutate(label_bottom_log2 = ifelse(singleton_count_log2_ratio_norm < -5,
                                      formatC(round(singleton_count_log2_ratio_norm, 1), format='f', digits=1), NA))
  
  
  p_singleton_log2 = ggplot(singletons, aes(pattern, singleton_count_log2_ratio_norm)) +
    geom_hline(aes(yintercept = unique(singletons$tel_content_log2_ratio), color="expected log2 ratio (= telomere content log2 ratio)"), 
               linetype="dotted", size=1)+
    scale_color_manual(values="darkgrey") +
    theme(legend.position="top", legend.justification = "left") +
    theme(legend.title = element_blank()) +
    geom_bar(aes(fill=pattern), stat="identity") +
    scale_fill_manual(values=col_repeat_types, na.value="grey50", guide=FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title.x = element_blank()) +
    ylab("Normalized singleton count tumor/control (log2)") +
    theme(legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) + # rotate x-axis labels
    coord_cartesian(ylim=c(-5, 5)) +
    geom_text(y=5, aes(x=pattern, label=label_top_log2), na.rm=TRUE) +
    geom_text(y=-5, aes(x=pattern, label=label_bottom_log2), na.rm=TRUE)
  
  #-------------------------------------------------------------
  # plot distance
  #------------------------------------------------------------- 
  
  #set NA to 0 to prevent warning
  singletons = singletons %>%
    mutate(distance_to_expected_singleton_log2_ratio = ifelse(is.na(distance_to_expected_singleton_log2_ratio),
                                                              0,
                                                              distance_to_expected_singleton_log2_ratio))
  
  # if a absolute log2 ratio is larger than 5, add the number to the plot (the plot is cut off at 5)
  singletons = singletons %>%
    mutate(label_top_dist = ifelse(distance_to_expected_singleton_log2_ratio > 5,
                                   formatC(round(distance_to_expected_singleton_log2_ratio, 1), format='f', digits=1), NA)) %>%
    mutate(label_bottom_dist = ifelse(distance_to_expected_singleton_log2_ratio < -5,
                                      formatC(round(distance_to_expected_singleton_log2_ratio, 1), format='f', digits=1), NA))
  
  p_distance = ggplot(singletons, aes(pattern, distance_to_expected_singleton_log2_ratio, fill=pattern)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=col_repeat_types, na.value="grey50", guide=FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title.x = element_blank()) +
    ylab("Distance to expected singleton count\n(norm. singleton log2 ratio - tel. content log2 ratio)") +
    theme(legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) + # rotate x-axis labels
    coord_cartesian(ylim=c(-5, 5)) +
    geom_text(y=5, aes(x=pattern, label=label_top_dist), na.rm=TRUE) +
    geom_text(y=-5, aes(x=pattern, label=label_bottom_dist), na.rm=TRUE)
  
  #-------------------------------------------------------------
  # combine plots
  #------------------------------------------------------------- 
  
  p_grid = plot_grid(p_counts + theme(axis.title.y = element_text(margin = margin(r = -30))), 
                     p_singletons_norm, 
                     p_singleton_log2 + theme(axis.title.y = element_text(margin = margin(r = -30))), 
                     p_distance, nrow=2, align="vh", axis="tb", scale=0.93)

  p_grid_title = plot_grid(title, p_grid, ncol=1, rel_heights=c(0.1, 1)) 
  
  for (plot_type in plot_file_format){
    ggsave(paste0(plot_file_prefix, ".", plot_type), p_grid_title, width=10, height=10)
  }
  
  save(singletons, p_singleton_log2, file = paste0(temp_dir, "/", pid, "_singletons.rds"))
}



remove = file.remove("Rplots.pdf")




