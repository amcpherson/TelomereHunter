# Usage: R --no-save --slave --args <PID> <SPECTRUM_SUMMARY_DIR> <REPEAT_THRESHOLD> <CONSECUTIVE_FLAG> <MAPQ> <T_TYPE> <C_TYPE> <G_TYPE> <J_TYPE> <PLOT_FILE_FORMAT> < plot_repeat_frequency_intratelomeric.R
# Description: Makes a histograms of the telomere repeats per intratelomeric read in the tumor and control sample

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
pid = commandArgs[5]
spectrum_dir = paste(commandArgs[6], "/", pid, "/", sep="")
repeat_threshold = commandArgs[7]
consecutive_flag = commandArgs[8]
mapq_threshold = commandArgs[9]
repeat_types = commandArgs[10]
plot_file_format = commandArgs[11]

library(ggplot2, quietly=TRUE, warn.conflicts=FALSE)
library(grid, quietly=TRUE, warn.conflicts=FALSE)
library(gridExtra, quietly=TRUE, warn.conflicts=FALSE)
library(cowplot, quietly=TRUE, warn.conflicts=FALSE)

if (consecutive_flag == "True"){
  count_type = "consecutive"
}else{
  count_type = "non-consecutive"
}

if (plot_file_format=="all"){
  plot_file_format=c("pdf", "png", "svg")
}


frequency_table_file_T = paste0(spectrum_dir, "/tumor_TelomerCnt_", pid, "/", pid, "_repeat_frequency_per_intratelomeric_read.tsv")


colors = c() 

if (file.exists(frequency_table_file_T)){
  colors = c(colors, "blue") 
  frequency_table_T = read.table(frequency_table_file_T, header=TRUE)
  frequency_table_T$sample="tumor"
  frequency_table_T$percent = frequency_table_T$count / sum(frequency_table_T$count) * 100
  frequency_table_T$percent_cumulative = cumsum(frequency_table_T$percent)
  #frequency_table_T$percent_cumulative = rev(cumsum(rev(frequency_table_T$percent)))
}else{
  frequency_table_T = matrix(, nrow = 0, ncol = 5)
}


frequency_table_file_C = paste0(spectrum_dir, "/control_TelomerCnt_", pid, "/", pid, "_repeat_frequency_per_intratelomeric_read.tsv")

if (file.exists(frequency_table_file_C)){
  colors = c(colors, "darkorange") 
  frequency_table_C = read.table(frequency_table_file_C, header=TRUE)
  frequency_table_C$sample="control"
  frequency_table_C$percent = frequency_table_C$count / sum(frequency_table_C$count) * 100
  frequency_table_C$percent_cumulative = cumsum(frequency_table_C$percent)
  #frequency_table_C$percent_cumulative = rev(cumsum(rev(frequency_table_C$percent)))
}else{
  frequency_table_C = matrix(, nrow = 0, ncol = 5)
}


df = data.frame(rbind(as.matrix(frequency_table_T),as.matrix(frequency_table_C)))
df$sample = factor(df$sample, levels=c("tumor", "control"))
df$number_repeats = as.numeric(levels(df$number_repeats))[df$number_repeats]
df$percent = as.numeric(levels(df$percent))[df$percent]
df$percent_cumulative = as.numeric(levels(df$percent_cumulative))[df$percent_cumulative]


# make histogram
p_hist = ggplot(df, aes(x=number_repeats, y=percent)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sample) +
  xlab("Number of repeats") +
  ylab("Percent of intratelomeric reads")


#make cumulative line plot
df$group = factor(df$sample, levels=c("control", "tumor"))  #order so that control is plotted first

p_cum = ggplot(df, aes(x=number_repeats, y=percent_cumulative, group=sample)) +
  geom_line(aes(group = group, colour = sample)) +
  #ggtitle(paste0(pid, ": Frequency of telomere repeat occurrences in intratelomeric reads")) +
  theme(plot.title = element_text(face="bold", hjust=0.5)) +
  scale_color_manual(values=colors) +
  theme(legend.title=element_blank()) +
  theme(legend.position="top") +
  xlab("Number of repeats") +
  ylab("Percent (cumulative)") +
  scale_x_reverse()


#merge
p = plot_grid(p_hist, p_cum, 
                     rel_widths = c(1.7, 1))

title = paste0(pid, ": Frequency of telomere repeats in intratelomeric reads")

filtering_criteria = paste0("Filtering criteria: ",
                            repeat_threshold, " ", count_type, " repeats",
                            ", mapq threshold = ", mapq_threshold,
                            "\nrepeat types = ", repeat_types)

p_title = plot_grid(ggdraw() + draw_label(title, fontface='bold'), 
              p,
              ggdraw() + draw_label(filtering_criteria), 
              ncol=1, rel_heights=c(0.1, 1, 0.15))


#save plot
plot_dir = paste0(spectrum_dir, "plots/")

if (!(file.exists(plot_dir))){
  dir.create(plot_dir, recursive=TRUE)
}

plot_file = paste0(plot_dir, pid, "_hist_telomere_repeats_per_intratelomeric_read")

for (plot_type in plot_file_format){
  ggsave(paste0(plot_file, ".", plot_type), p_title, width=10, height=5)
}

#save histogram plot as R object
temp_dir = paste0(plot_dir, "/temp")
if (!(file.exists(temp_dir))){dir.create(temp_dir, recursive=TRUE)}
save(p_hist, file = paste0(temp_dir, "/", pid, "_hist_telomere_repeats_per_intratelomeric_read.rds"))

remove = file.remove("Rplots.pdf")