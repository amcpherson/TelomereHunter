
# usage: R --no-save --slave --args <PIPELINE_DIR> <MAIN_PATH> <PID> <PLOT_FILE_FORMAT> < singleton_plot.R
# description: merges the most important output plots of TelomereHunter

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

library(cowplot, quietly=TRUE, warn.conflicts=FALSE)


# get commandline arguments
commandArgs = commandArgs()
main_path = commandArgs[5]
pid = commandArgs[6]
plot_file_format = commandArgs[7]

if (plot_file_format=="all"){
  plot_file_format=c("pdf", "png", "svg")
}

plot_dir = paste0(main_path, "/plots")
temp_dir = paste0(plot_dir, "/temp")


#-----------------------------------------------
# load plot objects
#-----------------------------------------------

#telomere content
tel_content_file = paste0(temp_dir, "/", pid, "_plot_tel_content.rds")
load(tel_content_file)

#quality control
qc_file = paste0(temp_dir, "/", pid, "_hist_telomere_repeats_per_intratelomeric_read.rds")
load(qc_file)

#TVRs
TVR_file = paste0(temp_dir, "/", pid, "_TVR_plots.rds")
load(TVR_file)


if(exists("p_TVR_log2")){
  p_TVR = p_TVR_log2
}else{
  p_TVR = p_TVR_counts
}


#singletons
singletons_file = paste0(temp_dir, "/", pid, "_singletons.rds")
load(singletons_file)

if(exists("p_singleton_log2")){
  p_singletons = p_singleton_log2
}else{
  p_singletons = p_singletons_norm
}


#-----------------------------------------------
# remove R object files
#-----------------------------------------------

remove = file.remove(tel_content_file, qc_file, TVR_file, singletons_file)
remove = file.remove(temp_dir, recursive=TRUE)

#-----------------------------------------------
# combine plots
#-----------------------------------------------

p = plot_grid(p_tel_content + 
                theme(legend.position="top") +
                theme(legend.title=element_blank()) +
                guides(fill = guide_legend(nrow=2)) + 
                ggtitle("GC corrected telomere content") +
                ylab("Telomere content (intratelomeric reads per\nmillion reads with GC content of 48-52%)"), #+
                #theme(axis.title.y = element_text(margin = margin(r = -30))),
          p_hist + ggtitle("Telomere repeats in intratelomeric reads") +
            theme(axis.title.x = element_text(margin = margin(t = -60))), # +
            # theme(axis.title.y = element_text(margin = margin(r = -30))),
          p_TVR + ggtitle("Telomere variant repeats (arbitrary context)"), # +
            # theme(axis.title.y = element_text(margin = margin(r = -30))),
          p_singletons + ggtitle("Singleton telomere variant repeats"),
          axis = "bl",
          align = "hv",
          nrow=2)


p_title = plot_grid(ggdraw() + draw_label(pid, fontface='bold'), 
                         p, ncol=1, rel_heights=c(0.05, 1)) 


#add TTAGG TVR plot
if(!exists("p_TVR_log2")){
  p = ggdraw() +
    draw_plot(p_title) +
    draw_plot(p_TVR_counts_TTAGGG, x = 0.35, y = 0.32, width = 0.12, height = 0.15) 
}


for (plot_type in plot_file_format){
  ggsave(paste0(main_path, "/", pid, "_telomerehunter_summary_plot.", plot_type), p_title, width=10, height=10)
}

remove = file.remove("Rplots.pdf")

