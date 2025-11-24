# Usage: R --no-save --slave --args <FUNCTION_DIR> <PID> <OUTPUT_DIR> <PLOT_FILE_FORMAT> < plot_gc_content_simple.R
# Description: Makes a summary plot of the GC content per read for all reads and for intratelomeric reads in the input sample (tumor or control)

# get commandline arguments
commandArgs = commandArgs()
pipeline_dir = commandArgs[5]
pid = commandArgs[6]
out_dir = paste0(commandArgs[7], "/", pid, "/")
plot_file_format = commandArgs[8]

if (plot_file_format=="both"){
  plot_file_format=c("pdf", "png")
}

source(file.path(pipeline_dir, "functions_for_plots.R"))

# get gc content tables
gc_content_file_T = paste0(out_dir, "/tumor_TelomerCnt_", pid, "/", pid, "_tumor_gc_content.tsv")
gc_content_file_C = paste0(out_dir, "/control_TelomerCnt_", pid, "/", pid, "_control_gc_content.tsv")

if (file.exists(gc_content_file_T)){
  sample = "tumor"
  color = "maroon2"
  gc_content_file_intratel = paste0(out_dir, "/tumor_TelomerCnt_", pid, "/", pid, "_intratelomeric_tumor_gc_content.tsv")

  gc_content = read.table(gc_content_file_T, header=TRUE)
  gc_content_intratel = read.table(gc_content_file_intratel, header=TRUE)
  
  summary_table=read.table(paste0(out_dir, "/tumor_TelomerCnt_", pid, "/", pid, "_tumor_summary.tsv"), header=TRUE, sep="\t")
  
}else if(file.exists(gc_content_file_C)){
  sample = "control"
  color = "turquoise3"
  gc_content_file_intratel = paste0(out_dir, "/control_TelomerCnt_", pid, "/", pid, "_intratelomeric_control_gc_content.tsv")
  
  gc_content = read.table(gc_content_file_C, header=TRUE)
  gc_content_intratel = read.table(gc_content_file_intratel, header=TRUE)
  
  summary_table=read.table(paste0(out_dir, "/control_TelomerCnt_", pid, "/", pid, "_control_summary.tsv"), header=TRUE, sep="\t")
}


# add percentage of reads to gc content tables
gc_content$fraction_of_reads = gc_content$read_count/sum(as.numeric(gc_content$read_count))
gc_content_intratel$fraction_of_reads = gc_content_intratel$read_count/sum(as.numeric(gc_content_intratel$read_count))

# get gc bins used for correction
gc_bins=as.character(summary_table$gc_bins_for_correction)
gc_bins=as.integer(strsplit(gc_bins, ", ")[[1]])



#make plot
plot_dir = file.path(out_dir, "plots")
plot_file_prefix = paste0(plot_dir, "/",  pid, "_gc_content")

if (!(file.exists(plot_dir))){
  dir.create(plot_dir, recursive=TRUE)
}


dfm_all = data.frame(gc_content_percent=gc_content$gc_content_percent, fraction_of_reads=gc_content$fraction_of_reads, sample=sample)

plot_all = ggplot(dfm_all, aes(gc_content_percent, fraction_of_reads, group=sample, color = sample)) +
  annotate("rect", xmin = gc_bins-0.5, xmax = gc_bins+0.5, ymin = 0, ymax = max(dfm_all$fraction_of_reads),alpha = 0.2, fill = rep(color, length(gc_bins))) +
  geom_line() +
  theme(legend.title=element_blank()) +
  ggtitle("All Reads") +
  xlab("GC Content [%]") +
  ylab("Fraction of Reads") +
  scale_color_manual(values=color) 

if (length(gc_bins)==5 && all.equal(gc_bins,c(48, 49, 50, 51, 52))){
  plot_all = plot_all +
    annotate("text", x = max(gc_bins), y = max(dfm_all$fraction_of_reads), label = "Bins used for\nGC correction", color="grey40", size=4, hjust = 0, vjust=2) 
}else{
  plot_all = plot_all +
    annotate("text", x = max(gc_bins), y = max(dfm_all$fraction_of_reads), label = "Bins used for adaptive\nGC correction", color="grey40", size=4, hjust = 0, vjust=2)
}

dfm_intratel = data.frame(gc_content_percent=gc_content_intratel$gc_content_percent, fraction_of_reads=gc_content_intratel$fraction_of_reads, sample=sample)

plot_intratel = ggplot(dfm_intratel, aes(gc_content_percent, fraction_of_reads, group=sample, color = sample)) +
  geom_line() +
  ggtitle("Intratelomeric Reads") +
  xlab("GC Content [%]") +
  ylab("Fraction of Reads")  +
  scale_color_manual(values=color) #+
#  geom_segment(aes(x = 0, y = 0.01, xend = 100, yend = 0.01), color="darkgrey")  #add line at 1%


grid_arrange_shared_legend(plot_list=list(plot_all, plot_intratel), plot_file_format=plot_file_format, title=paste0(pid, ": GC Content"), plot_file_prefix=plot_file_prefix)
