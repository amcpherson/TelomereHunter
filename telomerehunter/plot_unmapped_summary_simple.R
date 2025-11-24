# Usage: R --no-save --slave --args <FUNCTION_DIR> <PID> <SPECTRUM_SUMMARY_DIR> <REPEAT_THRESHOLD> <CONSECUTIVE_FLAG> <MAPQ> <PLOT_FILE_FORMAT> < plot_spectrum_simple.R
# Description: makes a bar plot of the number of intratelomeric reads (per million reads) in the tumor and control sample.

# get commandline arguments
commandArgs = commandArgs()
pipeline_dir = commandArgs[5]
pid = commandArgs[6]
spectrum_dir = paste0(commandArgs[7], "/", pid, "/")
repeat_threshold = commandArgs[8]
consecutive_flag = commandArgs[9]
mapq_threshold = commandArgs[10]
plot_file_format = commandArgs[11]

if (consecutive_flag == "True"){
  count_type = "Consecutive"
}else{
  count_type = "Non-consecutive"
}

if (plot_file_format=="both"){
  plot_file_format=c("pdf", "png")
}

source(file.path(pipeline_dir, "/functions_for_plots.R"))



# get spectrum
spectrum_tumor_file = paste0(spectrum_dir, "/tumor_TelomerCnt_", pid,"/", pid, ".spectrum")
spectrum_control_file = paste0(spectrum_dir, "/control_TelomerCnt_", pid,"/", pid, ".spectrum")

if (file.exists(spectrum_tumor_file)){
  sample="Tumor"
  spectrum = read.table(spectrum_tumor_file, header=TRUE, comment.char="")
  read_count_file = paste0(spectrum_dir, "/tumor_TelomerCnt_", pid,"/", pid, "_readcount.tsv")
}else{
  sample="Control"
  spectrum = read.table(spectrum_control_file, header=TRUE, comment.char="")
  read_count_file = paste0(spectrum_dir, "/control_TelomerCnt_", pid,"/", pid, "_readcount.tsv")
}

spectrum = spectrum[spectrum[,"chr"]=="unmapped",]

# get total number of reads 
readcount = read.table(read_count_file, header=TRUE, comment.char="")
total_reads = sum(as.numeric(readcount[ ,"reads"]))


#normalize
spectrum[ ,4:ncol(spectrum)] = spectrum[ ,4:ncol(spectrum)] * (spectrum[ ,"reads_with_pattern"] / sum(spectrum[ ,4:ncol(spectrum)]))
spectrum[ ,4:ncol(spectrum)] = spectrum[ ,4:ncol(spectrum)]*1000000 / total_reads


height = t(spectrum[ ,4:ncol(spectrum)])
colnames(height) = NULL

plot_dir = file.path(spectrum_dir, "plots")

if (!(file.exists(plot_dir))){
  dir.create(plot_dir, recursive=TRUE)
}


if (nchar(pid)<9){
  main = paste0(pid,": Telomere Repeat Types in Intratelomeric Reads")
}else{
  main = paste0(pid,":\nTelomere Repeat Types in Intratelomeric Reads")
}


plot_file_prefix = paste0(plot_dir,"/", pid, "_unmapped.pdf")

barplot_repeattype(height=height,
                   width=14,
                   plot_file_prefix=plot_file_prefix,
                   plot_file_format=plot_file_format, 
                   main=main,
                   cex.main=1,
                   repeat_threshold=repeat_threshold, count_type=count_type, mapq_threshold=mapq_threshold,
                   axis=FALSE, axis_simple=TRUE, labels=c(sample), xlas=1, tick=FALSE, cex.axis=0.85, cex.lab=0.9, inset_legend=c(-0.63,0), cex.names=0.9) 

