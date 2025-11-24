# Usage: R --no-save --slave --args <FUNCTION_DIR> <PID> <SPECTRUM_SUMMARY_DIR> <REPEATHRESHOLD> <CONSECUTIVE_FLAG> <MAPQ> <PLOT_FILE_FORMAT>< plot_spectrum_summary_simple.R
# Description: Makes a bar plot of the number of telomere reads (per million reads) in the different fractions of a sample (only tumor or control).

# get commandline arguments
commandArgs = commandArgs()
pipeline_dir = commandArgs[5]
pid = commandArgs[6]
spectrum_dir = paste0(commandArgs[7], "/", pid, "/")
repeat_threshold = commandArgs[8]
consecutive_flag = commandArgs[9]
mapq_threshold = commandArgs[10]
plot_file_format = commandArgs[11]


if (consecutive_flag == "yes"){
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



plot_dir = file.path(spectrum_dir, "plots")

if (!(file.exists(plot_dir))){
  dir.create(plot_dir, recursive=TRUE)
}

# make empty spectrum summary tables
spectrum_summary = matrix(data=0, nrow=4, ncol=(ncol(spectrum)-2), byrow=TRUE)
row.names(spectrum_summary)=c("intra_chromosomal", "subtelomeric", "junction_spanning", "intra_telomeric")
colnames(spectrum_summary) = colnames(spectrum)[3:ncol(spectrum)]
spectrum_summary = data.frame(spectrum_summary)


# get intra-telomeric reads (= unmapped reads)
spectrum_summary["intra_telomeric",] = spectrum[spectrum[,"chr"]=="unmapped",3:ncol(spectrum)]


for (chr in c(1:22, "X", "Y")){
  spectrum_chr = spectrum[spectrum[,"chr"]==chr,3:ncol(spectrum)]
  
  # get junction spanning reads
  spectrum_summary["junction_spanning", ] = spectrum_summary["junction_spanning", ] + spectrum_chr[1,] + spectrum_chr[nrow(spectrum_chr),]
  
  # get subtelomeric reads (= first and last band)
  spectrum_summary["subtelomeric", ] = spectrum_summary["subtelomeric", ] + spectrum_chr[2,] + spectrum_chr[nrow(spectrum_chr)-1,]
  
  # get intra-chromosomal reads (= all other bands) 
  spectrum_summary["intra_chromosomal", ] = spectrum_summary["intra_chromosomal", ] + apply(spectrum_chr[3:(nrow(spectrum_chr)-2),], 2, sum)   
  
}


# get total number of reads 
readcount = read.table(read_count_file, header=TRUE, comment.char="")
total_reads = sum(as.numeric(readcount[ ,"reads"]))


#normalize
spectrum_summary_norm = spectrum_summary[,2:ncol(spectrum_summary)] * (spectrum_summary[,"reads_with_pattern"] / apply(spectrum_summary[ ,2:ncol(spectrum_summary)], 1,sum))
spectrum_summary_norm = spectrum_summary_norm * (1000000 / total_reads)


height = t(spectrum_summary_norm)
colnames(height) = NULL


if (nchar(pid)<11){
  main = paste0(pid,": Telomere Repeat Types in ", sample, " Sample")
}else{
  main = paste0(pid,":\nTelomere Repeat Types in ", sample, " Sample")
}


plot_file_prefix = paste0(plot_dir,"/", pid, "_summary")
barplot_repeattype(height=height,
                   plot_file_prefix=plot_file_prefix,
                   plot_file_format=plot_file_format,
                   width=20,
                   main=main,
                   repeat_threshold=repeat_threshold, count_type=count_type, mapq_threshold=mapq_threshold,
                   axis=FALSE, axis_simple=TRUE, labels=c("Intrachromosomal", "Subtelomeric", "Junction Spanning", "Intratelomeric"), xlas=1, tick=FALSE, cex.axis=0.85,
                   cex.lab=0.9) 

