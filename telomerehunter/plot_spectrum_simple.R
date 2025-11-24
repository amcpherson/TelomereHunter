# Usage: R --no-save --slave --args <FUNCTION_DIR> <PID> <SPECTRUM_SUMMARY_DIR> <REPEAT_THRESHOLD> <CONSECUTIVE_FLAG> <MAPQ> <BANDING_FILE> <PLOT_FILE_FORMAT> < plot_spectrum_simple.R
# Description: Makes a bar plot for each chromosome. The chromosome bands are plotted on the x-axis and the specific telomere reads per million bases 
#              and per billion reads are plotted on the y-axis. Attention: the "junction bands" are not divided by band length (for these bands
#              the y-axis is "specific telomere reads per billion reads"), so the scaling is not correct in comparison to "normal" bands!!! 


# get commandline arguments
commandArgs = commandArgs()
pipeline_dir = commandArgs[5]
pid = commandArgs[6]
spectrum_dir = paste0(commandArgs[7], "/", pid, "/")
repeat_threshold = commandArgs[8]
consecutive_flag = commandArgs[9]
mapq_threshold = commandArgs[10]
banding_file = commandArgs[11]
plot_file_format = commandArgs[12]

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


# get total number of reads 
readcount = read.table(read_count_file, header=TRUE, comment.char="")
total_reads = sum(as.numeric(readcount[ ,"reads"]))


#normalize
spectrum[ ,4:ncol(spectrum)] = spectrum[ ,4:ncol(spectrum)] * (spectrum[ ,"reads_with_pattern"] / apply(spectrum[ ,4:ncol(spectrum)], 1,sum))
spectrum[ ,4:ncol(spectrum)] = spectrum[ ,4:ncol(spectrum)]*1000000000 / total_reads

# get band lengths
band_info = read.table(banding_file)
colnames(band_info) = c("chr", "start", "end", "band_name", "stain")
band_info[,"chr"] = gsub("chr", "", band_info[,1])
band_info[,"length"] = band_info[,"end"] - band_info[,"start"]

plot_dir = file.path(spectrum_dir, "plots")

if (!(file.exists(plot_dir))){
  dir.create(plot_dir, recursive=TRUE)
}

for (chr in c(1:22, "X", "Y")){
  
  spectrum_chr = spectrum[spectrum[,"chr"]==chr,]
  bands=spectrum_chr[ , "band"]
  
  #normalize by band length
  band_info_chr = band_info[band_info[,"chr"]==chr,]
  spectrum_chr[2:(nrow(spectrum_chr)-1),3:ncol(spectrum_chr)] = spectrum_chr[2:(nrow(spectrum_chr)-1),3:ncol(spectrum_chr)] *(1000000/ band_info_chr[,"length"])
  
  height = t(spectrum_chr[ ,4:ncol(spectrum_chr)])

  colnames(height)=NULL
  
  main = paste0(pid,": Telomere Repeat Types in Chr", chr, " (", sample, " Sample)")
  plot_file_prefix = paste0(plot_dir,"/", pid, "_", chr)

  barplot_repeattype(height=height,
                     plot_file_prefix=plot_file_prefix, plot_file_format=plot_file_format, width=28, mar=c(5.1, 5.1, 5.1, 9.3),
                     main=main, ylab="Telomere Reads (per Million Bases and per Billion Reads)",
                     repeat_threshold=repeat_threshold, count_type=count_type, mapq_threshold=mapq_threshold,
                     inset_legend=c(-0.18,0),
                     axis=FALSE, axis_simple=TRUE, labels=bands, cex.axis=0.75, cex.lab=0.9) 
 
}










