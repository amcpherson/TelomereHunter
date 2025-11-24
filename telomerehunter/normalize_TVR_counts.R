# usage: R --no-save --slave --args <MAIN_PATH> <PID> <TVRs_for_summary> < normalize_TVR_counts.R

# description: -uses the output tables of TVR_screen.py 
#              - normalizes counts of variant repeats by dividing by the total number of intratelomeric reads (result: mean number of variant repeats per intratelomeric read)
#              - also calcualtes other normalized values: 
#                       - by dividing by all reads in the sample
#                       - by dividing by the summed intratelomeric read length (to get count per 100 bp of intratelomeric read)
#              - calculates log2 ratio
#              - output: table containing normalized variant counts for tumor and control sample and the log2 ratio
#              - adds counts normalized by all reads of selected patterns to the summary table

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
TVRs_for_summary = commandArgs[7]

TVRs_for_summary = unlist(strsplit(TVRs_for_summary, split=","))

#####################################################################################################################################################

pattern_table_list = list()

for (sample in c("tumor", "control")){
  
  sample_dir = paste0(main_path, "/", sample, "_TelomerCnt_", pid, "/")

  #read TVR table
  pattern_table_path = paste0(sample_dir, "/TVRs/", pid, "_", sample, "_TVRs.txt")    
  if(!file.exists(pattern_table_path)){
    pattern_table_list[[sample]] = data.frame(Pattern="TTA", Count_norm_by_all_reads=NA, Count_norm_by_intratel_reads=NA, Count_per_100_bp_intratel_read=NA)
    next}    
  pattern_table = read.table(pattern_table_path, header=TRUE, comment.char = "")
  
  #get total number of reads
  readcount_table = read.table(file.path(sample_dir, paste0(pid, "_readcount.tsv")), header=TRUE)
  total_reads = sum(as.numeric(readcount_table$reads))  
  
  #normalize pattern count by total number of reads
  pattern_table$Count_norm_by_all_reads = pattern_table$Count / total_reads

  #get total number of intratelomeric reads
  total_intratelomeric = as.numeric(system(paste("samtools view -c", paste0(sample_dir, pid, "_filtered_intratelomeric.bam")), intern=TRUE))

  #normalize pattern count by intratelomeric reads (= mean number of repeat occurences per intratelomeric read)
  pattern_table$Count_norm_by_intratel_reads = pattern_table$Count / total_intratelomeric
  
  #get summed read length
  summed_read_length_table = read.table(paste0(sample_dir, "/TVRs/", pid, "_", sample, "_summed_read_length.tsv"), 
                                      header=TRUE, sep="\t")

  #normalize pattern count by summed intratelomeric read length
  pattern_table$Count_per_100_bp_intratel_read = pattern_table$Count*100 / summed_read_length_table[1, "summed_intratel_read_length"]
  
  pattern_table_list[[sample]] = pattern_table[, c("Pattern", "Count_norm_by_all_reads", "Count_norm_by_intratel_reads", "Count_per_100_bp_intratel_read")]
  
}

#merge tables
table_merged = merge(pattern_table_list[["tumor"]], pattern_table_list[["control"]], by="Pattern", all=TRUE, sort=FALSE, suffixes = c("_T","_C"))

#if pattern tables for both samples exist: set missing counts to zero
if(!(is.na(pattern_table_list[["tumor"]][1,"Count_norm_by_intratel_reads"]) || is.na(pattern_table_list[["control"]][1,"Count_norm_by_intratel_reads"]))){
  table_merged[is.na(table_merged)] = 0
}

#order according to normalized pattern count
#table_merged = table_merged[order(pmax(table_merged$Count_norm_by_intratel_reads_T,table_merged$Count_norm_by_intratel_reads_C, na.rm = TRUE), decreasing = TRUE),]  
table_merged = table_merged[order(pmax(table_merged$Count_norm_by_intratel_reads_C, na.rm = TRUE), decreasing = TRUE),]  
table_merged = table_merged[order(pmax(table_merged$Count_norm_by_intratel_reads_T, na.rm = TRUE), decreasing = TRUE),]  


#calculate log2 ratio
table_merged$log2_ratio_count_norm_by_intratel_reads = log2((table_merged$Count_norm_by_intratel_reads_T)/ (table_merged$Count_norm_by_intratel_reads_C))
table_merged$log2_ratio_count_per_100_bp_intratel_read = log2((table_merged$Count_per_100_bp_intratel_read_T)/ (table_merged$Count_per_100_bp_intratel_read_C))

#write table
write.table(table_merged, file=paste0(main_path, "/", pid, "_normalized_TVR_counts.tsv"), quote=FALSE, row.names = FALSE, sep="\t")


#----------------------------------------
# add results to summary table
#----------------------------------------

#get summary table
summary_file = paste0(main_path, "/", pid, "_summary.tsv")
summary = read.table(summary_file, header=TRUE, sep="\t")

#add column names
TVR_colnames = paste0(TVRs_for_summary, "_arbitrary_context_norm_by_intratel_reads")
summary[, TVR_colnames] = NA

#TVR normalized counts
row.names(table_merged) = table_merged$Pattern

summary[summary$sample == "tumor", TVR_colnames] = table_merged[TVRs_for_summary, "Count_norm_by_intratel_reads_T"]
summary[summary$sample == "control", TVR_colnames] = table_merged[TVRs_for_summary, "Count_norm_by_intratel_reads_C"]

#replace NA with 0
summary[is.na(summary)] = 0

#save extended summary table
write.table(summary, file=summary_file, quote=FALSE, row.names = FALSE, sep="\t")


