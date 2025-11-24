# usage: R --no-save --slave --args <MAIN_PATH> <PID> <CONTEXT> <TVRs_for_context> < TVR_context_summary_tables.R

# Copyright 2018 Lina Sieverling, Lars Feuerbach

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


#################################################################################################
###
### makes tables with top contexts for each selected TVR 
### and with the number of TVR singletons (i.e. TVRs in t-type context)
###
#################################################################################################

# get commandline arguments
commandArgs = commandArgs()
main_path = commandArgs[5]
pid = commandArgs[6]
context = as.numeric(commandArgs[7])
TVRs_for_context = commandArgs[8]


library('reshape2', quietly=TRUE, warn.conflicts=FALSE)

TVRs_for_context = unlist(strsplit(TVRs_for_context, split=","))

hexamer_repeats = context/6

#----------------------------------------------------------------------------------------------------
# get telomere content
#----------------------------------------------------------------------------------------------------

summary_file = paste0(main_path, "/", pid, "_summary.tsv")

summary = read.table(file=summary_file,
                                 header=TRUE, sep="\t", stringsAsFactors=FALSE)

# which samples were run?
samples = summary$sample

#----------------------------------------------------------------------------------------------------
# get pattern neighborhood results (top context and t-type context X bp before and after pattern)
#----------------------------------------------------------------------------------------------------


df_list_all=list()


for (pattern in TVRs_for_context){

  top_contexts = data.frame()
  
  for (sample in samples){

    TVR_context_dir = paste0(main_path, "/", sample, "_TelomerCnt_", pid, "/TVR_context/") 
    
    neighborhood_table_file = paste0(TVR_context_dir, "/", pid, "_", sample, "_", pattern, "_", context, "bp_", context, "bp_neighborhood.tsv")
    
    #if(!file.exists(neighborhood_table_file)){next}
    
    neighborhood_table = read.table(file=neighborhood_table_file,
                                    header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    #get top contexts
    row_name = paste(pid,sample,context, sep="_")      
    top_contexts[row_name,"PID"] = pid
    top_contexts[row_name,"Sample"] = sample
    top_contexts[row_name,"Context_bp"] = context
    top_contexts[row_name, "Bases"] = neighborhood_table[1, "Bases"]
    top_contexts[row_name, "Count"] = neighborhood_table[1, "Count"]
    top_contexts[row_name, "Percent"] = neighborhood_table[1, "Percent"]   
    
    #get singleton counts
    pattern_t_type_context = paste0(paste0(rep("TTAGGG", hexamer_repeats), collapse = ''), "-", pattern, "-", paste0(rep("TTAGGG", hexamer_repeats), collapse = ''))
    
    if (sum(neighborhood_table$Bases==pattern_t_type_context)!=0){
      top_contexts[row_name, "Bases_t_type"] = neighborhood_table[neighborhood_table$Bases==pattern_t_type_context, "Bases"]
      top_contexts[row_name, "Count_t_type"] = neighborhood_table[neighborhood_table$Bases==pattern_t_type_context, "Count"]
      top_contexts[row_name, "Percent_t_type"] = neighborhood_table[neighborhood_table$Bases==pattern_t_type_context, "Percent"]         
    } else{
      top_contexts[row_name, "Bases_t_type"] = NA
      top_contexts[row_name, "Count_t_type"] = NA
      top_contexts[row_name, "Percent_t_type"] = NA
    }             
  }
  
  # add telomere content
  if("tumor" %in% samples & "control" %in% samples){
    tel_content_tumor = summary[summary$sample=="tumor", "tel_content"]
    tel_content_control = summary[summary$sample=="control", "tel_content"]

    top_contexts$tel_content_log2_ratio = log2(tel_content_tumor/tel_content_control)
  }else{
    top_contexts$tel_content_log2_ratio = NA
  }

  #save
  df_list_all[[pattern]] = top_contexts
  
}

pattern_contexts = data.frame()

for(pattern in TVRs_for_context){
  df_pattern = df_list_all[[pattern]]
  df_pattern$pattern = pattern
  pattern_contexts = rbind(pattern_contexts, df_pattern)
}


#----------------------------------------------------------------------------------------------------
# get counts of patterns in t-type context
#----------------------------------------------------------------------------------------------------

# get t-type count per sample
singleton_table = pattern_contexts[, c("PID", "pattern", "tel_content_log2_ratio", "Sample", "Count_t_type")]
singleton_table = dcast(data = singleton_table,formula = PID + pattern + tel_content_log2_ratio ~Sample, value.var = "Count_t_type")

if("control" %in% samples){
  colnames(singleton_table)[colnames(singleton_table)=="control"] = "singleton_count_control"
}else{
  singleton_table$singleton_count_control = NA
}

if("tumor" %in% samples){
  colnames(singleton_table)[colnames(singleton_table)=="tumor"] = "singleton_count_tumor"
}else{
  singleton_table$singleton_count_tumor = NA
}


# get log2 ratio of counts
singleton_table$singleton_count_log2_ratio = log2(singleton_table$singleton_count_tumor/singleton_table$singleton_count_control)


#----------------------------------------------------------------------------------------------------
# normalize by total number of reads
#----------------------------------------------------------------------------------------------------

if("control" %in% samples){

  total_reads_control = summary[summary$sample=="control", "total_reads"]
  singleton_table$singleton_count_control_norm = singleton_table$singleton_count_control/total_reads_control

}else{

 singleton_table$singleton_count_control_norm = NA

}


if("tumor" %in% samples){

  total_reads_tumor = summary[summary$sample=="tumor", "total_reads"]
  singleton_table$singleton_count_tumor_norm = singleton_table$singleton_count_tumor/total_reads_tumor

}else{

 singleton_table$singleton_count_tumor_norm = NA

}


# get log2 ratio of normalized counts
singleton_table$singleton_count_log2_ratio_norm = log2(singleton_table$singleton_count_tumor_norm/singleton_table$singleton_count_control_norm)


#----------------------------------------------------------------------------------------------------
# difference to expected pattern occurrence (=diagonal)
#----------------------------------------------------------------------------------------------------

singleton_table$distance_to_expected_singleton_log2_ratio =  singleton_table$singleton_count_log2_ratio_norm - singleton_table$tel_content_log2_ratio


#----------------------------------------------------------------------------------------------------
# save tables
#----------------------------------------------------------------------------------------------------

# set missing counts to zero
if("control" %in% samples){
  singleton_table[is.na(singleton_table$singleton_count_control), "singleton_count_control"] = 0
  singleton_table[is.na(singleton_table$singleton_count_control_norm), "singleton_count_control_norm"] = 0
}

if("tumor" %in% samples){
  singleton_table[is.na(singleton_table$singleton_count_tumor), "singleton_count_tumor"] = 0
  singleton_table[is.na(singleton_table$singleton_count_tumor_norm), "singleton_count_tumor_norm"] = 0
}

#arrange columns
singleton_table = singleton_table[ , c("PID", "pattern", "singleton_count_tumor", "singleton_count_control", "singleton_count_log2_ratio", 
                                        "singleton_count_tumor_norm", "singleton_count_control_norm", "singleton_count_log2_ratio_norm", 
                                       "tel_content_log2_ratio", "distance_to_expected_singleton_log2_ratio")]


#arrange top context columns
patterns_top_contexts = pattern_contexts[ , c("PID", "Sample", "pattern", "Context_bp", "Bases", "Count", "Percent")]


#save
write.table(patterns_top_contexts, paste0(main_path, "/", pid, "_TVR_top_contexts.tsv"), 
            sep="\t", quote=FALSE, row.names=FALSE)

write.table(singleton_table, paste0(main_path, "/", pid, "_singletons.tsv"), 
            sep="\t", quote=FALSE, row.names=FALSE)




#----------------------------------------
# add results to summary table
#----------------------------------------

#add column names
singleton_colnames = paste0(TVRs_for_context, "_singletons_norm_by_all_reads")
summary[, singleton_colnames] = NA

#TVR normalized counts
row.names(singleton_table) = singleton_table$pattern

if("tumor" %in% samples){
  summary[summary$sample == "tumor", singleton_colnames] = singleton_table[TVRs_for_context, "singleton_count_tumor_norm"]
}

if("control" %in% samples){
  summary[summary$sample == "control", singleton_colnames] = singleton_table[TVRs_for_context, "singleton_count_control_norm"]
}

#save extended summary table
write.table(summary, file=summary_file, quote=FALSE, row.names = FALSE, sep="\t")

