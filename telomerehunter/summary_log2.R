# usage: R --no-save --slave --args <MAIN_PATH> <PID> < summary_log2.R

# description: - calculates the log2 ratio for the summary table


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

#################################################################################################

#get summary table
summary_file = paste0(main_path, "/", pid, "_summary.tsv")

summary = read.table(file=summary_file,
                     header=TRUE, sep="\t", stringsAsFactors=FALSE)

# check which samples were run
samples = summary$sample


# if tumor and control were run, calculate summary 
if("tumor" %in% samples & "control" %in% samples){
  summary[3, 11:dim(summary)[2]] = log2(summary[1,11:dim(summary)[2]] /summary[2,11:dim(summary)[2]])
  summary[3, "PID"] = pid
  summary[3, "sample"] = "log2(tumor/control)"
  
  summary[3, 3:10] = ""
  
  #save extended summary table
  write.table(summary, file=summary_file, quote=FALSE, row.names = FALSE, sep="\t")
}
