#!/usr/bin/python

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

import os
import sys
import pysam

#########################################################################
### gives the sum of all intratelomeric read lengths (base pairs)     ###
### secondary and supplementary alignments are skipped                ###
#########################################################################

def summed_intratelomeric_read_length(main_path, pid, sample):

	outfile = open(main_path + "/TVRs/" + pid + "_" + sample + "_summed_read_length.tsv", "w")

	# header
	outfile.write("PID\tsample\tsummed_intratel_read_length\n")

	# open input bam_file for reading
	bam_file = main_path + "/" + pid + "_filtered_intratelomeric.bam"
	bamfile = pysam.Samfile( bam_file, "rb" )

	summed_read_length = 0

	for read in bamfile.fetch(until_eof=True):   

	  if read.is_secondary:        #skips all secondary alignments
	      continue

	  if read.flag >= 2048:        # skip supplementary alignments
	      continue
	    
	  summed_read_length += len(read.seq)

	outfile.write(pid + "\t" + sample + "\t" + str(summed_read_length) + "\n")

	outfile.close()