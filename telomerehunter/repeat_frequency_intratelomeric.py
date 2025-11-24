#!/usr/bin/python

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

import os
import sys
import re
import pysam
from telomerehunter.filter_telomere_reads import getReverseComplement

###############################################################################################
### get the distribution of telomeric repeats per intratelomeric read in the input BAM file ###
###############################################################################################

def repeat_frequency_intratelomeric(input_path, out_dir, pid, repeats):
	################################################
	### get patterns and make regular expression ###
	################################################

	patterns_regex_forward = ""
	patterns_regex_reverse = ""
		
	for repeat in repeats:
		patterns_regex_forward += repeat + "|"
		patterns_regex_reverse += getReverseComplement(repeat) + "|"

	# remove last '|' from regex
	patterns_regex_forward=patterns_regex_forward[:-1]
	patterns_regex_reverse=patterns_regex_reverse[:-1]

		
	#########################
	### open file handles ###
	#########################

	# open input bam_file for reading
	bamfile = pysam.Samfile( input_path + "/" + pid + "_filtered_intratelomeric.bam", "rb" )


	##################################
	### initialize frequency table ###
	##################################

	frequency_table = { repeats:0 for repeats in range(0,2+1) }

	######################################
	### loop through filtered BAM file ###
	######################################

	for read in bamfile.fetch(until_eof=True):   

		sequence = read.seq

		number_repeats_forward = len(re.findall(patterns_regex_forward, sequence))
		number_repeats_reverse = len(re.findall(patterns_regex_reverse, sequence))
		
		if number_repeats_forward > number_repeats_reverse:
			number_repeats = number_repeats_forward
		else:
			number_repeats = number_repeats_reverse

		try:
		      frequency_table[number_repeats] += 1
		except:
		      frequency_table[number_repeats] = 1       # if key does not exist: add to frequency_table


	##################################
	### write frequency table file ###
	##################################

	# open frequency table file for writing
	frequency_table_file = open(out_dir + "/" + pid + "_repeat_frequency_per_intratelomeric_read.tsv", "w")

	# header
	frequency_table_file.write( "number_repeats\tcount\n")

	#write line for each chromosome band
	for frequency in frequency_table:

		frequency_table_file.write( "%i\t%i\n" % (frequency, frequency_table[frequency]))


	##########################
	### close file handles ###
	##########################

	bamfile.close()
	frequency_table_file.close()


