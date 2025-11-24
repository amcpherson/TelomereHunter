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


#############################################################################################
### Returns the repeat threshold for the input BAM file:                                  ###
### threshold = floor(read length*repeat_threshold_per_100_bp/100).                       ###
### If the read lengths of the first 100 reads differ, it returns "n".                    ###
#############################################################################################

import os
import sys
import pysam



# ----------------------------------------------------------------
# get all read lengths occurring in BAM file (first 1000 reads)
# (without secondary or supplementary alignments)
# ----------------------------------------------------------------
def get_read_lengths(bam_file, reads_to_parse=1000):

    # open input bam_file for reading
    bamfile = pysam.Samfile( bam_file, "rb" )

    # print unique read lengths of the first N non-supplementary or secondary alignments
    cntr=0
    read_lengths=[]

    for read in bamfile.fetch(until_eof=True):   

        if read.is_secondary:        #skips all secondary alignments
            continue

        if read.flag >= 2048:        # skip supplementary alignments
            continue

        read_lengths.append(len(read.seq))

        cntr+=1
        if cntr == reads_to_parse:    
            break
  
    read_lengths = sorted(list(set(read_lengths)))

    read_lengths_str = ",".join(str(i) for i in read_lengths)

    return read_lengths_str





def get_repeat_threshold(read_lengths_str, repeat_threshold_per_100_bp):

    read_lengths = read_lengths_str.split(",")

    read_lengths = [int(i) for i in read_lengths]

    repeat_thresholds = [int(round(float(i)*repeat_threshold_per_100_bp/100)) for i in read_lengths]

    repeat_thresholds = sorted(list(set(repeat_thresholds)))

    if len(repeat_thresholds)==1:
        repeat_threshold = int(repeat_thresholds[0])
        repeat_threshold_str = repeat_threshold
        print("Read length is " + read_lengths_str + ". Repeat threshold is set to " + str(repeat_threshold_str) + ".")


    elif len(repeat_thresholds)>1:
        print "Read lengths in sample differ: repeat threshold will be set individually for each read."
        repeat_threshold_str = ",".join(str(i) for i in repeat_thresholds)
        repeat_threshold = "n"

    else:
        print("Error in calculating the repeat threshold.")

    return [repeat_threshold, repeat_threshold_str]

