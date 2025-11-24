#!/usr/bin/python

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


import os
import sys
import re
import pysam

import telomerehunter


#####################################################################################
###
### get the context of telomere variant repeats (TVRs)
###
#####################################################################################



#####################################################################################################

def invertDictionary(dictionary):
    #Inverts keys and values in dictionary
    dicInvert = {}
    keys = dictionary.keys()
    for key in keys:
       dicInvert[(dictionary[key],key)] = key
    return dicInvert

def Dictionary2Table(dictionary, outfile_path, cutoff):
    # write sorted dictionary as a table to file

    dictionaryInv = invertDictionary(dictionary)

    table = "\t".join(["Bases", "Count", "Percent"]) + "\n"

    keys = dictionaryInv.keys()
    keys.sort()
    keys.reverse()

    total_counts = sum(dictionary.values())

    for key in keys:

        count = key[0]
        percent = float(count)/total_counts*100

        if count>=cutoff:
            table += key[1] + "\t" + str(count) + "\t" + str(round(percent, 2)) + "\n"


    outfile = open( outfile_path, "w")  
    outfile.write(table)        
    outfile.close


#####################################################################################################

#main_path: Path to TelomereHunter results
#pattern: Pattern for which to get context
#min_base_quality: Minimum base quality required for pattern
#context_before: Number of bases before start of pattern
#context_after: Number of bases after end of pattern
#telomere_pattern: Pattern with which to identify G-rich reads
#cutoff: Count cutoff for diplaying neighborhood in output table
#tel_file: the telomere file in which to search for TVRs (filtered, filtered_intratelomeric, ...)

def TVR_context(main_path,           
                pid, 
                sample, 
                pattern, 
                min_base_quality=20, 
                context_before=18, 
                context_after=18, 
                telomere_pattern="GGG", 
                cutoff=0,
                tel_file="filtered_intratelomeric"):

    outdir = main_path + "/TVR_context"

    if not os.path.exists(outdir):
        os.makedirs(outdir)


    bam_file = main_path + "/" + pid + "_" + tel_file + ".bam"

    bamfile = pysam.Samfile( bam_file, "rb" )


    #directories for counting of patterns in neighborhood
    neighborhood_before = {}       
    neighborhood_after = {} 
    neighborhood = {} 

    for read in bamfile.fetch(until_eof=True):
            
        seq = read.seq
        qual = read.qual

        indices_telomeric = [m.start() for m in re.finditer(telomere_pattern, seq)]


        #get reverse complement of sequence if the telomeric pattern was not found often enough
        if len(indices_telomeric)<6:
             seq = telomerehunter.TVR_screen.getReverseComplement(seq)
             qual = read.qual[::-1]

        indices_pattern = [m.start() for m in re.finditer(pattern, seq)]

       
        for i in indices_pattern:

            bases_before=''
            bases_after=''

            # get base qualities of pattern
            q1 = telomerehunter.TVR_screen.convertAsciiToPhred(qual[i])
            q2 = telomerehunter.TVR_screen.convertAsciiToPhred(qual[i+1])
            q3 = telomerehunter.TVR_screen.convertAsciiToPhred(qual[i+2])
            q4 = telomerehunter.TVR_screen.convertAsciiToPhred(qual[i+3])
            q5 = telomerehunter.TVR_screen.convertAsciiToPhred(qual[i+4])
            q6 = telomerehunter.TVR_screen.convertAsciiToPhred(qual[i+5])

            # skip if one of the positions has a base quality lower than the threshold
            if sum(i < min_base_quality for i in [q1, q2, q3, q4, q5, q6])>0:       
                continue

            if i-context_before >= 0:
                bases_before = seq[i-context_before:i]

                try:
                    neighborhood_before[bases_before] += 1
                except KeyError:
                    neighborhood_before[bases_before] = 1


            if i+len(pattern)+context_after <= len(seq):
                bases_after = seq[i+len(pattern):i+len(pattern)+context_after]

                try:
                    neighborhood_after[bases_after] += 1
                except KeyError:
                    neighborhood_after[bases_after] = 1


            if bases_before and bases_after:
                bases_all = bases_before + "-" + pattern + "-" + bases_after

                try:
                    neighborhood[bases_all] += 1
                except KeyError:
                    neighborhood[bases_all] = 1


    Dictionary2Table(neighborhood_before, outdir + "/" + pid + "_" + sample + "_" + pattern + "_" + str(context_before) + "bp_neighborhood_before.tsv", cutoff=cutoff)
    Dictionary2Table(neighborhood_after, outdir + "/" + pid + "_" + sample + "_" + pattern + "_" + str(context_after) + "bp_neighborhood_after.tsv", cutoff=cutoff)
    Dictionary2Table(neighborhood, outdir + "/" + pid + "_" + sample + "_" + pattern + "_" + str(context_before) + "bp_" + str(context_after) + "bp_neighborhood.tsv", cutoff=cutoff)

