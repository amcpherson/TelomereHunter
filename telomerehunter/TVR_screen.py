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


##################################################################################################################################################################
### script loops through filtered BAM file containing intratelomeric reads and searches for patterns of the type XXXGGG (and reverse complement)               ###
### only counts patterns if the base qualities are all greater than 20                                                                                         ###
### pattern counts are written to output tables (with the frequency of the pattern and the average base quality at each position)                              ###
##################################################################################################################################################################


# define functions
def assure_path_exists(path):
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)


def screenBamFile(bam_file, patterns, quals, qual_threshold=20):

    bamfile = pysam.Samfile( bam_file, "rb" )

    pattern = "GGG"
    offset = -3

    for read in bamfile.fetch(until_eof=True):
      seq = read.seq
      indices = [m.start() for m in re.finditer(pattern, seq)]
      
      #get reverse complement of sequence (and quality) if the pattern was not found often enough
      if len(indices)<4:
         seq = getReverseComplement(seq)
         qual = read.qual[::-1]
         indices = [m.start() for m in re.finditer(pattern, seq)]
      else:
         qual = read.qual
      
      for i in indices:
        
        if i+offset < 0:    # skip if pattern does not have 3 bases in front of it
          continue
              
        p = seq[i+offset:i]
        
        if p.find('N')!=-1:    # skip if pattern contains an "N"
          continue
        
        # quality of all positions
        q1 = convertAsciiToPhred(qual[i+offset])
        q2 = convertAsciiToPhred(qual[i+1+offset])
        q3 = convertAsciiToPhred(qual[i+2+offset])
        q4 = convertAsciiToPhred(qual[i+3+offset])
        q5 = convertAsciiToPhred(qual[i+4+offset])
        q6 = convertAsciiToPhred(qual[i+5+offset])
        #print p,qual[i+offset:i], q1, q2, q3  # for debugging
        
        # skip if one of the positions has a base quality lower than the threshold
        if sum(i < qual_threshold for i in [q1, q2, q3, q4, q5, q6])>0:       
          continue
        
        
        try:
          patterns[p] += 1
          quals[p][0].append(q1) # record base qualities
          quals[p][1].append(q2) 
          quals[p][2].append(q3) 
          quals[p][3].append(q4)
          quals[p][4].append(q5)
          quals[p][5].append(q6)
        except KeyError:
          patterns[p] = 1
          quals[p] = ([q1],[q2],[q3],[q4],[q5],[q6]) # create new list for base qualities



def convertAsciiToPhred(asciiSymbol):
    return ord(asciiSymbol)-33 # base quality is encodead as Ascii Symbol



def getReverseComplement(sequence):
    sequence_temp = sequence.replace("A", "1")
    sequence_temp = sequence_temp.replace("C", "2")
    sequence_temp = sequence_temp.replace("G", "3")
    sequence_temp = sequence_temp.replace("T", "4")

    sequence_temp2 = sequence_temp.replace("1", "T")
    sequence_temp2 = sequence_temp2.replace("2", "G")
    sequence_temp2 = sequence_temp2.replace("3", "C")
    sequence_temp2 = sequence_temp2.replace("4", "A")

    sequence_reverse_complement = sequence_temp2[::-1]  # reverses a string

    return sequence_reverse_complement


    
def invertPatterns(patterns):
    """Inverts keys and values in dictionary"""
    p2 = {}
    keys = patterns.keys()
    for key in keys:
       p2[(patterns[key],key)] = key
    return p2

def showPatterns(patterns,quals):
  """Sorts and displays patterns and qualites"""
#  p2 = {}
  keys = patterns.keys()
  keys.sort()
  keys.reverse()
  output = "\t".join(["Pattern", "Count", "Frequency_in_Percent", "Avg_Qual_pos1", "Avg_Qual_pos2", "Avg_Qual_pos3", "Avg_Qual_pos4", "Avg_Qual_pos5", "Avg_Qual_pos6"]) + "\n"
  counts = 0
  for key in keys:
    counts += key[0]
  for key in keys:
#    if patterns[key].find('N')==-1:         # don't output patterns with N
      (qs1,qs2,qs3,qs4,qs5,qs6) = quals[patterns[key]] # lists of qualites for position 1, 2 and 3
      q1_mean = sum(qs1)/len(qs1)
      q2_mean = sum(qs2)/len(qs2)
      q3_mean = sum(qs3)/len(qs3)
      q4_mean = sum(qs4)/len(qs4)
      q5_mean = sum(qs5)/len(qs5)
      q6_mean = sum(qs6)/len(qs6)
      output += "\t".join( [key[1] + "GGG" ,str(key[0]), str(key[0]*100.0/counts), str(q1_mean), str(q2_mean), str(q3_mean), str(q4_mean), str(q5_mean), str(q6_mean)] ) +"\n"
  return output



def TVR_screen(main_path, pid, sample, min_base_quality):
  
  path_bam_intratelomeric = main_path + "/" + pid + "_filtered_intratelomeric.bam"
  
  patterns = {}     #directory for counting of patterns
  quals = {}    #directory: patterns are keys, for each pattern 3 lists containing qualities at position 1, 2 and 3
  screenBamFile(path_bam_intratelomeric, patterns, quals, qual_threshold=min_base_quality)
  p2 = invertPatterns(patterns) # invert keys and values

  outfile_path = main_path + "/TVRs/" + pid + "_" + sample + "_TVRs.txt"
  assure_path_exists(outfile_path)
  
  #print showPatterns(p2,quals)
  outfile = open( outfile_path, "w")  
  outfile.write(showPatterns(p2,quals))        # display sorted list
  outfile.close









