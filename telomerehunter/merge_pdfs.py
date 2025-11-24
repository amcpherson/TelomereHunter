#!/usr/bin/env python

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
from PyPDF2 import PdfFileReader, PdfFileMerger


###############################################################################
### merge all PDF plots produced as TelomereHunter output into one PDF file ###
###############################################################################

def mergeTelomereHunterPDFs(pid, outdir):

    # all possible pdf file names in correct order for merged pdf
    chromosomes = [ str(i) for i in range(1,22+1) ] + ["X","Y"]
    possible_file_names = [pid + '_telomere_content',pid + '_sorted_telomere_read_counts', pid + '_hist_telomere_repeats_per_intratelomeric_read', 
                           pid + '_gc_content', pid + '_TVR_barplot', pid + '_TVR_scatterplot', pid + '_singletons'] + [ pid + "_" + i for i in chromosomes] 
    possible_pdf_names = [ i + '.pdf' for i in possible_file_names ]


    # check which of the possible pdf files exist and sort
    files_dir = os.path.join(outdir, pid, "plots")
    pdf_files = [f for f in os.listdir(files_dir) if f.endswith("pdf")]
    pdf_files_ordered =  [pdf for pdf in possible_pdf_names if pdf in pdf_files]


    # merge files
    merger = PdfFileMerger()

    for filename in pdf_files_ordered:
        merger.append(PdfFileReader(os.path.join(files_dir, filename), "rb"))

    merger.write(os.path.join(outdir, pid, pid + "_all_plots_merged.pdf"))


