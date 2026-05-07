#!/usr/bin/env python3

# Single-cell adaptation of TelomereHunter filter_telomere_reads
# Processes a barcoded BAM (CB tag) and tracks per-cell statistics

import os
import sys
import subprocess
import re
import pysam
from collections import defaultdict
from telomerehunter.filter_telomere_reads import getReverseComplement


def filter_telomere_reads_sc(bam_file, band_file, out_dir, pid, repeat_threshold_calc,
                             repeat_threshold_set, mapq_threshold, repeats, consecutive_flag,
                             remove_duplicates, barcode_tag='CB', min_reads=500):
    """
    Single-pass filter of a barcoded BAM file.
    
    For each read:
      - Extracts cell barcode from CB tag
      - Tracks per-cell total read count and GC histogram
      - Writes telomeric reads (with CB preserved) to filtered BAM
    
    Returns:
      per_cell_gc: dict {barcode: {0..100: count}}
      per_cell_total_reads: dict {barcode: int}
      per_cell_tel_reads: dict {barcode: int}  (reads passing telomere threshold)
    """

    ################################################
    ### get patterns and make regular expression ###
    ################################################

    patterns_regex_forward = ""
    patterns_regex_reverse = ""

    for repeat in repeats:
        patterns_regex_forward += repeat + "|"
        patterns_regex_reverse += getReverseComplement(repeat) + "|"

    patterns_regex_forward = patterns_regex_forward[:-1]
    patterns_regex_reverse = patterns_regex_reverse[:-1]


    #########################
    ### open file handles ###
    #########################

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    filtered_file_path = os.path.join(out_dir, pid + "_filtered.bam")
    filtered_file = pysam.AlignmentFile(filtered_file_path, "wb", template=bamfile)


    ############################
    ### make chromosome list ###
    ############################

    references = bamfile.references

    if len(references) > 0 and references[0][0:3] == "chr":
        bam_chr_prefix = "chr"
    else:
        bam_chr_prefix = ''

    chromosome_list = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chromosome_list_with_prefix = [bam_chr_prefix + c for c in chromosome_list]


    ################################################
    ### make band list (for read count tracking) ###
    ################################################

    bands_list = {c: {"band_name": [], "end": []} for c in chromosome_list + ["unmapped"]}
    spectrum_list = {c: {} for c in chromosome_list}

    for line in open(band_file, "r"):
        try:
            parts = line.rstrip().split()
            end = parts[2]
            band_name = parts[3]

            if parts[0][:3] == "chr":
                chrom_name = parts[0][3:]
            else:
                chrom_name = parts[0]

            bands_list[chrom_name]["band_name"].append(band_name)
            bands_list[chrom_name]["end"].append(int(end))
            spectrum_list[chrom_name][band_name] = {"reads": 0}
        except:
            pass

    spectrum_list["unmapped"] = {"unmapped": {"reads": 0}}
    bands_list["unmapped"]["band_name"].append("unmapped")
    bands_list["unmapped"]["end"].append(0)


    ###############################
    ### per-cell data structures ###
    ###############################

    per_cell_gc = defaultdict(lambda: defaultdict(int))        # {barcode: {gc_pct: count}}
    per_cell_total_reads = defaultdict(int)                    # {barcode: total_read_count}
    per_cell_tel_reads = defaultdict(int)                      # {barcode: telomeric_read_count}

    # Bulk GC (for fallback normalization)
    bulk_gc = defaultdict(int)

    # For band-based tracking
    chr_offset = len(bam_chr_prefix)
    lastChromosome = ''
    chromosomeLsEnd = None
    chromosomeLsBand = None
    i = 0


    #############################
    ### loop through BAM file ###
    #############################

    reads_processed = 0

    for read in bamfile.fetch(until_eof=True):

        if read.is_secondary:
            continue
        if remove_duplicates and read.is_duplicate:
            continue
        if read.flag >= 2048:  # skip supplementary alignments
            continue

        sequence = read.seq
        try:
            read_length = len(sequence)
        except TypeError:
            continue

        # Get cell barcode
        try:
            barcode = read.get_tag(barcode_tag)
        except KeyError:
            continue  # skip reads without barcode

        per_cell_total_reads[barcode] += 1
        reads_processed += 1

        # GC content
        n_count = sequence.count('N')
        if float(n_count) / float(read_length) <= 0.2:
            gc_content = int(round(float(sequence.count('C') + sequence.count('G')) / float(read_length - n_count) * 100))
            per_cell_gc[barcode][gc_content] += 1
            bulk_gc[gc_content] += 1

        # Track per-band read counts (bulk, for spectrum file compatibility)
        tid = read.tid
        ref_name = ''
        if tid != -1:
            ref_name = references[tid]

        if read.is_unmapped or ref_name not in chromosome_list_with_prefix or read.mapq < mapq_threshold:
            spectrum_list["unmapped"]["unmapped"]["reads"] += 1
        else:
            chromosome = ref_name[chr_offset:]
            if chromosome != lastChromosome:
                chromosomeLsEnd = bands_list[chromosome]["end"]
                chromosomeLsBand = bands_list[chromosome]["band_name"]
                lastChromosome = chromosome
                i = 0

            read_start_pos = read.pos
            while i < len(chromosomeLsEnd) - 1 and read_start_pos > chromosomeLsEnd[i]:
                i += 1

            band = chromosomeLsBand[i]
            spectrum_list[chromosome][band]["reads"] += 1

        # Check telomere repeat threshold
        if repeat_threshold_calc == 'n':
            repeat_threshold = int(round(float(read_length) * repeat_threshold_set / 100))
        else:
            repeat_threshold = repeat_threshold_calc

        if consecutive_flag:
            is_telomeric = (re.search("(" + patterns_regex_forward + "){" + str(repeat_threshold) + "}", sequence) or
                           re.search("(" + patterns_regex_reverse + "){" + str(repeat_threshold) + "}", sequence))
        else:
            is_telomeric = (len(re.findall(patterns_regex_forward, sequence)) >= repeat_threshold or
                           len(re.findall(patterns_regex_reverse, sequence)) >= repeat_threshold)

        if is_telomeric:
            per_cell_tel_reads[barcode] += 1
            filtered_file.write(read)


    #############################
    ### write read count file ###
    #############################

    readcount_file_path = os.path.join(out_dir, pid + "_readcount.tsv")
    with open(readcount_file_path, "w") as readcount_file:
        readcount_file.write("chr\tband\treads\n")
        for chromosome in chromosome_list + ["unmapped"]:
            for band in bands_list[chromosome]["band_name"]:
                readcount_file.write("%s\t%s\t%i\n" % (chromosome, band, spectrum_list[chromosome][band]["reads"]))


    ##########################
    ### close file handles ###
    ##########################

    bamfile.close()
    filtered_file.close()

    ############################
    ### index filtered file  ###
    ############################

    pysam.index(filtered_file_path)

    ##################################
    ### sort filtered file by name ###
    ##################################

    name_sorted_path = os.path.join(out_dir, pid + "_filtered_name_sorted.bam")
    try:
        subprocess.check_call("samtools sort -n " + filtered_file_path + " -o " + name_sorted_path,
                             shell=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        subprocess.call("samtools sort -n " + filtered_file_path + " -o " + name_sorted_path, shell=True)

    sys.stdout.write("Processed %d reads from %d cells\n" % (reads_processed, len(per_cell_total_reads)))

    return per_cell_gc, per_cell_total_reads, per_cell_tel_reads, bulk_gc
