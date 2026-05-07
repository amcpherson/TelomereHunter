#!/usr/bin/env python3

# Single-cell adaptation of TelomereHunter sort_telomere_reads
# Classifies telomeric reads into fractions and tracks per-cell counts

import os
import sys
import re
import pysam
from collections import defaultdict
from telomerehunter.filter_telomere_reads import getReverseComplement


def sort_telomere_reads_sc(input_dir, band_file, out_dir, pid, mapq_threshold, repeats,
                           barcode_tag='CB'):
    """
    Sort filtered telomeric reads into fractions (intratelomeric, junction-spanning,
    subtelomeric, intrachromosomal) while tracking per-cell counts.
    
    Returns:
      per_cell_fractions: dict {barcode: {'intratelomeric': n, 'junctionspanning': n,
                                          'subtelomeric': n, 'intrachromosomal': n}}
    """

    #####################
    ### get patterns  ###
    #####################

    patterns = []
    for repeat in repeats:
        patterns.append(repeat)
        patterns.append(getReverseComplement(repeat))

    # open input bam_file for reading
    bamfile = pysam.AlignmentFile(os.path.join(input_dir, pid + "_filtered_name_sorted.bam"), "rb")

    # open filtered files for writing
    intratelomeric_file = pysam.AlignmentFile(os.path.join(out_dir, pid + "_filtered_intratelomeric.bam"), "wb", template=bamfile)
    junctionspanning_file = pysam.AlignmentFile(os.path.join(out_dir, pid + "_filtered_junctionspanning.bam"), "wb", template=bamfile)
    subtelomeric_file = pysam.AlignmentFile(os.path.join(out_dir, pid + "_filtered_subtelomeric.bam"), "wb", template=bamfile)
    intrachromosomal_file = pysam.AlignmentFile(os.path.join(out_dir, pid + "_filtered_intrachromosomal.bam"), "wb", template=bamfile)


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


    ####################################
    ### make band and spectrum lists ###
    ####################################

    bands_list = {c: {"band_name": [], "end": []} for c in chromosome_list + ["unmapped"]}
    spectrum_list = {c: {} for c in chromosome_list}

    for c in chromosome_list:
        spectrum_list[c]["junction1"] = {p: 0 for p in patterns}
        spectrum_list[c]["junction1"]["other"] = 0.0
        spectrum_list[c]["junction1"]["reads_with_pattern"] = 0
        spectrum_list[c]["junction2"] = {p: 0 for p in patterns}
        spectrum_list[c]["junction2"]["reads_with_pattern"] = 0
        spectrum_list[c]["junction2"]["other"] = 0.0

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

            spectrum_list[chrom_name][band_name] = {p: 0 for p in patterns}
            spectrum_list[chrom_name][band_name]["other"] = 0.0
            spectrum_list[chrom_name][band_name]["reads_with_pattern"] = 0
        except:
            pass

    spectrum_list["unmapped"] = {"unmapped": {p: 0 for p in patterns}}
    spectrum_list["unmapped"]["unmapped"]["other"] = 0.0
    spectrum_list["unmapped"]["unmapped"]["reads_with_pattern"] = 0

    bands_list["unmapped"]["band_name"].append("unmapped")
    bands_list["unmapped"]["end"].append(0)


    ###############################
    ### per-cell data structures ###
    ###############################

    per_cell_fractions = defaultdict(lambda: {'intratelomeric': 0, 'junctionspanning': 0,
                                              'subtelomeric': 0, 'intrachromosomal': 0})

    # Per-cell GC distribution for intratelomeric reads (for GC bias plot)
    per_cell_intratel_gc = defaultdict(lambda: defaultdict(int))  # {barcode: {gc_pct: count}}


    #################################################################
    ### go through name sorted BAM file containing telomere reads ###
    #################################################################

    break_flag = False

    while True:
        try:
            read1 = next(bamfile)
        except StopIteration:
            break

        try:
            read2 = next(bamfile)
        except StopIteration:
            _sort_single_read_sc(read1, references, bands_list, chromosome_list, mapq_threshold,
                                spectrum_list, patterns, per_cell_fractions,
                                intratelomeric_file, subtelomeric_file, intrachromosomal_file,
                                barcode_tag, per_cell_intratel_gc)
            break

        # Handle unpaired reads (different query names)
        read1_first_flag = True
        while read1.qname != read2.qname:
            if read1_first_flag:
                _sort_single_read_sc(read1, references, bands_list, chromosome_list, mapq_threshold,
                                    spectrum_list, patterns, per_cell_fractions,
                                    intratelomeric_file, subtelomeric_file, intrachromosomal_file,
                                    barcode_tag, per_cell_intratel_gc)
                try:
                    read1 = next(bamfile)
                except StopIteration:
                    _sort_single_read_sc(read2, references, bands_list, chromosome_list, mapq_threshold,
                                        spectrum_list, patterns, per_cell_fractions,
                                        intratelomeric_file, subtelomeric_file, intrachromosomal_file,
                                        barcode_tag, per_cell_intratel_gc)
                    break_flag = True
                    break
                read1_first_flag = False
            else:
                _sort_single_read_sc(read2, references, bands_list, chromosome_list, mapq_threshold,
                                    spectrum_list, patterns, per_cell_fractions,
                                    intratelomeric_file, subtelomeric_file, intrachromosomal_file,
                                    barcode_tag, per_cell_intratel_gc)
                try:
                    read2 = next(bamfile)
                except StopIteration:
                    _sort_single_read_sc(read1, references, bands_list, chromosome_list, mapq_threshold,
                                        spectrum_list, patterns, per_cell_fractions,
                                        intratelomeric_file, subtelomeric_file, intrachromosomal_file,
                                        barcode_tag, per_cell_intratel_gc)
                    break_flag = True
                    break
                read1_first_flag = True

        if break_flag:
            break

        ### reads with mates - classify the pair
        r1_chrom, r1_band, r1_unmapped, r1_junction_flag, r1_junction = _read_check(
            read1, references, bands_list, chromosome_list, mapq_threshold)
        r2_chrom, r2_band, r2_unmapped, r2_junction_flag, r2_junction = _read_check(
            read2, references, bands_list, chromosome_list, mapq_threshold)

        # Get barcodes
        try:
            bc1 = read1.get_tag(barcode_tag)
        except KeyError:
            bc1 = None
        try:
            bc2 = read2.get_tag(barcode_tag)
        except KeyError:
            bc2 = None

        barcode = bc1 or bc2  # use whichever is available

        # INTRATELOMERIC: both mates unmapped
        if r1_unmapped and r2_unmapped:
            _write_and_count(read1, intratelomeric_file, "unmapped", "unmapped", spectrum_list, patterns)
            _write_and_count(read2, intratelomeric_file, "unmapped", "unmapped", spectrum_list, patterns)
            if barcode:
                per_cell_fractions[barcode]['intratelomeric'] += 2
                _track_gc(read1, barcode, per_cell_intratel_gc)
                _track_gc(read2, barcode, per_cell_intratel_gc)

        # JUNCTION SPANNING: one unmapped, other in terminal band
        elif (r1_unmapped and r2_junction_flag) or (r2_unmapped and r1_junction_flag):
            if r1_unmapped:
                chromosome, band = r2_chrom, r2_junction
            else:
                chromosome, band = r1_chrom, r1_junction
            _write_and_count(read1, junctionspanning_file, chromosome, band, spectrum_list, patterns)
            _write_and_count(read2, junctionspanning_file, chromosome, band, spectrum_list, patterns)
            if barcode:
                per_cell_fractions[barcode]['junctionspanning'] += 2

        # SUBTELOMERIC: both in terminal bands
        elif r1_junction_flag and r2_junction_flag:
            _write_and_count(read1, subtelomeric_file, r1_chrom, r1_band, spectrum_list, patterns)
            _write_and_count(read2, subtelomeric_file, r2_chrom, r2_band, spectrum_list, patterns)
            if barcode:
                per_cell_fractions[barcode]['subtelomeric'] += 2

        # SUBTELOMERIC/INTRACHROMOSOMAL mixed
        elif r1_junction_flag or r2_junction_flag:
            if r1_junction_flag:
                r1_file, r2_file = subtelomeric_file, intrachromosomal_file
                r1_frac, r2_frac = 'subtelomeric', 'intrachromosomal'
            else:
                r1_file, r2_file = intrachromosomal_file, subtelomeric_file
                r1_frac, r2_frac = 'intrachromosomal', 'subtelomeric'
            _write_and_count(read1, r1_file, r1_chrom, r1_band, spectrum_list, patterns)
            _write_and_count(read2, r2_file, r2_chrom, r2_band, spectrum_list, patterns)
            if barcode:
                per_cell_fractions[barcode][r1_frac] += 1
                per_cell_fractions[barcode][r2_frac] += 1

        # INTRACHROMOSOMAL: one mapped, one unmapped (not in terminal band)
        elif r1_unmapped and not r2_junction_flag or r2_unmapped and not r1_junction_flag:
            if r1_unmapped:
                chromosome, band = r2_chrom, r2_band
            else:
                chromosome, band = r1_chrom, r1_band
            _write_and_count(read1, intrachromosomal_file, chromosome, band, spectrum_list, patterns)
            _write_and_count(read2, intrachromosomal_file, chromosome, band, spectrum_list, patterns)
            if barcode:
                per_cell_fractions[barcode]['intrachromosomal'] += 2

        # INTRACHROMOSOMAL: both mapped non-terminal
        else:
            _write_and_count(read1, intrachromosomal_file, r1_chrom, r1_band, spectrum_list, patterns)
            _write_and_count(read2, intrachromosomal_file, r2_chrom, r2_band, spectrum_list, patterns)
            if barcode:
                per_cell_fractions[barcode]['intrachromosomal'] += 2


    ###########################
    ### write spectrum file ###
    ###########################

    spectrum_file_path = os.path.join(out_dir, pid + "_spectrum.tsv")
    with open(spectrum_file_path, "w") as spectrum_file:
        spectrum_file.write("chr\tband\treads_with_pattern")
        for pattern in patterns:
            spectrum_file.write("\t" + pattern)
        spectrum_file.write("\tother\n")

        for chromosome in chromosome_list:
            for band in ["junction1"] + bands_list[chromosome]["band_name"] + ["junction2"]:
                spectrum_file.write("%s\t%s\t%i" % (chromosome, band, spectrum_list[chromosome][band]["reads_with_pattern"]))
                for pattern in patterns:
                    spectrum_file.write("\t" + str(spectrum_list[chromosome][band][pattern]))
                spectrum_file.write("\t" + str(int(round(spectrum_list[chromosome][band]["other"]))) + "\n")

        spectrum_file.write("unmapped\tunmapped\t%i" % (spectrum_list["unmapped"]["unmapped"]["reads_with_pattern"]))
        for pattern in patterns:
            spectrum_file.write("\t" + str(spectrum_list["unmapped"]["unmapped"][pattern]))
        spectrum_file.write("\t" + str(int(round(spectrum_list["unmapped"]["unmapped"]["other"]))) + "\n")


    ##########################
    ### close file handles ###
    ##########################

    bamfile.close()
    intratelomeric_file.close()
    junctionspanning_file.close()
    subtelomeric_file.close()
    intrachromosomal_file.close()

    return per_cell_fractions, per_cell_intratel_gc


def _read_check(read, references, bands_list, chromosome_list, mapq_threshold):
    """Check read mapping status and band position."""
    tid = read.tid
    ref_name = ''
    if tid != -1:
        ref_name = references[tid]
        ref_name = ref_name.replace("chr", "")

    read_pos = read.pos
    read_junctionspanning = False
    read_junction = ""

    if read.is_unmapped or ref_name not in chromosome_list or read.mapq < mapq_threshold:
        return ("unmapped", "unmapped", True, False, "")
    else:
        chromosome = ref_name
        i = 0
        while i < len(bands_list[chromosome]["end"]) - 1 and read_pos > bands_list[chromosome]["end"][i]:
            i += 1
        band = bands_list[chromosome]["band_name"][i]

        if i == 0:
            read_junctionspanning = True
            read_junction = "junction1"
        elif i == len(bands_list[chromosome]["end"]) - 1:
            read_junctionspanning = True
            read_junction = "junction2"

        return (chromosome, band, False, read_junctionspanning, read_junction)


def _write_and_count(read, fraction_file, chromosome, band, spectrum_list, patterns):
    """Write read to fraction BAM and update spectrum counts."""
    fraction_file.write(read)
    spectrum_list[chromosome][band]["reads_with_pattern"] += 1

    read_total_pattern_count = 0
    seq = read.seq
    for pattern in patterns:
        count = seq.count(pattern)
        spectrum_list[chromosome][band][pattern] += count
        read_total_pattern_count += count

    spectrum_list[chromosome][band]["other"] += float(len(seq)) / 6 - read_total_pattern_count


def _track_gc(read, barcode, per_cell_intratel_gc):
    """Track GC content of a read for per-cell intratelomeric GC distribution."""
    seq = read.seq
    if seq is None:
        return
    read_length = len(seq)
    n_count = seq.count('N')
    if float(n_count) / float(read_length) <= 0.2:
        gc_content = int(round(float(seq.count('C') + seq.count('G')) / float(read_length - n_count) * 100))
        per_cell_intratel_gc[barcode][gc_content] += 1


def _sort_single_read_sc(read, references, bands_list, chromosome_list, mapq_threshold,
                         spectrum_list, patterns, per_cell_fractions,
                         intratelomeric_file, subtelomeric_file, intrachromosomal_file,
                         barcode_tag, per_cell_intratel_gc=None):
    """Sort a read without a mate into the correct fraction."""
    chromosome, band, is_unmapped, is_junction, junction = _read_check(
        read, references, bands_list, chromosome_list, mapq_threshold)

    try:
        barcode = read.get_tag(barcode_tag)
    except KeyError:
        barcode = None

    if is_unmapped:
        _write_and_count(read, intratelomeric_file, "unmapped", "unmapped", spectrum_list, patterns)
        if barcode:
            per_cell_fractions[barcode]['intratelomeric'] += 1
            if per_cell_intratel_gc is not None:
                _track_gc(read, barcode, per_cell_intratel_gc)
    elif is_junction:
        _write_and_count(read, subtelomeric_file, chromosome, band, spectrum_list, patterns)
        if barcode:
            per_cell_fractions[barcode]['subtelomeric'] += 1
    else:
        _write_and_count(read, intrachromosomal_file, chromosome, band, spectrum_list, patterns)
        if barcode:
            per_cell_fractions[barcode]['intrachromosomal'] += 1
