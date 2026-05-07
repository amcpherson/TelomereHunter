#!/usr/bin/env python3

# Single-cell adaptation of TelomereHunter estimate_telomere_content
# Computes per-cell telomere content with GC correction

import os
import sys
import pysam
from collections import defaultdict


def estimate_telomere_content_sc(per_cell_gc, per_cell_total_reads, per_cell_tel_reads,
                                 per_cell_fractions, per_cell_intratel_gc, bulk_gc,
                                 out_dir, pid,
                                 read_length, repeat_threshold_set, per_read_length,
                                 repeat_threshold_calc, gc_lower, gc_upper, min_reads=500,
                                 min_gc_reads=100):
    """
    Estimate per-cell telomere content using per-cell GC-corrected normalization.
    
    For each cell:
      tel_content = intratelomeric_reads / reads_in_GC_bins * 1,000,000
    
    Uses per-cell GC distribution for correction. Falls back to bulk GC
    (scaled to cell's total reads) only if cell has < min_gc_reads in window.
    
    Also writes per-cell GC distribution files for all reads and intratelomeric
    reads, for downstream GC bias plotting.
    """

    gc_bins = range(gc_lower, gc_upper + 1)

    # Compute bulk GC total in correction window
    bulk_gc_in_window = sum(bulk_gc.get(gc, 0) for gc in gc_bins)
    bulk_total_reads = sum(bulk_gc.values())

    per_cell_results = []

    for barcode in sorted(per_cell_total_reads.keys()):
        total_reads = per_cell_total_reads[barcode]

        # Skip cells below minimum read threshold
        if total_reads < min_reads:
            continue

        tel_reads = per_cell_tel_reads.get(barcode, 0)
        fractions = per_cell_fractions.get(barcode, {})
        intratel_reads = fractions.get('intratelomeric', 0)
        junction_reads = fractions.get('junctionspanning', 0)
        subtel_reads = fractions.get('subtelomeric', 0)
        intrachrom_reads = fractions.get('intrachromosomal', 0)

        # Per-cell GC reads in correction window
        cell_gc = per_cell_gc.get(barcode, {})
        cell_gc_in_window = sum(cell_gc.get(gc, 0) for gc in gc_bins)

        # Decide normalization strategy
        if cell_gc_in_window >= min_gc_reads:
            # Use per-cell GC correction
            gc_denominator = cell_gc_in_window
            gc_method = 'per_cell'
        else:
            # Fall back to bulk GC, scaled to cell's total reads
            if bulk_total_reads > 0:
                gc_denominator = int(round(float(bulk_gc_in_window) / float(bulk_total_reads) * total_reads))
            else:
                gc_denominator = 0
            gc_method = 'bulk'

        # Compute telomere content
        if gc_denominator > 0:
            tel_content = float(intratel_reads) / float(gc_denominator) * 1000000
        else:
            tel_content = 0.0

        per_cell_results.append({
            'barcode': barcode,
            'total_reads': total_reads,
            'tel_reads': tel_reads,
            'intratel_reads': intratel_reads,
            'junction_reads': junction_reads,
            'subtel_reads': subtel_reads,
            'intrachrom_reads': intrachrom_reads,
            'gc_reads_in_window': cell_gc_in_window,
            'gc_denominator_used': gc_denominator,
            'gc_correction_method': gc_method,
            'tel_content': tel_content,
        })


    ##################################
    ### write per-cell results TSV ###
    ##################################

    output_path = os.path.join(out_dir, pid + "_per_cell_telomere_content.tsv")
    with open(output_path, "w") as f:
        header = ["cell_barcode", "total_reads", "tel_reads", "intratel_reads",
                  "junction_reads", "subtel_reads", "intrachrom_reads",
                  "gc_reads_in_window", "gc_denominator_used",
                  "gc_correction_method", "tel_content"]
        f.write("\t".join(header) + "\n")

        for cell in per_cell_results:
            f.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%.6f\n" % (
                cell['barcode'],
                cell['total_reads'],
                cell['tel_reads'],
                cell['intratel_reads'],
                cell['junction_reads'],
                cell['subtel_reads'],
                cell['intrachrom_reads'],
                cell['gc_reads_in_window'],
                cell['gc_denominator_used'],
                cell['gc_correction_method'],
                cell['tel_content'],
            ))

    sys.stdout.write("Wrote per-cell telomere content for %d cells to %s\n" % (len(per_cell_results), output_path))


    ################################################
    ### write per-cell GC distribution TSVs      ###
    ### (for GC bias plotting)                   ###
    ################################################

    # All reads: per-cell GC distribution
    gc_all_path = os.path.join(out_dir, pid + "_per_cell_gc_content.tsv")
    with open(gc_all_path, "w") as f:
        f.write("cell_barcode\tgc_content_percent\tread_count\n")
        for barcode in sorted(per_cell_total_reads.keys()):
            if per_cell_total_reads[barcode] < min_reads:
                continue
            cell_gc = per_cell_gc.get(barcode, {})
            for gc in range(0, 101):
                count = cell_gc.get(gc, 0)
                if count > 0:
                    f.write("%s\t%d\t%d\n" % (barcode, gc, count))

    # Intratelomeric reads: per-cell GC distribution
    gc_intratel_path = os.path.join(out_dir, pid + "_per_cell_intratel_gc_content.tsv")
    with open(gc_intratel_path, "w") as f:
        f.write("cell_barcode\tgc_content_percent\tread_count\n")
        for barcode in sorted(per_cell_total_reads.keys()):
            if per_cell_total_reads[barcode] < min_reads:
                continue
            cell_gc = per_cell_intratel_gc.get(barcode, {})
            for gc in range(0, 101):
                count = cell_gc.get(gc, 0)
                if count > 0:
                    f.write("%s\t%d\t%d\n" % (barcode, gc, count))

    # Bulk GC (aggregated across all cells) - for reference line in plot
    gc_bulk_path = os.path.join(out_dir, pid + "_bulk_gc_content.tsv")
    with open(gc_bulk_path, "w") as f:
        f.write("gc_content_percent\tread_count\n")
        for gc in range(0, 101):
            f.write("%d\t%d\n" % (gc, bulk_gc.get(gc, 0)))


    ##################################
    ### write bulk summary file    ###
    ##################################

    summary_path = os.path.join(out_dir, pid + "_summary.tsv")
    total_all = sum(per_cell_total_reads.values())
    tel_all = sum(per_cell_tel_reads.values())
    intratel_all = sum(f.get('intratelomeric', 0) for f in per_cell_fractions.values())

    if bulk_gc_in_window > 0:
        bulk_tel_content = float(intratel_all) / float(bulk_gc_in_window) * 1000000
    else:
        bulk_tel_content = 0.0

    if repeat_threshold_calc == "n":
        repeat_threshold_calc = "heterogeneous"
    if per_read_length:
        repeat_threshold_display = str(repeat_threshold_set) + " per 100 bp"
    else:
        repeat_threshold_display = str(repeat_threshold_set)

    gc_bins_str = str(gc_lower) + "-" + str(gc_upper)

    with open(summary_path, "w") as f:
        f.write("PID\tsample\ttotal_reads\tread_length\trepeat_threshold_set\trepeat_threshold_used\t"
                "tel_reads\tintratel_reads\tgc_bins_for_correction\ttotal_reads_with_tel_gc\t"
                "tel_content\tnum_cells\n")
        f.write("%s\t%s\t%i\t%s\t%s\t%s\t%i\t%i\t%s\t%i\t%f\t%i\n" % (
            pid, "single_cell", total_all, read_length, repeat_threshold_display,
            repeat_threshold_calc, tel_all, intratel_all, gc_bins_str,
            bulk_gc_in_window, bulk_tel_content, len(per_cell_results)))

    return per_cell_results


def tvr_screen_sc(input_dir, pid, barcode_tag='CB', min_base_quality=20):
    """
    Screen intratelomeric reads for TVR patterns, tracking per-cell counts.
    
    Returns per-cell TVR count table.
    """
    import re

    bam_path = os.path.join(input_dir, pid + "_filtered_intratelomeric.bam")
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    per_cell_tvr = defaultdict(lambda: defaultdict(int))  # {barcode: {pattern: count}}

    pattern_str = "GGG"
    offset = -3

    for read in bamfile.fetch(until_eof=True):
        try:
            barcode = read.get_tag(barcode_tag)
        except KeyError:
            continue

        seq = read.seq
        indices_fwd = [m.start() for m in re.finditer(pattern_str, seq)]

        # Check both orientations, use whichever has more GGG hits
        seq_rc = _getReverseComplement(seq)
        indices_rc = [m.start() for m in re.finditer(pattern_str, seq_rc)]

        if len(indices_rc) > len(indices_fwd):
            seq = seq_rc
            qual = read.qual[::-1]
            indices = indices_rc
        else:
            qual = read.qual
            indices = indices_fwd

        for i in indices:
            if i + offset < 0:
                continue

            p = seq[i + offset:i]
            if 'N' in p:
                continue

            # Check base quality at all 6 positions
            quals_ok = True
            for pos in range(6):
                if ord(qual[i + offset + pos]) - 33 < min_base_quality:
                    quals_ok = False
                    break

            if quals_ok:
                per_cell_tvr[barcode][p + "GGG"] += 1

    bamfile.close()

    return per_cell_tvr


def write_tvr_table_sc(per_cell_tvr, out_dir, pid, min_reads=500, per_cell_total_reads=None):
    """Write per-cell TVR counts to a TSV file."""

    # Collect all observed patterns
    all_patterns = set()
    for cell_patterns in per_cell_tvr.values():
        all_patterns.update(cell_patterns.keys())
    all_patterns = sorted(all_patterns)

    output_path = os.path.join(out_dir, pid + "_per_cell_TVR_counts.tsv")
    with open(output_path, "w") as f:
        header = ["cell_barcode"] + all_patterns
        f.write("\t".join(header) + "\n")

        for barcode in sorted(per_cell_tvr.keys()):
            # Optionally filter by min reads
            if per_cell_total_reads and per_cell_total_reads.get(barcode, 0) < min_reads:
                continue
            row = [barcode]
            for pattern in all_patterns:
                row.append(str(per_cell_tvr[barcode].get(pattern, 0)))
            f.write("\t".join(row) + "\n")

    sys.stdout.write("Wrote per-cell TVR counts for %d cells to %s\n" % (len(per_cell_tvr), output_path))


def _getReverseComplement(sequence):
    """Get the reverse complement of a DNA sequence."""
    trans = str.maketrans("ACGT", "TGCA")
    return sequence.translate(trans)[::-1]
