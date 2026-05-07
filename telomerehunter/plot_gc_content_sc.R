# Usage: R --no-save --slave --args <OUT_DIR> <PID> <PLOT_FILE_FORMAT> <GC_LOWER_LIMIT> <GC_UPPER_LIMIT> < plot_gc_content_sc.R
# Description: Makes a per-cell GC content plot showing GC bias for all reads and intratelomeric reads

library(ggplot2, quietly=TRUE, warn.conflicts=FALSE)
library(cowplot, quietly=TRUE, warn.conflicts=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE)

# get commandline arguments
commandArgs = commandArgs()
out_dir = commandArgs[5]
pid = commandArgs[6]
plot_file_format = commandArgs[7]
gc_lower_limit = as.numeric(commandArgs[8])
gc_upper_limit = as.numeric(commandArgs[9])

if (plot_file_format == "all") {
  plot_file_format = c("pdf", "png", "svg")
}

pid_dir = file.path(out_dir, pid)
plot_dir = file.path(pid_dir, "plots")
if (!file.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }


##################################
### Read per-cell GC data      ###
##################################

gc_all_file = file.path(pid_dir, paste0(pid, "_per_cell_gc_content.tsv"))
gc_intratel_file = file.path(pid_dir, paste0(pid, "_per_cell_intratel_gc_content.tsv"))

gc_all = read.table(gc_all_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gc_intratel = read.table(gc_intratel_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Compute per-cell fraction of reads at each GC bin
gc_all_totals = aggregate(read_count ~ cell_barcode, data = gc_all, FUN = sum)
names(gc_all_totals)[2] = "total"
gc_all = merge(gc_all, gc_all_totals, by = "cell_barcode")
gc_all$fraction_of_reads = gc_all$read_count / gc_all$total

gc_intratel_totals = aggregate(read_count ~ cell_barcode, data = gc_intratel, FUN = sum)
names(gc_intratel_totals)[2] = "total"
gc_intratel = merge(gc_intratel, gc_intratel_totals, by = "cell_barcode")
gc_intratel$fraction_of_reads = gc_intratel$read_count / gc_intratel$total

# Label read types
gc_all$read_type = "All reads"
gc_intratel$read_type = "Intratelomeric reads"

dfm = rbind(gc_all[, c("cell_barcode", "gc_content_percent", "fraction_of_reads", "read_type")],
            gc_intratel[, c("cell_barcode", "gc_content_percent", "fraction_of_reads", "read_type")])

gc_bins = c(gc_lower_limit:gc_upper_limit)


#################
### make plot ###
#################

cells = sort(unique(gc_all$cell_barcode))
n_cells = length(cells)

# For PDF: one page per cell using pdf() device
pdf_path = file.path(plot_dir, paste0(pid, "_gc_content_per_cell.pdf"))
pdf(pdf_path, width = 10, height = 5)

for (cell in cells) {

  df_cell_all = gc_all[gc_all$cell_barcode == cell, ]
  df_cell_intratel = gc_intratel[gc_intratel$cell_barcode == cell, ]

  df_cell_all$read_type = "All reads"
  df_cell_intratel$read_type = "Intratelomeric reads"

  df_cell = rbind(df_cell_all[, c("gc_content_percent", "fraction_of_reads", "read_type")],
                  df_cell_intratel[, c("gc_content_percent", "fraction_of_reads", "read_type")])

  ymax_cell = max(df_cell_all$fraction_of_reads, na.rm = TRUE)

  p = ggplot(df_cell, aes(x = gc_content_percent, y = fraction_of_reads)) +
    geom_rect(data = data.frame(read_type = factor(c("All reads"))),
              aes(xmin = min(gc_bins) - 0.5, xmax = max(gc_bins) + 0.5,
                  ymin = 0, ymax = ymax_cell),
              inherit.aes = FALSE, alpha = 0.3, fill = "grey") +
    geom_line(color = "blue") +
    theme_cowplot() +
    ggtitle(paste0(pid, ": ", cell)) +
    xlab("GC content [%]") +
    ylab("Fraction of reads") +
    geom_text(data = data.frame(read_type = factor(c("All reads"))),
              aes(x = min(gc_upper_limit, 88), y = ymax_cell,
                  label = "Bins used for\nGC correction"),
              inherit.aes = FALSE, color = "grey40", size = 3.5, hjust = 0, vjust = 2) +
    xlim(0, 100) +
    facet_wrap(~read_type, scales = "free")

  print(p)
}

dev.off()

# Also save individual png/svg if requested
other_formats = plot_file_format[plot_file_format != "pdf"]
if (length(other_formats) > 0) {
  for (cell in cells) {

    df_cell_all = gc_all[gc_all$cell_barcode == cell, ]
    df_cell_intratel = gc_intratel[gc_intratel$cell_barcode == cell, ]

    df_cell_all$read_type = "All reads"
    df_cell_intratel$read_type = "Intratelomeric reads"

    df_cell = rbind(df_cell_all[, c("gc_content_percent", "fraction_of_reads", "read_type")],
                    df_cell_intratel[, c("gc_content_percent", "fraction_of_reads", "read_type")])

    ymax_cell = max(df_cell_all$fraction_of_reads, na.rm = TRUE)

    p = ggplot(df_cell, aes(x = gc_content_percent, y = fraction_of_reads)) +
      geom_rect(data = data.frame(read_type = factor(c("All reads"))),
                aes(xmin = min(gc_bins) - 0.5, xmax = max(gc_bins) + 0.5,
                    ymin = 0, ymax = ymax_cell),
                inherit.aes = FALSE, alpha = 0.3, fill = "grey") +
      geom_line(color = "blue") +
      theme_cowplot() +
      ggtitle(paste0(pid, ": ", cell)) +
      xlab("GC content [%]") +
      ylab("Fraction of reads") +
      geom_text(data = data.frame(read_type = factor(c("All reads"))),
                aes(x = min(gc_upper_limit, 88), y = ymax_cell,
                    label = "Bins used for\nGC correction"),
                inherit.aes = FALSE, color = "grey40", size = 3.5, hjust = 0, vjust = 2) +
      xlim(0, 100) +
      facet_wrap(~read_type, scales = "free")

    # Sanitize barcode for filename
    cell_safe = gsub("[^A-Za-z0-9_-]", "_", cell)

    for (fmt in other_formats) {
      ggsave(file.path(plot_dir, paste0(pid, "_gc_content_", cell_safe, ".", fmt)), p, width = 10, height = 5)
    }
  }
}

# Suppress Rplots.pdf
if (file.exists("Rplots.pdf")) { file.remove("Rplots.pdf") }
