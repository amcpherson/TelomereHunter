#!/bin/bash

#runs all r scripts for making plots


if [ "$SIMPLE_PLOT" == 'True' ];then

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR $REPEAT_THRESHOLD $CONSECUTIVE $MAPQ $BANDING_FILE < $PIPELINE_DIR/plot_spectrum_simple.R

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR $REPEAT_THRESHOLD $CONSECUTIVE $MAPQ $PID_INFO < $PIPELINE_DIR/plot_spectrum_summary_simple.R

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR $REPEAT_THRESHOLD $CONSECUTIVE $MAPQ < $PIPELINE_DIR/plot_unmapped_summary_simple.R

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR < $PIPELINE_DIR/plot_gc_content_simple.R

else

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR $REPEAT_THRESHOLD $CONSECUTIVE $MAPQ $BANDING_FILE < $PIPELINE_DIR/plot_spectrum.R

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR $REPEAT_THRESHOLD $CONSECUTIVE $MAPQ $PID_INFO < $PIPELINE_DIR/plot_spectrum_summary.R

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR $REPEAT_THRESHOLD $CONSECUTIVE $MAPQ < $PIPELINE_DIR/plot_unmapped_summary.R

	R --no-save --slave --args $PIPELINE_DIR $PID $OUTDIR < $PIPELINE_DIR/plot_gc_content.R

fi


R --no-save --slave --args $PID $OUTDIR $REPEAT_THRESHOLD $CONSECUTIVE $MAPQ $T_TYPE $C_TYPE $G_TYPE $J_TYPE < $PIPELINE_DIR/plot_repeat_frequency_intratelomeric.R
