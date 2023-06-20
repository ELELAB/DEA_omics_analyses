#!/bin/bash
NTHREADS=4
export OMP_NUM_THREADS=$NTHREADS
tsp -L yi -N $NTHREADS Rscript omics_wrangle.R
