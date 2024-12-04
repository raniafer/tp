#!/bin/bash

echo "compiling main program"
gfortran diffusion2.f90 -o diffusion2 -fopenmp
echo "done"

echo "running program"
export OMP_NUM_THREADS=14

echo "running with $OMP_NUM_THREADS threads"
./diffusion2

gnuplot << EOF

set term pdfcairo font "cmr10,14"
set output 'result.pdf'

set key outside

set grid


plot "initial.txt" using 1:2 with lines title "initiale",\
     "solution.txt" using 1:2 with lines title "finale"
     
EOF

gnuplot << EOF

set term pdfcairo font "cmr10,14"
set output 'performance.pdf'

set key outside

set xlabel "number of threads"
set ylabel "run time"

plot "time.txt" using 1:2 with lines

EOF
