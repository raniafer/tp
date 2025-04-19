#!/bin/bash

gnuplot << EOF
set terminal pdfcairo font "cmr10,14"

set size 2,1

set output "resultat.pdf"

set xrange [0. : 1.]
set grid
set xlabel 'x'
set ylabel 'h(x)'

set multiplot layout 1,2

set title "sequentiel"
	plot 'finalf.dat' using 1:2 title "" with lines
set title "parallel"
	plot 'finalfMPI.dat' using 1:2 title "" with lines

EOF
