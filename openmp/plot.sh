#!/bin/bash

gnuplot << EOF

set term pdfcairo font "cmr10,14"
set output 'plot.pdf'

set xlabel 'nombre de coeurs'
set ylabel 'erreur'

set key outside

set grid


plot "resultat.dat" with lines
     
EOF
