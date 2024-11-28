#!/bin/bash

gnuplot << EOF

set term pdfcairo font "cmr10,14"
set output 'plot.pdf'

set key outside

set grid


plot "result.txt" using 1:2 with lines title "initiale",\
     "result.txt" using 1:2 with lines title "finale"
     
EOF
