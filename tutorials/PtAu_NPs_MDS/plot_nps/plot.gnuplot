set term pngcairo size 1200,480

set output "mds_np.png"

set multiplot layout 1,3

set size ratio -1
set format x ""
set format y ""

unset xtics
unset ytics

set key right box

set title "PtAu NP according to cluster"
plot "../mds.dat" u 1:2:3 lc var pt 7 not, \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2 pt 7 ps 3 lc "green" t "Medoids", \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2:(sprintf("%i", $3+1)) w labels not

set title "PtAu NP according to composition"
set cblabel "x in Pt_{x}Au_{1-x}"
plot "../mds.dat" u 1:2:7 palette pt 7 not, \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2 pt 6 ps 1.4 lw 4 lc "green" t "Medoids"

set title "PtAu NP according to energy"
set cblabel "Total energy (eV/atom)"
plot "../mds.dat" u 1:2:($6/$5) palette pt 7 not, \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2 pt 6 ps 1.4 lw 4 lc "green" t "Medoids"
