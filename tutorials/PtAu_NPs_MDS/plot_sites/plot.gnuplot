set term pngcairo size 1200,480

set output "mds_sites.png"

set multiplot layout 1,3

set size ratio -1
set format x ""
set format y ""

unset tics

set key right box opaque

set title "Atomic sites according to cluster"
plot "../mds_sites.dat" u 1:2:3 lc var pt 7 not, \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2 pt 7 ps 3 lc "green" t "Medoids", \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2:(sprintf("%i", $3+1)) w labels not

set title "Sites according to composition and surface/bulk"
plot "../mds_sites.dat" u 1:(stringcolumn(7) eq "Pt" && stringcolumn(5) eq "True" ? $2 : 1/0) pt 7 t "Pt (surface)", \
     "../mds_sites.dat" u 1:(stringcolumn(7) eq "Pt" && stringcolumn(5) eq "False" ? $2 : 1/0) pt 7 t "Pt (bulk)", \
     "../mds_sites.dat" u 1:(stringcolumn(7) eq "Au" && stringcolumn(5) eq "True" ? $2 : 1/0) pt 7 t "Au (surface)", \
     "../mds_sites.dat" u 1:(stringcolumn(7) eq "Au" && stringcolumn(5) eq "False" ? $2 : 1/0) pt 7 t "Au (bulk)", \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2 pt 6 ps 1.4 lw 4 lc "green" t "Medoids"

set tics

set title "Local site energy"
set cblabel "Local GAP energy (eV/atom)"
plot "../mds_sites.dat" u 1:2:6 palette pt 7 not, \
     "" u (stringcolumn(4) eq "True" ? $1 : 1/0):2 pt 6 ps 1.4 lw 4 lc "green" t "Medoids"
