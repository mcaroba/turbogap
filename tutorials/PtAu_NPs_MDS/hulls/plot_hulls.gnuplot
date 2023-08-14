set term pngcairo size 1200,480

set output "hulls.png"

set multiplot layout 1,6

set xlabel "x in Pt_{x}Au_{1-x}"
set ylabel "Energy above hull (eV/atom)"

set yrange [-0.02:0.12]
set xtics 0.5
set mxtics 1

do for [i=25:50:5] {
set title "N = ".i
plot "".i.".dat" index 1 pt 7 not, "" index 0 pt 6 ps 1.2 lw 4 w lp not
}
