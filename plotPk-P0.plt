reset
set grid
#set yrange [0:1000]
set yrange [0.000000000001:100000]
#set xrange [5000:100000]
set ylabel '|P_k(t)-P_k(0)|'
set logscale y
set term png
###
#set xrange [0:1.2]
set xlabel 'n'
set output sprintf('plotting-Pk-P0-test.png')
plot "P_K-data.txt" using ($1):($2) title '|P_0(t)-P_0(0)|' lt rgb 'blue', "P_K-data.txt" using ($1):($3) title '|P_1(t)-P_1(0)|' lt rgb 'red', "P_K-data.txt" using ($1):($4) title '|P_2(t)-P_2(0)|' lt rgb 'green'#
#plot "Korrekturen/P_K-data-N2.txt" using ($1):($2) title '|P_0(t)-P_0(0)|' lt rgb 'blue', "Korrekturen/P_K-data-N2.txt" using ($1):($3) title '|P_1(t)-P_1(0)|' lt rgb 'red', "Korrekturen/P_K-data-N2.txt" using ($1):($4) title '|P_2(t)-P_2(0)|' lt rgb 'green'
set term qt
replot
