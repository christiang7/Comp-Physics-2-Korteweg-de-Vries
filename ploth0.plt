reset
set grid
set yrange [0:20]
set ylabel 'A'
set term png
###
set xrange [0:1.2]
set xlabel 'h_0'
set output sprintf('plotting-A-h0-.png')
plot "h0-daten.txt" using ($1):($4) title 'Amplitude A(h_0)' lt rgb 'blue'
A_1(x) = a*x + b
A_2(x) = c/(x+e)+ d
fit [0:0.9] [0:20] A_1(x) 'h0-daten.txt' using ($1):($5) via a,b
fit [0.4:1] [0:20] A_2(x) 'h0-daten.txt' using ($1):($5) via c,d,e#
set output sprintf('plotting-A-h0-fit-1.png')
plot "h0-daten.txt" using ($1):($5) title 'Amplitude A(h_0)' lt rgb 'blue', A_1(x) lt rgb 'red'#, A_2(x) lt rgb 'brown'
set output sprintf('plotting-A-h0-fit-2.png')
plot [0.3:1.2][0:20] "h0-daten.txt" using ($1):($5) title 'Amplitude A(h_0)' lt rgb 'blue', A_2(x) lt rgb 'brown'
###
#set xrange [-0.2:1]
#set xlabel '1-h_0'
#set output sprintf('plotting-A-(1-h0).png')
#plot "h0-daten.txt" using (1-$1):($4)
###
set term qt
replot
