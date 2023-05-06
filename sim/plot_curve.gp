set term png enhanced
set output 'sim_curve.png'
set xrange [4:11]
set logscale y
set xlabel 'Eb/N0 (dB)'
set ylabel 'BER'
set grid 
plot 'uncoded_res.csv' u 2:3 w lp t 'uncoded', 'cc_res.csv' u 2:3 with lp t 'CC', 0.5*erfc(sqrt(10**(0.1*x)))
