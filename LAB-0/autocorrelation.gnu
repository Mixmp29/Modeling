set terminal pngcairo size 1000, 600
set output "Autocorrelation.png"

set xlabel "Временной сдвиг τ"
set ylabel "Значение a(τ)"
set title "Автокорреляция 10000 случайных чисел от 0 до 1"
set grid
set zeroaxis
set xtics 1
set ytics 0.1

set xrange [0:14]
set yrange [-1:1]

plot "autocorrelation.txt" u 1:2 w lp title "автокорреляция"
