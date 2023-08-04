set terminal pngcairo size 1000, 600
set output "Hi2.png"

set xlabel "Случайные числа"
set ylabel "Значение ni/n"
set title "10000 случайных чисел от 0 до 1"
set grid
set key left
set zeroaxis
set xtics 0.1
set ytics 0.01

set xrange [0:1]
set yrange [0:0.1]

plot "Hi2.txt" u 1:2 w lp title 'k = 14'
