set terminal pngcairo size 1000, 600
set output "res.png"

set xlabel "Значения  a + γ1(b − a)  "
set ylabel "Значение γ2"
set title "10000 случайных чисел от 0 до 2 сгенерированных по заданному распределению"
set grid
set key left
set zeroaxis


set xrange [0:2]
set yrange [0:2]

plot "DistFunc.txt" u 1:2 w p title 'случайная величина',\
     [0:2] 2-x w line title 'y = 2 - x',\
     [0:2] sin(x) w line title 'y = sin(x)'
