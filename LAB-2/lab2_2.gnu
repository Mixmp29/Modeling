set terminal pngcairo size 1000, 600
set output "res2.png"

set xlabel "Случайная величина Xi "
set ylabel "Вероятность Pi"
set title "Соотноешение теоретических и эмпирических вероятностей с возвратом"
set grid
#set key left
set zeroaxis

set xrange [-0.5:7.5]
set yrange [0:40000]

plot 'Results.txt' usi ($1-0.1):4:(0.2) w boxes fs pattern 2 title "эмпирическая p(X)",\
     'Results.txt' usi ($1+0.1):5:(0.2) w boxes fs solid 0.7 title "теоретическая p^'(X)"