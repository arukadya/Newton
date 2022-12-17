set terminal jpeg
set output "arrow2_2.jpg"
plot "arrow2_2.dat" using 1:2:3:4 w vector
