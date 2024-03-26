set terminal wxt size 1024,512
set xrange[221:273]
set yrange[:]
set key font ",16"
set xtics font ", 16"
set ytics font ", 16"
#set multiplot layout 3,1 rowsfirst
plot "ecg205.txt" with steps title "ECG 205" lt rgb "black" lw 3
pause -1