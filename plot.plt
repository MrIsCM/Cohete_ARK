set terminal gif animate delay 3
set output 'Plots/CL3.gif'

posC = 'Data/Pos.dat'
posL = 'Data/PosLuna.dat'

set xrange [-1:1.5]
set yrange [-1.5:1.5]


# set object circle at first 0,0 radius char 0.01 \
#     fillcolor rgb 'blue' fillstyle solid noborder

do for [bb = 0:400]{
	set title sprintf('Frame:%03.0f',bb)
	plot posC every ::::bb u 2:3 title 'Cohete' ,\
		posL every ::::bb+40 u 2:3 title 'Luna'
}