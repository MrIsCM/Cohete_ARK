set terminal gif animate delay 3
set output 'Plots/CL6.gif'

posC = 'Data2/Pos.dat'
posL = 'Data2/PosLuna.dat'

set xrange [-2:2.5]
set yrange [-2.5:2.5]


# set object circle at first 0,0 radius char 0.01 \
#     fillcolor rgb 'blue' fillstyle solid noborder

do for [bb = 0:250]{
	set title sprintf('Frame:%03.0f',bb)
	plot posC every ::::bb u 2:3 w lines title 'Cohete' ,\
		posL every ::::bb u 2:3 title 'Luna'
}