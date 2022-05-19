set terminal gif animate delay 2
set output 'Gifs/Cohete_Luna.gif'

posC = 'Data/Pos.dat'
posL = 'Data/PosLuna.dat'

set xrange [-1.5:1.5]
set yrange [-1.5:1.5]


# set object circle at first 0,0 radius char 1 \
#     fillcolor rgb 'blue' fillstyle solid noborder

do for [bb = 0:1000]{
	plot posC every ::::bb u 2:3 w lines title sprintf('Frame:%03.0f',bb),\
		posL every ::::bb notitle
}