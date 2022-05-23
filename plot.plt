set terminal gif animate delay 2
set output 'Gifs/Cohete_Luna3.gif'

posC = 'Data/Pos.dat'
posL = 'Data/PosLuna.dat'

set xrange [-1.5:1.5]
set yrange [-1.5:1.5]


# set object circle at first 0,0 radius char 1 \
#     fillcolor rgb 'blue' fillstyle solid noborder

do for [bb = 0:200]{
	plot posC every ::::1000*bb u 2:3 notitle ,\
		posL every ::::1000*bb notitle
}
# plot posC u 2:3, \
# posL 