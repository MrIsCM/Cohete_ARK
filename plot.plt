set terminal gif animate delay 3
set output 'Gifs/Cohete_Luna7.gif'

posC = 'Data/Pos.dat'
posL = 'Data/PosLuna.dat'

set xrange [-1.5:3]
set yrange [-2:1.5]


# set object circle at first 0,0 radius char 0.01 \
#     fillcolor rgb 'blue' fillstyle solid noborder

do for [bb = 0:250]{
	set title sprintf('Frame:%03.0f',bb)
	set object circle at first 0,0 radius char 1 \
    fillcolor rgb 'blue' fillstyle solid noborder
	# set label at 0,0 'Tierra'
	plot posC every ::::bb u 2:3 title 'Cohete' ,\
		posL every ::::bb u 2:3 title 'Luna'
}
# plot posC u 2:3, \
# posL 