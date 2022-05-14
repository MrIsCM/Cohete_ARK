# set terminal pngcairo
# set output 'prueba.png'

posC = 'Data/Pos.dat'
posL = 'Data/PosLuna.dat'



set object circle at first 0,0 radius char 1 \
    fillcolor rgb 'blue' fillstyle solid noborder

plot posC u 2:3 w lines,\
	posL