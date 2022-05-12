!=======================================
! 		Ismael Charpentier
!
! 			12/05/2022
!
! 	  Algoritmo de Runge-Kutta
! 	Problema de los tres cuerpos
!
!=======================================
program cohete
	implicit none

	real :: h 		! Parametro de discretizacion 
	integer :: t 	! Parametro de iteración temporal

	integer, parameter :: t_max = 1000
	real, parameter :: m = 20000 		! Masa del cohete [Kg]
	
	! Datos del problema
	double precision :: Mt, Ml, m, G, dtl, w, Rt, Rl  

	!Factores de escala
	double precision :: fr, fpr, fpphi


	!========================
	!	Data: link = https://ergodic.ugr.es/cphys/index.php?id=lec_cohete
	!========================

	G = 6.67E-11 		! Constante de grav. univ. [N m**2 Kg**-2]
	Mt = 5.9736E24		! Masa Tierra [Kg]
	Ml = 0.07349E24		! Masa Luna [Kg]
	dtl = 3.844E8 		! Distancia Tierra-Luna (unidad astro.) [m] 
	w = 2.6617E-6 		! Vel. Angular Luna alrededor Tierra [s**-1]
	Rt = 6.378160E6		! Radio Tierra [m] 
	Rl = 1.7374E6		! Radio Luna [m] 

	!-------------------------
	! 	Factores de escala
	!-------------------------
	fr = 1/dtl 
	fpr = 1/(m*dtl)
	fpphi = 1/(m*dtl**2)



	! Archivo para guardar la posición del cohete
	open(10, file='Data/Pos.dat', status='Unknown')

	! Bucle temporal
	do t = 0, t_max 	

	
	end do

	close(10)
	
end program cohete