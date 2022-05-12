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

	m = 10000 			! Masa del cohete [Kg]

	!-------------------------
	! 	Factores de escala
	!-------------------------
	fr = 1/dtl 
	fpr = 1/(m*dtl)
	fpphi = 1/(m*dtl**2)

	!----------------------
	!	Parametro h: 
	!----------------------
	! En primer lugar lo escojo "arbitrariamente" 
	h = 60
	! En segundo lugar, se utiliza el reajuste de h descrito en el guión



	! Archivo para guardar la posición del cohete
	open(10, file='Data/Pos.dat', status='Unknown')

	! Bucle temporal
	do t = 0, t_max 	

	end do

	close(10)
	
end program cohete

subroutine alg_RK(G, Mt, Ml, dtl, w, Rt, Rl, m, x_pos, y_pos)
	double precision, intent(in) :: G, Mt, Ml, dtl, w, Rt, Rl
	double precision, intent(inout) ::  x_pos, y_pos


	
end subroutine alg_RK

!----------------------
! 	(1) r/dt = ...
!----------------------
subroutine r_dt(pr, result)
	real, intent(in) :: pr 
	real, intent(out) :: result

	result = pr 
	
end subroutine r_dt

!----------------------
! 	(2) phi/dt = ...
!----------------------
subroutine phi_dt(pphi, result)
	real, intent(in) :: pphi
	real, intent(out) :: result

	result = pphi 
	
end subroutine phi_dt


!----------------------
! 	(3) dp_r/dt = ...
!----------------------
subroutine dp_dt(p_phi, delta, mu, r, r_p, phi, w, t, result)
	real, intent(in) :: p_phi, delta, mu, r, r_p, phi, w, t
	real, intent(out) :: result

	result = p_phi**2/r**3 - delta*(1/r**2 + mu/r_p**3 *(r - cos(phi-w*t)))
	
end subroutine dp_dt


!----------------------
! 	(4) dp_phi/dt = ...
!----------------------
subroutine dpo_dt(delta, mu, r, r_p, phi, w, result)
	real, intent(in) :: delta, mu, r, r_p, phi, w
	real, intent(out) :: result

	result = - delta*mu*r*sin(phi-w*t)/r_p**3 
	
end subroutine dpo_dt


!---------------------------------------------------------------
!	Función genérica que permite ampliar el algoritmo
!	 para EDOs de cualquier tipo, definiendo f en este apartado
!---------------------------------------------------------------
! subroutine f(x, y, result)
! 	real, intent(in) :: x, y
! 	real, intent(out) :: result

! end subroutine f