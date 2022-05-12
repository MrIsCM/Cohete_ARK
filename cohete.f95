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
	real :: Mt, Ml, m, G, dtl, w, Rt, Rl

	! Coordenadas
	real :: r, phi

	!Factores de escala
	real :: fr, fpr, fpphi

	!Variables utilidad. Simplificar expresiones
	real :: delta, mu, r_p


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

	! Asignacion de condiciones iniciales
	! r = 
	! phi = 
	! pr =
	! pphi = 


	! Bucle temporal
	do t = 0, t_max
		delta = G*Mt/dtl**3
		mu = Ml/Mt
		r_p = (1+r**2-2*r*cos(phi-w*t))



	end do

	close(10)
	
end program cohete

subroutine alg_RK_rdot(pr, h, result)
	implicit none
	real, intent(in) :: h
	real, intent(out) :: pr, result

	real, dimension(4) :: K
	real :: fy, y0

	call r_dt(pr, fy)
	y0 = fy
	K(1) = h*fy

	call r_dt(pr+K(1)/2, fy)
	K(2) = h*fy

	call r_dt(pr+K(2)/2, fy)
	K(3)=h*fy 

	call r_dt(pr+K(3), fy)
	K(4)=h*fy

	result = y0 + (K(1)+2*K(2)+2*K(3)+K(4))/6

end subroutine alg_RK_rdot

subroutine alg_RK_phidot(pphi, h, result)
	implicit none
	real, intent(in) :: h
	real, intent(out) :: pphi, result

	real, dimension(4) :: K
	real :: fy, y0

	call phi_dt(pphi, fy)
	y0 = fy
	K(1) = h*fy

	call phi_dt(pphi+K(1)/2, fy)
	K(2) = h*fy

	call phi_dt(pphi+K(2)/2, fy)
	K(3)=h*fy 

	call phi_dt(pphi+K(3), fy)
	K(4)=h*fy

	result = y0 + (K(1)+2*K(2)+2*K(3)+K(4))/6

end subroutine alg_RK_phidot

subroutine alg_RK_prdot(h, pphi, delta, mu, r, r_p, phi, w, t, result)
	implicit none
	integer, intent(in) :: t
	real, intent(in) :: h, pphi, delta, mu, r, r_p, phi, w
	real, intent(out) :: result

	real, dimension(4) :: K
	real :: fy, y0

	call dpr_dt(pphi, delta, mu, r, r_p, phi, w, t, fy)
	y0 = fy
	K(1) = h*fy

	call dpr_dt(pphi+K(1)/2, delta, mu, r, r_p, phi, w, t+h/2, fy)
	K(2) = h*fy

	call dpr_dt(pphi+K(2)/2, delta, mu, r, r_p, phi, w, t+h/2, fy)
	K(3)=h*fy 

	call dpr_dt(pphi+K(3), delta, mu, r, r_p, phi, w, t+h, fy)
	K(4)=h*fy

	result = y0 + (K(1)+2*K(2)+2*K(3)+K(4))/6

end subroutine alg_RK_prdot

subroutine alg_RK_pphidot(h, delta, mu, r, r_p, phi, w, t, result)
	implicit none
	integer, intent(in) :: t
	real, intent(in) :: h, delta, mu, r, r_p, phi, w
	real, intent(out) :: result

	real, dimension(4) :: K
	real :: fy, y0

	call dpo_dt(t, delta, mu, r, r_p, phi, w, fy)
	y0 = fy
	K(1) = h*fy

	call dpo_dt(pphi+K(1)/2, delta, mu, r, r_p, phi, w, t+h/2, fy)
	K(2) = h*fy

	call dpo_dt(pphi+K(2)/2, delta, mu, r, r_p, phi, w, t+h/2, fy)
	K(3)=h*fy 

	call dpo_dt(pphi+K(3), delta, mu, r, r_p, phi, w, t+h, fy)
	K(4)=h*fy

	result = y0 + (K(1)+2*K(2)+2*K(3)+K(4))/6

end subroutine alg_RK_pphidot

!----------------------
! 	(1) r/dt = ...
!----------------------
subroutine r_dt(pr, result)
	implicit none
	real, intent(in) :: pr 
	real, intent(out) :: result

	result = pr 
	
end subroutine r_dt

!----------------------
! 	(2) phi/dt = ...
!----------------------
subroutine phi_dt(pphi, result)
	implicit none
	real, intent(in) :: pphi
	real, intent(out) :: result

	result = pphi 
	
end subroutine phi_dt


!----------------------
! 	(3) dp_r/dt = ...
!----------------------
subroutine dpr_dt(p_phi, delta, mu, r, r_p, phi, w, t, result)
	implicit none
	real, intent(in) :: p_phi, delta, mu, r, r_p, phi, w, t
	real, intent(out) :: result

	result = p_phi**2/r**3 - delta*(1/r**2 + mu/r_p**3 *(r - cos(phi-w*t)))
	
end subroutine dpr_dt


!----------------------
! 	(4) dp_phi/dt = ...
!----------------------
subroutine dpo_dt(t, delta, mu, r, r_p, phi, w, result)
	implicit none
	integer, intent(in) :: t 
	real, intent(in) :: delta, mu, r, r_p, phi, w
	real, intent(out) :: result

	result = - delta*mu*r*sin(phi-w*t)/r_p**3 
	
end subroutine dpo_dt


!---------------------------------------------------------------
!	Función genérica que permite ampliar el algoritmo
!	 para EDOs de cualquier tipo, definiendo f en este apartado
!---------------------------------------------------------------
! subroutine f(x, y, result)
!	implicit none
! 	real, intent(in) :: x, y
! 	real, intent(out) :: result

! end subroutine f