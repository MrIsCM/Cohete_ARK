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

	! Parametro de iteración temporal
	! Parametro de discretizacion
	real :: t, h 		 

	integer, parameter :: t_max = 6000
	
	! Datos del problema
	real :: Mt, Ml, m, G, dtl, w, Rt, Rl

	! Coordenadas
	real :: r, phi, pr, pphi
	real :: y(4)

	!Factores de escala
	real :: fr, fpr, fpphi

	!Variables utilidad. Simplificar expresiones
	real :: delta, mu, r_p

	! Variable de control del bucle temporal
	logical :: run 


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



	! Asignacion de condiciones iniciales
	r = Rt*fr 			! Despega de la Tierra 
	phi = 0 		! *Parece* arbitrario
	pr = 11*fr 		! ~ Vel. escape [m/s] 
	pphi = 2

	y = (/r, phi, pr, pphi/)

	! Constantes de utilidad (simplificacion de expresiones)
	delta = G*Mt/dtl**3
	mu = Ml/Mt

	! Archivo para guardar la posición del cohete
	open(12, file='Data/Pos.dat', status='Unknown')
	open(13, file='Data/PosLuna.dat', status='Unknown')


	! Bucle temporal
	t=0
	run = .true.
	do while (t<=t_max)
		! Calculo distancia Luna-Nave (reescalado)
		r_p = (1+r**2-2*r*cos(phi-w*t)) 
		if (r_p <= Rl) then
			run = .False.
		else 
		end if 
			call alg_RK(delta, mu, r_p, w, h, y, t)
			write(12,*) t, (y(1))*cos(phi), (y(1))*sin(phi)
			write(13,*) dtl*fr*cos(w*t), dtl*fr*sin(w*t)
			t = t+h 
	end do 

	close(13)
	
end program cohete

subroutine alg_RK(delta, mu, r_p, w, h, y, t)
	implicit none
	real, intent(in) :: h

	! Constantes presentes en las ecs. mov
	real, intent(in) :: delta, mu, r_p, w	

	! Coordenadas
	real, intent(inout) :: t 		! t
	real, intent(inout) :: y(4) 	! y_n(t) 	[n=1,2,3,4]

	integer :: i, j
	real :: K(4,4), f_ynt

	do i = 1, 4
		call fn(y(1), y(2), y(3), y(4), t, delta, mu, r_p, w, i, f_ynt)
		K(i,1) = h*f_ynt

		call fn(y(1) + K(1,1)/2, y(2)+ K(2,1)/2, y(3)+ K(3,1)/2, y(4)+ K(4,1)/2, t+h/2, delta, mu, r_p, w, i, f_ynt)
		K(i,1) = h*f_ynt

		call fn(y(1) + K(1,2)/2, y(2)+ K(2,2)/2, y(3)+ K(3,2)/2, y(4)+ K(4,2)/2, t+h/2, delta, mu, r_p, w, i, f_ynt)
		K(i,1) = h*f_ynt

		call fn(y(1) + K(1,3), y(2)+ K(2,3), y(3)+ K(3,3), y(4)+ K(4,3), t+h/2, delta, mu, r_p, w, i, f_ynt)
		K(i,1) = h*f_ynt
	end do

	do i = 1, 4
		y(i) = y(i) + (K(i,1) + 2*K(i,2) + 2*K(i,3) + K(i,4))/6
	end do

	! open(10, file='Data/k.dat', status='unknown')
	! open(11, file='Data/y.dat', status='unknown')
	! do i = 1, 4
	! 	write(10,*) (K(i, j), j=1,4)
	! 	write(11,*) t, y(i)
	! end do
	! close(10)
	! close(11)
	
end subroutine alg_RK


!----------------------
! 	(1) r/dt = ...
!----------------------
subroutine dr_dt(pr, result)
	implicit none
	real, intent(in) :: pr 
	real, intent(out) :: result

	result = pr 
	
end subroutine dr_dt

!----------------------
! 	(2) phi/dt = ...
!----------------------
subroutine dphi_dt(pphi, result)
	implicit none
	real, intent(in) :: pphi
	real, intent(out) :: result

	result = pphi 
	
end subroutine dphi_dt


!----------------------
! 	(3) dp_r/dt = ...
!----------------------
subroutine dpr_dt(r, phi, p_phi, t, delta, mu, r_p, w, result)
	implicit none
	real, intent(in) :: r, phi, p_phi, t, delta, mu, r_p, w
	real, intent(out) :: result

	result = p_phi**2/r**3 - delta*(1/r**2 + mu/r_p**3 *(r - cos(phi-w*t)))
	
end subroutine dpr_dt


!----------------------
! 	(4) dp_phi/dt = ...
!----------------------
subroutine dpo_dt(r, phi, t, delta, mu, r_p, w, result)
	implicit none
	real, intent(in) :: delta, mu, r_p, w, r, phi, t
	real, intent(out) :: result

	result = - delta*mu*r*sin(phi-w*t)/r_p**3 
	
end subroutine dpo_dt

!=============================================================
!	Subrutina que unifica todas mis funciones en una
!	
!	La idea es que me proporcione claridad al resto 
!	del codigo. Concretamente al llamarlo en la sub
!	del algortimo Runge-Kutta
!=============================================================
subroutine fn(r, phi, pr, pphi, t, delta, mu, r_p, w, i, f_yt)
	implicit none 
	integer, intent(in) :: i 
	real, intent(in) :: r, phi, pr, pphi, t, delta, mu, r_p, w

	real, intent(out) :: f_yt

	if (i==1) then
		call dr_dt(pr, f_yt)
	else if (i==2) then
		call dphi_dt(pphi, f_yt)
	else if (i==3) then 
		call dpr_dt(r, phi, pphi, t, delta, mu, r_p, w, f_yt)
	else if (i==4) then
		call dpo_dt(r, phi, t, delta, mu, r_p, w, f_yt)
	end if

	
end subroutine fn

!---------------------------------------------------------------
!	Función genérica que permite ampliar el algoritmo
!	 para EDOs de cualquier tipo, definiendo f en este apartado
!---------------------------------------------------------------
! subroutine f(x, y, result)
!	implicit none
! 	real, intent(in) :: x, y
! 	real, intent(out) :: result

! end subroutine f