!=======================================
! 				ICM
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
	real*8 :: t, h

	integer :: i 

	!integer, parameter :: t_max = 60000
	real*8, parameter :: pi = 3.14159265359
	
	! Datos del problema
	real*8 :: Mt, Ml, m, G, dtl, w, Rt, Rl

	! Coordenadas
	real*8 :: r, phi, pr, pphi
	real*8 :: y(4)

	!Factores de escala
	real*8 :: fr, fpr, fpphi

	!Variables utilidad. Simplificar expresiones
	real*8 :: delta, mu, r_p


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
	h = 0.01
	! En segundo lugar, se utiliza el reajuste de h descrito en el guión



	! Asignacion de condiciones iniciales
	r = Rt*fr 					! Despega de la Tierra 
	phi = 0.3 				! *Parece* arbitrario
	pr = 0.1*fpr  			! ~ Vel. escape [m/s] 
	pphi = 0.001*fpphi

	y = (/r, phi, pr, pphi/)

	! Constantes de utilidad (simplificacion de expresiones)
	delta = G*Mt/dtl**3
	mu = Ml/Mt

	! Archivo para guardar la posición del cohete
	open(12, file='Data/Pos.dat', status='Unknown')
	open(13, file='Data/PosLuna.dat', status='Unknown')
	open(14, file='Data/DistL_C.dat', status='Unknown')
	open(15, file='Data/PosPolar1.dat', status='Unknown')

	! Bucle temporal
	!run = .true.
	do i = 0, 100000
		t = i*h 
		r_p = (1+y(1)**2-2*y(1)*cos(phi-w*t)) 

		if (mod(i,100) == 0) then
			! Cohete: tiempo, x, y 
			write(12,*) t, y(1)*cos(y(2)), y(1)*sin(y(2))
			! Luna: tiempo, x_luna, y_luna
			write(13,*) t , dtl*fr*cos(w*t), dtl*fr*sin(w*t)
			! Distancia Luna-Cohete: tiempo, distancia
			write(14,*) t, r_p
			! Cohete Polares: tiempo, radio, angulo
			write(15,*) t, y(1), y(2)
		end if
		! Calculo distancia Luna-Nave (reescalado)
		call alg_RK(delta, mu, r_p, w, h, y, t)
			
	end do 

	close(12)
	close(13)
	close(14)
	close(15)
	
end program cohete

subroutine alg_RK(delta, mu, r_p, w, h, y, t)
	implicit none
	real*8, intent(in) :: h

	! Constantes presentes en las ecs. mov
	real*8, intent(in) :: delta, mu, r_p, w	

	! Coordenadas
	real*8, intent(in) :: t 		! t
	real*8, intent(inout) :: y(4) 	! y_n(t) 	[n=1,2,3,4]

	integer :: i, j
	real*8 :: K(4,4), f_ynt

	do i = 1, 4
		call fn(y(1), y(2), y(3), y(4), t, delta, mu, r_p, w, i, f_ynt)
		K(i,1) = h*f_ynt
	end do
		
	do i=1,4
		call fn(y(1) + K(1,1)/2, y(2)+ K(2,1)/2, y(3)+ K(3,1)/2, y(4)+ K(4,1)/2, t+h/2, delta, mu, r_p, w, i, f_ynt)
		K(i,2) = h*f_ynt
	end do 

	do i=1,4
		call fn(y(1) + K(1,2)/2, y(2)+ K(2,2)/2, y(3)+ K(3,2)/2, y(4)+ K(4,2)/2, t+h/2, delta, mu, r_p, w, i, f_ynt)
		K(i,3) = h*f_ynt
	end do 

	do i=1,4
		call fn(y(1) + K(1,3), y(2)+ K(2,3), y(3)+ K(3,3), y(4)+ K(4,3), t+h/2, delta, mu, r_p, w, i, f_ynt)
		K(i,4) = h*f_ynt
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
	real*8, intent(in) :: pr 
	real*8, intent(out) :: result

	result = pr 
	
end subroutine dr_dt

!----------------------
! 	(2) phi/dt = ...
!----------------------
subroutine dphi_dt(pphi, result)
	implicit none
	real*8, intent(in) :: pphi
	real*8, intent(out) :: result

	result = pphi 
	
end subroutine dphi_dt


!----------------------
! 	(3) dp_r/dt = ...
!----------------------
subroutine dpr_dt(r, phi, p_phi, t, delta, mu, r_p, w, result)
	implicit none
	real*8, intent(in) :: r, phi, p_phi, t, delta, mu, r_p, w
	real*8, intent(out) :: result

	result = p_phi**2/r**3 - delta*(1/r**2 + mu/r_p**3 *(r - cos(phi-w*t)))
	
end subroutine dpr_dt


!----------------------
! 	(4) dp_phi/dt = ...
!----------------------
subroutine dpo_dt(r, phi, t, delta, mu, r_p, w, result)
	implicit none
	real*8, intent(in) :: delta, mu, r_p, w, r, phi, t
	real*8, intent(out) :: result

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
	real*8, intent(in) :: r, phi, pr, pphi, t, delta, mu, r_p, w

	real*8, intent(inout) :: f_yt

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
! 	real*8, intent(in) :: x, y
! 	real*8, intent(out) :: result

! end subroutine f