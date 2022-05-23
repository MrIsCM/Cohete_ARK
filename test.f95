!=======================================
! 				ICM
!
! 	  Algoritmo de Runge-Kutta
! 	Simplificacion: Tierra y Cohete
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
	real*8 :: Mt, m, G, Rt

	! Coordenadas
	real*8 :: r_0, phi_0, pr_0, pphi_0
	real*8 :: y(4)

	! !Factores de escala
	! real*8 :: fr, fpr, fpphi



	!========================
	!	Data: link = https://ergodic.ugr.es/cphys/index.php?id=lec_cohete
	!========================

	G = 6.67E-11 		! Constante de grav. univ. [N m**2 Kg**-2]
	Mt = 5.9736E24		! Masa Tierra [Kg]
	Rt = 6.378160E6		! Radio Tierra [m] 

	m = 1000 			! Masa del cohete [Kg]

	! !-------------------------
	! ! 	Factores de escala
	! !-------------------------
	! fr = 1/dtl 
	! fpr = 1/(m*dtl)
	! fpphi = 1/(m*dtl**2)

	!----------------------
	!	Parametro h: 
	!----------------------
	! En primer lugar lo escojo "arbitrariamente" 
	h = 20
	! En segundo lugar, se utiliza el reajuste de h descrito en el guión



	! Asignacion de condiciones iniciales
	r_0 = Rt					! Despega de la Tierra 
	phi_0 = 0.3 				! *Parece* arbitrario
	pr_0 = m*Rt*0.0014 			! ~ Vel. escape [m/s] 
	pphi_0 = m*Rt**2 *0.001

	y = (/r_0, phi_0, pr_0, pphi_0/)

	open(10, file='DataTest/Pos.dat', status='Unknown')
	open(11, file='DataTest/PosPolar.dat', status='Unknown')
	open(12, file='DataTest/Momentos.dat', status='Unknown')



	do i=0, 4*10000
		t = h*i
		if (mod(i,100) == 0) then 
			write(10,*) t, y(1)*cos(y(2))/Rt, y(1)*sin(y(2))/Rt
			write(11,*) t, y(1)/Rt, y(2)
			write(12,*) t, y(3)/(m*Rt), y(4)/(m*Rt**2)
		end if 

		call alg_RK(h, m, y, t)
	end do

	close(10)
	close(11)
	close(12)

	
end program cohete

subroutine alg_RK(h, m, y, t)
	implicit none
	real*8, intent(in) :: h, m 

	! Coordenadas
	real*8, intent(in) :: t 		! t
	real*8, intent(inout) :: y(1:4)

	integer :: i
	real*8 :: K(4,4), f_ynt

	do i = 1, 4
		call fn(m, y(1), y(2), y(3), y(4), t, i, f_ynt)
		K(i,1) = h*f_ynt
	end do 

	do i=1,4
		call fn(m, y(1)+K(1,1)/2, y(2)+k(2,1)/2, y(3)+K(3,1)/2, y(4)+K(4,1)/2, t+h/2, i, f_ynt)
		K(i,2) = h*f_ynt
	end do 

	do i=1,4
		call fn(m, y(1)+K(1,2)/2, y(2)+k(2,2)/2, y(3)+K(3,2)/2, y(4)+K(4,2)/2, t+h/2, i, f_ynt)
		K(i,3) = h*f_ynt
	end do 

	do i=1,4 
		call fn(m, y(1)+K(1,3), y(2)+k(2,3), y(3)+K(3,3), y(4)+K(4,3), t+h, i, f_ynt)
		K(i,4) = h*f_ynt
	end do

	do i = 1, 4
		y(i) = y(i) + (K(i,1) + 2*K(i,2) + 2*K(i,3) + K(i,4))/6
	end do
	
end subroutine alg_RK


!----------------------
! 	(1) r/dt = ...
!----------------------
subroutine dr_dt(m, r, phi, p_r, p_phi, t, result)
	implicit none
	real*8, intent(in) :: r, phi, p_r, p_phi, m, t  
	real*8, intent(out) :: result

	result = p_r/m  
	
end subroutine dr_dt

!----------------------
! 	(2) phi/dt = ...
!----------------------
subroutine dphi_dt(m, r, phi, p_r, p_phi, t, result)
	implicit none
	real*8, intent(in) ::r, phi, p_r, p_phi, m, t 
	real*8, intent(out) :: result

	result = p_phi/(m*r**2)
	
end subroutine dphi_dt


!----------------------
! 	(3) dp_r/dt = ...
!----------------------
subroutine dpr_dt(m, r, phi, p_r, p_phi, t, result)
	implicit none
	real*8, intent(in) :: r, phi, p_r, p_phi, m, t 
	real*8, intent(out) :: result 

	real*8, parameter :: G = 6.67E-11 
	real*8, parameter :: Mt = 5.9736E24 

	result = p_phi**2/(m*r**3) - G*Mt*m/r**2
	
end subroutine dpr_dt


!----------------------
! 	(4) dp_phi/dt = ...
!----------------------
subroutine dpphi_dt(m, r, phi, p_r, p_phi, t, result)
	implicit none
	real*8, intent(in) :: r, phi, p_r, p_phi, m, t
	real*8, intent(out) :: result

	result = 0
	
end subroutine dpphi_dt

!=============================================================
!	Subrutina que unifica todas mis funciones en una
!	
!	La idea es que me proporcione claridad al resto 
!	del codigo. Concretamente al llamarlo en la sub
!	del algortimo Runge-Kutta
!=============================================================
subroutine fn(m, r, phi, p_r, p_phi, t, i, result)
	implicit none 
	integer, intent(in) :: i 
	real*8, intent(in) :: r, phi, p_r, p_phi, t, m 

	real*8, intent(inout) :: result

	if (i==1) then
		call dr_dt(m, r, phi, p_r, p_phi, t, result)
	else if (i==2) then
		call dphi_dt(m, r, phi, p_r, p_phi, t, result)
	else if (i==3) then 
		call dpr_dt(m, r, phi, p_r, p_phi, t, result)
	else if (i==4) then
		call dpphi_dt(m, r, phi, p_r, p_phi, t, result)
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