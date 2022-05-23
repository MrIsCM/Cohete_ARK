!==========================================
! 					ICM
!
! 	  	Algoritmo de Runge-Kutta
! 	Sobre-Simplificacion: Tierra y Cohete 
! 			 (1 sola coord.)
!
!==========================================
program simple
	implicit none
	
	integer :: i 
	real*8 :: t, r, pr 

	real*8, parameter :: m = 1000
	real*8, parameter :: h = 1
	real*8, parameter :: Rt = 6.378160E6

	open(10, file='DataTest/SimplePos.dat', status='Unknown')

	! Condiciones iniciales
	r = Rt
	pr = m*Rt*0.0015 
	
	do i = 0, 10000
		t = t+h 
		if (mod(i,10)==0) then 
			write(10, *) t, r/Rt, pr/(m*Rt) 
		end if 
		call Runge_Kutta(r, pr, t, m, h)
	end do

	close(10)

	
end program simple

subroutine Runge_Kutta(r, pr, t, m, h)
	real*8, intent(in) :: t, m, h 
	real*8, intent(inout) ::  r, pr

	real*8 :: K_r(4), K_pr(4)
	real*8 :: result 

	! --------------------------------------
	! 					K^(1)
	! --------------------------------------
	call dr_dt(r, pr, t, m, result)
	K_r(1) = h*result

	call dpr_dt(r, pr, t, m, result)
	K_pr(1) = h*result 


	! --------------------------------------
	! 					K^(2)
	! --------------------------------------
	call dr_dt(r+K_r(1)/2, pr+K_pr(1)/2, t+h/2, m, result)
	K_r(2) = h*result

	call dpr_dt(r+K_r(1)/2, pr+K_pr(1)/2, t+h/2, m, result)
	K_pr(2) = h*result


	! --------------------------------------
	! 					K^(3)
	! --------------------------------------
	call dr_dt(r+K_r(2)/2, pr+K_pr(2)/2, t+h/2, m, result)
	K_r(3) = h*result

	call dpr_dt(r+K_r(2)/2, pr+K_pr(2)/2, t+h/2, m, result)
	K_pr(3) = h*result


	! --------------------------------------
	! 					K^(4)
	! --------------------------------------
	call dr_dt(r+K_r(3), pr+K_pr(3), t+h, m, result)
	K_r(4) = h*result

	call dpr_dt(r+K_r(3), pr+K_pr(3), t+h, m, result)
	K_pr(4) = h*result


	! CALCULO
	r = r + (K_r(1)+2*K_r(2)+2*K_r(3)+K_r(4))/6
	pr = pr + (K_pr(1)+2*K_pr(2)+2*K_pr(3)+K_pr(4))/6

	
end subroutine Runge_Kutta

subroutine dr_dt(r, pr, t, m, result)
	implicit none
	real*8, intent(in) :: r, pr, t, m 
	real*8, intent(out) :: result

	result = pr/m
	
end subroutine dr_dt

subroutine dpr_dt(r, pr, t, m, result)
	implicit none
	real*8, intent(in) :: r, pr, t, m 
	real*8, intent(out) :: result

	real*8, parameter :: G = 6.67E-11
	real*8, parameter :: Mt = 5.9736E24

	result = -G*Mt*m/r**2
	
end subroutine dpr_dt