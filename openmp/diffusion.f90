module m
	use OMP_LIB
	contains

	subroutine init(n, dx, f0)
		implicit none
		integer :: i
		integer, parameter :: n
		real, dimension(n) :: f0
		real, parameter :: pi=acos(-1.)
		real :: a, b, dx, x

		! condition initiale
		a=0; b=2
		dx=(b-a)/(n)

		do i=1, n
			x=i*dx
			f0(i)=2+cos(pi*x)
		end do

	end subroutine

	function smb(f0) result(F)
		implicit none
		integer :: n, i
		real :: dx
		real, dimension (n) :: f0
		real, dimension (n-1) :: df1
		real, dimension (n-2) :: df2, F
		real :: D, V

		call init(n, dx, f0)

		do i=2, n
			df1(i)=(f0(i)-f0(i-1))/dx
		end do
		
		do i=2, n-1
			df2(i)=(f0(i+1)-2*f0(i)+f0(i-1))/dx**2
		end do	
		
		do i=1, n
			F(i)=-V*df1(i)+D*df2(i)
		end do


	end function
	
	subroutine integre(n, dt, D, V, f0, F)
		implicit none 
		integer :: n, i
		real, dimension(n-1) :: df1
		real, dimension(n-2) :: df2
		real, dimension(n-2) :: F, f0
		real :: dt, D, V
	
		!call smb(n, df1, df2, f0)
		

	!	do i=1, n-2
	!		F(i)=-V*df1(i)+D*df2(i)
	!		f0(i)=f0(i)+dt*F(i)
	!	end do
	end subroutine
end module

program diffusion
	use m
	implicit none 
	
	integer, parameter :: n=20
	real, dimension(n-2) :: f0, F
	integer :: i, Nt=500
	real :: D, V, dt=0.1
	
	
!	call init(n, dx, f0)
!	call smb(n, df1, df2)
	
	call integre(n, dt, D, V, f0, F)

	do i=1, Nt
		f0(i)=f0(i)+dt*F(i)	
	end do

	D=0.1; V=15

	do i=1, n-2
		print*, f0(i)
	end do

end program
