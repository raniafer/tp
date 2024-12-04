module m
	use OMP_LIB
	use iso_fortran_env
	contains

	subroutine init(n, dx, f0)
		implicit none
		integer :: i
		integer :: n
		real, dimension(n) :: f0
		real, parameter :: pi=acos(-1.)
		real :: a, b
		real, intent(out) :: dx

		! condition initiale
		a=0.; b=2.
		dx=(b-a)/(n)

		do i=1, n
			f0(i)=2+cos(pi*i*dx)
			!print*, f0(i)
		end do

	end subroutine 	

	subroutine smb(n, f0, F, dx, D, V)
                 implicit none
                 integer :: n, i
                 real :: dx, D, V
                 real, dimension (n), intent(in) :: f0
                 real, dimension (n) :: df1, df2
		 real, dimension (n), intent(out) :: F                

 		 !print*, "df1"
                 do i=2, n-1
                         df1(i)=(f0(i)-f0(i-1))/dx
			! print*, df1(i)
                 end do
 
		 !print*, "df2"
                 do i=2, n-1
                         df2(i)=(f0(i+1)-2*f0(i)+f0(i-1))/(dx**2)
			! print*, df2(i)
                 end do  

		! ghost cells
		df1(1)=df1(2)
		df1(n-1)=df1(n)
		df2(1)=df2(2)
		df2(n-1)=df2(n)

                 do i=1, n
                         F(i)=-V*df1(i)+D*df2(i)
                        ! print*, F(i)
                 end do

         end subroutine

	subroutine integre(n, f0, F, dt)
		implicit none

		integer :: n, j
		real, dimension(n) :: F, f0
		real :: dt

                 do j=1, n
                         f0(j)=f0(j)+dt*F(j)
                         !print*, "f0", f0(j)
                         !print*, "F", F(j)
                 end do

	end subroutine

end module

program diffusion
	use m
	implicit none 
	
	integer :: n=50
	real, dimension(50) :: f0 ! solution
	real, dimension(50) :: F !, df1, df2
	integer :: i, j, Nt
	real :: D=0.06, V=1.5, dt, dx
	real :: t1, t2

	open(unit=20, file="initial.txt")
	open(unit=30, file="result.txt")

	call init(50, dx, f0)

	dt=0.001*dx/V; Nt=2./dt
	print*, dx, dt, Nt

	do i=1, n
		write(*, *) i, f0(i)
	end do

	do i=1, Nt
		call smb(n, f0, F, dx, D, V)
		call integre(n, f0, F, dt)
		!print*, "__________"
		!do j=1,n
		!	print*, f0(j)
		!end do 
	end do


	do i=1, n
		write(*, *) i, f0(i)
	end do

	close(20)
	close(30)
end program
