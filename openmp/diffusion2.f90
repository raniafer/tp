module m
	use OMP_LIB
	contains

	subroutine init(n, dx, f0)
		implicit none
		integer :: i
		integer :: n
		real, dimension(n) :: f0
		real, parameter :: pi=acos(-1.)
		real :: a, b, dx, x

		! condition initiale
		a=0; b=2
		dx=(b-a)/(n)

		do i=1, n
			x=i*dx
			f0(i)=2+cos(pi*x)
			print*, f0(i)
		end do

	end subroutine 	

	subroutine smb(n, df1, df2, f0, D, V, F)
                 implicit none
                 integer :: n, i
                 real :: dx
                 real, dimension (n) :: f0
                 real, dimension (n) :: df1, df2
                 real, dimension (n) :: F
                 real :: D, V

                 call init(n, dx, f0)

 		 print*, "df1"
                 do i=2, n-1
			df1(1)=0
                         df1(i)=(f0(i)-f0(i-1))/dx
			df1(n)=df1(n-1)
			 print*, df1(i)
                 end do
 
		 print*, "df2"
                 do i=2, n-1
			df2(1)=0
                         df2(i)=(f0(i+1)-2*f0(i)+f0(i-1))/(dx**2)
			df2(n)=df2(n-1)
			 print*, df2(i)
                 end do  

 		 print*, "F(i)"
                 do i=1, n
                         F(i)=-V*df1(i)+D*df2(i)
			 print*, F(i)
                 end do

         end subroutine

end module

program diffusion
	use m
	implicit none 
	
	integer :: n=20
	real, dimension(20) :: f0
	real, dimension(20) :: df1, df2, F
	integer :: i, Nt=500
	real :: D=1, V=0.01, dt=0.1, dx
	n=20
	
	call smb(20, df1, df2, f0, D, V, F)

	open(unit=20, file="result.txt")

	do i=1, n
		write(20,*) i, f0(i), F(i)
	end do

	close(20)
end program
