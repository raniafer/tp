module m
	use OMP_LIB
	use iso_fortran_env

	contains

	subroutine init(n, dx, f0, i)
		implicit none
		integer :: i
		integer :: n
		real, dimension(n), intent(out) :: f0
		real, parameter :: pi=acos(-1.)
		real :: a, b, dx, x

		! initial condition
		a=-1.; b=1.
		dx=(b-a)/(n)

		do i=1, n
			x=-1.+real(i-1)*dx
			f0(i)=2.+cos(pi*x)
		end do
	end subroutine 	

end module

program diffusion
	use m
	implicit none 
	
	integer :: n=100
	real, dimension(100) :: f0=0., df1, df2
	real, dimension(100) :: F=0.
	integer :: io, i, j, k, Nt, i1
	real :: D=0.5, V=5., dt, dx
	real :: t1, t2
	integer :: n_threads=6
	
	open(unit=20, file="initial.txt", iostat=io)
	open(unit=30, file="solution.txt", iostat=io)
	open(unit=40, file="time.txt", iostat=io)

	! initialize
	call init(n, dx, f0, i1)

	! stability condition
	dt=0.0001*dx/V; Nt=0.5/dt
	print*, "dx=", dx, "dt=", dt, "Nt=", Nt

	do i=1, n
		write(20, *) i, f0(i)
		print*, f0(i)
	end do


	call omp_set_num_threads(n_threads)

	t1=omp_get_wtime()		
	!$omp parallel private(i, i1, j) 
		do i=1, Nt

			!$omp do schedule(runtime)
			do i1=2, n-1
                         df1(i1)=(f0(i1)-f0(i1-1))/dx
			 df2(i1)=(f0(i1+1)-2*f0(i1)+f0(i1-1))/(dx**2)
			 F(i1)=-V*df1(i1)+D*df2(i1)
                 	end do
			!$omp end do

!$omp single
		! periodic boundary conditions
			df1(1)=df1(n-1)
			df1(n)=df1(2)
			df2(1)=df2(n-1)
			df2(n)=df2(2)
			F(1)=-V*df1(1)+D*df2(1)
			F(n)=-V*df1(n)+D*df2(n)
		! time stepping loop, not parallelized
!$omp end single			
			!$omp do schedule(runtime)
			do j=1, n
				f0(j)=f0(j)+dt*F(j)
			end do
			!$omp end do

		end do
	!$omp end parallel		
	t2=omp_get_wtime()			

		print*, "apres parallelisation", (t2-t1), "s pour ", n_threads, "coeurs"
		!write(40,*) (t2-t1), n_threads
	!end do

	do i=1, n
		write(*, *) i, f0(i)
	end do

	close(20)
	close(30)
	close(40)
end program
