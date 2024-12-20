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

	subroutine smb(n, f0, D, V, F, dx, i, df1, df2)
                 implicit none
                 integer :: n, i
                 real, intent(in) :: dx, D, V
                 real, dimension (n), intent(in) :: f0
                 real, dimension (n) :: df1, df2
                 real, dimension (n), intent(out) :: F

 		 !print*, "df1"
		! first derivative
!$omp do
                 do i=2, n-1
                         df1(i)=(f0(i)-f0(i-1))/dx
			 df2(i)=(f0(i+1)-2*f0(i)+f0(i-1))/(dx**2)
			 F(i)=-V*df1(i)+D*df2(i)
                 end do
!$omp end do

		! periodic boundary conditions
		df1(1)=df1(n-1)
		df1(n)=df1(2)
		df2(1)=df2(n-1)
		df2(n)=df2(2)
		F(1)=-V*df1(1)+D*df2(1)
		F(n)=-V*df1(n)+D*df2(n)
		
 		 !print*, "F(i)"
                ! do i=1, n
                !         F(i)=-V*df1(i)+D*df2(i)
			! print*, F(i)
                ! end do

         end subroutine

end module

program diffusion
	use m
	implicit none 
	
	integer :: n=50
	real, dimension(50) :: f0=0., f0_new=0., df1, df2
	real, dimension(50) :: F=0.
	integer :: io, i, j, k, Nt, i1
	real :: D=0.5, V=5., dt, dx
	real :: t1, t2
	integer :: n_threads=4
	
	open(unit=20, file="initial.txt", iostat=io)
	open(unit=30, file="solution.txt", iostat=io)
	open(unit=40, file="time.txt", iostat=io)

	! initialize
	call init(n, dx, f0, i1)

	! stability condition
	dt=0.0001*dx/V; Nt=1.5/dt
	print*, "dx=", dx, "dt=", dt, "Nt=", Nt

	do i=1, n
		write(20, *) i, f0(i)
		print*, f0(i)
	end do


	call omp_set_num_threads(n_threads)
	t1=omp_get_wtime()
	
	!$omp parallel private(i1) shared(n, df1, df2)
	!do k=1, n_threads
		do i=1, Nt
			
			call smb(n, f0, D, V, F, dx, i1, df1, df2)

			do j=1, n
				f0(j)=f0(j)+dt*F(j)
			end do
			
		end do

	t2=omp_get_wtime()
	!$omp end parallel

		print*, "apres parallelisation", (t2-t1), "s pour ", n_threads, "processeurs"
		!write(40,*) (t2-t1), n_threads
	!end do

	do i=1, n
		write(*, *) i, f0(i)
	end do

	close(20)
	close(30)
	close(40)
end program
