program p
	implicit none
	include 'mpif.h'
	
	integer :: nc, n=500000, i
	real :: x, y
	! mpi
	integer :: erreur 

	call mpi_init(erreur)

	nc=0
	do i=1,n
		call random_number(x)
		call random_number(y)
		if (sqrt(x**2+y**2) .le. 1) then
			nc=nc+1
		end if
		!print*, nc !x, y, sqrt(x**2+y**2)
	end do

	print*, 4*nc/n

	call mpi_finalize(erreur)
end program
