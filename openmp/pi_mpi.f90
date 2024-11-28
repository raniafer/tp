program p
	use iso_fortran_env
	implicit none
	include 'mpif.h'
	
	integer :: nc, n=5000000, i, ncpu=2
	real :: x, y, pi
	! mpi
	integer :: erreur, monCPU, totalCPU

	call mpi_init(erreur)
	call mpi_comm_size(mpi_comm_world, totalCPU, erreur)
	call mpi_comm_rank(mpi_comm_world, monCPU, erreur)

!	print*, monCPU

	nc=0

	!if (monCPU .eq. 0) then 
	!	call mpi_send(nc, 1, mpi_real, 1, 100, mpi_comm_world, erreur)
	!else if (monCPU .eq. 1) then
	!	call mpi_recv(nc, 1, mpi_real, 0, 100, mpi_comm_world, erreur)
	!end if

if (monCPU .eq. 0) then

	do i=1, int(n/ncpu)
		call random_number(x)
		call random_number(y)
		if (sqrt(x**2+y**2) .le. 1) then
			nc=nc+1
		end if
	end do

else if (monCPU .eq. 1) then

	do i=1, int(n/ncpu)
		call random_number(x)
		call random_number(y)
		if (sqrt(x**2+y**2) .le. 1) then
			nc=nc+1
		end if
	end do

endif


		!print*, nc !x, y, sqrt(x**2+y**2)

	pi=4*nc/n
	print*, pi
	call mpi_finalize(erreur)
end program
