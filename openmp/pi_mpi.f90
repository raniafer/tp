program p
	use iso_fortran_env
	implicit none
	include 'mpif.h'
	
	integer :: nc, n=5000, i, tmp
	real :: x, y, pi
	! mpi
	integer :: erreur, monCPU, totalCPU

	call mpi_init(erreur)
	call mpi_comm_size(mpi_comm_world, totalCPU, erreur)
	call mpi_comm_rank(mpi_comm_world, monCPU, erreur)

	nc=0

	do i=1, int(n/totalCPU)
		call random_number(x)
		call random_number(y)
		if ((x**2+y**2) .le. 1) then
			nc=nc+1
		end if
	end do

	! rassemblement des resultats dans le coeurs 0
	if (monCPU .eq. 0) then
		do i=1, totalCPU-1
			call mpi_recv(tmp, 1, mpi_integer, i, 100+i, mpi_comm_world, erreur)
			nc=nc+tmp
		end do
		pi=4.*nc/n
		print*, "pi = ", pi
	else
		call mpi_send(nc, 1, mpi_integer, 0, 100+monCPU, mpi_comm_world, erreur)
	end if

	call mpi_finalize(erreur)
end program
