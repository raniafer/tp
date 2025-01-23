program p
        use iso_fortran_env
        implicit none
        include 'mpif.h'

	integer, parameter :: n=10
	real, dimension(:), allocatable :: x
	real, dimension(:), allocatable :: y
	real :: tmp
	integer :: i, np	
        integer :: erreur, monCPU, totalCPU ! monCPU commence par 0

        call mpi_init(erreur)
        call mpi_comm_size(mpi_comm_world, totalCPU, erreur)
        call mpi_comm_rank(mpi_comm_world, monCPU, erreur)
	
	np=int(n/totalCPU)
	allocate(x(np))

	do i=1, np
		x(i)=1.
	end do

!	if(monCPU .eq. 0) then
!		call mpi_send(x(n*nCPU), 1, mpi_real, 1, 100, mpi_comm_world, erreur)
!	else if (monCPU.eq.1) then 
!		call mpi_recv(x(n*nCPU), 1, mpi_real, 0, 100, mpi_comm_world, erreur)
!	end if


end program
