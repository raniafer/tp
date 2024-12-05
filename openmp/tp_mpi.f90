program p
        use iso_fortran_env
        implicit none
        include 'mpif.h'

	integer, parameter :: n=10, nCPU=2
	real, dimension(n*nCPU) :: x=1.
	real, dimension(n*nCPU) :: y
	integer :: i	
        integer :: erreur, monCPU, totalCPU

        call mpi_init(erreur)
        call mpi_comm_size(mpi_comm_world, totalCPU, erreur)
        call mpi_comm_rank(mpi_comm_world, monCPU, erreur)



if(monCPU .eq. 0) then

	do i=2,n-1
		y(i)=x(i)+x(i-1)
	end do

else if (monCPU .eq. 1) then

	do i=n+2, n*nCPU
		y(i)=x(i)+x(i-1)
	end do	

end if

	

end program
