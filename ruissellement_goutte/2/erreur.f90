program average
	use iso_fortran_env

	implicit none

	integer :: i, n, nmax
	integer :: io
	character(len=5) :: file
	real, dimension(:), allocatable :: h
	real :: line

	!print*, "enter number of cells"
	!read*, file

	open(unit=20, file="125.txt")

	n=0
print*, "reading"
	       	do 
        		read(20,*,iostat=io) line
        		if (io==iostat_end) exit
    			n=n+1
			print*, line
        	end do
	rewind(20)

	allocate(h(n))

end program
