module m
	use OMP_LIB
	implicit none

	contains
	subroutine moy_tab(n)
		implicit none
		integer :: i, n
		real, dimension(n) :: a, b

		do i=1, n
			a(i)=real(i)
		end do

		print*, "b= "

		!$OMP PARALLEL PRIVATE(i) SHARED(n, a, b)
		!$OMP DO
		do i=2, n-1

			!if (i==1) then
			!	b(i)= 0.5*(a(i)+a(i+1))
			!else if (i==n) then
			!	b(i)= 0.5*(a(i)+a(i-1))
			!else
				b(i)=0.5*a(i)+0.25*a(i-1)+0.25*a(i+1)
			!end if
		end do
		!$OMP END DO
		!$OMP END PARALLEL

		b(1)=0.5*(a(2)+a(1))
		b(n)=0.5*(a(n)+a(n-1))

		!do i=1, n
		!	print*, b(i)
		!end do		
	end subroutine

end module

program tab
	use m
	implicit none 
	integer :: k
	real :: t1, t2

open(unit=10, file="resultat.dat")

	do k=2, 16, 2
                t1=omp_get_wtime()      

		print*, "k= ", k

		call omp_set_num_threads(k)
		call moy_tab(5000000)

                t2=omp_get_wtime()      
                print*, "number of threads= ", k, "t= ", t2-t1
		write (10,*) k, t2-t1
	end do

close(10)
end program
