module m
	use OMP_LIB

	contains
	function f(x) result(y)
		implicit none
		real :: x, y

		y=4/(1+x**2)
	end function

	subroutine integ(a, b, n, intg)
		implicit none
		integer :: i, n
		real :: a, b, h, intg, som, x
	
		som=0
		h=(b-a)/n

		!$OMP PARALLEL do reduction(+:som)
		do i=0,n
			x=a+i*h
			som=som+f(x)
		end do
		!$OMP END PARALLEL do

		intg=h/2*(f(a)+2*som+f(b))
	end subroutine
end module
program valeur_pi
	use m
	implicit none
	real :: a=0, b=1, intg
	integer :: n=10000000, k
	real :: t1, t2
	
	do k=2, 16, 2
		call omp_set_num_threads(k)

		t1=omp_get_wtime()	
		call integ(a, b, n, intg)
		t2=omp_get_wtime()

		print*, "pi= ", intg, ", number of threads= ", k, "t= ", t2-t1
	end do

end program
