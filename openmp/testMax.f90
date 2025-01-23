program para
	implicit none 
	include "mpif.h"
	
         integer :: N2, i
         integer :: erreur, NCPU, myID, statut(MPI_STATUS_SIZE)
         real :: cm, cmMax

         ! initialisation mpi
         call mpi_init(erreur)
         call mpi_comm_size(mpi_comm_world, NCPU, erreur)
         call mpi_comm_rank(mpi_comm_world, myID, erreur)

	do i=0,NCPU-1
		cm=i-5
	end do 

         if (myID.gt.0) then
                 cmMax=cm
                 call MPI_SEND(cmMax, 1, MPI_real, 0, 100+myID, mpi_comm_world, erreur)
         else 
                 do i=1, NCPU-1
                         call MPI_RECV(cmMax, 1, MPI_real, i, 100+i, mpi_comm_world, statut, erreur)
                         print*, "cm envoye ", cm
                         if (cm.gt.cmMax) cmMax=cm
                         print*, "cm maximale ", cmMax
                 end do
         end if
         
         call mpi_finalize(erreur)


end program para
