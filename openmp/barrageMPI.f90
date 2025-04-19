module dpres
  integer, parameter :: dp=kind(1.D0) ! double precision
end module dpres

module FVHLL ! module for finite volume scheme
contains
subroutine integre(N,dt,dx,h,hu,fh,fu) ! subroutine d'integration (Godunov Scheme)
  use dpres
  implicit none
  integer :: i,N
  real(kind=dp) :: dt,dx
  real(kind=dp), dimension(0:N+1) :: h,hu,fh,fu
  ! schema volume finis u[i]^(n+1) = u[i]^n + dt( F[i+1/2] - F[i-1/2] )
  do i=1,N
    h(i) = h(i)  - dt*(fh(i) - fh(i-1))/dx
    hu(i)= hu(i) - dt*(fu(i) - fu(i-1))/dx
  enddo
end subroutine integre

  ! routine de calcul des flux, schema HLL : flux solution au travers des surfaces
  ! => resolution d'un problème de Riemann

subroutine flux(N,h,hu,fh,fu,cmax)
  use dpres
  implicit none
  integer :: i,N
  real(kind=dp) :: g=9.81D0,cm,cmax,hg,hd,hug,hud,ug,ud,cg,cd,c1,c2
  real(kind=dp), dimension(0:N+1) :: h,hu,fh,fu
  cmax=0.D0
  do i=0,N
    hg = h(i)     ; hd = h(i+1);
    ug = hu(i)/hg ; ud = hu(i+1)/hd; !   ___
    ! calcul des vitesses d'ondes c  =  Vg.h
    cg = dsqrt(g*hg) ; cd = dsqrt(g*hd)
    ! calcul des vitesses d'ondes
    c1 = dmin1(ug-cg,ud-cd);
    c2 = dmax1(ug+cg,ud+cd);
    if(c1.ge.0.D0) then ! toutes les ondes traversent à droite
      fh(i) = hg*ug
      fu(i) = hg*ug*ug + 0.5D0*g*hg*hg
      cm    = dabs(c2)
    else if(c2.le.0.D0) then ! toutes les ondes traversent à gauche
      fh(i) = hd*ud
      fu(i) = hd*ud*ud + 0.5D0*g*hd*hd
      cm    = dabs(c1)
    else ! cas d'un problème de Riemann
         ! Flux HLL, F = (c2*fg - c1*fd)/(c2-c1)  + c1*c2*(Ud - Ug)/(c2-c1) 
      fh(i) = ( c2*hg*ug - c1*hd*ud )/(c2-c1) + c1*c2*(hd-hg)/(c2-c1)
      fu(i) = ( c2*( hg*ug*ug + 0.5d0*g*hg*hg ) - c1*( hd*ud*ud + 0.5d0*g*hd*hd ))/(c2-c1) + c1*c2*(hd*ud-hg*ug)/(c2-c1)
      cm    = dmax1(dabs(c1),dabs(c2)) ! onde la plus rapide
    endif
    cmax = dmax1(cmax,cm)      
  enddo
end subroutine flux
end module FVHLL

module In_Out ! module for initialisation and output
contains
subroutine init(N,dx,x,h,hu,x0) ! initialisation : dam break problem
  use dpres                    ! h : discontinuous, u = 0
  implicit none
  integer :: i,N
  real(kind=dp), dimension(0:N+1) :: h,hu ,x
  real(kind=dp) :: dx, x0
  do i=1,N
    x(i) = x0 + ( dfloat(i-1) + 0.5D0) *dx
    if(x(i).lt.0.5D0) h(i) = 10.D0
    if(x(i).ge.0.5D0) h(i) = 1.D0
    hu(i) = 0.D0
  enddo
end subroutine init
end module In_Out

module para ! le module regroupe les routines utilisées pour la parallelization
   include "mpif.h"

   contains
   subroutine hEchange(myID,NCPU,N,u)
   use dpres
   implicit none
   integer :: myID,NCPU,N
   integer :: source, dest, erreur,statut(MPI_STATUS_SIZE)
   real(kind=dp), dimension(0:N+1) :: u

        ! echange des cellules fantome
        dest = myID + 1
        if (myID .eq. NCPU-1) dest = mpi_proc_null	! n'envoie a personne
        call MPI_SEND(u(N),1,MPI_DOUBLE_PRECISION,dest,100+dest,MPI_COMM_WORLD,erreur)
        source = myID - 1
        if (myID .eq. 0) source = mpi_proc_null         ! ne recoit de personne
        call MPI_RECV(u(0),1,MPI_DOUBLE_PRECISION,source,100+myID,MPI_COMM_WORLD,statut,erreur)


        dest = myID - 1
        if(myID .eq. 0) dest = mpi_proc_null            ! n'envoie a personne
        call MPI_SEND(u(1),1,MPI_DOUBLE_PRECISION,dest,200+dest,MPI_COMM_WORLD,erreur)
        source = myID + 1
        if (myID .ge. NCPU-1) source = mpi_proc_null    ! ne recoit de personne
        call MPI_RECV(u(N+1),1,MPI_DOUBLE_PRECISION,source,200+myID,MPI_COMM_WORLD,statut,erreur)

        ! echange des conditions aux limites
        if(myID.eq.0) u(1)=u(2)
        if(myID.eq.NCPU-1) u(N)=u(N-1)
   end subroutine hEchange

   subroutine huEchange(myID,NCPU,N,u)
   use dpres
   implicit none
   integer :: myID,NCPU,N
   integer :: source, dest, erreur,statut(MPI_STATUS_SIZE)
   real(kind=dp), dimension(0:N+1) :: u

        ! echange des cellules fantome
        dest = myID + 1
        if (myID .eq. NCPU-1) dest = mpi_proc_null	! n'envoie a personne
        call MPI_SEND(u(N),1,MPI_DOUBLE_PRECISION,dest,100+dest,MPI_COMM_WORLD,erreur)
        source = myID - 1
        if (myID .eq. 0) source = mpi_proc_null         ! ne recoit de personne
        call MPI_RECV(u(0),1,MPI_DOUBLE_PRECISION,source,100+myID,MPI_COMM_WORLD,statut,erreur)

	dest = myID - 1
        if(myID .eq. 0) dest = mpi_proc_null            ! n'envoie a personne
        call MPI_SEND(u(1),1,MPI_DOUBLE_PRECISION,dest,200+dest,MPI_COMM_WORLD,erreur)
        source = myID + 1
        if (myID .ge. NCPU-1) source = mpi_proc_null    ! ne recoit de personne
        call MPI_RECV(u(N+1),1,MPI_DOUBLE_PRECISION,source,200+myID,MPI_COMM_WORLD,statut,erreur)

        ! conditions aux limites
        if(myID.eq.0) u(1)=u(2)
        if(myID.eq.NCPU-1) u(N)=-u(N-1)
   end subroutine huEchange

   subroutine decoupe(iCPU,NCPU,N,dx,x0,NP)
   use dpres
   implicit none
   integer :: N,NP,iCPU,NCPU,reste
   real(kind=dp) :: dx,x0
   NP = N/NCPU
   reste = mod(N,NCPU)
   if(iCPU<reste) then
     NP = NP + 1
     x0 = dfloat(NP*iCPU)*dx
   else
     x0 = dfloat((NP+1)*reste+NP*(iCPU-reste))*dx
   endif
   end subroutine decoupe

   subroutine ecritMPIf(myID,NCPU,N,x,h,hu)
   use dpres
   implicit none
   integer :: N,i,id,myID,NCPU,erreur
   real(kind=dp) :: x0,dx
   real(kind=dp), dimension(0:N+1) :: x, h, hu

   do id=0,NCPU-1
     if(myID.eq.id) then
       if(myID.eq.0) then
           open(unit=50,file="finalfMPI.dat", status="unknown")
       else
           open(unit=50,file="finalfMPI.dat",status="old",access="append")
       end if
       do i=1,N
         write(50,*) x(i), h(i), hu(i)/h(i)
       enddo
       close(50)
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,erreur)
   enddo
   end subroutine ecritMPIf


   subroutine ecritMPI0(myID,NCPU,N,x,h,hu)
   use dpres
   implicit none
   integer :: N,i,id,myID,NCPU,erreur
   real(kind=dp) :: x0,dx
   real(kind=dp), dimension(0:N+1) :: x, h, hu

   do id=0,NCPU-1
     if(myID.eq.id) then
       if(myID.eq.0) then
           open(unit=51,file="initialfMPI.dat", status="unknown")
       else
           open(unit=51,file="initialfMPI.dat",status="old",access="append")
       end if
       do i=1,N
         write(51,*) x(i), h(i), hu(i)/h(i)
       enddo
       close(51)
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,erreur)
   enddo
   end subroutine ecritMPI0
end module


program ShallowWater
  use dpres
  use para ! module parallelisme
  use In_Out
  use FVHLL
  implicit none
  integer :: i,j,N,Nt,it
  real(kind=dp) :: dx,dt,cm, cmMax
  real(kind=dp), dimension(:), allocatable :: h,hu,fh,fu ,x

  ! variables mpi
  integer :: N2
  integer :: erreur, NCPU, myID, statut(MPI_STATUS_SIZE)
  real(kind=dp) :: x0
  real(kind=dp) :: startTime, endTime

  ! initialisation mpi
  call mpi_init(erreur)
  call mpi_comm_size(mpi_comm_world, NCPU, erreur)
  call mpi_comm_rank(mpi_comm_world, myID, erreur)

  open(unit=23, file="param.dat", status="old") ! fichier parametre

! lecture des donnes 
  read(23,*) N 
  read(23,*) Nt
  close(23)
  write(*,*) N,Nt
  dx = 1.D0/dfloat(N);
  dt = 1.e-8

  ! divise le domaine en NCPU-1 sous-domaines à traiter chacun par un thread
  call decoupe(myID,NCPU,N,dx,x0,N2)

  ! allocation des champs position (x), hauteur (h), impulsion (hu) et flux (fh,fu)
  allocate(x(0:N2+1),h(0:N2+1),hu(0:N2+1),fh(0:N2+1),fu(0:N2+1))

  call init(N2,dx,x,h,hu,x0)  ! initialisation et sauvegarde
  call ecritMPI0(myID,NCPU,N2,x,h,hu) ! écriture de la condition initiale
  close(24)
	
  !call MPI_BARRIER(MPI_COMM_WORLD, erreur) ! synchronisation entre les threads
  startTime=MPI_Wtime()

  ! boucle en temps
  do j = 1,Nt
  
    ! conditions aux limites
    call hEchange(myID, NCPU, N2, h)
    call huEchange(myID, NCPU, N2, hu)
    
    call flux(N2,h,hu,fh,fu,cm);      ! calcul des flux HLL et vitesse d'onde (cm)

    ! stock dans cmMax la valeur maximale sur l'ensemble des threads
    call MPI_ALLREDUCE(cm, cmMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, erreur)

    dt = 0.75D0*dx/cmMax;               ! pas de temps admissible (CFL=0.75)

    call integre(N2,dt,dx,h,hu,fh,fu) ! integration sur 1 pas de temps (Godunov Scheme)

  enddo
  endTime=MPI_Wtime()

  !print*, "elapsed time", endTime-startTime

  call ecritMPIf(myID,NCPU,N2,x,h,hu) ! sauvegarde du resultat

  deallocate(x);deallocate(h);deallocate(hu);deallocate(fh);deallocate(fu)
  call MPI_Finalize(erreur)
end program ShallowWater
