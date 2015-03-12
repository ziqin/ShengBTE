!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2013 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2013 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2013 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2013 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Main program.
program ShengBTE
  use config
  use input
  use processes
  use wedgetc
  use iterations
  use conductivity
  use scaling
  use phonon_routines
  use gruneisen
  use integrals
  implicit none

  include "mpif.h"

  real(kind=8) :: omega,sigma,kappa_sg(3,3),kappa_old(3,3),relchange
  integer(kind=4) :: i,j,ii,jj,kk,ll,mm,nn
  real(kind=8),allocatable :: energy(:,:),q0(:,:),velocity(:,:,:),velocity_z(:,:)
  complex(kind=8),allocatable :: eigenvect(:,:,:)

  real(kind=8),allocatable :: grun(:,:)

  real(kind=8),allocatable :: rate_scatt(:,:),rate_scatt_reduce(:,:)
  real(kind=8),allocatable :: tau_zero(:,:),tau(:,:),tau_b(:,:),tau2(:,:),tau_b2(:,:)
  real(kind=8),allocatable :: dos(:,:),pdos(:,:,:),rate_scatt_isotope(:,:)
  real(kind=8),allocatable :: F_n(:,:,:),F_n_0(:,:,:),F_n_aux(:,:)
  real(kind=8),allocatable :: ThConductivity(:,:,:)
  real(kind=8),allocatable :: ThConductivityMode(:,:,:,:)
  real(kind=8),allocatable :: ticks(:),cumulative_kappa(:,:,:,:)

  integer(kind=4) :: Ntri
  ! IJK are vectors of product(ngrid) points in the lattice coordinate
  ! system, with each component ranging from 0 to Ngrid(:)-1.  Index_N
  ! is a mapping of 3 indices for an individual q-point into one index.
  integer(kind=4),allocatable :: Index_i(:),Index_j(:),Index_k(:),IJK(:,:),Index_N(:,:,:)
  real(kind=8),allocatable :: Phi(:,:,:,:),R_j(:,:),R_k(:,:)

  integer(kind=4) :: nlist,Ntotal_plus,Ntotal_minus
  integer(kind=4),allocatable :: nequi(:),list(:),AllEquiList(:,:),TypeofSymmetry(:,:),eqclasses(:)
  integer(kind=4),allocatable :: N_plus(:),N_minus(:),N_plus_reduce(:),N_minus_reduce(:)
  integer(kind=4),allocatable :: Naccum_plus(:),Naccum_minus(:)
  integer(kind=4),allocatable :: Indof2ndPhonon_plus(:),Indof3rdPhonon_plus(:)
  integer(kind=4),allocatable :: Indof2ndPhonon_minus(:),Indof3rdPhonon_minus(:)
  integer(kind=4),allocatable :: Indof2ndPhonon_plus_reduce(:),Indof3rdPhonon_plus_reduce(:)
  integer(kind=4),allocatable :: Indof2ndPhonon_minus_reduce(:),Indof3rdPhonon_minus_reduce(:)
  real(kind=8) :: radnw,kappa_or_old
  real(kind=8),allocatable :: Gamma_plus(:),Gamma_minus(:)
  real(kind=8),allocatable :: Pspace_plus_partial(:,:),Pspace_plus_total(:,:),Pspace_plus_tmp(:)
  real(kind=8),allocatable :: Pspace_minus_partial(:,:),Pspace_minus_total(:,:),Pspace_minus_tmp(:)
  real(kind=8),allocatable :: Gamma_plus_reduce(:),Gamma_minus_reduce(:)
  real(kind=8),allocatable :: ffunc(:,:),radnw_range(:),v_or(:,:),F_or(:,:)
  real(kind=8),allocatable :: kappa_or(:),kappa_wires(:,:),kappa_wires_reduce(:,:)

  integer(kind=4) :: iorient,ierr
  character(len=4) :: aux
  character(len=128) :: sorientation

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  ! Read CONTROL and initialize variables such as the
  ! direct/reciprocal lattice vectors, the symmetry operations, the
  ! Pearson deviation coefficients related to isotopic disorder and so
  ! on and so forth.
  call read_config()

  allocate(energy(nptk,nbands),eigenvect(nptk,nbands,nbands),grun(nptk,nbands))
  allocate(q0(nptk,3),velocity(nptk,nbands,3),velocity_z(nptk,nbands))
  allocate(IJK(3,nptk),Index_N(0:(ngrid(1)-1),0:(ngrid(2)-1),0:(ngrid(3)-1)))
  allocate(F_n(nbands,nptk,3),F_n_0(nbands,nptk,3),F_n_aux(nbands,nptk))
  allocate(ThConductivity(nbands,3,3),ThConductivityMode(nptk,nbands,3,3),kappa_wires(nbands,nwires),&
       kappa_wires_reduce(nbands,nwires))
  allocate(Nequi(nptk),list(nptk),AllEquiList(nsymm,nptk),TypeOfSymmetry(nsymm,nptk),&
       eqclasses(nptk))
  allocate(ffunc(nptk,nbands),v_or(nptk,nbands),F_or(nbands,nptk),kappa_or(nbands))

  ! Obtain the q-point equivalence classes defined by symmetry
  ! operations.
  call Id2Ind(IJK)
  call wedge(Nlist,Nequi,List,ALLEquiList,TypeofSymmetry)
  do ll=1,Nlist
     do kk=1,Nequi(ll)
        eqclasses(ALLEquiList(kk,ll))=List(ll)
     end do
  end do
  

  do ii=1,Ngrid(1)
     do jj=1,Ngrid(2)
        do kk=1,Ngrid(3)
           ll=((kk-1)*Ngrid(2)+(jj-1))*Ngrid(1)+ii
           q0(ll,:)=rlattvec(:,1)*(ii-1.0)/ngrid(1)+&
                rlattvec(:,2)*(jj-1.0)/ngrid(2)+&
                rlattvec(:,3)*(kk-1.0)/ngrid(3)
        end do
     end do
  end do

  if(myid.eq.0) then
     open(1,file="BTE.qpoints",status="replace")
     do ll=1,Nlist
        write(1,"(I9,x,I9,x,3(E20.10,x))") List(ll),Nequi(ll),q0(List(ll),:)
     end do
     close(1)
     open(1,file="BTE.qpoints_full",status="replace")
     do ll=1,nptk
        write(1,"(I9,x,I9,x,3(E20.10,x))") ll,eqclasses(ll),q0(ll,:)
     end do
     close(1)
  end if

  ! Obtain the phonon spectrum from the 2nd-order force constants (in
  ! either of the two formats suported) and the dielectric parameters
  ! in CONTROL.
  if(myid.eq.0) then
     write(*,*) "Info: about to obtain the spectrum"
     if(espresso) then
        write(*,*) "Info: expecting Quantum Espresso 2nd-order format"
     else
        write(*,*) "Info: expecting Phonopy 2nd-order format"
     end if
  end if
  call eigenDM(energy,eigenvect,velocity)
  if(myid.eq.0)write(*,*) "Info: spectrum calculation finished"

  ! Compute the harmonic integrals: lattice specific heat and
  ! small-grain-limit reduced thermal conductivity. Write out this
  ! information, as well as the spectrum itself.
  if(myid.eq.0) then
     open(1,file="BTE.cv",status="replace")
     write(1,*) cv(energy)
     close(1)
  end if

  call kappasg(energy,velocity,kappa_sg)
  if(myid.eq.0) then
     open(1,file="BTE.kappa_sg",status="replace")
     write(1,"(9E20.10)") kappa_sg
     close(1)
  end if

  write(aux,"(I0)") 3*Nbands
  if(myid.eq.0) then
     open(1,file="BTE.v",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") velocity(list(ll),:,:)
     end do
     close(1)
     open(1,file="BTE.v_full",status="replace")
     do ll=1,nptk
        write(1,"("//trim(adjustl(aux))//"E20.10)") velocity(ll,:,:)
     end do
  end if
  write(aux,"(I0)") Nbands
  if(myid.eq.0) then
     open(1,file="BTE.omega",status="replace")
     open(2,file="BTE.omega_full",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") energy(list(ll),:)
     end do
     do ll=1,nptk
        write(2,"("//trim(adjustl(aux))//"E20.10)") energy(ll,:)
     end do
     close(1)
     close(2)
  end if

  ! Locally adaptive estimates of the total and projected densities of states.
  allocate(dos(Nbands,Nlist),pdos(Nbands,Nlist,natoms),&
       rate_scatt_isotope(Nbands,Nlist))
  dos=0.d0
  pdos=0.d0
  rate_scatt_isotope=0.d0

  do mm=1,Nlist
     do nn=1,Nbands
        omega=energy(list(mm),nn)
        do ii=1,nptk
           do jj=1,Nbands
              sigma=base_sigma(velocity(ii,jj,:))
              if(abs(omega-Energy(ii,jj)).lt.2.5*sigma) then
                 dos(nn,mm)=dos(nn,mm)+exp(-(omega-Energy(ii,jj))**2/(sigma**2))/sigma/sqrt(pi)
              end if
           end do
        end do
     end do
  end do
  dos=dos/float(nptk)

  if(myid.eq.0) then
     open(1,file="BTE.dos",status="replace")
     do mm=1,Nlist
        do nn=1,Nbands
           write(1,"(2E25.17)") energy(list(mm),nn),dos(nn,mm)
        end do
     end do
     close(1)
  end if

  do mm=1,Nlist
     do nn=1,Nbands
        omega=energy(list(mm),nn)
        do ii=1,nptk
           do jj=1,Nbands
              sigma=base_sigma(velocity(ii,jj,:))
              if(abs(omega-Energy(ii,jj)).lt.2.5*sigma) then
                 do kk=1,natoms
                    pdos(nn,mm,kk)=pdos(nn,mm,kk)+&
                         exp(-(omega-Energy(ii,jj))**2/(sigma**2))/sigma/sqrt(Pi)*&
                         (abs(dot_product(eigenvect(list(mm),nn,((kk-1)*3+1):((kk-1)*3+3)),&
                         eigenvect(ii,jj,((kk-1)*3+1):((kk-1)*3+3)))))**2
                 end do
              end if
           end do
        end do
     end do
  end do
  pdos=pdos/float(nptk)

  write(aux,"(I0)") Natoms
  if(myid.eq.0) then
     open(1,file="BTE.pdos",status="replace")
     do mm=1,Nlist
        do nn=1,Nbands
           write(1,"(E25.17,"//trim(adjustl(aux))//"E25.17)") energy(list(mm),nn),pdos(nn,mm,:)
        end do
     end do
     close(1)
  end if

  ! Isotopic scattering, closely related to the projected DOS.
  if(isotopes) then
     do mm=1,Nlist
        do nn=1,Nbands
           omega=energy(list(mm),nn)
           do ii=1,nptk
              do jj=1,Nbands
                 sigma=base_sigma(velocity(ii,jj,:))
                 if(abs(omega-Energy(ii,jj)).lt.2.5*sigma) then
                    do kk=1,natoms
                       rate_scatt_isotope(nn,mm)=rate_scatt_isotope(nn,mm)+&
                            exp(-(omega-Energy(ii,jj))**2/(sigma**2))/sigma/sqrt(Pi)*&
                            (abs(dot_product(eigenvect(list(mm),nn,((kk-1)*3+1):((kk-1)*3+3)),&
                            eigenvect(ii,jj,((kk-1)*3+1):((kk-1)*3+3))))**2*gfactors(types(kk)))
                    end do
                 end if
              end do
           end do
           rate_scatt_isotope(nn,mm)=rate_scatt_isotope(nn,mm)/(2.d0*nptk)*pi*omega**2
        end do
     end do

     write(aux,"(I0)") Nbands
     if(myid.eq.0) then
        open(1,file="BTE.w_isotopic",status="replace")
        do ll=1,Nlist
           write(1,"("//trim(adjustl(aux))//"E20.10)") rate_scatt_isotope(:,ll)
        end do
        close(1)
     end if
  end if

  ! N_plus for absorption processes and N_minus for emission processes.
  allocate(N_plus(Nlist*Nbands),N_minus(Nlist*Nbands))
  N_plus=0
  N_minus=0
  allocate(N_plus_reduce(Nlist*Nbands),N_minus_reduce(Nlist*Nbands))
  N_plus_reduce=0
  N_minus_reduce=0
  allocate(Naccum_plus(Nlist*Nbands),Naccum_minus(Nlist*Nbands))
  allocate(rate_scatt_reduce(Nbands,Nlist))
  rate_scatt_reduce=0.d0
  allocate(rate_scatt(Nbands,Nlist),tau_zero(Nbands,Nlist),tau(Nbands,Nlist),tau_b(Nbands,Nlist),tau2(Nbands,nptk),tau_b2(Nbands,nptk))
  rate_scatt=0.d0
  allocate(radnw_range(nwires))
  do ii=1,nwires
     radnw_range(ii)=rmin+(ii-1.0)*dr
  end do

  ! Compute and print the number of allowed absorption and emission processes.
  do mm=myid+1,Nbands*Nlist,numprocs
     call Nprocesses(mm,ii,jj,energy,velocity,Nlist,List(1:Nlist),IJK)
     N_plus_reduce(mm)=ii
     N_minus_reduce(mm)=jj
  end do

  call MPI_ALLREDUCE(N_plus_reduce,N_plus,Nbands*Nlist,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(N_minus_reduce,N_minus,Nbands*Nlist,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  deallocate(N_plus_reduce,N_minus_reduce)

  Ntotal_plus=0
  Ntotal_minus=0
  do mm=1,Nbands*Nlist
     Ntotal_plus=Ntotal_plus+N_plus(mm)
     Ntotal_minus=Ntotal_minus+N_minus(mm)
  end do
  if(myid.eq.0)write(*,*) "Info: Ntotal_plus =",Ntotal_plus
  if(myid.eq.0)write(*,*) "Info: Ntotal_minus =",Ntotal_minus

  ! Obtain the phase space volume per mode and their sum.
  allocate(Pspace_plus_total(Nbands,Nlist),Pspace_plus_partial(Nbands,Nlist),&
       Pspace_plus_tmp(maxval(N_plus)))
  allocate(Pspace_minus_total(Nbands,Nlist),Pspace_minus_partial(Nbands,Nlist),&
       Pspace_minus_tmp(maxval(N_minus)))
  Pspace_plus_total=0.
  Pspace_plus_partial=0.
  Pspace_minus_total=0.
  Pspace_minus_partial=0.
  do mm=myid+1,Nbands*Nlist,numprocs
     Pspace_plus_tmp=0.
     if(N_plus(mm).ne.0) then
        ii=modulo(mm-1,Nbands)+1
        jj=int((mm-1)/Nbands)+1
        call D_plus(mm,N_plus(mm),energy,velocity,Nlist,List,IJK,&
             Pspace_plus_tmp)
        Pspace_plus_partial(ii,jj)=Pspace_plus_partial(ii,jj)+sum(Pspace_plus_tmp)
     end if
     Pspace_minus_tmp=0.
     if(N_minus(mm).ne.0) then
        ii=modulo(mm-1,Nbands)+1
        jj=int((mm-1)/Nbands)+1
        call D_minus(mm,N_minus(mm),energy,velocity,Nlist,List,IJK,&
             Pspace_minus_tmp)
        Pspace_minus_partial(ii,jj)=Pspace_minus_partial(ii,jj)+sum(Pspace_minus_tmp)
     end if
  end do

  call MPI_ALLREDUCE(Pspace_plus_partial,Pspace_plus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
       MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Pspace_minus_partial,Pspace_minus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
       MPI_SUM,MPI_COMM_WORLD,ierr)

  if(myid.eq.0) then
     open(1,file="BTE.P3_plus",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_plus_total(:,ll)
     end do
     close(1)
     open(1,file="BTE.P3_minus",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_minus_total(:,ll)
     end do
     close(1)
     open(1,file="BTE.P3",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") 2.*(Pspace_plus_total(:,ll)+&
             Pspace_minus_total(:,ll)/2.)/3.
     end do
     close(1)
     do ii=1,Nlist
        Pspace_plus_total(:,ii)=Pspace_plus_total(:,ii)*Nequi(ii)
     end do
     open(1,file="BTE.P3_plus_total",status="replace")
     write(1,*) sum(Pspace_plus_total)
     close(1)
     do ii=1,Nlist
        Pspace_minus_total(:,ii)=Pspace_minus_total(:,ii)*Nequi(ii)
     end do
     open(1,file="BTE.P3_minus_total",status="replace")
     write(1,*) sum(Pspace_minus_total)
     close(1)
     open(1,file="BTE.P3_total",status="replace")
     write(1,*) 2.*(sum(Pspace_plus_total)+sum(Pspace_minus_total)/2.)/3.
     close(1)
  end if

  deallocate(Pspace_plus_total,Pspace_plus_partial,Pspace_plus_tmp)
  deallocate(Pspace_minus_total,Pspace_minus_partial,Pspace_minus_tmp)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Up to this point, no anharmonic information is used.
  if(onlyharmonic) then
     write(*,*) "Info: onlyharmonic=.true., stopping here"
     call MPI_FINALIZE(ierr)
     stop
  end if

  ! Load the anharmonic IFCs from FORCE_CONSTANTS_3RD.
  call read3fc(Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)

  call mode_grun(energy,eigenvect,Ntri,Phi,R_j,R_k,&
       Index_i,Index_j,Index_k,grun)
  write(aux,"(I0)") Nbands
  if(myid.eq.0) then
     open(1,file="BTE.gruneisen",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") grun(list(ll),:)
     end do
     close(1)
     open(1,file="BTE.gruneisen_total",status="replace")
     write(1,*) total_grun(energy,grun)
     close(1)
  end if
  deallocate(grun)

  Naccum_plus(1)=0
  Naccum_minus(1)=0
  do mm=2,Nbands*Nlist
     Naccum_plus(mm)=Naccum_plus(mm-1)+N_plus(mm-1)
     Naccum_minus(mm)=Naccum_minus(mm-1)+N_minus(mm-1)
  end do

  ! Here begins the most expensive part of the calculation: obtaining the
  ! three-phonon scattering amplitudes for all processes.
  allocate(Indof2ndPhonon_plus_reduce(Ntotal_plus),Indof3rdPhonon_plus_reduce(Ntotal_plus),&
       Indof2ndPhonon_minus_reduce(Ntotal_minus),Indof3rdPhonon_minus_reduce(Ntotal_minus),&
       Gamma_plus_reduce(Ntotal_plus),Gamma_minus_reduce(Ntotal_minus))
  Indof2ndPhonon_plus_reduce=0
  Indof3rdPhonon_plus_reduce=0
  Indof2ndPhonon_minus_reduce=0
  Indof3rdPhonon_minus_reduce=0
  Gamma_plus_reduce=0.d0
  Gamma_minus_reduce=0.d0
  do mm=myid+1,Nbands*NList,numprocs
     i=modulo(mm-1,Nbands)+1
     ll=int((mm-1)/Nbands)+1
     if((N_plus(mm).ne.0)) then
        allocate(Indof2ndPhonon_plus(N_plus(mm)),Indof3rdPhonon_plus(N_plus(mm)))
        allocate(Gamma_plus(N_plus(mm)))
        call Ind_plus(mm,N_plus(mm),energy,velocity,eigenvect,Nlist,List,&
             Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
             Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus)
        Indof2ndPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=Indof2ndPhonon_plus
        Indof3rdPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=Indof3rdPhonon_plus
        Gamma_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=Gamma_plus
        rate_scatt_reduce(i,ll)=rate_scatt_reduce(i,ll)+sum(Gamma_plus)
        deallocate(Indof2ndPhonon_plus,Indof3rdPhonon_plus)
        deallocate(Gamma_plus)
     end if
     if((N_minus(mm).ne.0)) then
        allocate(Indof2ndPhonon_minus(N_minus(mm)),Indof3rdPhonon_minus(N_minus(mm)))
        allocate(Gamma_minus(N_minus(mm)))
        call Ind_minus(mm,N_minus(mm),energy,velocity,eigenvect,Nlist,List,&
             Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
             Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus)
        Indof2ndPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=Indof2ndPhonon_minus
        Indof3rdPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=Indof3rdPhonon_minus
        Gamma_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=Gamma_minus
        rate_scatt_reduce(i,ll)=rate_scatt_reduce(i,ll)+sum(Gamma_minus)*5.D-1
        deallocate(Indof2ndPhonon_minus,Indof3rdPhonon_minus)
        deallocate(Gamma_minus)
     end if
  end do

  call MPI_ALLREDUCE(rate_scatt_reduce,rate_scatt,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
       MPI_SUM,MPI_COMM_WORLD,ierr)

  write(aux,"(I0)") Nbands
  if(myid.eq.0) then
     open(1,file="BTE.w_anharmonic",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") rate_scatt(:,ll)
     end do
     close(1)
  end if

  ! Obtain the total scattering rates in the relaxation time approximation.
  rate_scatt=rate_scatt+rate_scatt_isotope
  if(myid.eq.0) then
     open(1,file="BTE.w",status="replace")
     do ll = 1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") rate_scatt(:,ll)
     end do
     close(1)
  end if

  deallocate(rate_scatt_reduce)
  allocate(Indof2ndPhonon_plus(Ntotal_plus),Indof3rdPhonon_plus(Ntotal_plus),&
       Indof2ndPhonon_minus(Ntotal_minus),Indof3rdPhonon_minus(Ntotal_minus))
  allocate(Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus))

  call MPI_ALLREDUCE(Indof2ndPhonon_plus_reduce,Indof2ndPhonon_plus,Ntotal_plus,&
       MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Indof3rdPhonon_plus_reduce,Indof3rdPhonon_plus,Ntotal_plus,&
       MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Indof2ndPhonon_minus_reduce,Indof2ndPhonon_minus,Ntotal_minus,&
       MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Indof3rdPhonon_minus_reduce,Indof3rdPhonon_minus,Ntotal_minus,&
       MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Gamma_plus_reduce,Gamma_plus,Ntotal_plus,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Gamma_minus_reduce,Gamma_minus,Ntotal_minus,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  deallocate(Indof2ndPhonon_plus_reduce,Indof3rdPhonon_plus_reduce,&
       Indof2ndPhonon_minus_reduce,Indof3rdPhonon_minus_reduce,&
       Gamma_plus_reduce,Gamma_minus_reduce)

  tau_zero=0.d0
  do ll = 1,Nlist
     do i=1,Nbands
        if(rate_scatt(i,ll).ne.0) then
           tau_zero(i,ll)=1./rate_scatt(i,ll)
        end if
     end do
  end do

  ! Set up everything to start the iterative process.
  call iteration0(Nlist,Nequi,ALLEquiList,energy,velocity,tau_zero,F_n)
  F_n_0=F_n
  if(myid.eq.0) then
     ! Open all output files.
     open(2001,file="BTE.kappa",status="replace")
     open(2002,file="BTE.kappa_tensor",status="replace")
     open(2003,file="BTE.kappa_scalar",status="replace")
     call TConduct(energy,velocity,F_n,ThConductivity,ThConductivityMode)
     write(aux,"(I0)") 9*Nbands
     write(2001,"(I9,"//trim(adjustl(aux))//"E20.10)") 0,ThConductivity
     write(2002,"(I9,9E20.10,E20.10)") 0,sum(ThConductivity,dim=1)
     write(2003,"(I9,E20.10,E20.10)") 0,&
          sum(sum(ThConductivity,dim=1),reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.

     open(2004,file="BTE.kappa_mode",status="replace")
     do ll = 1,nptk
        write(2004,"("//trim(adjustl(aux))//"E20.10)") ThConductivityMode(ll,:,:,:)
     end do
     close(2004)

     ! Iterate to convergence if desired.
     if(convergence) then
        do ii=1,maxiter
           kappa_old=sum(ThConductivity,dim=1)
           call iteration(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
                Ntotal_plus,Ntotal_minus,Indof2ndPhonon_plus,Indof3rdPhonon_plus,&
                Indof2ndPhonon_minus,Indof3rdPhonon_minus,energy,velocity,&
                Gamma_plus,Gamma_minus,tau_zero,F_n)
           call TConduct(energy,velocity,F_n,ThConductivity,ThConductivityMode)
           write(2001,"(I9,"//trim(adjustl(aux))//"E20.10)") ii,ThConductivity
           write(2002,"(I9,9E20.10)") ii,sum(ThConductivity,dim=1)
           write(2003,"(I9,E20.10)") ii,&
                sum(sum(ThConductivity,dim=1),reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.
           relchange=twonorm3x3(sum(ThConductivity,dim=1)-kappa_old)/&
                twonorm3x3(kappa_old)
           write(*,*) "Info: Iteration",ii
           write(*,*) "Info:","Relative change","=",relchange
           if(relchange.lt.eps)exit
        end do

        open(2004,file="BTE.kappa_mode",status="replace")       ! output the converged values
        do ll = 1,nptk
           write(2004,"("//trim(adjustl(aux))//"E20.10)") ThConductivityMode(ll,:,:,:)
        end do
        close(2004)

     end if
     close(2001)
     close(2002)
     close(2003)

     ! Write out the converged scattering rates.
     do ll=1,Nlist
        do ii=1,Nbands
           tau(ii,ll)=dot_product(F_n(ii,List(ll),:),velocity(List(ll),ii,:))/&
                (dot_product(velocity(List(ll),ii,:),velocity(List(ll),ii,:))*energy(List(ll),ii))
           tau_b(ii,ll)=1/dot_product(velocity(List(ll),ii,:),velocity(List(ll),ii,:))
        end do
     end do
     write(aux,"(I0)") Nbands
     open(1,file="BTE.w_final",status="replace")
     open(2,file="BTE.w_boundary",status="replace")
     do ll = 1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") 1./tau(:,ll)
        write(2,"("//trim(adjustl(aux))//"E20.10)") 1./tau_b(:,ll)
     end do
     close(1)
     close(2)
     do ll=1,nptk
        do ii=1,Nbands
           tau2(ii,ll)=dot_product(F_n(ii,ll,:),velocity(ll,ii,:))/&
                (dot_product(velocity(ll,ii,:),velocity(ll,ii,:))*energy(ll,ii))
           tau_b2(ii,ll)=1/dot_product(velocity(ll,ii,:),velocity(ll,ii,:))
        end do
     end do
     write(aux,"(I0)") Nbands
     open(1,file="BTE.w_final_full",status="replace")
     open(2,file="BTE.w_boundary_full",status="replace")
     do ll = 1,nptk
        write(1,"("//trim(adjustl(aux))//"E20.10)") 1./tau2(:,ll)
        write(2,"("//trim(adjustl(aux))//"E20.10)") 1./tau_b2(:,ll)
     end do
     close(1)
     close(2)

     ! If results for nanowires have been requested, obtain a lower bound
     ! for the thermal conductivity along each crystallographic orientation
     ! by using the bulk RTA results.
     if(nanowires) then
        do iorient=1,norientations
           write(sorientation,"(I128)") iorient
           open(3001,file="BTE.kappa_nw_"//trim(adjustl(sorientation))//"_lower",status="replace")
           do ii=1,nptk
              do jj=1,Nbands
                 v_or(ii,jj)=dot_product(velocity(ii,jj,:),uorientations(:,iorient))
                 F_or(jj,ii)=dot_product(F_n_0(jj,ii,:),uorientations(:,iorient))
              end do
           end do
           do mm=1,Nwires
              radnw=radnw_range(mm)
              call ScalingOfTau(Nlist,Nequi,ALLEquiList,v_or,velocity,tau_zero,radnw,ffunc)
              do ii=1,nptk
                 do jj=1,Nbands
                    F_n_aux(jj,ii)=F_or(jj,ii)*ffunc(ii,jj)
                 end do
              end do
              call TConductScalar(energy,v_or,F_n_aux,kappa_or)
              write(3001,"(E30.20,"//trim(adjustl(aux))//"E20.10,E20.10)") 2.d0*radnw,&
                   kappa_or,sum(kappa_or)
           end do
           close(3001)
        end do
     end if
     allocate(ticks(nticks),cumulative_kappa(nbands,3,3,nticks))
     ! Cumulative thermal conductivity.
     call CumulativeTConduct(energy,velocity,F_n,ticks,cumulative_kappa)
     write(aux,"(I0)") 9*nbands+1
     open(2001,file="BTE.cumulative_kappa",status="replace")
     open(2002,file="BTE.cumulative_kappa_tensor",status="replace")
     open(2003,file="BTE.cumulative_kappa_scalar",status="replace")
     do ii=1,nticks
        write(2001,"("//trim(adjustl(aux))//"E20.10)") ticks(ii),cumulative_kappa(:,:,:,ii)
        write(2002,"(10E20.10)") ticks(ii),&
             sum(cumulative_kappa(:,:,:,ii),dim=1)
        write(2003,"(2E20.10)") ticks(ii),&
             sum(sum(cumulative_kappa(:,:,:,ii),dim=1),&
             reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.
     end do
     close(2001)
     close(2002)
     close(2003)
     deallocate(ticks,cumulative_kappa)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! If results for nanowires have been requested, repeat the iterative process
  ! for each desired orientation, introducing the appropriate scaling of
  ! relaxation times.
  if(nanowires) then
     kappa_wires=0.d00
     kappa_wires_reduce=0.d00
     kk=ceiling(float(nwires)/numprocs)
     do iorient=1,norientations
        write(*,"(A,I0,A,3(x,I0))") "Info: nanowires with orientation ",&
             iorient,":",orientations(:,iorient)
        write(sorientation,"(I128)") iorient
        do ii=1,nptk
           do jj=1,Nbands
              v_or(ii,jj)=dot_product(velocity(ii,jj,:),uorientations(:,iorient))
              F_or(jj,ii)=dot_product(F_n_0(jj,ii,:),uorientations(:,iorient))
           end do
        end do
        do mm=myid+1,nwires,numprocs
           radnw=radnw_range(mm)
           call ScalingOfTau(Nlist,Nequi,ALLEquiList,v_or,velocity,tau_zero,radnw,ffunc)
           do ii=1,nptk
              do jj=1,Nbands
                 F_n_aux(jj,ii)=F_or(jj,ii)*ffunc(ii,jj)
              end do
           end do
           call TConductScalar(energy,v_or,F_n_aux,kappa_or)
           if(convergence) then
              do ii=1,maxiter
                 kappa_or_old=sum(kappa_or)
                 call iteration_scalar(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
                      Ntotal_plus,Ntotal_minus,Indof2ndPhonon_plus,Indof3rdPhonon_plus,&
                      Indof2ndPhonon_minus,Indof3rdPhonon_minus,energy,v_or,&
                      Gamma_plus,Gamma_minus,tau_zero,F_n_aux)
                 do ll=1,nptk
                    do jj=1,nbands
                       F_n_aux(jj,ll)=F_n_aux(jj,ll)*ffunc(ll,jj)
                    end do
                 end do
                 call TConductScalar(energy,v_or,F_n_aux,kappa_or)
                 relchange=abs((sum(kappa_or)-kappa_or_old)/kappa_or_old)
                 if(relchange.lt.eps)exit
              end do
           end if
           kappa_wires_reduce(:,mm)=kappa_or
        end do
        call MPI_ALLREDUCE(kappa_wires_reduce,kappa_wires,Nbands*Nwires,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        write(aux,"(I0)") 3*nbands
        open(3001,file="BTE.kappa_nw_"//trim(adjustl(sorientation)),status="replace")
        do ii=1,Nwires
           radnw=radnw_range(ii)
           write(3001,"(E30.20,"//trim(adjustl(aux))//"E20.10,E20.10)") 2.d0*radnw,&
                kappa_wires(:,ii),sum(kappa_wires(:,ii))
        end do
        close(3001)
     end do
  end if

  if(myid.eq.0)write(*,*) "Info: normal exit"

  call free_config()
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
end program ShengBTE
