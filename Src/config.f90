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

! Variables and routines used to read and store the configuration.
module config
  use iso_fortran_env
  use misc
  use data
  use symmetry
  implicit none

  integer(kind=4) :: nelements,natoms,ngrid(3),norientations
  namelist /allocations/ nelements,natoms,ngrid,norientations
  real(kind=8) :: lfactor,lattvec(3,3),epsilon(3,3)
  character(len=3),allocatable :: elements(:)
  integer(kind=4),allocatable :: types(:),orientations(:,:)
  integer(kind=4) :: scell(3)
  real(kind=8),allocatable :: positions(:,:),masses(:),gfactors(:),born(:,:,:)
  namelist /crystal/ lfactor,lattvec,elements,types,positions,masses,gfactors,&
       epsilon,born,scell,orientations
  integer(kind=4) :: maxiter,nticks
  real(kind=8) :: T,scalebroad,rmin,rmax,dr,eps
  namelist /parameters/ T,scalebroad,rmin,rmax,dr,maxiter,nticks,eps
  logical :: nonanalytic,convergence,isotopes,autoisotopes,nanowires,onlyharmonic,espresso
  namelist /flags/ nonanalytic,convergence,isotopes,autoisotopes,&
       nanowires,onlyharmonic,espresso

  integer(kind=4) :: nbands,nptk,nwires
  real(kind=8) :: cgrid,V,rV,rlattvec(3,3),slattvec(3,3)
  real(kind=8),allocatable :: cartesian(:,:),uorientations(:,:)

  integer(kind=4) :: nsymm
  integer(kind=4),allocatable :: rotations(:,:,:)
  real(kind=8),allocatable :: crotations(:,:,:),qrotations(:,:,:)
  character(len=10) :: international

  ! MPI variables, assigned in ShengBTE.f90.
  integer(kind=4) :: myid,numprocs

contains

  subroutine read_config()
    integer(kind=4) :: i,j,k,ii,jj,kk,ll,info
    integer(kind=4) :: P(3)
    integer(kind=4),allocatable :: rtmp(:,:,:),ID_equi(:,:)
    logical,allocatable :: valid(:)
    real(kind=8) :: dnrm2,tmp1(3,3),tmp2(3,3),tmp3(3,3)
    real(kind=8),allocatable :: crtmp(:,:,:),qrtmp(:,:,:)
    real(kind=8),allocatable :: translations(:,:),ctranslations(:,:)

    ! Set the defaults and read the namelist.
    nelements=0
    natoms=0
    ngrid=0
    norientations=0
    open(1,file="CONTROL",status="old")
    read(1,nml=allocations)
    if(norientations.lt.0) then
       if(myid.eq.0)write(error_unit,*) "Error: norientations must be >=0"
       stop 1
    end if
    if(nelements.lt.1.or.natoms.lt.1.or.natoms.lt.nelements) then
       if(myid.eq.0)write(error_unit,*) "Error: nelements,natoms must be >0, natoms must be >=nelements"
       stop 1
    end if
    if(any(ngrid.lt.1)) then
       if(myid.eq.0)write(error_unit,*) "Error: all components of ngrid must be must be >0"
       stop 1
    end if
    allocate(elements(nelements),types(natoms),positions(3,natoms),&
         masses(nelements),gfactors(nelements),born(3,3,natoms),&
         cartesian(3,natoms))
    if(norientations.ne.0) then
       allocate(orientations(3,norientations))
       allocate(uorientations(3,norientations))
       orientations=0
    end if
    lfactor=1.
    scell=-1
    epsilon=0.
    epsilon(1,1)=1.
    epsilon(2,2)=1.
    epsilon(3,3)=1.
    born=0.
    gfactors=0.
    types=0
    read(1,nml=crystal)
    if(.not.all(scell.gt.0)) then
       if(myid.eq.0)write(error_unit,*) "Error: all supercell sizes must be >0"
       stop 1
    end if
    if(.not.all(types.ne.0)) then
       if(myid.eq.0)write(error_unit,*) "Error: atom types must be initialized correctly"
       stop 1
    end if
    do i=1,norientations
       if(all(orientations(:,i).eq.0)) then
          if(myid.eq.0)write(error_unit,*) "Error: orientation uninitialized or zero"
          stop 1
       end if
    end do
    T=-1.
    scalebroad=1.0
    rmin=5.0
    rmax=505.0
    dr=100.0
    maxiter=1000
    nticks=100
    eps=1e-5
    read(1,nml=parameters)
    if(T.le.0.) then
       if(myid.eq.0)write(error_unit,*) "Error: T must be >0 K"
       stop 1
    end if
    if(rmin.le.0.or.rmax.le.rmin.or.dr.le.0) then
       if(myid.eq.0)write(error_unit,*) "Error: rmin and dr must be >0, and rmax must be > rmin"
       stop 1
    end if
    if(maxiter.le.0) then
       if(myid.eq.0)write(error_unit,*) "Error: maxiter must be >0"
       stop 1
    end if
    nonanalytic=.true.
    convergence=.true.
    isotopes=.true.
    autoisotopes=.true.
    nanowires=.false.
    onlyharmonic=.false.
    espresso=.false.
    read(1,nml=flags)
    if(nanowires.and.norientations.eq.0) then
       if(myid.eq.0)write(error_unit,*) "Error: nanowires=.TRUE. but norientations=0"
       stop 1
    end if
    close(1)

    nptk=product(ngrid)
    cgrid=nptk**(1./3.)
    nbands=3*natoms
    nwires=ceiling((rmax-rmin)/dr)

    lattvec=lfactor*lattvec

    do i=1,3
       j=mod(i,3)+1
       k=mod(j,3)+1
       call cross_product(lattvec(:,j),lattvec(:,k),rlattvec(:,i))
    end do

    V=abs(dot_product(lattvec(:,1),rlattvec(:,1)))
    rV=2.*pi/V
    rlattvec=rV*rlattvec

    cartesian=matmul(lattvec,positions)

    if(nanowires) then
       uorientations=matmul(lattvec,orientations)
       do i=1,norientations
          uorientations(:,i)=uorientations(:,i)/&
               dnrm2(3,uorientations(:,i),1)
       end do
    end if

    ! Compute the average masses and g-factors automatically, if
    ! requested.
    if(autoisotopes) then
       call data_fill_isotopes()
       do i=1,nelements
          call data_calc_mandg(elements(i),masses(i),gfactors(i))
       end do
       call data_free_isotopes()
    end if

    ! Find out the symmetries of the system.
    nsymm=get_num_operations(lattvec,natoms,types,positions)
    allocate(translations(3,nsymm),rotations(3,3,nsymm),&
         ctranslations(3,nsymm),crotations(3,3,nsymm),qrotations(3,3,nsymm))
    call get_operations(lattvec,natoms,types,positions,nsymm,&
         rotations,translations,international)
    if(myid.eq.0)write(*,*) "Info: symmetry group ",trim(international)," detected"
    if(myid.eq.0)write(*,*) "Info: ",nsymm," symmetry operations"
    call get_cartesian_operations(lattvec,nsymm,rotations,translations,&
         crotations,ctranslations)
    deallocate(translations,ctranslations)
    ! Transform the rotation matrices to the reciprocal-space basis.
    do i=1,nsymm
       tmp1=matmul(transpose(lattvec),lattvec)
       tmp2=transpose(rotations(:,:,i))
       tmp3=tmp1
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       qrotations(:,:,i)=transpose(matmul(tmp2,tmp3))
    end do
    ! Find rotations that are either duplicated or incompatible with
    ! the q-point grid.
    allocate(ID_Equi(nsymm,nptk),valid(nsymm))
    valid=.TRUE.
    ll=0
    do ii=2,nsymm
       do jj=1,ii-1
          if(.not.valid(jj))cycle
          if(all(rotations(:,:,ii).eq.rotations(:,:,jj))) then
             valid(ii)=.FALSE.
             ll=ll+1
             exit
          end if
       end do
    end do
    if(myid.eq.0.and.ll.ne.0)write(*,*) "Info:",ll,"duplicated rotations will be discarded"
    call symmetry_map(ID_equi)
    jj=0
    do ii=1,nsymm
       if(valid(ii).and.any(ID_equi(ii,:).eq.-1)) then
          valid(ii)=.FALSE.
          jj=jj+1
       end if
    end do
    if(myid.eq.0.and.jj.ne.0)write(*,*) "Info:",jj,"rotations are incompatible with the q-point grid and will be discarded"
    ! Filter out those rotations through a series of move_alloc calls.
    ! Arrays to take into account: rotations,crotations,qrotations.
    if(ll+jj.ne.0) then
       call move_alloc(rotations,rtmp)
       call move_alloc(crotations,crtmp)
       call move_alloc(qrotations,qrtmp)
       allocate(rotations(3,3,nsymm-ll-jj),crotations(3,3,nsymm-ll-jj),&
            qrotations(3,3,nsymm-ll-jj))
       kk=0
       do ii=1,nsymm
          if(valid(ii)) then
             kk=kk+1
             rotations(:,:,kk)=rtmp(:,:,ii)
             crotations(:,:,kk)=crtmp(:,:,ii)
             qrotations(:,:,kk)=qrtmp(:,:,ii)
          end if
       end do
       nsymm=nsymm-ll-jj
       deallocate(rtmp,crtmp,qrtmp)
    end if
    deallocate(ID_Equi,valid)
  end subroutine read_config

  subroutine free_config()
    deallocate(elements,types,positions,masses,gfactors,born,cartesian,&
         rotations,crotations,qrotations)
    if(nanowires) then
       deallocate(orientations,uorientations)
    end if
  end subroutine free_config

  
  ! Compute all images of a reciprocal-space point using the rotational part of the
  ! symmetry operations. Everything is performed in lattice coordinates.
  subroutine symm(r_in,r_out)
    implicit none
    integer(kind=4),intent(in) :: r_in(3)
    real(kind=8),intent(out) :: r_out(3,nsymm)

    integer(kind=4) :: ii

    do ii=1,nsymm
       r_out(:,ii)=ngrid*matmul(qrotations(:,:,ii),dble(r_in)/ngrid)
    end do
  end subroutine symm

  ! Find the equivalences among points.
  subroutine symmetry_map(ID_equi)
    implicit none
    integer(kind=4),intent(out) :: ID_equi(nsymm,nptk)

    integer(kind=4) :: Ind_cell(3,nptk)
    integer(kind=4) :: i,isym,ivec(3)
    real(kind=8) :: vec(3),vec_symm(3,nsymm),dnrm2

    call Id2Ind(Ind_cell)
    do i=1,nptk
       call symm(Ind_cell(:,i),vec_symm)
       do isym=1,Nsymm
          vec=vec_symm(:,isym)
          ivec=nint(vec)
          if(dnrm2(3,abs(vec-dble(ivec)),1).gt.1e-6) then
             ID_equi(isym,i)=-1
          else
             ID_equi(isym,i)=Ind2Id(modulo(ivec,ngrid))
          end if
       end do
    end do

  end subroutine symmetry_map

  ! Create a table that can be used to demultiplex cell indices.
  subroutine Id2ind(Ind_cell)
    implicit none
    integer(kind=4),intent(out) :: Ind_cell(3,nptk)

    integer(kind=4) :: ii

    do ii=1,nptk
       Ind_cell(3,ii)=int((ii-1)/(Ngrid(1)*Ngrid(2)))
       Ind_cell(2,ii)=int(modulo(ii-1,Ngrid(1)*Ngrid(2))/Ngrid(1))
       Ind_cell(1,ii)=modulo(ii-1,Ngrid(1))
    end do
  end subroutine Id2ind

  ! Multiplex three cell indices into one.
  function Ind2Id(Ind_cell)
    implicit none
    integer(kind=4),intent(in) :: Ind_cell(3)

    integer(kind=4) :: Ind2Id

    Ind2Id=1+Ind_cell(1)+(Ind_cell(2)+Ind_cell(3)*Ngrid(2))*Ngrid(1)
  end function Ind2Id

  ! Return the base broadening (without prefactor) for a mode.
  function base_sigma(v)
    implicit none
    real(kind=8),intent(in) :: v(3)

    real(kind=8) :: base_sigma

    integer(kind=4) :: nu

    base_sigma=0.
    do nu=1,3
       base_sigma=base_sigma+(dot_product(rlattvec(:,nu),v)/ngrid(nu))**2
    end do

    base_sigma=sqrt(base_sigma/6.)
  end function base_sigma
end module config
