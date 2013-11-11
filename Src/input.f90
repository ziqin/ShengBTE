!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2013 Wu Li <wu.li@cea.fr>
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

! Routines to read some input formats.
module input
  use iso_fortran_env
  use config
  implicit none

  character(len=*),parameter :: filename_2fc="FORCE_CONSTANTS_2ND"
  character(len=*),parameter :: filename_3fc="FORCE_CONSTANTS_3RD"
  real(kind=8),parameter :: unitfactor=9648.53336213 ! from eV/(A^2*amu) to THz^2

contains

  subroutine read2fc(fc)
    implicit none
    real(kind=8),allocatable,intent(out) :: fc(:,:,:,:,:,:,:)
    
    integer(kind=4) :: ntot,atom1,atom2,i,j,ip
    integer(kind=4) :: ix1,iy1,iz1,ix2,iy2,iz2,iatom1,iatom2
    real(kind=8) :: mm(natoms,natoms)

    do i=1,natoms
       mm(i,i)=masses(types(i))
       do j=i+1,natoms
          mm(i,j)=sqrt(masses(types(i))*masses(types(j)))
          mm(j,i)=mm(i,j)
       end do
    end do

    allocate(fc(natoms,3,scell(1),scell(2),scell(3),natoms,3))

    open(1,file=filename_2fc,status="old")
    read(1,*) ntot
    if(ntot.ne.scell(1)*scell(2)*scell(3)*natoms) then
       if(myid.eq.0)write(error_unit,*) "Error: wrong number of force constants for the specified scell"
       stop 1
    end if
    do i=1,ntot
       do j=1,ntot
          read(1,*) atom1,atom2
          call split_index(atom1,scell(1),scell(2),scell(3),&
               ix1,iy1,iz1,iatom1)
          call split_index(atom2,scell(1),scell(2),scell(3),&
               ix2,iy2,iz2,iatom2)
          if(ix1.eq.1.and.iy1.eq.1.and.iz1.eq.1) then
             do ip=1,3
                read(1,*) fc(iatom1,ip,ix2,iy2,iz2,iatom2,:)
             end do
          else
             do ip=1,3
                read(1,*)
             end do
          end if
       end do
    end do
    close(1)

    do iatom1=1,natoms
       do iatom2=1,natoms
          fc(iatom1,:,:,:,:,iatom2,:)=&
               fc(iatom1,:,:,:,:,iatom2,:)/mm(iatom1,iatom2)
       end do
    end do
    fc=unitfactor*fc
  end subroutine read2fc

  subroutine split_index(index,nx,ny,nz,ix,iy,iz,iatom)
    implicit none

    integer(kind=4),intent(in) :: index,nx,ny,nz
    integer(kind=4),intent(out) :: ix,iy,iz,iatom

    integer(kind=4) :: tmp1,tmp2

    call divmod(index-1,nx,tmp1,ix)
    call divmod(tmp1,ny,tmp2,iy)
    call divmod(tmp2,nz,iatom,iz)

    ix=ix+1
    iy=iy+1
    iz=iz+1
    iatom=iatom+1
  end subroutine split_index

  subroutine read3fc(Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
    implicit none
    integer(kind=4),intent(out) :: Ntri
    integer(kind=4),allocatable,intent(out) :: Index_i(:),Index_j(:),Index_k(:)
    real(kind=8),allocatable,intent(out) :: Phi(:,:,:,:),R_j(:,:),R_k(:,:)

    real(kind=8) :: tmp(3,3)
    integer(kind=4) :: ii,jj,ll,mm,nn,ltem,mtem,ntem,info,P(3)

    open(1,file=filename_3fc,status="old")
    read(1,*) Ntri
    allocate(Index_i(Ntri),Index_j(Ntri),Index_k(Ntri))
    allocate(Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri))
    do ii=1,Ntri
       read(1,*) jj
       read(1,*) R_j(:,ii)
       read(1,*) R_k(:,ii)
       read(1,*) Index_i(ii),Index_j(ii),Index_k(ii)
       do ll=1,3
          do mm=1,3
             do nn=1,3
                read(1,*) ltem,mtem,ntem,Phi(ll,mm,nn,ii)
             end do
          end do
       end do
    end do
    close(1)
    tmp=lattvec
    call dgesv(3,Ntri,tmp,3,P,R_j,3,info)
    R_j=matmul(lattvec,anint(R_j/10.))
    tmp=lattvec
    call dgesv(3,Ntri,tmp,3,P,R_k,3,info)
    R_k=matmul(lattvec,anint(R_k/10.))
  end subroutine read3fc
end module input
