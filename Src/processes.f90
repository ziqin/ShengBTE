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

! Compute the number of allowed three-phonon processes and their
! scattering amplitudes.

module processes
  use iso_fortran_env
  use data
  use config
  implicit none

  real(kind=8),parameter :: hbarp=hbar*1e22

contains

  subroutine Nprocesses(mm,N_plus,N_minus,energy,velocity,Nlist,List,IJK)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
    integer(kind=4),intent(out) :: N_plus,N_minus
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k ! Each element
                                                       ! of q or
                                                       ! qprime ranges
                                                       ! from 0 to
                                                       ! Ngrid(:)-1.
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll
    real(kind=8) :: sigma ! Gaussian broadening parameter, in THz.
    real(kind=8) :: omega,omegap,omegadp ! omega of the first, second
                                         ! and third phonon

    real(kind=8) :: dnrm2

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus=0
    N_minus=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,List(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    if (omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if (abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                      N_plus=N_plus+1
                   endif
                end if
             end do ! k
             !--------END absorption process-------------
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if (abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                      N_minus=N_minus+1
                   endif
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
    end if
  end subroutine Nprocesses

  subroutine Ind_plus(mm,N_plus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_plus,Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(out) :: Indof2ndPhonon_plus(N_plus),Indof3rdPhonon_plus(N_plus)
    real(kind=8),intent(out) :: Gamma_plus(N_plus)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_plus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,rr,ss,tt
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) :: omega,omegap,omegadp
    real(kind=8) :: realqprime(3),realqdprime(3)
    complex(kind=8) :: Vp,Vp0,prefactor ! The expression of Vp can be found in Natalio's book chapter.

    real(kind=8) :: dnrm2

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                      N_plus_count=N_plus_count+1
                      Indof2ndPhonon_plus(N_plus_count)=(index_N(qprime(1),qprime(2),qprime(3))-1)*Nbands+j
                      Indof3rdPhonon_plus(N_plus_count)=(index_N(qdprime(1),qdprime(2),qdprime(3))-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      !--------BEGIN calculation of Vp-----------
                      Vp=0.
                      do ll=1,Ntri
                         prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*&
                              masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
                              exp(iunit*dot_product(realqprime/ngrid,R_j(:,ll)))*&
                              exp(-iunit*dot_product(realqdprime/ngrid,R_k(:,ll)))
                         Vp0=0.
                         do rr=1,3
                            do ss=1,3
                               do tt=1,3
                                  Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                                       eigenvect(index_N(q(1),q(2),q(3)),i,tt+3*(Index_i(ll)-1))*&
                                       eigenvect(index_N(qprime(1),qprime(2),qprime(3)),j,ss+3*(Index_j(ll)-1))*&
                                       conjg(eigenvect(index_N(qdprime(1),qdprime(2),qdprime(3)),k,rr+3*(Index_k(ll)-1)))
                               end do
                            end do
                         end do
                         Vp=Vp+prefactor*Vp0
                      end do
                      !--------END calculation of Vp-------------
                      Gamma_plus(N_plus_count)=hbarp*pi/4.d0*(fBEprime-fBEdprime)*&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)*abs(Vp)**2
                      ! At this point, Gamma's units are
                      ! (1.d-34J*s)*(1.d12/s)^(-4)*1amu^(-3)*(ev/angstrom**3)^2,
                      ! that is, 5.60626442*1.d8 THz
                      Gamma_plus(N_plus_count)=Gamma_plus(N_plus_count)*5.60626442*1.d8/nptk ! THz

                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
    end if
    if(N_plus_count.ne.N_plus) write(error_unit,*) "Error: in Ind_plus, N_plus_count!=N_plus"
  end subroutine Ind_plus

  subroutine Ind_minus(mm,N_minus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_minus,Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(out) :: Indof2ndPhonon_minus(N_minus),Indof3rdPhonon_minus(N_minus)
    real(kind=8),intent(out) :: Gamma_minus(N_minus)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_minus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,rr,ss,tt
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) ::  omega,omegap,omegadp
    real(kind=8) :: realqprime(3),realqdprime(3)
    complex(kind=8) :: Vp,Vp0,prefactor ! The expression of Vp can be found in Natalio's book chapter.

    real(kind=8) :: dnrm2

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if (abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                      N_minus_count=N_minus_count+1
                      Indof2ndPhonon_minus(N_minus_count)=(index_N(qprime(1),qprime(2),qprime(3))-1)*Nbands+j
                      Indof3rdPhonon_minus(N_minus_count)=(index_N(qdprime(1),qdprime(2),qdprime(3))-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      !--------BEGIN calculation of Vp-----------
                      Vp=0.
                      do ll=1,Ntri
                         prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*&
                              masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
                              exp(-iunit*dot_product(realqprime/ngrid,R_j(:,ll)))*&
                              exp(-iunit*dot_product(realqdprime/ngrid,R_k(:,ll)))
                         Vp0=0.
                         do rr=1,3
                            do ss=1,3
                               do tt=1,3
                                  Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                                       eigenvect(index_N(q(1),q(2),q(3)),i,tt+3*(Index_i(ll)-1))*&
                                       conjg(eigenvect(index_N(qprime(1),qprime(2),qprime(3)),j,ss+3*(Index_j(ll)-1)))*&

                                       conjg(eigenvect(index_N(qdprime(1),qdprime(2),qdprime(3)),k,rr+3*(Index_k(ll)-1)))
                               end do
                            end do
                         end do
                         Vp=Vp+prefactor*Vp0
                      end do
                      !--------END calculation of Vp-------------
                      Gamma_minus(N_minus_count)=hbarp*pi/4.d0*(fBEprime+fBEdprime+1)*&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)*abs(Vp)**2
                      Gamma_minus(N_minus_count)=Gamma_minus(N_minus_count)*5.60626442*1.d8/nptk
                   end if
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
    end if
    if(N_minus_count.ne.N_minus) write(error_unit,*) "Error: in Ind_minus, N_minus_count!=N_minus"
  end subroutine Ind_minus

  subroutine D_plus(mm,N_plus,energy,velocity,Nlist,List,IJK,P3_plus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_plus
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(out) :: P3_plus(N_plus)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_plus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll
    real(kind=8) :: sigma
    real(kind=8) :: omega,omegap,omegadp

    real(kind=8) :: dnrm2

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(index_N(q(1),q(2),q(3)),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             omegap=energy(index_N(qprime(1),qprime(2),qprime(3)),j)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                omegadp=energy(index_N(qdprime(1),qdprime(2),qdprime(3)),k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(index_N(qprime(1),qprime(2),qprime(3)),j,:)-&
                        velocity(index_N(qdprime(1),qdprime(2),qdprime(3)),k,:))
                   if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                      N_plus_count=N_plus_count+1
                      P3_plus(N_plus_count)=exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
    end if
    if(N_plus_count.ne.N_plus) write(error_unit,*) "Error: in D_plus, N_plus_count!=N_plus"
  end subroutine D_plus

end module processes
