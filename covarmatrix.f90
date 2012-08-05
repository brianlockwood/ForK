!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file covarmatrix.f90
!> \brief This module contains the subroutines used to compute the covariance matrix. For efficiency, different subroutines are used for the different covariance matrices needed in the code. Although the naming may be confusing, this module contains the subroutines required to calculate the block covariance matrix needed for a gradient-enhanced Kriging model. 

module covarmatrix_mod
    implicit none
  interface covarmatrix
     module procedure covarmatrix_training, covarmatrix_training_grad, covarmatrix_gen, covarmatrix_gen_grad
  end interface covarmatrix
contains
!> \brief Subroutine to create the covariance matrix between training points for a function-only Kriging model
subroutine covarmatrix_training(ndim,ntot,X,theta,Kmat)
  use covars
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  
  real(8), intent(in) :: theta(ndim)

  real(8), intent(out) :: Kmat(ntot,ntot)

  integer i,j,k,l
  real(8) kr
   
  do i=1,ntot
     do j=1,i
        call covarfunc(ndim,X(:,i),X(:,j),theta,kr)
        Kmat(i,j)=kr
     end do
  end do

  return
end subroutine covarmatrix_training
!> \brief Subroutine to create covariance matrix between training points for gradient-enhanced Kriging model. This matrix is a block matrix where the diagonals represent covariance between function values and derivative values and the off-diagonals represent the covaraince between function and derivative values.
subroutine covarmatrix_training_grad(ndim,ntot,gtot,pts,dims,X,theta,Kmat)
  use covars
  implicit none
  integer, intent(in) :: ndim, ntot,gtot
  integer, intent(in) :: pts(gtot),dims(gtot)
  real(8), intent(in) :: X(ndim,ntot)
  
  real(8), intent(in) :: theta(ndim)

  real(8), intent(out) :: Kmat(ntot+gtot,ntot+gtot)

  integer i,j,k,l
  real(8) kr,dkr(ndim),d2kr(ndim,ndim)

  integer nder(ntot)
  integer npoint(ntot+1)

  nder=0
  do i=1,gtot
     nder(pts(i))=nder(pts(i))+1
  end do
  
  npoint(1)=1
  do i=1,ntot
     npoint(i+1)=npoint(i)+nder(i)
  end do
    
  do i=1,ntot
     do j=1,i
        call covarfunc(ndim,X(:,i),X(:,j),theta,kr)
        Kmat(i,j)=kr
     end do
  end do
  
  do i=1,ntot
     do j=1,ntot
        call dcovarfunc(ndim,X(:,i),X(:,j),theta,dkr)
        do k=npoint(i),npoint(i+1)-1
           l=dims(k)
           Kmat(ntot+k,j)=dkr(l)
        end do
     end do
  end do

  do i=1,ntot
     do j=1,i-1
        call d2covarfunc(ndim,X(:,i),X(:,j),theta,d2kr)
        do k=npoint(i),npoint(i+1)-1
           do l=npoint(j),npoint(j+1)-1
              Kmat(ntot+k,ntot+l)=d2kr(dims(k),dims(l))
           end do
        end do
     end do
     call d2covarfunc(ndim,X(:,i),X(:,i),theta,d2kr)
     do k=npoint(i),npoint(i+1)-1
        do l=npoint(i),k
           Kmat(ntot+k,ntot+l)=d2kr(dims(k),dims(l))
        end do
     end do
  end do

  return
end subroutine covarmatrix_training_grad
!> \brief Builds a general covariance matrix between two sets of points for a function-only Kriging model. It is used to construct the covariance between training points and test points. 
subroutine covarmatrix_gen(ndim,ntot,X,mtot,Xm,theta,Kmat)
  use covars
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  
  integer, intent(in) :: mtot
  real(8), intent(in) :: Xm(ndim,mtot)

  real(8), intent(in) :: theta(ndim)

  real(8), intent(out) :: Kmat(ntot,mtot)

  integer i,j,k,l
  real(8) kr

   
  do i=1,ntot
     do j=1,mtot
        call covarfunc(ndim,X(:,i),Xm(:,j),theta,kr)
        Kmat(i,j)=kr
     end do
  end do
  
  return
end subroutine covarmatrix_gen
!> \brief Builds a general covariance matrix between two sets of points for a gradient-enhanced Kriging model. Includes the covariance between function and derivative values. 
subroutine covarmatrix_gen_grad(ndim,ntot,gtot,pts,dims,X,mtot,htot,ptsm,dimsm,Xm,theta,Kmat)
  use covars
  implicit none
  integer, intent(in) :: ndim, ntot,gtot
  integer, intent(in) :: pts(gtot),dims(gtot)
  real(8), intent(in) :: X(ndim,ntot)
  
  integer, intent(in) :: mtot,htot
  integer, intent(in) :: ptsm(htot),dimsm(htot)
  real(8), intent(in) :: Xm(ndim,mtot)

  real(8), intent(in) :: theta(ndim)

  real(8), intent(out) :: Kmat(ntot+gtot,mtot+htot)

  integer i,j,k,l
  real(8) kr,dkr(ndim),d2kr(ndim,ndim)

  integer nder(ntot), nderm(mtot)
  integer npoint(ntot+1),npointm(mtot+1)

  nder=0
  do i=1,gtot
     nder(pts(i))=nder(pts(i))+1
  end do
  nderm=0
  do i=1,htot
     nderm(ptsm(i))=nderm(ptsm(i))+1
  end do
  
  npoint(1)=1
  do i=1,ntot
     npoint(i+1)=npoint(i)+nder(i)
  end do
  npointm(1)=1
  do i=1,mtot
     npointm(i+1)=npointm(i)+nderm(i)
  end do
    
  do i=1,ntot
     do j=1,mtot
        call covarfunc(ndim,X(:,i),Xm(:,j),theta,kr)
        Kmat(i,j)=kr
     end do
  end do
  
  do i=1,ntot
     do j=1,mtot
        call dcovarfunc(ndim,X(:,i),Xm(:,j),theta,dkr)
        do k=npointm(j),npointm(j+1)-1
           l=dimsm(k)
           Kmat(i,mtot+k)=-dkr(l)
        end do
     end do
  end do

  do i=1,ntot
     do j=1,mtot
        call dcovarfunc(ndim,X(:,i),Xm(:,j),theta,dkr)
        do k=npoint(i),npoint(i+1)-1
           l=dims(k)
           Kmat(ntot+k,j)=dkr(l)
        end do
     end do
  end do

  do i=1,ntot
     do j=1,mtot
        call d2covarfunc(ndim,X(:,i),Xm(:,j),theta,d2kr)
        do k=npoint(i),npoint(i+1)-1
           do l=npointm(j),npointm(j+1)-1
              Kmat(ntot+k,mtot+l)=d2kr(dims(k),dimsm(l))
           end do
        end do
     end do
  end do

  return
end subroutine covarmatrix_gen_grad
end module covarmatrix_mod
