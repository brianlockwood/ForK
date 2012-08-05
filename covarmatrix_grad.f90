!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file covarmatrix_grad.f90
!> \brief This module contains the subroutines required to compute the covariance matrix between derivative values and function values. For a function-only Kriging model, the matrix consists of the covariance between derivative and function values with elements corresponding to the derivative of the covariance matrix. For a gradient-enhanced Kriging model, the matrix consists of the covariance between derivative and function values as well as the covariance with other derivative values. These subroutines are required when the derivative predictions are made from the Kriging surface. This module does NOT contain the subroutine required to compute the covariance matrix required in the construction of a gradient-enhanced Kriging model (that subroutine is in the module <a href=covarmatrix_8f90.html>covarmatrix</a>).

module covarmatrix_grad_mod
  implicit none
  interface covarmatrix_grad
     module procedure covarmatrix_grad_gek, covarmatrix_grad_func
  end interface covarmatrix_grad
contains
!> \brief Subroutine to create the covariance matrix between a set of function values and a set of derivative values (typically located at the test points).
  subroutine covarmatrix_grad_func(ndim,ntot,X,mtot,htot,ptsm,dimsm,Xm,theta,Kmat)
    use covars
    implicit none
    integer, intent(in) :: ndim, ntot
    real(8), intent(in) :: X(ndim,ntot)

    integer, intent(in) :: mtot,htot
    integer, intent(in) :: ptsm(htot),dimsm(htot)
    real(8), intent(in) :: Xm(ndim,mtot)

    real(8), intent(in) :: theta(ndim)

    real(8), intent(out) :: Kmat(ntot,htot)

    integer i,j,k,l
    real(8) dkr(ndim)

    integer nderm(mtot)
    integer npointm(mtot+1)

    nderm=0
    do i=1,htot
       nderm(ptsm(i))=nderm(ptsm(i))+1
    end do

    npointm(1)=1
    do i=1,mtot
       npointm(i+1)=npointm(i)+nderm(i)
    end do

    do i=1,ntot
       do j=1,mtot
          call dcovarfunc(ndim,X(:,i),Xm(:,j),theta,dkr)
          do k=npointm(j),npointm(j+1)-1
             l=dimsm(k)
             Kmat(i,k)=-dkr(l)
          end do
       end do
    end do

    return
  end subroutine covarmatrix_grad_func
!> \brief Subroutine to compute the covariance between a set of training data and a set of derivative values for a gradient-enhanced model. For a gradient-enhanced model, the set of training data includes function and derivative values evaluated at the training points. 
  subroutine covarmatrix_grad_gek(ndim,ntot,gtot,pts,dims,X,mtot,htot,ptsm,dimsm,Xm,theta,Kmat)
    use covars
    implicit none
    integer, intent(in) :: ndim, ntot,gtot
    integer, intent(in) :: pts(gtot),dims(gtot)
    real(8), intent(in) :: X(ndim,ntot)

    integer, intent(in) :: mtot,htot
    integer, intent(in) :: ptsm(htot),dimsm(htot)
    real(8), intent(in) :: Xm(ndim,mtot)

    real(8), intent(in) :: theta(ndim)

    real(8), intent(out) :: Kmat(ntot+gtot,htot)

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
          call dcovarfunc(ndim,X(:,i),Xm(:,j),theta,dkr)
          do k=npointm(j),npointm(j+1)-1
             l=dimsm(k)
             Kmat(i,k)=-dkr(l)
          end do
       end do
    end do

    do i=1,ntot
       do j=1,mtot
          call d2covarfunc(ndim,X(:,i),Xm(:,j),theta,d2kr)
          do k=npoint(i),npoint(i+1)-1
             do l=npointm(j),npointm(j+1)-1
                Kmat(ntot+k,l)=d2kr(dims(k),dimsm(l))
             end do
          end do
       end do
    end do

    return
  end subroutine covarmatrix_grad_gek
end module covarmatrix_grad_mod
