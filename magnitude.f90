!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file magnitude.f90
!> \brief This module contains the subroutine that can be used to calculate the optimal covariance magnitude for a function-only and a gradient-enhanced Kriging model. 
module magnitude_mod
    interface magnitude
     module procedure magnitude_func,magnitude_grad
  end interface magnitude
contains

!> \brief Subroutine for computing the optimal magnitude associated with a function-only Kriging model. It is used when the length scale and noise ratio are determined and the magnitude is needed. This code is repeated from the <a href=likelihood_mle_8f90.html>likelihood_mle</a> subroutine. If the magnitude formula changes, both this and <a href=likelihood_mle_8f90.html>likelihood_mle</a> need to be modified.

!> \detail The optimal magnitude given lenght scales and noise ratio is given as: 
!>\f[\sigma^{2} = \frac{(Y-H^{T} \beta)^{T} \hat{K}^{-1} (Y-H^{T} \beta )}{N} \f]
!> where \f$\hat{K}\f$ is the unit magnitude covariance matrix, the elements of which are given as:
!> \f[ \hat{K}_{i,j} = k(\vec{x}_{i}, \vec{x}_{j},\theta) + \Gamma \delta_{i,j}\f]

!> This subroutine is NOT used within the <a href=likelihood_mle_8f90.html>likelihood_mle</a> subroutine and is used only when the length scale is known and the magnitude is needed.

!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
!>         University of Wyoming
!> \date   May 18, 2012
!> \param(in) <b> ndim </b>: The dimension of the problem
!> \param(in) <b> ntot </b>: The number of training points
!> \param(in) <b> X </b>:  Location of the training points (size=[ndimxntot])
!> \param(in) <b>Y </b>:  Funcation values for the training points (size=[ntot])
!> \param(in) <b>stot </b>: Number of Terms in the regression 
!> \param(in) <b>H</b>: The collocation matrix for the regression (size=[stotxntot]) 
!> \param(in) <b> theta </b>: Length scale in each direction (size=[ndim])
!> \param(in) <b> gamma </b>: Ratio of Noise Magnitude to covariance matrix (added to diagonal of covariance matrix with magnitude of 1)
!> \param(out) <b> sigma </b>: Optimal covariance magnitude (\f$\sigma\f$)

subroutine magnitude_func(ndim,ntot,X,Y,stot,H,theta,gamma,sigma)
  use choleskymod
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot) 
  real(8), intent(in) :: theta(ndim),gamma

  real(8), intent(out) :: sigma

  integer i,j,k,l

  real(8) W(ntot),B(ntot),Col(ntot)

  real(8) beta(stot)

  real(8) Kmat(ntot,ntot)
  real(8) Kmati(ntot,ntot)
  real(8) A(stot,stot)

  call covarmatrix(ndim,ntot,X,theta,Kmat)
  
  !  Add in Noise and magnitude  
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+gamma**2
  end do

!  Cholesky and Invert
  call cholesky(ntot,Kmat)

!  Invert Covariance Matrix
  call invertcholesky(ntot,Kmat,Kmati)
 
!  Construct Regression Matrix
  do i=1,stot
     do j=1,stot
        A(i,j)=0.D0
        do k=1,ntot
           do l=1,ntot
              A(i,j)=A(i,j)+H(i,k)*Kmati(k,l)*H(j,l)
           end do
        end do
     end do
  end do

!  Cholesky Regression Matrix and solve for Beta
  call cholesky(stot,A)

  do i=1,ntot
     W(i)=0.D0
     do j=1,ntot
        W(i)=W(i)+Kmati(i,j)*Y(j)
     end do
  end do

  do i=1,stot
     Col(i)=0.D0
     do j=1,ntot
        Col(i)=Col(i)+H(i,j)*W(j)
     end do
  end do
  
  call choleskysolve(stot,A,col,Beta)

 !  Compute processed training data
  do k=1,ntot
     B(k)=0.D0
     do i=1,ntot
        do j=1,stot
           B(k)=B(k)+Kmati(k,i)*H(j,i)*Beta(j)
        end do
     end do
  end do
  
  sigma=0.D0
  do i=1,ntot
     sigma=sigma+Y(i)*(W(i)-B(i))
  end do
  sigma=sqrt(sigma/(ntot-stot))
  
  return
end subroutine magnitude_func
!> \brief Subroutine for computing the optimal magnitude associated with a gradient-enhanced Kriging model. It is used when the length scale and noise ratio are determined and the magnitude is needed. This code is repeated from the <a href=likelihood_mle_8f90.html>likelihood_mle</a> subroutine. If the magnitude formula changes, both this and <a href=likelihood_mle_8f90.html>likelihood_mle</a> need to be modified.
subroutine magnitude_grad(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,gamma,gammaG,sigma)
  use choleskymod
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: gtot
  integer, intent(in) :: pts(gtot), dims(gtot)
  real(8), intent(in) :: dY(gtot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot+gtot) 
  real(8), intent(in) :: theta(ndim),gamma,gammaG

  real(8), intent(out) :: sigma

  real(8) Z(ntot+gtot)

  integer i,j,k,l

  real(8) W(ntot+gtot),B(ntot+gtot),Col(ntot+gtot)

  real(8) beta(stot)

  real(8) Kmat(ntot+gtot,ntot+gtot)
  real(8) Kmati(ntot+gtot,ntot+gtot)
  real(8) A(stot,stot)

  !  Training Data
  do i=1,ntot
     Z(i)=Y(i)
  end do
  do i=ntot+1,ntot+gtot
     Z(i)=dY(i-ntot)
  end do

  call covarmatrix(ndim,ntot,gtot,pts,dims,X,theta,Kmat)
  
  !  Add in Noise and magnitude  
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+gamma**2
  end do
  do i=1,gtot
     Kmat(ntot+i,ntot+i)=Kmat(ntot+i,ntot+i)+gammaG**2
  end do

!  Cholesky and Invert
  call cholesky(ntot+gtot,Kmat)

!  Invert Covariance Matrix
  call invertcholesky(ntot+gtot,Kmat,Kmati)
 
!  Construct Regression Matrix
  do i=1,stot
     do j=1,stot
        A(i,j)=0.D0
        do k=1,ntot+gtot
           do l=1,ntot+gtot
              A(i,j)=A(i,j)+H(i,k)*Kmati(k,l)*H(j,l)
           end do
        end do
     end do
  end do

!  Cholesky Regression Matrix and solve for Beta
  call cholesky(stot,A)

  do i=1,ntot+gtot
     W(i)=0.D0
     do j=1,ntot+gtot
        W(i)=W(i)+Kmati(i,j)*Z(j)
     end do
  end do

  do i=1,stot
     Col(i)=0.D0
     do j=1,ntot+gtot
        Col(i)=Col(i)+H(i,j)*W(j)
     end do
  end do
  
  call choleskysolve(stot,A,col,Beta)

 !  Compute processed training data
  do k=1,ntot+gtot
     B(k)=0.D0
     do i=1,ntot+gtot
        do j=1,stot
           B(k)=B(k)+Kmati(k,i)*H(j,i)*Beta(j)
        end do
     end do
  end do
  
  sigma=0.D0
  do i=1,ntot+gtot
     sigma=sigma+Z(i)*(W(i)-B(i))
  end do
  sigma=sqrt(sigma/(ntot+gtot-stot))
  
  return
end subroutine magnitude_grad
end module magnitude_mod
