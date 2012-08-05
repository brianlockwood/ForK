!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file likelihood.f90
!> \brief Module containing the likelihood calculation for function-only and gradient-enhanced Kriging models. This likelihood formula treats all covariance parameters (length scales, magnitude and noise magnitudes) as hyperparameters that can be varied independently. The regression parameters are still determined based on an optimality condition.
module likelihood_mod
    implicit none
  interface likelihood
     module procedure likelihood_func, likelihood_grad
  end interface likelihood
contains
!>  \brief Subroutine for computing the likelihood formula for a function only Kriging model where noise level and covariance magnitude can be varied independently

!>  \detail The likelihood formula for a Gaussian process with a regression mean function is given by:
!> \f[
!>  \log p(y | X; \theta) = -\frac{1}{2} Y^{T} K^{-1} (Y - H^{T} \beta) - \frac{1}{2} \log \vert K \vert - \frac{1}{2} \log \vert A \vert - \frac{N-S}{2} \log 2\pi \f]
!>          where \f$N\f$ is the number of training points (denoted by the argument \a ntot) and \f$S\f$ is the number of terms in the regression (given by the argument \a stot). Here a vague prior on the regression parameters, \f$ \beta \f$ has been assumed (meaning that the distribution of these parameters is not known). The optimal regression parameters are given by:
!>   \f[ \beta = A^{-1} H Y \f]
!>   where the regression matrix \f$  A \f$ is defined as:
!>   \f[ A = H K^{-1} H^{T} \f].
!>   Both the covariance matrix and regression matrix are symmetric positive definite, allowing cholesky to be performed on both matrices. This cholesky factorization is used to invert these matrices as well as compute the determinant.

!>   For a cholesky-factorized matrix (applicable for the regression and covariance matrix), the logarithm of the determinant is given as:
!>    \f[ \log \vert C \vert = \log \vert L L^{T} \vert = 2 \sum_{k} L_{k,k} \f]

!>   The subroutine below stores the inverse of the covariance in place. In the case that the cholesky of the inverse of the matrix is available, the determinant of the original matrix is found by applying a negative sign.
!>    \f[ \log \vert C \vert = \log \vert (L^{-1} L^{-T})^{-1} \vert = - 2 \sum_{k} L^{-1}_{k,k} \f]

!>  Further details of the likelihood formula used in this subroutine can be found in Chapter 2 of <a href=http://www.gaussianprocess.org/gpml/><i> Gaussian Processes for Machine Learning </i> </a> by Rasmussen and Williams.

!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
!>         University of Wyoming
!> \date   May 17, 2012
!> \param(in) <b> ndim </b>: The dimension of the problem
!> \param(in) <b> ntot </b>: The number of training points
!> \param(in) <b> X </b>:  Location of the training points (size=[ndimxntot])
!> \param(in) <b>Y </b>:  Funcation values for the training points (size=[ntot])
!> \param(in) <b>stot </b>: Number of Terms in the regression 
!> \param(in) <b>H</b>: The collocation matrix for the regression (size=[stotxntot]) 
!> \param(in) <b> theta </b>: Length scale in each direction (size=[ndim])
!> \param(in) <b> sigma </b>: Covariance Magnitude
!> \param(in) <b> sigmaN </b>: Noise Magnitude (added to diagonal of covariance matrix)
!> \param(out) <b> logpy </b>: Natural Logarithm of the likelihood (\f$ p(y | X,\theta, \sigma, \sigma_{N})\f$

  subroutine likelihood_func(ndim,ntot,X,Y,stot,H,theta,sigma,sigmaN,logpy)
  use choleskymod
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot) 
  real(8), intent(in) :: theta(ndim),sigma,sigmaN

  real(8), intent(out) :: logpy

  integer i,j,k,l

  real(8) W(ntot),B(ntot),M(ntot)

  real(8) beta(stot)

  real(8) Kmat(ntot,ntot)
  real(8) A(stot,stot)

  real(8) K_det,A_det

  real(8) PI

  real(8) C(stot,ntot)
  real(8) Q(stot)

  PI=4.D0*ATAN(1.D0)

  call covarmatrix(ndim,ntot,X,theta,Kmat)
  
  !  Add in Noise and magnitude  
  Kmat=sigma**2*Kmat
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+sigmaN**2
  end do

!  Cholesky and Invert
  call invertsymmetric(ntot,Kmat)
 
!  Construct Regression Matrix
  call symmatrixmulttrans(ntot,Kmat,stot,H,C)
  call matrixmulttrans(stot,ntot,H,stot,C,A)

  call symmatvec(ntot,Kmat,Y,W)
  call matvec(stot,ntot,H,W,Q)
  Beta=Q
  call symmetricsolve(stot,A,Beta)

 !  Compute processed training data
  call matvectrans(stot,ntot,C,Beta,M)

  ! Calculate Determinants
  K_det=0.D0
  do i=1,ntot
     K_det=K_det+log(Kmat(i,i))
  end do
  K_det=-2.*K_det

  A_det=0.D0
  do i=1,stot
     A_det=A_det+log(A(i,i))
  end do
  A_det=-2.*A_det

  logpy=0.D0
  do i=1,ntot
     logpy=logpy+Y(i)*(W(i)-M(i))
  end do
  logpy=-0.5*logpy-0.5*K_det-0.5*A_det-0.5*(ntot-stot)*log(2.D0*PI)
 
  return
end subroutine likelihood_func
!> \brief Subroutine used for computing the likelihood formula for a gradient-enhanced Kriging model where funtion noise level, derivative noise level and covariance magnitude can be varied independently.
subroutine likelihood_grad(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,sigma,sigmaN,sigmaNG,logpy)
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
  real(8), intent(in) :: theta(ndim),sigma,sigmaN,sigmaNG

  real(8), intent(out) :: logpy

  real(8) Z(ntot+gtot)

  integer i,j,k,l

  real(8) W(ntot+gtot),Q(stot),C(stot,ntot+gtot),M(ntot+gtot)

  real(8) beta(stot)

  real(8) Kmat(ntot+gtot,ntot+gtot)
  real(8) A(stot,stot)
  
  real(8) K_det,A_det

  real(8) PI

  real(8) maxEig, minEig

  PI=4.D0*ATAN(1.D0)

  !  Training Data
  do i=1,ntot
     Z(i)=Y(i)
  end do
  do i=ntot+1,ntot+gtot
     Z(i)=dY(i-ntot)
  end do

  call covarmatrix(ndim,ntot,gtot,pts,dims,X,theta,Kmat)
  
  !  Add in Noise and magnitude  
  Kmat=sigma**2*Kmat
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+sigmaN**2
  end do
  do i=1,gtot
     Kmat(ntot+i,ntot+i)=Kmat(ntot+i,ntot+i)+sigmaNG**2
  end do

!  Cholesky and Invert
  call invertsymmetric(ntot+gtot,Kmat)
 
!  Construct Regression Matrix
  call symmatrixmulttrans(ntot+gtot,Kmat,stot,H,C)
  call matrixmulttrans(stot,ntot+gtot,H,stot,C,A)

!  Cholesky Regression Matrix and solve for Beta

  call symmatvec(ntot+gtot,Kmat,Z,W)  
  call matvec(stot,ntot+gtot,H,W,Q)

  Beta=Q
  call symmetricsolve(stot,A,Beta)

 !  Compute processed training data
  call matvectrans(stot,ntot+gtot,C,Beta,M)
  
  K_det=0.D0
  do i=1,ntot+gtot
     K_det=K_det+log(Kmat(i,i))
  end do
  K_det=-2.*K_det

  A_det=0.D0
  do i=1,stot
     A_det=A_det+log(A(i,i))
  end do
  A_det=2.*A_det

  logpy=0.D0
  do i=1,ntot+gtot
     logpy=logpy+Z(i)*(W(i)-M(i))
  end do
  logpy=-0.5*logpy-0.5*K_det-0.5*A_det-0.5*(ntot+gtot-stot)*log(2.D0*PI)
  
  return
end subroutine likelihood_grad
end module likelihood_mod
