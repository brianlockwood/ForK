!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file likelihood_mle.f90
!> \brief Module containing the likelihood calculation for function-only and gradient-enhanced Kriging models. This likelihood formula treats the length scale in each dimension of the covariance function as hyperparameters that can be varied independently. The regression parameters and covariance magnitude are determined based on an optimality condition. The ratio of the noise level to the covariance magnitude is fixed value and is introduced to ensure proper conditioning of the covariance matrix. 
module likelihood_mle_mod
  implicit none
  interface likelihood_mle
     module procedure likelihood_mle_func, likelihood_mle_grad
  end interface likelihood_mle
contains
!>  \brief Subroutine for computing the likelihood formula for a function only Kriging model using an explicit mean function and covariance magnitude determined based on an optimality condition.

!>  \detail The likelihood formula for a Gaussian process with explicit mean function is given by:
!>          \f[ \log p(y | X; \theta) = -\frac{1}{2} (Y - H^{T} \beta)^{T} K^{-1} (Y - H^{T} \beta) - \frac{1}{2} \log \vert K \vert - \frac{N}{2} \log 2\pi \f]
!>          where \f$N\f$ is the number of training points (denoted by the argument \a ntot)
!>          The optimal regression parameters are determined by differentiation of the above equation and solving for the parameters \f$\beta\f$ which set this derivative to zero. Based on this process, the optimal regression parameters are given by:
!>          \f[ \beta = (H K^{-1} H^{T})^{-1} H K^{-1} Y \f]
!>          The optimal covariance magnitude can be derived in a similar manner as the regression parameters (differentiating the likelihood formula and determining the magnitude that sets this derivative to zero). To determine this optimal magnitude, the covariance matrix should be decomposed into the magnitude and a covariance matrix with magnitude of 1.
!>          \f[ K = \sigma^{2} \hat{K} \f]
!>          The magnitude 1 covariance matrix includes the noise term added to the diagonal. The elements of the covariance matrix are given as:
!>          \f[ \hat{K}_{i,j} = k(\vec{x}_{i}, \vec{x}_{j},\theta) + \Gamma \delta_{i,j}\f]
!>          where \f$\Gamma\f$ is the parameter specified in the argument \a sigmaN and is the assumed level of noise in the magnitude 1 covariance matrix. The parameter is hard coded and is used to ensure \f$ \hat{K} \f$ is invertable.
!>          Using this decomposition, the magnitude 1 covariance matrix can be used directly to determine the optimal regression parameters,
!>          \f[ \beta = (H \hat{K}^{-1} H^{T})^{-1} H \hat{K}^{-1} Y \f]
!>          as the inverse of the covariance matrix can be written as:
!>          \f[ K^{-1} = \frac{1}{\sigma^{2}} \hat{K}^{-1} \f]
!>          Using these definitions, the optimal covariance magnitude is given as:
!>          \f[\sigma^{2} = \frac{(Y-H^{T} \beta)^{T} \hat{K}^{-1} (Y-H^{T} \beta )}{N} \f].
!>          Using this optimal covariance, the final likelihood formula can be written as:
!>          \f[ \log p(y | X; \theta) = - \frac{N}{2} \log \sigma^{2} -\frac{1}{2} \log \vert \hat{K} \vert - \frac{N}{2} - \frac{N}{2} \log  2\pi \f]

!>          This final formula is a function of the length scales, \f$\theta\f$, and the noise level ratio \f$\Gamma\f$. The optimization is performed only over the length scales and the noise level ratio is fixed throughout the optimization. 

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
!> \param(in) <b> sigmaN </b>: Ratio of Noise Magnitude to covariance matrix (added to diagonal of covariance matrix with magnitude of 1)
!> \param(out) <b> logpy </b>: Natural Logarithm of the likelihood (\f$ p(y | X,\theta, \sigma, \sigma_{N})\f$
  subroutine likelihood_mle_func(ndim,ntot,X,Y,stot,H,theta,sigmaN,logpy)
  use choleskymod
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot) 
  real(8), intent(in) :: theta(ndim),sigmaN

  real(8), intent(out) :: logpy

  integer i,j,k,l

  real(8) W(ntot),B(ntot),M(ntot)

  real(8) beta(stot)

  real(8) Kmat(ntot,ntot)
  real(8) A(stot,stot)

  real(8) K_det,A_det

  real(8) PI

  real(8) sigma

  real(8) C(stot,ntot)
  
  real(8) Q(stot)

  PI=4.D0*ATAN(1.D0)

  call covarmatrix(ndim,ntot,X,theta,Kmat)
    
  !  Add in Noise and magnitude  
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+sigmaN**2
  end do

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

 !  Compute processed training data
  call matvectrans(stot,ntot,H,Beta,B)

  sigma=0.D0
  do i=1,ntot
     sigma=sigma+(Y(i)-B(i))*(W(i)-M(i))
  end do
  sigma=(sigma/(ntot))

  K_det=0.D0
  do i=1,ntot
     K_det=K_det+log(Kmat(i,i))
  end do
  K_det=-2.*K_det

  logpy=-(ntot)*log(sigma)-K_det
    return
  end subroutine likelihood_mle_func
!>  \brief Subroutine for computing the likelihood formula for a gradient-enhanced Kriging model using an explicit mean function and covariance magnitude determined based on an optimality condition.
  subroutine likelihood_mle_grad(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,sigmaN,sigmaNG,logpy)
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
  real(8), intent(in) :: theta(ndim),sigmaN,sigmaNG

  real(8), intent(out) :: logpy

  real(8) Z(ntot+gtot)

  integer i,j,k,l

  real(8) W(ntot+gtot),Q(stot),M(ntot+gtot),C(stot,ntot+gtot), B(ntot+gtot)

  real(8) beta(stot)

  real(8) Kmat(ntot+gtot,ntot+gtot)
  real(8) A(stot,stot)

  real(8) K_det,A_det

  real(8) PI

  real(8) sigma

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
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+sigmaN**2
  end do
  do i=1,gtot
     Kmat(ntot+i,ntot+i)=Kmat(ntot+i,ntot+i)+sigmaNG**2
  end do

  call invertsymmetric(ntot+gtot,Kmat)
 
!  Construct Regression Matrix
  call symmatrixmulttrans(ntot+gtot,Kmat,stot,H,C)
  call matrixmulttrans(stot,ntot+gtot,H,stot,C,A)

  call symmatvec(ntot+gtot,Kmat,Z,W)
  call matvec(stot,ntot+gtot,H,W,Q)
  Beta=Q
  call symmetricsolve(stot,A,Beta)

 !  Compute processed training data
  call matvectrans(stot,ntot+gtot,C,Beta,M)

 !  Compute processed training data
  call matvectrans(stot,ntot+gtot,H,Beta,B)

  sigma=0.D0
  do i=1,ntot+gtot
     sigma=sigma+(Z(i)-B(i))*(W(i)-M(i))
  end do
  sigma=(sigma/(ntot+gtot))
  
  K_det=0.D0
  do i=1,ntot+gtot
     K_det=K_det+log(Kmat(i,i))
  end do
  K_det=-2.*K_det

  logpy=-(ntot+gtot)*log(sigma)-K_det


    return
  end subroutine likelihood_mle_grad
end module likelihood_mle_mod
