!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file kriginggradvariance.f90
!> \brief Subroutine to compute the variance of derivative values from Kriging Surface (requires  <a href=buildkriging_8f90.html>buildkriging</a> to be called first)

!> \detail In order to provide a confidence measure for the derivative predictions associated with the Kriging surface, the variance of the derivative prediction can be calculated based on the covariance structure of the Kriging model. To predict the variance of derivative values, the variance between derivative and function values as well as derivative and derivative values must be determined. Due to the linearity of differentiation, these covariance values are determined by differentiation of the covariance function. As seen previously, the covariance between a derivative value and function value is given as:

!> \f[ cov(\frac{\partial y_{i}}{\partial x_{k}}, y_{j}) = \frac{\partial}{\partial x_{k}} k(\vec{x}_{i},\vec{x}_{j}) \f]

!> while the correlation between derivative values is given by:

!> \f[ cov(\frac{\partial y_{i}}{\partial x_{k}}, \frac{\partial y_{j}}{\partial x_{l}}) = \frac{\partial^{2}}{\partial x_{k} \partial x_{l}} k(\vec{x}_{i},\vec{x}_{j}) \f].

!> Using these definitions, the variance between points can be calculated as:

!>\f[
!> V\left[\frac{\partial y(\vec{x}_{*})}{\partial x_{k}} \right]= \frac{\partial^{2}}{\partial x_{k}^{2}} k(\vec{x}_{*},\vec{x}_{*})-l_{*,k}^T K^{-1} l_{*,k}+ S_{k}(\vec{x}_*) A^{-1} S_{k}(\vec{x}_*)^{T} \f]

!> where \f$l_{*,k}^{T}\f$ is a vector with elements representing the covariance between the derivative value with respect to \f$x_{k}\f$ and the training point function values, given as: 

!> \f[  l_{*,k}^{T} = \left[\frac{\partial}{\partial x_{k}} k(\vec{x}_{*},\vec{X}_{1}) \dots \frac{\partial}{\partial x_{k}} k(\vec{x}_{*},\vec{X}_{N}) \right] \f]

!> and \f$S_{k}(\vec{x}_{*})\f$ represents the derivative of the mean function defined as:
!> \f[
!>     S_{k}(\vec{x}_{*}) = g_{k}(\vec{x}_{*}) - l_{*,k}^T K^{-1} H \f]
!>  where \f$g_{k}\f$ is a vector containing the derivative of the regression basis functions with respect to \f$x_{k}\f$ evaluated at the test point \f$\vec{x}_{*}\f$ and \f$G_{k}\f$ is a matrix containing the derivative of the regression basis evaluated at the training points \f$\vec{X}\f$. Using this variance prediction, confidence bounds can be placed on the derivative predictions from the Kriging Model. <br>

!> Because the variance of the derivative prediction requires the covariance function to be twice differentiable, only the Matern functions with \f$ \nu \geq 3/2 \f$ can be used.


!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 2, 2012
!> \param(in) <b>ndim </b>: The dimension of the problem
!> \param(in) <b>ntot </b>: The number of Training points
!> \param(in) <b>X </b>:  The location of the training points (size=[ndimxntot])
!> \param(in) <b>stot </b>: Number of Terms in the regression 
!> \param(in) <b>H</b>: The collocation matrix for the regression (size=[stotxntot]) 
!> \param(in) <b> beta</b>: Regression coefficients based on the optimal estimate for the Kriging model (size=[stot]) <br>
!>                          Supplied by <a href=buildkriging_8f90.html>buildkriging</a> subroutine
!> \param(in) <b>hyper</b>: Hyperparameters for the Kriging Model (size=[ndim+2]) <br>
!>                          Supplied by <a href=buildkriging_8f90.html>buildkriging</a> subroutine
!> \param(in) <b>mtot </b>: The number of test points, the places where function prediction are desired 
!> \param(in) <b>Xm </b>:  The location of the test points (size=[ndimxmtot])
!> \param(in) <b>gtotm</b>: The total number of derivatives to be predicted (ndim*mtot if the gradient at every test point is desired)
!> \param(in) <b> ptsm</b>: A vector identifying the test point where a particular derivative is specified (size=[gtotm] with values ranging from 1 to mtot)
!> \param(in) <b> dimsm</b>: A vector identifying the dimension the derivative is taken with respect to (size=[gtotm] with values ranging from 1 to ndim)
!> \param(in) <b>Gm </b>:  The derivative of the regression basis evaluated at the test points (size=[stotxgtotm])
!> \param(out) <b>Var</b>:  The variance of the derivative values (size=[gtotm]) <br>
!>                          Because the variance prediction requires the covariance matrix of the training data, it is more expensive than the mean prediction <br>
!>                          If possible, predict the variance associated with multiple points using a single subroutine call
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul>
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. When the variance of the derivative predictions is desired, a twice differentiable covariance function is required. This requires  \f$\nu \geq 3/2 \f$.
!>              Must supply the same covariance flag as used in <a href=buildkriging_8f90.html>buildkriging</a>.

subroutine kriginggradvariance(ndim,ntot,X,stot,H,hyper,mtot,Xm,gtotm,ptsm,dimsm,Gm,Var,covarflagi)
  use choleskymod
  use covars
  use covarmatrix_mod
  use covarmatrix_grad_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot)
  real(8), intent(in) :: hyper(ndim+2)

  integer, intent(in) :: mtot
  real(8), intent(in) :: Xm(ndim,mtot)
  integer, intent(in) :: gtotm
  integer, intent(in) :: ptsm(gtotm), dimsm(gtotm)
  real(8), intent(in) :: Gm(stot,gtotm)

  real(8), intent(out) :: Var(gtotm)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN
  real(8) Kmat(ntot,gtotm)
  real(8) Kcov(ntot,ntot)
 
  integer i,j,k,l

  real(8) S(stot,gtotm)

  real(8) A(stot,stot)

  real(8) Hs(ndim,ndim)

  real(8) W(ntot), Q(stot)

  real(8) C(stot,ntot)
  
  covarflag=covarflagi
  
!  Find Hyper Parameters
  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)

  call covarmatrix(ndim,ntot,X,theta,Kcov)

  Kcov=sigma**2*Kcov
  do i=1,ntot
     Kcov(i,i)=Kcov(i,i)+sigmaN**2
  end do

!  Cholesky and Invert
  call invertsymmetric(ntot,Kcov)

!  Calculate Covariance Matrix
  call covarmatrix_grad(ndim,ntot,X,mtot,gtotm,ptsm,dimsm,Xm,theta,Kmat)

!  Add magnitude  
  Kmat=sigma**2*Kmat

!  Construct Regression Matrix
  call symmatrixmulttrans(ntot,Kcov,stot,H,C)
  call matrixmulttrans(stot,ntot,H,stot,C,A)

  call invertsymmetric(stot,A)

!  Regression
  do j=1,stot
     do i=1,gtotm
        S(j,i)=Gm(j,i)
        do k=1,ntot
           S(j,i)=S(j,i)-Kmat(k,i)*C(j,k)
        end do
     end do
  end do

  !  Give Prediction
  do i=1,gtotm
     j=ptsm(i)
     k=dimsm(i)
     call d2covarfunc(ndim,Xm(:,j),Xm(:,j),theta(:),Hs)
     Var(i)=Hs(k,k)*sigma**2
     
     call symmatvec(ntot,Kcov,Kmat(:,i),W)

     do j=1,ntot
        Var(i)=Var(i)-Kmat(j,i)*W(j)
     end do

     call symmatvec(stot,A,S(:,i),Q)

     do j=1,stot
        Var(i)=Var(i)+S(j,i)*Q(j)
     end do

  end do

  return
end subroutine kriginggradvariance
