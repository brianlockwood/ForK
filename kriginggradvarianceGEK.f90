!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file kriginggradvarianceGEK.f90
!> \brief Subroutine to predict the variance of derivative values from Gradient-Enhanced Kriging Surface (requires  <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> to be called first)

!> \detail This subroutine calculates the variance associated with derivative predictions from a gradient-enhanced Kriging model. The details of this variance prediction are given in <a href=kriginggradvariance_8f90.html>kriginggradvariance</a>. The form of these variance predictions are the same for a gradient-enhanced and function-only model. The variance predictions for derivative values are given as:

!>\f[
!> V\left[\frac{\partial y(\vec{x}_{*})}{\partial x_{k}} \right]= \frac{\partial^{2}}{\partial x_{k}^{2}} k(\vec{x}_{*},\vec{x}_{*})-\underline{l}_{*,k}^T \underline{K}^{-1} \underline{l}_{*,k}+ \underline{S}_{k}(\vec{x}_*) \underline{A}^{-1} \underline{S}_{k}(\vec{x}_*)^{T} \f]

!> See the subroutines <a href = krigingfuncpredictGEK_8f90.html> krigingfuncpredictGEK</a>, <a href = krigingfuncvarianceGEK_8f90.html> krigingfuncvarianceGEK</a> and <a href = kriginggradpredictGEK_8f90.html> kriginggradpredictGEK</a> for the definitions of the gradient-enhanced versions of the covariance matrices, collocation matrix, regression matrix and training data. The only previously undefined term is \f$\underline{S}_{k}(\vec{x}_*)\f$ which is given as:

!> \f[ \underline{S}_{k}(\vec{x}_*) = g_{k}(\vec{x}_{*}) - \underline{l}_{*,k}^T \underline{K}^{-1} \underline{H} \f]

!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
!>         University of Wyoming
!> \date   May 2, 2012
!> \param(in) <b>ndim </b>: The dimension of the problem
!> \param(in) <b>ntot </b>: The number of Training points
!> \param(in) <b>X </b>:  The location of the training points (size=[ndimxntot])
!> \param(in) <b>gtot</b>: Number of derivative values included in training data (ndim*ntot if all derivatives are included at the training points)
!> \param(in) <b>pts</b>:  List identifying what point the derivative value is enforced at (size=[gtot] with values ranging from 1 to ntot)
!> \param(in) <b>dims</b>:  List identifying the dimension the derivative is taken with respect to (size=[gtot] with values ranging from 1 to ndim)
!> \param(in) <b>stot </b>: Number of Terms in the regression 
!> \param(in) <b>H</b>: The collocation matrix for the regression including derivative values. (size=[stotxntot+gtot]) <br>
!>                      Columns 1:ntot are the basis evaluated at the training points <br>
!>                      Columns ntot+1:ntot+gtot are the derivative of the basis evaluated at the training points 
!> \param(in) <b> beta</b>: Regression coefficients based on the optimal estimate for the Kriging model (size=[stot]) <br>
!>                          Supplied by <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> subroutine
!> \param(in) <b>hyper</b>: Hyperparameters for the Kriging Model (size=[ndim+3]) <br>
!>                          Supplied by <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> subroutine
!> \param(in) <b>mtot </b>: The number of test points, the places where function prediction are desired 
!> \param(in) <b>Xm </b>:  The location of the test points (size=[ndimxmtot])
!> \param(in) <b>gtotm</b>: The total number of derivatives to be predicted (ndim*mtot if the gradient at every test point is desired)
!> \param(in) <b> ptsm</b>: A vector identifying the test point where a particular derivative is specified (size=[gtotm] with values ranging from 1 to mtot)
!> \param(in) <b> dimsm</b>: A vector identifying the dimension the derivative is taken with respect to (size=[gtotm] with values ranging from 1 to ndim)
!> \param(in) <b>Gm </b>:  The derivative of the regression basis evaluated at the test points (size=[stotxgtotm])
!> \param(out) <b>Var</b>:  The predicted function values (size=[gtotm]) <br>
!>                          Because the variance prediction requires the covariance matrix of the training data, it is more expensive than the mean prediction <br>
!>                          If possible, predict the variance associated with multiple points using a single subroutine call
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul>
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. <br>
!>              When using gradient values, \f$\nu=1/2\f$ is not differentiable enough so \f$\nu \geq 3/2\f$ must be used<br>
!>              Must supply the same covariance flag as used in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>.

subroutine kriginggradvarianceGEK(ndim,ntot,X,gtot,pts,dims,stot,H,hyper,mtot,Xm,gtotm,ptsm,dimsm,Gm,Var,covarflagi)
  use choleskymod
  use covars
  use covarmatrix_mod
  use covarmatrix_grad_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: gtot
  integer, intent(in) :: pts(gtot), dims(gtot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot+gtot)
  real(8), intent(in) :: hyper(ndim+3)

  integer, intent(in) :: mtot
  real(8), intent(in) :: Xm(ndim,mtot)
  integer, intent(in) :: gtotm
  integer, intent(in) :: ptsm(gtotm), dimsm(gtotm)
  real(8), intent(in) :: Gm(stot,gtotm)

  real(8), intent(out) :: Var(gtotm)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN,sigmaNG
  real(8) Kmat(ntot+gtot,gtotm)
  real(8) Kcov(ntot+gtot,ntot+gtot)
 
  integer i,j,k,l

  real(8) S(stot,gtotm)

  real(8) A(stot,stot)

  real(8) Hs(ndim,ndim)

  real(8) C(stot,ntot+gtot),Q(stot),W(ntot+gtot)
  
  covarflag=covarflagi

!  Find Hyper Parameters
  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)
  sigmaNG=hyper(ndim+3)

  call covarmatrix(ndim,ntot,gtot,pts,dims,X,theta,Kcov)

  Kcov=sigma**2*Kcov
  do i=1,ntot
     Kcov(i,i)=Kcov(i,i)+sigmaN**2
  end do
  do i=1,gtot
     Kcov(ntot+i,ntot+i)=Kcov(ntot+i,ntot+i)+sigmaNG**2
  end do

!  Cholesky and Invert
  call invertsymmetric(ntot+gtot,Kcov)

!  Calculate Covariance Matrix
  call covarmatrix_grad(ndim,ntot,gtot,pts,dims,X,mtot,gtotm,ptsm,dimsm,Xm,theta,Kmat)

!  Add magnitude  
  Kmat=sigma**2*Kmat

!  Construct Regression Matrix
  call symmatrixmulttrans(ntot+gtot,Kcov,stot,H,C)
  call matrixmulttrans(stot,ntot+gtot,H,stot,C,A)

  call invertsymmetric(stot,A)

!  Regression
  do j=1,stot
     do i=1,gtotm
        S(j,i)=Gm(j,i)
        do k=1,ntot+gtot
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

     call symmatvec(ntot+gtot,Kcov,Kmat(:,i),W)

     do j=1,ntot+gtot
        Var(i)=Var(i)-Kmat(j,i)*W(j)
     end do

     call symmatvec(stot,A,S(:,i),Q)

     do j=1,stot
        Var(i)=Var(i)+S(j,i)*Q(j)
     end do

  end do

  return
end subroutine kriginggradvarianceGEK
