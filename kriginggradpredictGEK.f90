!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file kriginggradpredictGEK.f90
!> \brief Subroutine to predict mean derivative values from Gradient-Enhanced Kriging Surface (requires  <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> to be called first)

!> \detail Predicting derivative values from a gradient-enhanced Kriging model follows the same procedure as predicting derivative values from function-only Kriging models (see <a href=kriginggradpredict_8f90.html>kriginggradpredict</a> for details). For a gradient-enhanced model, derivative values are predicted as:

!> \f[ \frac{\partial y(\vec{x}_{*})}{\partial x_{k}} = g_{k}(\vec{x}_{*}) \beta + \underline{l}_{*,k}^{T} \underline{K}^{-1} (\underline{Y} - \underline{H}^{T} \beta) \f]

!> where the underline denotes the block definition required to account for derivative values. The block definitions of the training data, collocation matrix and covariance matrix are given in the documentation of <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>. The matrix \f$ \underline{l}_{*,k}^{T} \f$ is a block matrix whose elements represent covariance between derivative values and function values as well as derivative values with other derivative values:

!> \f[ \underline{l}_{*,k}^{T} = \left[ cov( \frac{\partial y^{*}}{\partial x_{k}}, Y) cov(\frac{\partial y^{*}}{\partial x_{k}}, \delta Y) \right]

!> where the elements of \f$cov( \frac{\partial y^{*}}{\partial x_{k}}, Y)\f$ are given by:

!> \f[ cov(\frac{\partial y^{*}}{\partial x_{k}}, Y_{j}) = \frac{\partial}{\partial x_{k}} k(\vec{x}_{*},X_{j}) \f]

!> and the elements of \f$cov(\frac{\partial y^{*}}{\partial x_{k}}, \delta Y)\f$ are given as:

!> \f[ cov(\frac{\partial y^{*}}{\partial x_{k}}, \delta Y_{j,l}) = \frac{\partial^{2}}{\partial x_{k} \partial x_{l}} k(\vec{x}_{*},X_{j}) \f]

!>  Using the processed data produced by <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>, derivative values can be calculated as:

!> \f[ \frac{\partial y(\vec{x}_{*})}{\partial x_{k}} = g_{k}(\vec{x}_{*}) \beta + \underline{l}_{*,k}^{T} \underline{V} \f]

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
!> \param(in) <b>V</b>:     Processed Training Data (size=[ntot+gtot]) <br>
!>                          Supplied by <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> subroutine
!> \param(in) <b>hyper</b>: Hyperparameters for the Kriging Model (size=[ndim+3]) <br>
!>                          Supplied by <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> subroutine
!> \param(in) <b>mtot </b>: The number of test points, the places where function prediction are desired 
!> \param(in) <b>Xm </b>:  The location of the test points (size=[ndimxmtot])
!> \param(in) <b>gtotm</b>: The total number of derivatives to be predicted (ndim*mtot if the gradient at every test point is desired)
!> \param(in) <b> ptsm</b>: A vector identifying the test point where a particular derivative is specified (size=[gtotm] with values ranging from 1 to mtot)
!> \param(in) <b> dimsm</b>: A vector identifying the dimension the derivative is taken with respect to (size=[gtotm] with values ranging from 1 to ndim)
!> \param(in) <b>Gm </b>:  The derivative of the regression basis evaluated at the test points (size=[stotxgtotm])
!> \param(out) <b>dYm</b>:  The predicted function values (size=[gtotm]) <br>
!>                         Using the processed data V, predicting derivative values is essentially linear with respect to the number of test points so mtot can be set to one and this subroutine can be called multiple times if desired. <br>
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul>
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. <br>
!>              When using gradient values, \f$\nu=1/2\f$ is not differentiable enough so \f$\nu \geq 3/2\f$ must be used<br>
!>              Must supply the same covariance flag as used in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>.

subroutine kriginggradpredictGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,V,hyper,mtot,Xm,gtotm,ptsm,dimsm,Gm,dYm,covarflagi)
  use covars, only: covarflag
  use covarmatrix_grad_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: gtot
  integer, intent(in) :: pts(gtot), dims(gtot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot+gtot)
  real(8), intent(in) :: V(ntot+gtot),hyper(ndim+3),Beta(stot)

  integer, intent(in) :: mtot
  real(8), intent(in) :: Xm(ndim,mtot)
  integer, intent(in) :: gtotm
  integer, intent(in) :: ptsm(gtotm), dimsm(gtotm)
  real(8), intent(in) :: Gm(stot,gtotm)

  real(8), intent(out) :: dYm(gtotm)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN,sigmaNG
  real(8) Kmat(ntot+gtot,gtotm)

  integer i,j,k,l

  covarflag=covarflagi

!  Find Hyper Parameters
  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)
  sigmaNG=hyper(ndim+3)

!  Calculate Covariance Matrix
  call covarmatrix_grad(ndim,ntot,gtot,pts,dims,X,mtot,gtotm,ptsm,dimsm,Xm,theta,Kmat)

!  Add magnitude  
  Kmat=sigma**2*Kmat

!  Give Prediction
  do i=1,gtotm
     dYm(i)=0.d0
     do j=1,stot
        dYm(i)=dYm(i)+Gm(j,i)*Beta(j)
     end do
     do j=1,ntot+gtot
        dYm(i)=dYm(i)+Kmat(j,i)*V(j)
     end do
  end do

  return
end subroutine kriginggradpredictGEK
