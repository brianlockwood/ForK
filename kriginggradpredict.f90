!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file kriginggradpredict.f90
!> \brief Subroutine to predict mean derivative values from Kriging Surface (requires  <a href=buildkriging_8f90.html>buildkriging</a> to be called first)

!> \detail Using the Kriging model, the derivative values from the Kriging model can be calculated using the covariance between derivative and function values.  Due to the linearity of differentiation, the covariance between the derivative and function values is given as:
!> \f[ cov(\frac{\partial y_{i}}{\partial x_{k}}, y_{j}) = \frac{\partial}{\partial x_{k}} k(\vec{x}_{i},\vec{x}_{j}) \f]

!> Additionally, the derivative of the regression basis can be used to predict the mean value of the derivative. Let the elements of the vector \f$g_{k}(\vec{x})\f$ represent the derivatives of basis functions with respect to \f$x_{k}\f$. Using this definition, the derivative predictions from the Kriging surface are given by:

!> \f[ \frac{\partial y(\vec{x}_{*})}{\partial x_{k}} = g_{k}(\vec{x}_{*}) \beta + l_{*,k}^{T} K^{-1} (Y - H^{T} \beta) \f]

!> where \f$l_{*,k}^{T}\f$ is a vector with elements representing the covariance between the derivative value with respect to \f$x_{k}\f$ and the training point function values. 

!> \f[  l_{*,k}^{T} = \left[\frac{\partial}{\partial x_{k}} k(\vec{x}_{*},\vec{X}_{1}) \dots \frac{\partial}{\partial x_{k}} k(\vec{x}_{*},\vec{X}_{N}) \right] \f]

!> Using the processed data produced by <a href=buildkriging_8f90.html>buildkriging</a> can be used to predict the derivative values inexpensively. Using the variable \f$ V \f$, the derivative values can be predicted as: 

!> \f[ \frac{\partial y(\vec{x}_{*})}{\partial x_{k}} = g_{k}(\vec{x}_{*}) \beta + l_{*,k}^{T} V \f]

!> To predict the gradient, the derivative of the covariance matrix is required. Hence, the derivative must be differentiable once. This allows all of the covariance functions to be used in this context. In addition to predicting the mean derivative value, the variance associated with the gradient can be predicted with <a href=kriginggradvariance_8f90.html>kriginggradvariance</a>. If the variance of the derivative prediction is desired, the covariance must be twice differentiable. This allows only the Matern functions with \f$ \nu \geq 3/2 \f$ to be used.

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
!> \param(in) <b>V</b>:     Processed Training Data (size=[ntot]) <br>
!>                          Supplied by <a href=buildkriging_8f90.html>buildkriging</a> subroutine
!> \param(in) <b>hyper</b>: Hyperparameters for the Kriging Model (size=[ndim+2]) <br>
!>                          Supplied by <a href=buildkriging_8f90.html>buildkriging</a> subroutine
!> \param(in) <b>mtot </b>: The number of test points, the places where function prediction are desired 
!> \param(in) <b>Xm </b>:  The location of the test points (size=[ndimxmtot])
!> \param(in) <b>gtotm</b>: The total number of derivatives to be predicted (ndim*mtot if the gradient at every test point is desired)
!> \param(in) <b> ptsm</b>: A vector identifying the test point where a particular derivative is specified (size=[gtotm] with values ranging from 1 to mtot)
!> \param(in) <b> dimsm</b>: A vector identifying the dimension the derivative is taken with respect to (size=[gtotm] with values ranging from 1 to ndim)
!> \param(in) <b>Gm </b>:  The derivative of the regression basis evaluated at the test points (size=[stotxgtotm])
!> \param(out) <b>dYm</b>:  The predicted function values (size=[gtotm]) <br>
!>                         Using the processed data V, predicting derivative values is essentially linear with respect to the number of test points so mtot can be set to one and this subroutine can be called multiple times if desired. <br>
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. For function only Kriging, all three options are available. <br>
!>              Must supply the same covariance flag as used in <a href=buildkriging_8f90.html>buildkriging</a>.

subroutine kriginggradpredict(ndim,ntot,X,stot,H,beta,V,hyper,mtot,Xm,gtotm,ptsm,dimsm,Gm,dYm,covarflagi)
  use covars, only: covarflag
  use covarmatrix_grad_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot)
  real(8), intent(in) :: V(ntot),hyper(ndim+2),Beta(stot)

  integer, intent(in) :: mtot
  real(8), intent(in) :: Xm(ndim,mtot)
  integer, intent(in) :: gtotm
  integer, intent(in) :: ptsm(gtotm), dimsm(gtotm)
  real(8), intent(in) :: Gm(stot,gtotm)

  real(8), intent(out) :: dYm(gtotm)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN
  real(8) Kmat(ntot,gtotm)

  integer i,j,k,l

  covarflag=covarflagi

!  Find Hyper Parameters
  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)

!  Calculate Covariance Matrix

  call covarmatrix_grad(ndim,ntot,X,mtot,gtotm,ptsm,dimsm,Xm,theta,Kmat)

!  Add magnitude  
  Kmat=sigma**2*Kmat

!  Give Prediction
  do i=1,gtotm
     dYm(i)=0.d0
     do j=1,stot
        dYm(i)=dYm(i)+Gm(j,i)*Beta(j)
     end do
     do j=1,ntot
        dYm(i)=dYm(i)+Kmat(j,i)*V(j)
     end do
  end do

  return
end subroutine kriginggradpredict
