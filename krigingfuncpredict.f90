!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingfuncpredict.f90
!> \brief Subroutine to Predict mean function values from Kriging Surface (requires  <a href=buildkriging_8f90.html>buildkriging</a> to be called first)

!> \detail  Model predictions throughout the domain are determined by sampling from the conditional distribution \f$y_* | \vec{X},Y\f$ using the covariance between points in the domain where \f$\vec{X},Y\f$ are the input and output training data.  The posterior mean predictions for an explicit mean are given by the formula:
!>\f[
!>	y(\vec{x}_{*}) | \vec{X},Y,m(x) = m(\vec{x}_{*}) + k_*^T K^{-1} (Y-m(\vec{x}_{*}))
!>\f]
!>where \f$k_{*}^{T}\f$ represents the covariance between the test point, \f$\vec{x}_{*}\f$, and the training points \f$\vec{X}\f$ (a row vector of length ntot). <br>
!> For a regression mean function, the function predictions take the form of:
!> \f[ y(\vec{x}_{*}) | \vec{X},Y,\beta = h^{T}(\vec{x}_{*}) \beta + k_*^T K^{-1} (Y-H^{T} \beta) \f]
!> The regression parameters \f$\beta\f$ and the hyperparamters in the covariance function are supplied by the subroutine <a href=buildkriging_8f90.html>buildkriging</a>. Using only this data, function predictions can be made; however, the construction and inverse of the covariance matrix can make the function predictions expensive. Because this matrix is inverted during the construction of the Kriging model, this work can be re-used for function predictions. Defining the processed data \f$V\f$ as:
!> \f[ V = K^{-1} (Y - H^{T} \beta) \f]
!> the function predictions are given by:
!> \f[ y(\vec{x}_{*}) | \vec{X},Y,\beta = h^{T}(\vec{x}_{*}) \beta + k_*^T V \f]

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
!> \param(in) <b>Hm </b>:  The collocation matrix evaluated at the test points (size=[stotxmtot])
!> \param(out) <b>Ym</b>:  The predicted function values (size=[mtot]) <br>
!>                         Using the processed data V, predicting function values is essentially linear with respect to the number of test points so mtot can be set to one and this subroutine can be called multiple times. <br>
!>                          Often the function values at a set of test points is required, hence the ability to make the predictions in a single function call.
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. For function only Kriging, all three options are available. <br>
!>              Must supply the same covariance flag as used in <a href=buildkriging_8f90.html>buildkriging</a>.
subroutine krigingfuncpredict(ndim,ntot,X,stot,H,beta,V,hyper,mtot,Xm,Hm,Ym,covarflagi)
  use covars, only: covarflag
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot)
  real(8), intent(in) :: V(ntot),hyper(ndim+2),Beta(stot)

  integer, intent(in) :: mtot
  real(8), intent(in) :: Xm(ndim,mtot)
  real(8), intent(in) :: Hm(stot,mtot)

  real(8), intent(out) :: Ym(mtot)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN
  real(8) Kmat(ntot,mtot)

  integer i,j,k,l

  real(8) tot

  covarflag=covarflagi

!  Find Hyper Parameters
  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)

!  Calculate Covariance Matrix
  call covarmatrix(ndim,ntot,X,mtot,Xm,theta,Kmat)

!  Add magnitude  
  Kmat=sigma**2*Kmat

!  Give Prediction
  do i=1,mtot
     Ym(i)=0.d0
     do j=1,stot
        Ym(i)=Ym(i)+Hm(j,i)*Beta(j)
     end do
     do j=1,ntot
        Ym(i)=Ym(i)+Kmat(j,i)*V(j)
     end do
  end do

  return
end subroutine krigingfuncpredict
