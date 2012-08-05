!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingfuncpredictGEK.f90
!> \brief Subroutine to Predict mean function values from Gradient-Enhanced Kriging Surface (requires  <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> to be called first)

!> \detail This subroutine is used to predict function values from a Gradient-enhanced Kriging model (built using <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>). The details of predicting function values from a Kriging model are detailed in the <a href=krigingfuncpredict_8f90.html>krigingfuncpredict</a> documentation. To make predictions from a gradient-enhanced Kriging model, the covariance matrix and training data must be redefined. Using a regression mean function, the mean function predictions at a test point are given as:

!>\f[
!>	y(\vec{x}_{*}) | \vec{X},Y,m(x) = h^{T}(\vec{x}_{*}) \beta + \underline{k}_*^T K^{-1} (\underline{Y}-\underline{H}^{T} \beta)
!>\f]
!> where the covariance matrix between test and training points has been redefined to include the correlation between derivative values. Additionally, the collocation matrix and  training data are extend to include derivative values. The extension of the collocation matrix and training data for derivative values are found in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>. The covariance matrix between training data included derivative values is given as:

!> \f[
!>    \underline{K} = \left[ \begin{array}{cc}
!>                cov(Y,Y) & cov(Y,\delta Y) \\\
!>                cov(\delta Y,Y) & cov(\delta Y,\delta Y)
!>               \end{array} \right] \f]
!>  where \f$cov(Y,Y)\f$ is the covariance matrix between function values (same as the covariance matrix used in a function-only model), \f$cov(\delta Y, Y) \f$ is a rectangular matrix whose elements represent the covariance between the function values and derivative values at the training points and \f$cov(\delta Y, \delta Y) \f$ represents the covariance between derivative observations. The details of this matrix are also found in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>.

!> Finally, the covariance between the test points and training points must be extended to include derivatives (denoted as \f$\underline{k}_*^T\f$). To include derivatives, this covariance is a rectangular block matrix defined as:

!> \f[ \underline{k}_*^T = \left[cov(y^{*}, Y) cov(y^{*}, \delta Y) \right]

!> where \f$cov(y^{*}, Y)\f$ is a matrix (size = [mtot \times ntot] where mtot is the number of test points and ntot is the number of training points) whose elements represent the covariance between test function values and training function values. The elements of matrix \f$cov(y^{*}, \delta Y)\f$ (size = [mtot \times gtot] where gtot is the total number of derivative values in the training set) represent the covariance between test function values and training derivative values. 

!> The subroutine <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> produces a vector of processed data defined as:

!>  \f[ \underline{V} = \underline{K}^{-1} \left( \underline{Y} - \underline{H}^{T} \beta \right) \f] 

!>  Using this definition, function predictions can be calculated inexpensively using the formula:

!> \f[
!>	y(\vec{x}_{*}) | \vec{X},Y,m(x) = h^{T}(\vec{x}_{*}) \beta + \underline{k}_*^T \underline{V}
!>\f]

!> Because the covariance matrix between training data requires the second derivative of the covariance function, the Matern function with \f$\nu \geq 3/2 \f$ must be used for constructing and making predictions from the Gradient-enhanced Kriging model.

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
!> \param(in) <b>Hm </b>:  The collocation matrix evaluated at the test points. The derivative of the basis is NOT required for the test points to predict function values (size=[stotxmtot])
!> \param(out) <b>Ym</b>:  The predicted function values (size=[mtot]) <br>
!>                         Using the processed data V, predicting function values is essentially linear with respect to the number of test points so mtot can be set to one and this subroutine can be called multiple times. <br>
!>                          Often the function values at a set of test points is required, hence the ability to make the predictions in a single function call.
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul>
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. <br>
!>              When using gradient values, \f$\nu=1/2\f$ is not differentiable enough so \f$\nu \geq 3/2\f$ must be used<br>
!>              Must supply the same covariance flag as used in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>.
subroutine krigingfuncpredictGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,V,hyper,mtot,Xm,Hm,Ym,covarflagi)
  use covars, only: covarflag
  use covarmatrix_mod
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
  real(8), intent(in) :: Hm(stot,mtot)

  real(8), intent(out) :: Ym(mtot)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN,sigmaNG
  real(8) Kmat(ntot+gtot,mtot)

  integer i,j,k,l

  integer, parameter :: gtotm=0
  integer ptsm(gtotm),dimsm(gtotm)

  covarflag=covarflagi

!  Find Hyper Parameters
  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)
  sigmaNG=hyper(ndim+3)

!  Calculate Covariance Matrix
  call covarmatrix(ndim,ntot,gtot,pts,dims,X,mtot,gtotm,ptsm,dimsm,Xm,theta,Kmat)

!  Add magnitude  
  Kmat=sigma**2*Kmat

!  Give Prediction
  do i=1,mtot
     Ym(i)=0.d0
     do j=1,stot
        Ym(i)=Ym(i)+Hm(j,i)*Beta(j)
     end do
     do j=1,ntot+gtot
        Ym(i)=Ym(i)+Kmat(j,i)*V(j)
     end do
  end do

  return
end subroutine krigingfuncpredictGEK
