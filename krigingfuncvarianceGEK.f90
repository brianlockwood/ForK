!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingfuncvarianceGEK.f90
!> \brief Subroutine to Variance of function values from Gradient-Enhanced Kriging Surface (requires  <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> to be called first)

!> \detail This subroutine calculates the variance associated with function predictions from a Gradient-Enhanced Kriging model. The details of this variance prediction are given in <a href=krigingfuncvariance_8f90.html>krigingfuncvariance</a> for a function-only model. For a gradient-enhanced model, the variance predictions take the same form as a function-only model and are given by the formula:

!>\f[
!> V\left[y_*\right]= cov(\vec{x}_{*},\vec{x}_{*})-\underline{k}_*^T \underline{K}^{-1} \underline{k}_*+ \underline{R}(\vec{x}_*) \underline{A}^{-1} \underline{R}(\vec{x}_*)^{T} \f]

!> The gradient-enhanced versions of the covariance matrices can be found in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> and <a href=krigingfuncpredictGEK_8f90.html>krigingfuncpredictGEK</a>. The gradient-enhanced regression matrix is defined using the gradient-enhanced collocation matrix and covariance matrix, given as:

!> \f[ \underline{A} = \underline{H} \underline{K}^{-1} \underline{H}^{T} \f]

!> Finally, the term \f$\underline{R}\f$ is given as:

!> \f[ \underline{R}(\vec{x}_{*}) = h^{T} - \underline{k}^{T}_{*} \underline{K}^{-1} \underline{H}^{T} \f]


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
!> \param(in) <b>Hm </b>:  The collocation matrix evaluated at the test points. The derivative of the basis is NOT required for the test points to predict function values (size=[stotxmtot])
!> \param(out) <b>S</b>:  The variance of the function values (size=[mtot]) <br>
!>                        Predicting the variance values requires the construction of the covariance matrix for the training data so it is more expensive than function predictions <br>
!>                        Best to predict the variance associated with multiple points using a single function call if possible.
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul>
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. <br>
!>              When using gradient values, \f$\nu=1/2\f$ is not differentiable enough so \f$\nu \geq 3/2\f$ must be used<br>
!>              Must supply the same covariance flag as used in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>.

subroutine krigingfuncvarianceGEK(ndim,ntot,X,gtot,pts,dims,stot,H,hyper,mtot,Xm,Hm,S,covarflagi)
  use choleskymod
  use covars
  use covarmatrix_mod
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
  real(8), intent(in) :: Hm(stot,mtot)

  real(8), intent(out) :: S(mtot)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN,sigmaNG
  real(8) Kmat(ntot+gtot,mtot)
  real(8) Kcov(ntot+gtot,ntot+gtot)
 
  integer i,j,k,l

  real(8) R(stot,mtot)

  real(8) A(stot,stot)

  real(8) C(stot,ntot+gtot)

  real(8) W(ntot+gtot), Q(stot)

  integer, parameter :: gtotm=0
  integer ptsm(gtotm),dimsm(gtotm)

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
  call covarmatrix(ndim,ntot,gtot,pts,dims,X,mtot,gtotm,ptsm,dimsm,Xm,theta,Kmat)

!  Add magnitude  
  Kmat=sigma**2*Kmat

!  Construct Regression Matrix
  call symmatrixmulttrans(ntot,Kcov,stot,H,C)
  call matrixmulttrans(stot,ntot,H,stot,C,A)

  call invertsymmetric(stot,A)

!  Regression
  do j=1,stot
     do i=1,mtot
        R(j,i)=Hm(j,i)
        do k=1,ntot+gtot
           R(j,i)=R(j,i)-Kmat(k,i)*C(j,k)
        end do
     end do
  end do

  !  Give Prediction
  do i=1,mtot
     call covarfunc(ndim,Xm(:,i),Xm(:,i),theta(:),S(i))
     S(i)=S(i)*sigma**2

     call symmatvec(ntot+gtot,Kcov,Kmat(:,i),W)
     
     do j=1,ntot+gtot
        S(i)=S(i)-Kmat(j,i)*W(j)
     end do

     call symmatvec(stot,A,R(:,i),Q)
     
     do j=1,stot
        S(i)=S(i)+R(j,i)*Q(j)
     end do

  end do

  return
end subroutine krigingfuncvarianceGEK
