!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingfuncvariance.f90
!> \brief Subroutine to Predict variance of function values from Kriging Surface (requires  <a href=buildkriging_8f90.html>buildkriging</a> to be called first)

!> \detail Because of the stochastic nature of a Kriging model, the variance associated with mean function values can be predicted throughout the domain. The variance for function predictions are given as:
!>\f[
!> V\left[y_*\right]= cov(\vec{x}_{*},\vec{x}_{*})-k_*^T K^{-1} k_*+ R(\vec{x}_*) A^{-1} R(\vec{x}_*)^{T} \f]
!>
!> The first two terms, \f$cov(\vec{x}_{*},\vec{x}_{*})-k_*^T K^{-1} k_*\f$, account for the variance of associated with a zero mean gaussian process while the last term \f$R(\vec{x}_*) A^{-1} R(\vec{x}_*)^{T}\f$ accounts for the variance associated with the uncertainty associated with the regression parameters, where \f$ R(\vec{x}_*)\f$ is given as:

!> \f[ R(\vec{x}_{*}) = h^{T} - k_*^T K^{-1} H^{T} \f] <br>

!> To compute the function variance, the hyperparameters within the covariance function are required (supplied by <a href=buildkriging_8f90.html>buildkriging</a>). Hence, these subroutines must be used in conjunction with <a href=buildkriging_8f90.html>buildkriging</a> called first. <br>
!> Using the variance, a confidence interval can be established on the function predictions. For example, for a well training model, \f$ 95 \%\f$ of the true function values should lie within 2 standard deviations from the mean prediction, represented functionally as: 
!> \f{eqnarray*}{ y_{lower} &=& y_{*} - 2 \sqrt{V\left[y_*\right]} \\\
!>                y_{upper} &=& y_{*} + 2 \sqrt{V\left[y_*\right]}
!> \f}

!> Using this variance prediction, predictions from the model can be bounded. Additionally, the quality of the Kriging model can be tested using cross-validation procedures. 

!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 2, 2012
!> \param(in) <b>ndim </b>: The dimension of the problem
!> \param(in) <b>ntot </b>: The number of Training points
!> \param(in) <b>X </b>:  The location of the training points (size=[ndimxntot])
!> \param(in) <b>stot </b>: Number of Terms in the regression 
!> \param(in) <b>H</b>: The collocation matrix for the regression (size=[stotxntot]) 
!> \param(in) <b>hyper</b>: Hyperparameters for the Kriging Model (size=[ndim+2]) <br>
!>                          Supplied by <a href=buildkriging_8f90.html>buildkriging</a> subroutine
!> \param(in) <b>mtot </b>: The number of test points, the places where function prediction are desired 
!> \param(in) <b>Xm </b>:  The location of the test points (size=[ndimxmtot])
!> \param(in) <b>Hm </b>:  The collocation matrix evaluated at the test points (size=[stotxmtot])
!> \param(out) <b>S</b>:  The variance of the function values (size=[mtot]) <br>
!>                        Predicting the variance values requires the construction of the covariance matrix for the training data so it is more expensive than function predictions <br>
!>                        Best to predict the variance associated with multiple points using a single function call if possible.
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. For function only Kriging, all three options are available. <br>
!>              Must supply the same covariance flag as used in <a href=buildkriging_8f90.html>buildkriging</a>.
subroutine krigingfuncvariance(ndim,ntot,X,stot,H,hyper,mtot,Xm,Hm,S,covarflagi)
  use choleskymod
  use covars
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot)
  real(8), intent(in) :: hyper(ndim+2)

  integer, intent(in) :: mtot
  real(8), intent(in) :: Xm(ndim,mtot)
  real(8), intent(in) :: Hm(stot,mtot)

  real(8), intent(out) :: S(mtot)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN
  real(8) Kmat(ntot,mtot)
  real(8) Kcov(ntot,ntot)
 
  integer i,j,k,l

  real(8) R(stot,mtot)

  real(8) A(stot,stot)
  real(8) C(stot,ntot)

  real(8) Q(stot)
  real(8) W(ntot)

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
  call covarmatrix(ndim,ntot,X,mtot,Xm,theta,Kmat)

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
        do k=1,ntot
           R(j,i)=R(j,i)-Kmat(k,i)*C(j,k)
        end do
     end do
  end do

  !  Give Prediction
  do i=1,mtot
     call covarfunc(ndim,Xm(:,i),Xm(:,i),theta(:),S(i))
     S(i)=sigma**2*S(i)

     call symmatvec(ntot,Kcov,Kmat(:,i),W)
     
     do j=1,ntot
        S(i)=S(i)-Kmat(j,i)*W(j)
     end do

     call symmatvec(stot,A,R(:,i),Q)
     
     do j=1,stot
        S(i)=S(i)+R(j,i)*Q(j)
     end do
  end do

  return
end subroutine krigingfuncvariance
