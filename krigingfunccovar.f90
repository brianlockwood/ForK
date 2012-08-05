!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingfunccovar.f90
!> \brief This subroutine calculates the covariance matrix associated with the Kriging model's output distribution. This matrix is useful when the quantifying the uncertainty of general code predictions (such as a quantile prediction or statistic prediction). This subroutine is order \f$M^{2}\f$ in complexity, where \f$M\f$ is the number of test points. It is very slow, so use with caution.

!> \detail The covariance matrix for the Kriging output distribution can be useful when predicting the uncertainty (or confidence) in general predictions based on the Kriging surface. Although the variance can be used for quantifiying the uncertainty associated with point predictions, the uncertainty associated with global predictions such as moment statistics or quantiles requires sampling from the Kriging output distribution and computing multiple estimates of this global prediction. A distribution can then be constructed for these predictions, allowing a confidence interval to be associated with the predictions from the Kriging model. <br>

!> To sample from the Kriging model's output distribution, the output covariance matrix is constructed (\a S in this code) and Cholesky is performed (the result of which is denoted as \f$L\f$). Samples from the output distribution are given as:
!>\f[ Y_{p} = Y_{m} + L u\f]
!> where \f$Y_{m}\f$ is the mean of the Kriging function's output (the result of <a href=krigingfuncpredict_8f90.html>krigingfuncpredict</a>) and \f$u\f$ is a vector of normally distributioned values of length \f$M\f$ (normal distribution with zero mean and standard deviation of 1). The sampling itself is performed in the subroutine <a href=krigingfuncsample_8f90.html>krigingfuncsample</a> but it requires the covariance matrix to be constructed and cholesky performed.

!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 18, 2012
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
!> \param(out) <b>S</b>:  Covariance matrix between function values (size=[mtotxmtot]) <br>
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. For function only Kriging, all three options are available. <br>
!>              Must supply the same covariance flag as used in <a href=buildkriging_8f90.html>buildkriging</a>.
subroutine krigingfunccovar(ndim,ntot,X,stot,H,hyper,mtot,Xm,Hm,S,covarflagi)
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

  real(8), intent(out) :: S(mtot,mtot)

  integer, intent(in) :: covarflagi

  real(8) theta(ndim),sigma,sigmaN
  real(8) Kmat(ntot,mtot)
  real(8) Kcov(ntot,ntot)
 
  integer i,j,k,l

  real(8) R(stot,mtot)

  real(8) A(stot,stot)

  real(8) C(stot,ntot)

  real(8) W(ntot), Q(stot)

  real(8) Kmat2(ntot,mtot)
  real(8) R2(stot,mtot)

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
  call symmatrixmult(ntot,Kcov,mtot,Kmat,Kmat2)
  call symmatrixmult(stot,A,mtot,R,R2)

  do i=1,mtot
     do j=1,mtot
        WRITE(*,*) i,j
        call covarfunc(ndim,Xm(:,i),Xm(:,j),theta(:),S(i,j))
        S(i,j)=S(i,j)*sigma**2

        do l=1,ntot
           S(i,j)=S(i,j)-Kmat(l,i)*Kmat2(l,j)
        end do
        
        do l=1,stot
           S(i,j)=S(i,j)+R(l,i)*R2(l,j)
        end do
     end do
  end do

  do i=1,mtot
     S(i,i)=S(i,i)+sigmaN**2
  end do
  
  return
end subroutine krigingfunccovar
