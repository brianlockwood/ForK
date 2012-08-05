!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingeipredictGEK.f90
!> \brief Subroutine to Predict Expected improvement from Gradient-enhanced Kriging Surface based on a given minimum value (requires  <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> to be called first)

!> \detail See documentation for <a href=krigingeipredict_8f90.html>krigingeipredict</a> for details. This subroutine is the gradient-enhanced version and is functionally the same as the function-only subroutine.

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
subroutine krigingeipredictGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,V,hyper,mtot,Xm,Hm,EI,covarflagi,Ymin)
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

  real(8), intent(out) :: EI(mtot)

  integer, intent(in) :: covarflagi
  real(8), intent(in) :: Ymin

  real(8) Ym(mtot), Sm(mtot)
  real(8) s,chi,phi

  real(8) Pi

  integer i

  Pi=4.D0*ATAN(1.D0)

 call kriginfuncpredictGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,V,hyper,mtot,Xm,Hm,Ym,covarflagi)
 call kriginfuncvarianceGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,hyper,mtot,Xm,Hm,Sm,covarflagi)

  do i=1,mtot
     if (Sm(i)<=1.D-12) then
        EI(i)=0.D0
     else
        s=sqrt(Sm(i))
        chi=exp(-(Ymin-Ym(i))**2/(2.D0*s**2))/sqrt(2.D0*Pi)
        phi=0.5*(1.D0+erf((Ymin-Ym(i))/sqrt(2.D0*s**2)))
        EI(i)=(Ymin-Ym(i))*phi+s*chi
     end if
  end do

  return
end subroutine krigingeipredictGEK
