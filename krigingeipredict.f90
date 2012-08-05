!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingeipredict.f90
!> \brief Subroutine to Predict Expected improvement from Kriging Surface based on a given minimum value (requires  <a href=buildkriging_8f90.html>buildkriging</a> to be called first)

!> \detail Due to the stochastic nature of Kriging models, the variance and mean predictions can be used to predict the expected outcomes for various scenarios. A typical scenario is to determine the expected improvement of a current minimum value if a new sample point is generated at a particular location in the domain. This expected improvement is useful within the context of optimization, where the expected improvement is used as a criteria for determining new sample locations with a high probability of finding a new minimum. Because the Kriging surface is based on a gaussian process, the expected improvement can be computed analytically using the mean function predictions and the variance predictions from the Kriging model:

!>\f[	EI(x) = \begin{cases} (y_{min} - y^{*}(x)) \Phi \left( \frac{y_{min}-y^{*}(x)}{s(x)} \right) + s(x) \phi \left( \frac{y_{min}-y^{*}(x)}{s(x)} \right) & \text{if $s(x)>0$,} \\\
!>		0 & \text{if $s(x)=0$} \end{cases}\f]

!> where \Phi is the probability density function and \phi is the cummulative distribution function for the Normal distribution.


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
!> \param(out) <b>EI</b>:  The predicted expected improvement value (size=[mtot]) <br>
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. For function only Kriging, all three options are available. <br>
!>              Must supply the same covariance flag as used in <a href=buildkriging_8f90.html>buildkriging</a>.
!>   \param(in) <b> Ymin </b>: Minimum value to be improved upon
subroutine krigingeipredict(ndim,ntot,X,stot,H,beta,V,hyper,mtot,Xm,Hm,EI,covarflagi,Ymin)
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

  real(8), intent(out) :: EI(mtot)

  integer, intent(in) :: covarflagi
  real(8), intent(in) :: Ymin

  integer i,j,k,l

  real(8) Ym(mtot), Sm(mtot)
  real(8) s,chi,phi

  real(8) Pi

  Pi=4.D0*ATAN(1.D0)

  call kriginfuncpredict(ndim,ntot,X,stot,H,beta,V,hyper,mtot,Xm,Hm,Ym,covarflagi)
  call kriginfuncvariance(ndim,ntot,X,stot,H,beta,hyper,mtot,Xm,Hm,Sm,covarflagi)

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
end subroutine krigingeipredict
