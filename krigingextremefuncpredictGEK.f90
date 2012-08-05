!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingextremefuncpredictGEK.f90
!> \brief Subroutine to Predict mean function values from Gradient-Enhanced Kriging Surface (requires  <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a> to be called first) <br>
!>        Only Works for Ordinary Kriging (Constant Mean function)

!> \detail See documentation for <a href=krigingextremefuncpredict_8f90.html>krigingextremefuncpredict</a> for details. This subroutine is the gradient-enhanced version and is functionally the same as the function-only subroutine.


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
!> \param(out) <b> Xm </b>: Location of the extreme point (size=[ndim])
!> \param(out) <b> Ym </b>: Minimum or Maximum value
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul>
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. <br>
!>              When using gradient values, \f$\nu=1/2\f$ is not differentiable enough so \f$\nu \geq 3/2\f$ must be used<br>
!>              Must supply the same covariance flag as used in <a href=buildkrigingGEK_8f90.html>buildkrigingGEK</a>.
!>  \param(in) <b>optflagi</b>: Flag to govern the optimization algorithm used to determine minimum or maximum value <br>
!>                <ul><li>optflag=0  Simplex Search (Nelder-Mead Method) to determine extreme value <br>
!>                           Local optimization technique using simplex elements to determine search direction and line search to determine step sizes <br>
!>                           Fastest method but inherently sequential and may stall on local optimum
!>                <li>optflag=1  Pattern Search used to determine extreme value <br>
!>                           Global optimization method with fixed search directions (forward and backward in each direction) <br>
!>                           Each iteration requires \f$2*ndim\f$ function evaluations; however, these may be performed in parallel using openmp (Each thread performs own likelihood evaluation so OMP_STACK_SIZE must be set large enough) <br> </ul>
!> \param(in) <b> flag </b>  Determines whether minimum or maximum value is found <br>
!>                <ul><li> flag = 0 Subroutine finds the minimum value
!>                    <li> flag = 1 Subroutine finds the maximum value</ul>
!> \param(in) <b> Lb </b> Lower Bound for each variable (size = [ndim]) <br>
!>                        Should be the minimum possible value of Xm in each dimension
!> \param(in) <b> Ub </b> Upper Bound for each variable (size = [ndim]) <br>
!>                        Should be the maximum possible value of Xm in each dimension

subroutine krigingextremefuncpredictGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,V,hyper,Xm,Ym,covarflagi,optflag,flag,Lb,Ub)
  use argument
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: gtot
  integer, intent(in) :: pts(gtot), dims(gtot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot+gtot)
  real(8), intent(in) :: V(ntot+gtot),hyper(ndim+3),Beta(stot)

  real(8), intent(out) :: Xm(ndim)
  real(8), intent(out) :: Ym

  integer, intent(in) :: covarflagi,optflag,flag
  real(8), intent(in) :: Lb(ndim), Ub(ndim)

  type (generalarg) argwrap
  type (arg_kriging_grad) kriging

  kriging % descriptor = 'Kriging Grad'
  kriging % value = 'Function'
  kriging%ndim=ndim
  kriging%ntot=ntot
  allocate(kriging % X(ndim,ntot))
  kriging%X=X
  kriging%gtot=gtot
  allocate(kriging%pts(gtot))
  kriging%pts=pts
  allocate(kriging%dims(gtot))
  kriging%dims=dims
  kriging%stot=stot
  allocate(kriging%H(stot,ntot+gtot))
  kriging%H=H
  allocate(kriging%beta(stot))
  kriging%beta=beta
  allocate(kriging%V(ntot+gtot))
  kriging%V=V
  allocate(kriging%hyper(ndim+3))
  kriging%hyper=hyper
  kriging%covarflag=covarflagi
  

  argwrap % descriptor = kriging%descriptor
  argwrap % gradkrigingarg = kriging

  if (optflag==0) then
     call simplexsearch(ndim,Xm,Xm,Ym,Lb,Ub,flag,argwrap)
  else
     call patternsearch(ndim,Xm,Xm,Ym,Lb,Ub,flag,argwrap)
  end if


  return
end subroutine krigingextremefuncpredictGEK
