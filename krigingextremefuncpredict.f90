!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingextremefuncpredict.f90
!> \brief Subroutine used to determine the minimum or maximum value from a Kriging response surface (requires  <a href=buildkriging_8f90.html>buildkriging</a> to be called first) <br>
!>        Only Works for Ordinary Kriging (Constant Mean function)

!> \detail This subroutine predicts the global minimum or maximum function values from the Kriging surface (constructed using <a href=buildkriging_8f90.html>buildkriging</a>). Using an existing Kriging surface, the minimum or maximum value from the model is determined using either <a href=simplexsearch_8f90.html>simplexsearch</a> or <a href=patternsearch_8f90.html>patternsearch</a> (specified with the \a optflagi). Within these optimization subroutines, the code <a href=krigingfuncpredict_8f90.html>krigingfuncpredict </a> is invoked to predict the mean function value of the Kriging surface at a particular location. The optimization algorithms by design determine the minimum value of a function. To determine the maximum value (specified by \a flag=1), the objective is multiplied by \f$-1\f$ and minimization is performed. <br>
!> <br>
!> This optimization can only be performed on an Ordinary Kriging surface, meaning that the mean function is a constant value that is fitted during the construction of the mean function. From the universal Kriging framework used in <a href=buildkriging_8f90.html>buildkriging</a>, ordinary Kriging can be recovered by using a \f$p=0\f$ regression. This zeroth order regression requires a single term (meaning \f$stot=1\f$) and all the elements of the collocation matrix (\a H) are equal to one.

!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
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
!> \param(in out) <b> Xm </b>: Location of the extreme point (size=[ndim]) <br>
!>                             Starting point for the optimization is the input value for Xm. Final value overwrites this initial value.
!> \param(out) <b> Ym </b>: Minimum or Maximum value
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. For function only Kriging, all three options are available. <br>
!>              Good rule of thumb is to use the least smooth (\f$\nu=1/2\f$) unless there is a reason to assume more smoothness in the data <br>
!>              For small numbers of training points, higher \f$\nu\f$ can improve accuracy
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
subroutine krigingextremefuncpredict(ndim,ntot,X,stot,H,beta,V,hyper,Xm,Ym,covarflagi,optflag,flag,Lb,Ub)
  use argument
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot)
  real(8), intent(in) :: V(ntot),hyper(ndim+2),Beta(stot)

  real(8), intent(out) :: Xm(ndim)
  real(8), intent(out) :: Ym

  integer, intent(in) :: covarflagi,optflag,flag
  real(8), intent(in) :: Lb(ndim), Ub(ndim)

  type (generalarg) argwrap
  type (arg_kriging) kriging

  kriging % descriptor = 'Kriging'
  kriging % value = 'Function'
  kriging%ndim=ndim
  kriging%ntot=ntot
  allocate(kriging % X(ndim,ntot))
  kriging%X=X
  kriging%stot=stot
  allocate(kriging%H(stot,ntot))
  kriging%H=H
  allocate(kriging%beta(stot))
  kriging%beta=beta
  allocate(kriging%V(ntot))
  kriging%V=V
  allocate(kriging%hyper(ndim+2))
  kriging%hyper=hyper
  kriging%covarflag=covarflagi
  

  argwrap % descriptor = kriging%descriptor
  argwrap % krigingarg = kriging

  if (optflag==0) then
     call simplexsearch(ndim,Xm,Xm,Ym,Lb,Ub,flag,argwrap)
  else
     call patternsearch(ndim,Xm,Xm,Ym,Lb,Ub,flag,argwrap)
  end if


  return
end subroutine krigingextremefuncpredict
