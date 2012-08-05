!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file hyperparameters_all.f90
!> \brief This module contains the subroutines used to set-up the max likelihood optimization for the function-only and gradient-enhanced Kriging module. This likelihood formulation determines the regression parameters based on an optimality condition and fits all other covariance parameters through optimization(length, magnitude and noise level). The likelihood for this optimization is evaluated using <a href=likelihood_8f90.html>likelihood</a>. The bounds for these hyperparameters are set within these subroutines.

module hyperparameters_all_mod
    interface hyperparameters_all
     module procedure hyperparameters_all_func, hyperparameters_all_grad
  end interface hyperparameters_all
contains
!> \brief Subroutine used to set up the optimization required to determine the best hyperparameters for a function-only Kriging model using the Vague Prior mean function. Optimization is used to determine all hyperparameters, including the covariance magnitude and noise. Optimization bounds for each parameter are located here.
  subroutine hyperparameters_all_func(ndim,ntot,X,Y,stot,H,theta,sigma,sigmaN,logpy)
  !  Determine Optimimal Hyper parameters based on GP with Vague prior Mean function
  !  Length scale, Covariance Magnitude, Noise and Gradient Noise all included in the optimization
  !  Use patternsearch.f90 and likelihood.f90 
  use argument
  use opt,only: optflag
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot)
  
  real(8), intent(out) :: theta(ndim),sigma, sigmaN,logpy

  real(8) Lb(ndim+2),Ub(ndim+2)
  
  real(8) hyper(ndim+2)

  type (generalarg) argwrap
  type (arg_alllikelihood) like
  
!  Set bounds for Hyperparameters (Completely Heuristic)
 
  ! Length Scale bounds
  Lb(1:ndim)=1.D-12
  Ub(1:ndim)=1.D12

  ! Covariance Magnitude 
  Lb(ndim+1)=1.D-2
  Ub(ndim+1)=1.D12

  ! Noise on Function 
  Lb(ndim+2)=1.D-2
  Ub(ndim+2)=1.D12
  
  !  Starting point for Patternsearch (Best to use fixed starting point so parameters are deterministic)
  hyper(1:ndim)=1.D0
  hyper(ndim+1)=1.D0
  hyper(ndim+2)=0.01D0

  like % descriptor = 'AllLikelihood'
  like%ndim=ndim
  like%ntot=ntot
  allocate(like % X(ndim,ntot))
  like%X=X
  allocate(like % Y(ntot))
  like%Y=Y
  like%stot=stot
  allocate(like%H(stot,ntot))
  like%H=H

  argwrap % descriptor = like%descriptor
  argwrap % alllikearg = like 

  if (optflag==0) then
     call simplexsearch(ndim+2,hyper,hyper,logpy,Lb,Ub,1,argwrap)
  else
     call patternsearch(ndim+2,hyper,hyper,logpy,Lb,Ub,1,argwrap) 
  end if

  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)

  return
end subroutine hyperparameters_all_func
!> \brief Subroutine used to set up the optimization required to determine the best hyperparameters for a gradient-enhanced Kriging model using the Vague Prior mean function. Optimization is used to determine all hyperparameters, including the covariance magnitude, function noise level and derivative noise level. Optimization bounds for each parameter are located here.
subroutine hyperparameters_all_grad(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,sigma,sigmaN,sigmaNG,logpy)
  !  Determine Optimimal Hyper parameters based on GP with Vague prior Mean function
  !  Length scale, Covariance Magnitude, Noise and Gradient Noise all included in the optimization
  !  Use patternsearch.f90 and likelihood.f90 
  use argument
  use opt,only: optflag
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: gtot
  integer, intent(in) :: pts(gtot), dims(gtot)
  real(8), intent(in) :: dY(gtot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot+gtot)
  
  real(8), intent(out) :: theta(ndim),sigma, sigmaN,sigmaNG,logpy

  real(8) Lb(ndim+3),Ub(ndim+3)
  
  real(8) hyper(ndim+3)

  type (generalarg) argwrap
  type (arg_alllikelihood_grad) like
  
!  Set bounds for Hyperparameters (Completely Heuristic)
 
 ! Length Scale bounds
  Lb(1:ndim)=1.D-12
  Ub(1:ndim)=1.D12

  ! Covariance Magnitude 
  Lb(ndim+1)=1.D-3
  Ub(ndim+1)=1.D3

  ! Noise on Function 
  Lb(ndim+2)=1.D-2
  Ub(ndim+2)=1.D12
 
  ! Noise on Gradients
  Lb(ndim+3)=5.D-2
  Ub(ndim+3)=5.D12

  !  Starting point for Patternsearch (Best to use fixed starting point so parameters are deterministic)
  hyper(1:ndim)=1.D0
  hyper(ndim+1)=1.D0
  hyper(ndim+2)=0.01D0
  hyper(ndim+3)=0.05D0

  like % descriptor = 'AllLikelihood Grad' 
  like%ndim=ndim
  like%ntot=ntot
  like%gtot=gtot
  allocate(like % X(ndim,ntot))
  like%X=X
  allocate(like % Y(ntot))
  like%Y=Y
  allocate(like%pts(gtot))
  like%pts=pts
  allocate(like%dims(gtot))
  like%dims=dims
  allocate(like%dY(gtot))
  like%dY=dY
  like%stot=stot
  allocate(like%H(stot,ntot+gtot))
  like%H=H

  argwrap % descriptor = like%descriptor
  argwrap % gradalllikearg = like 

  if (optflag==0) then
     call simplexsearch(ndim+3,hyper,hyper,logpy,Lb,Ub,1,argwrap)
  else
     call patternsearch(ndim+3,hyper,hyper,logpy,Lb,Ub,1,argwrap) 
  end if

  theta(1:ndim)=hyper(1:ndim)
  sigma=hyper(ndim+1)
  sigmaN=hyper(ndim+2)
  sigmaNG=hyper(ndim+3)

  return
end subroutine hyperparameters_all_grad
end module hyperparameters_all_mod
