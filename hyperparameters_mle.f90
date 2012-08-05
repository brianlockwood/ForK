!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file hyperparameters_mle.f90
!> \brief This module contains the subroutines used to set-up the max likelihood optimization for the function-only and gradient-enhanced Kriging module using the MLE likelihood formulation. This MLE formulation determines the covariance magnitude and regression parameters based on the optimality condition and only performs optimization over the covariance length scales. The likelihood for this optimization is evaluated using <a href=likelihood_mle_8f90.html>likelihood_mle</a>. The noise level ratio is specified in these subroutines as well as the bounds for the length scale parameters. 

module hyperparameters_mle_mod
    interface hyperparameters_mle
     module procedure hyperparameters_mle_func, hyperparameters_mle_grad
  end interface hyperparameters_mle
contains
!> \brief Subroutine used to set up the optimization required to determine hyperparameters for a function-only Kriging model. Optimization is performed only for the length scales used in the covariance function and covariance magnitude is determined explicitly. Optimization bounds for the length scales are found here. The hard-coded noise level used to keep the covariance matrix positive-definite is also specified here.  
  subroutine hyperparameters_mle_func(ndim,ntot,X,Y,stot,H,theta,sigma,sigmaN,logpy)
    !  Determine Optimimal Hyper parameters based on MLE
    !  Length scale is the only thing included in the optimization
    use argument
    use opt, only: optflag
    use magnitude_mod
    implicit none
    integer, intent(in) :: ndim, ntot
    real(8), intent(in) :: X(ndim,ntot),Y(ntot)
    integer, intent(in) :: stot
    real(8), intent(in) :: H(stot,ntot)

    real(8), intent(out) :: theta(ndim),sigma, sigmaN,logpy

    real(8) Lb(ndim),Ub(ndim)

    type (generalarg) argwrap
    type (arg_likelihood) like

    !  Set bounds for Hyperparameters (Completely Heuristic)

    ! Length Scale bounds
    Lb(1:ndim)=1.D-1
    Ub(1:ndim)=1.D1

    !  Set Noise values (used only for conditioning)
    sigmaN=1.D-1
    
    theta(:)=1.D0
    like % descriptor = 'Likelihood'
    like%ndim=ndim
    like%ntot=ntot
    allocate(like % X(ndim,ntot))
    like%X=X
    allocate(like % Y(ntot))
    like%Y=Y
    like%stot=stot
    allocate(like%H(stot,ntot))
    like%H=H
    like%sigmaN=sigmaN

    argwrap % descriptor = like%descriptor
    argwrap % likearg = like 

    if (optflag==0) then
       call simplexsearch(ndim,theta,theta,logpy,Lb,Ub,1,argwrap)
    else
       call patternsearch(ndim,theta,theta,logpy,Lb,Ub,1,argwrap) 
    end if
    call magnitude(ndim,ntot,X,Y,stot,H,theta,sigmaN,sigma)
    sigmaN=sigma*sigmaN

    return
  end subroutine hyperparameters_mle_func
!> \brief Subroutine used to set up the optimization required to determine hyperparameters for a gradient-enhanced Kriging model. Optimization is performed only for the length scales used in the covariance function and covariance magnitude is determined explicitly. Optimization bounds for the length scales are found here. The hard-coded noise level for the function and gradient used to keep the covariance matrix positive-definite are also specified here.  
  subroutine hyperparameters_mle_grad(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,sigma,sigmaN,sigmaNG,logpy)
    !  Determine Optimimal Hyper parameters based on MLE
    !  Length scale is the only thing included in the optimization
    use argument
    use opt, only: optflag
    use magnitude_mod
    implicit none
    integer, intent(in) :: ndim, ntot
    real(8), intent(in) :: X(ndim,ntot),Y(ntot)
    integer, intent(in) :: gtot
    integer, intent(in) :: pts(gtot), dims(gtot)
    real(8), intent(in) :: dY(gtot)
    integer, intent(in) :: stot
    real(8), intent(in) :: H(stot,ntot+gtot)

    real(8), intent(out) :: theta(ndim),sigma, sigmaN,sigmaNG,logpy

    real(8) Lb(ndim),Ub(ndim)

    type (generalarg) argwrap
    type (arg_likelihood_grad) like

    !  Set bounds for Hyperparameters (Completely Heuristic)

    ! Length Scale bounds
    Lb(1:ndim)=1.D-1
    Ub(1:ndim)=10.D0

    !  Set Noise values (used only for conditioning)
    sigmaN=1.D-1
    sigmaNG=5.D-1

    theta(:)=1.D0
    like % descriptor = 'Likelihood Grad'
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
    like%sigmaN=sigmaN
    like%sigmaNG=sigmaNG

    argwrap % descriptor = like%descriptor
    argwrap % gradlikearg = like 

    if (optflag==0) then
       call simplexsearch(ndim,theta,theta,logpy,Lb,Ub,1,argwrap)
    else
       call patternsearch(ndim,theta,theta,logpy,Lb,Ub,1,argwrap) 
    end if
    call magnitude(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,sigmaN,sigmaNG,sigma)
    sigmaN=sigma*sigmaN
    sigmaNG=sigma*sigmaNG

    return
  end subroutine hyperparameters_mle_grad
end module hyperparameters_mle_mod
