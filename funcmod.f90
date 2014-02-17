!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file funcmod.f90
!> \brief This module contains the function calls that can be used in conjunction with the optimization algorithms <a href=patternsearch_8f90.html>patternsearch</a> and <a href=simplexsearch_8f90.html>simplexsearch</a>. The call to the subroutine "func" is made by the subroutine <a href=funcwrapper_8f90.html>funcwrapper</a>. The subroutine is overloaded based on the type of the \a arg argument. For each \a arg type, a different function is called (likelihood, likelihood_mle, or kriging). This module in addition to <a href=argument_8f90.html>argument</a> and <a href=funcwrapper_8f90.html>funcwrapper</a> must be modified if optimization is to be performed on different objectives.

module funcmod
  use argument
  implicit none

  interface func
     module procedure func_likelihood_grad,func_likelihood_all_grad,func_likelihood,func_likelihood_all, func_kriging, func_kriging_grad
  end interface func
contains
  subroutine func_likelihood(n,X,f,arg)
    use likelihood_mle_mod
    integer, intent(in) :: n
    real(8), intent(in) :: X(n)
    real(8), intent(out) :: f
    
    type (arg_likelihood), intent(in) :: arg

    call likelihood_mle(arg%ndim,arg%ntot,arg%X,arg%Y,arg%stot,arg%H,X(:),arg%sigmaN,f)

    return 
  end subroutine func_likelihood
  subroutine func_likelihood_grad(n,X,f,arg)
    use likelihood_mle_mod
    integer, intent(in) :: n
    real(8), intent(in) :: X(n)
    real(8), intent(out) :: f
    
    type (arg_likelihood_grad), intent(in) :: arg

    call likelihood_mle(arg%ndim,arg%ntot,arg%X,arg%Y,arg%gtot,arg%pts,arg%dims,arg%dY,arg%stot,arg%H,X(:),arg%sigmaN,arg%sigmaNG,f)

    return 
  end subroutine func_likelihood_grad
  subroutine func_likelihood_all(n,X,f,arg)
    use likelihood_mod
    integer, intent(in) :: n
    real(8), intent(in) :: X(n)
    real(8), intent(out) :: f
    
    type (arg_alllikelihood), intent(in) :: arg

    call likelihood(arg%ndim,arg%ntot,arg%X,arg%Y,arg%stot,arg%H,X(1:arg%ndim),X(arg%ndim+1),X(arg%ndim+2),f)
 
    return 
  end subroutine func_likelihood_all
  subroutine func_likelihood_all_grad(n,X,f,arg)
    use likelihood_mod
    integer, intent(in) :: n
    real(8), intent(in) :: X(n)
    real(8), intent(out) :: f
    
    type (arg_alllikelihood_grad), intent(in) :: arg

    call likelihood(arg%ndim,arg%ntot,arg%X,arg%Y,arg%gtot,arg%pts,arg%dims,arg%dY,arg%stot,arg%H,X(1:arg%ndim),X(arg%ndim+1),X(arg%ndim+2),X(arg%ndim+3),f)

    return 
  end subroutine func_likelihood_all_grad
  subroutine func_kriging(n,X,f,arg)
    integer, intent(in) :: n
    real(8), intent(in) :: X(n)
    real(8), intent(out) :: f
    
    type (arg_kriging), intent(in) :: arg

    real(8) Hm(arg%stot)
    real(8) Yt, Sm,s,chi,phi,Ymin
    real(8) Pi
    
    Pi=4.D0*ATAN(1.D0)

    Hm=1.D0
    if (arg%value=='Function') then
       call krigingfuncpredict(arg%ndim,arg%ntot,arg%X,arg%stot,arg%H,arg%beta,arg%V,arg%hyper,1,X,Hm,f,arg%covarflag)
    elseif (arg%value=='Variance') then
       call krigingfuncvariance(arg%ndim,arg%ntot,arg%X,arg%stot,arg%H,arg%beta,arg%hyper,1,X,Hm,f,arg%covarflag)
    elseif (arg%value=='EI') then
       Ymin=arg%ymin
       call krigingfuncpredict(arg%ndim,arg%ntot,arg%X,arg%stot,arg%H,arg%beta,arg%V,arg%hyper,1,X,Hm,Yt,arg%covarflag)
       call krigingfuncvariance(arg%ndim,arg%ntot,arg%X,arg%stot,arg%H,arg%beta,arg%hyper,1,X,Hm,Sm,arg%covarflag)

!       WRITE(*,*) Yt, Sm
     if (Sm<=1.D-12) then
        f=0.D0
     else
        s=sqrt(Sm)
        chi=exp(-(Ymin-Yt)**2/(2.D0*s**2))/sqrt(2.D0*Pi)
        phi=0.5*(1.D0+erf((Ymin-Yt)/sqrt(2.D0*s**2)))
        f=(Ymin-Yt)*phi+s*chi
!        WRITE(*,*) s, chi, phi, f
     end if
!     WRITE(*,*) f
!     pause
  end if
  
    return 
  end subroutine func_kriging
  subroutine func_kriging_grad(n,X,f,arg)
    integer, intent(in) :: n
    real(8), intent(in) :: X(n)
    real(8), intent(out) :: f
    
    type (arg_kriging_grad), intent(in) :: arg

    real(8) Hm(arg%stot)
    real(8) Yt, Sm,s,chi,phi,Ymin
    real(8) Pi
    
    Pi=4.D0*ATAN(1.D0)

    Hm=1.D0
    if (arg%value=='Function') then
       call krigingfuncpredictGEK(arg%ndim,arg%ntot,arg%X,arg%gtot,arg%pts,arg%dims,arg%stot,arg%H,arg%beta,arg%V,arg%hyper,1,X,Hm,f,arg%covarflag)
    elseif (arg%value=='Variance') then
       call krigingfuncvarianceGEK(arg%ndim,arg%ntot,arg%X,arg%gtot,arg%pts,arg%dims,arg%stot,arg%H,arg%beta,arg%hyper,1,X,Hm,f,arg%covarflag)
    elseif (arg%value=='EI') then
       Ymin=arg%ymin

       call krigingfuncpredictGEK(arg%ndim,arg%ntot,arg%X,arg%gtot,arg%pts,arg%dims,arg%stot,arg%H,arg%beta,arg%V,arg%hyper,1,X,Hm,Yt,arg%covarflag)
       call krigingfuncvarianceGEK(arg%ndim,arg%ntot,arg%X,arg%gtot,arg%pts,arg%dims,arg%stot,arg%H,arg%beta,arg%hyper,1,X,Hm,Sm,arg%covarflag)


       if (Sm<=1.D-12) then
          f=0.D0
       else
          s=sqrt(Sm)
          chi=exp(-(Ymin-Yt)**2/(2.D0*s**2))/sqrt(2.D0*Pi)
          phi=0.5*(1.D0+erf((Ymin-Yt)/sqrt(2.D0*s**2)))
          f=(Ymin-Yt)*phi+s*chi
       end if
    end if

    return 
  end subroutine func_kriging_grad
end module funcmod
