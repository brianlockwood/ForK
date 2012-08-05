!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file covars.f90
!> \brief Module containing the covariance functions required for the Kriging model. Multi-dimensional and the one dimensional used to build up multi-dimensional covariance functions are included. Derivatives of covariance functions also included here.
module covars
  implicit none
  integer covarflag 
  real(8), parameter :: root3=sqrt(3.d0) 
  real(8), parameter :: root5=sqrt(5.d0) 
  real(8), parameter :: third=1.D0/3.D0
contains
!>  \breif Function used to calculate the one dimensional covariance functions used to construct the multidimensional functions.

!>  \detail This function calculates the one dimensional Matern function with three possible values the parameter \f$ \nu \f$ (specified with the variable \a covarflag). The three possible functions are given as:
!>\f[	k(x,x') = \begin{cases}  e^{-r} & \text{for $\nu = 1/2$,} \\\
!>		 \left(1 + \sqrt{3} r \right) e^{-\sqrt{3} r} & \text{for $\nu = 3/2$,} \\\
!>               \left(1 + \sqrt{5} r+ \frac{5}{3} r^{2} \right) e^{-\sqrt{5} r} & \text{for $\nu = 5/2$,}  \end{cases}\f]
!>  where \f$r\f$ is defined as:
!>        \f[ r = \bigg \vert \frac{X - X'}{\theta} \bigg \vert \f]

!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 17, 2012
!> \param(in) <b> X </b>: Coordinate of the first point
!> \param(in) <b> Xt</b>: Coordinate of the second point
!> \param(in) <b> theta</b>: Length scale in the dimension
!> \param(in) <b>covarflag</b>: Flag to govern which covariance function is used (stored as global variable in this module) <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul>
!> \retval covar1d the 1-D Covariance value
    function covar1d(X,Xt,theta)
    implicit none
    real(8), intent(in) :: X,Xt,theta
    real(8) covar1d
    
    real(8) r
    
    r=abs((X-Xt)/theta)
    
    select case(covarflag)
    case(0)
       covar1d=exp(-r)
    case(1)
       covar1d=(1.D0+root3*r)*exp(-root3*r)
    case(2)
       covar1d=(1.D0+root5*r+5.D0*third*r**2)*exp(-root5*r)
    end select

  end function covar1d

!>  \breif Function used to calculate the derivative of the one dimensional covariance functions with respect to the first coordinate.

!>  \detail This function calculates the derivative of the one dimensional Matern function with three possible values the parameter \f$ \nu \f$ (specified with the variable \a covarflag). The derivative is performed with respect to the first coordinate (\a X). This derivative is required for the covariance between function and derivative values:
!>  \f[ cov(y, \frac{d y}{d x}) = \frac{d}{dx} k(X,X',\theta) \f].
 
!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 17, 2012
!> \param(in) <b> X </b>: Coordinate of the first point
!> \param(in) <b> Xt</b>: Coordinate of the second point
!> \param(in) <b> \theta</b>: Length scale in the dimension
!> \param(in) <b>covarflag</b>: Flag to govern which covariance function is used (stored as global variable in this module) <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul>
!> \retval dcovar1d the Derivative of the Covariance w.r.t. X

  function dcovar1d(X,Xt,theta)
    implicit none
    real(8), intent(in) :: X,Xt,theta
    real(8) dcovar1d
    
    real(8) r
    real(8) dr

    r=abs((X-Xt)/theta)
    dr=sign(1.D0,(X-Xt)/theta)*(1.D0/theta)
    
    select case(covarflag)
    case(0)
       dcovar1d=-exp(-r)*dr
    case(1)
       dcovar1d=-3.D0*r*exp(-root3*r)*dr
    case(2)
       dcovar1d=-5.D0*third*r*(1.D0+root5*r)*exp(-root5*r)*dr
    end select
    
  end function dcovar1d

!>  \breif Function used to calculate the second derivative of the one dimensional covariance functions with respect to the second coordinate.

!>  \detail This function calculates the second derivative of the one dimensional Matern function with three possible values the parameter \f$ \nu \f$ (specified with the variable \a covarflag). The second derivative of interest is the derivative of the result of \a dcovar1d with respect to the argument \a Xt. This derivative is required for the covariance between derivative values:
!>  \f[ cov(\frac{d y}{d x}, \frac{d y'}{d x}) = \frac{d^{2}}{dx dx'} k(X,X',\theta) \f]
 
!>  The practical result of differeniating with repsect to the second argument is that an extra negative sign appears. Hence, the second derivative of the Matern functions with respect to \f$r\f$ can be calculated and an additional negative sign can be applied to the result to yield the required derivative value. 

!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 17, 2012
!> \param(in) <b> X </b>: Coordinate of the first point
!> \param(in) <b> Xt</b>: Coordinate of the second point
!> \param(in) <b> \theta</b>: Length scale in the dimension
!> \param(in) <b>covarflag</b>: Flag to govern which covariance function is used (stored as global variable in this module) <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul>
!> \retval d2covar1d The second derivative of the Covariance w.r.t. Xt

  function d2covar1d(X,Xt,theta)
    implicit none
    real(8), intent(in) :: X,Xt,theta
    real(8) d2covar1d
    
    real(8) r
    real(8) dr1,dr2
    real(8) d2r


    r=abs((X-Xt)/theta)
    dr1=sign(1.D0,(X-Xt)/theta)*(1.D0/theta)
    dr2=sign(1.D0,(X-Xt)/theta)*(-1.D0/theta)
    d2r=0.D0
    
    select case(covarflag)
    case(0)
       WRITE(*,*) 'Covariance Function is not twice differentiable, try another'
       stop
    case(1)
       d2covar1d=-3.D0*((1.D0-root3*r)*dr1*dr2+r*d2r)*exp(-root3*r)
    case(2)
       d2covar1d=-5.D0*third*((1.D0+root5*r-5.D0*r**2)*dr2*dr1+r*(1.D0+root5*r)*d2r)*exp(-root5*r)
    end select
       
    
  end function d2covar1d

!> \brief Subroutine for computing the multidimensional covariance between two data points

!> \detail Multidimensional covariance is built up using a tensor product of 1 dimensional covariance functions: 
!> \f[
!>       k(\vec{X},\vec{X}'; \theta) = \prod_{i}^{d} k_{i}\left(\bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert\right)\f]
!>  where \f$D\f$ is the dimension of the problem. The 1-D covariance functions are defined in function covar1d

!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
!>         University of Wyoming
!> \date   May 3, 2012
!> \param(in) <b> ndim </b>: The dimension of the problem
!> \param(in) <b> X </b>:  Coordinates of the first point (size=[ndim])
!> \param(in) <b> Xt </b>:  Coordinates of the second point (size=[ndim])
!> \param(in) <b> theta </b>: Length scale in each dimension
!> \param(out) <b> k </b>: Covariance value
subroutine covarfunc(ndim,X,Xt,theta,k)
  implicit none
  integer, intent(in) :: ndim
  real(8), intent(in) :: X(ndim),Xt(ndim),theta(ndim)

  real(8), intent(out) :: k

  real(8) K1

  integer i

  if (ndim==1) then
     K=covar1d(X(1),Xt(1),theta(1))
  else    
     k=1.D0
     do i=1,ndim
        K1=covar1d(X(i),Xt(i),theta(i))
        K=K*K1
     end do
  end if

  return
end subroutine covarfunc

!> \brief Subroutine for computing the gradient of the multidimensional covariance function.

!> \detail The multidimensional covariance is built up using a tensor product of 1 dimensional covariance functions: 
!> \f[
!>       k(\vec{X},\vec{X}'; \theta) = \prod_{i}^{D} k_{i}\left(\bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right)\f]
!>  where \f$D\f$ is the dimension of the problem. The 1-D covariance functions are defined in function covar1d

!>  The derivative of this covariance with respect to the coordinate \f$ x_{j}\f$ can be written as:
!>  \f[ \frac{\partial}{\partial x_{j}} k(\vec{X},\vec{X}'; \theta) = \frac{d}{dx_{j}} k_{j} \left( \bigg \vert \frac{X_{j} - X'_{j}}{\theta_{j}} \bigg \vert \right) \prod_{i \ne j} k_{i} \left(\bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right)\f]

!>  Although conceptually easy, the above formula for the derivative of the covariance function is inefficient as the complexity is order \f$D^{2}\f$. Using memoization, the gradient can be computed in order \f$D\f$ time by computed and storing several intermediate products. Defining the lower product as:

!>  \f[ L_{j} = \prod_{i<j} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \f]
!>  and the upper product as:
!>  \f[ U_{j} = \prod_{i>j} k_{i}\left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \f], 
!>  the derivative can be computed as:
!>  \f[ \frac{\partial}{\partial x_{j}} = L_{j} \frac{d}{dx_{j}} k_{j} \left( \bigg \vert \frac{X_{j} - X'_{j}}{\theta_{j}} \bigg \vert \right) U_{j} \f]

!>  By storing the 1-D covariance function values, the upper and lower products can be built up recursively in order \f$D\f$ time. The use of these upper and lower products allows for the re-use of work for each component of the gradient, reducing the cost of constructing the entire gradient.

!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
!>         University of Wyoming
!> \date   May 3, 2012
!> \param(in) <b> ndim </b>: The dimension of the problem
!> \param(in) <b> X </b>:  Coordinates of the first point (size=[ndim])
!> \param(in) <b> Xt </b>:  Coordinates of the second point (size=[ndim])
!> \param(in) <b> theta </b>: Length scale in each dimension
!> \param(out) <b> dk </b>: Covariance value

subroutine dcovarfunc(ndim,X,Xt,theta,dk)
  implicit none
  integer, intent(in) :: ndim
  real(8), intent(in) :: X(ndim),Xt(ndim),theta(ndim)

  real(8), intent(out) :: dk(ndim)

  integer i

  real(8) L(ndim),U(ndim),Fs(ndim)

  real(8) dkr

  if (ndim==1) then
     dK(1)=dcovar1d(X(1),Xt(1),theta(1))
  else
     do i=1,ndim
        Fs(i)=covar1d(X(i),Xt(i),theta(i))
     end do

     L(1)=1.D0
     do i=1,ndim-1
        L(i+1)=L(i)*Fs(i)
     end do

     U(ndim)=1.d0
     do i=ndim,2,-1
        U(i-1)=U(i)*Fs(i)
     end do

     do i=1,ndim
        dKr=dcovar1d(X(i),Xt(i),theta(i))
        dK(i)=L(i)*dKr*U(i)
     end do
  end if

  return
end subroutine dcovarfunc

!> \brief Subroutine for computing the Hessian of the multidimensional covariance function. The second derivative is performed with respect to the second coordinate.

!> \detail The multidimensional covariance is built up using a tensor product of 1 dimensional covariance functions: 
!> \f[
!>       k(\vec{X},\vec{X}'; \theta) = \prod_{i}^{D} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \f]
!>  where \f$D\f$ is the dimension of the problem. The 1-D covariance functions are defined in function covar1d

!>  The derivative of this covariance with respect to the coordinate \f$ x_{j}\f$ can be written as:
!>  \f[ \frac{\partial}{\partial x_{j}} k(\vec{X},\vec{X}'; \theta) = \frac{d}{dx_{j}} k_{j} \left( \bigg \vert \frac{X_{j} - X'_{j}}{\theta_{j}} \bigg \vert \right) \prod_{i \ne j} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right)\f]

!>  The second derivative of the covariance with repect to the \f$l\f$ component of the second coordinate \f$ x'_{l}\f$ is given as:

!>  \f[ \frac{\partial^{2}}{\partial x_{j} \partial x'_{l}} k(\vec{X},\vec{X}'; \theta) = \begin{cases} \frac{d^{2}}{dx_{j} dx'_{j}} k_{j} \left( \bigg \vert \frac{X_{j} - X'_{j}}{\theta_{j}} \bigg \vert \right) \prod_{i \ne j} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \text{   for $l=j$, } \\\
!> -\frac{d}{dx_{j}} k_{j} \left( \bigg \vert \frac{X_{j} - X'_{j}}{\theta_{j}} \bigg \vert \right) \frac{d}{dx_{l}} k_{l} \left( \bigg \vert \frac{X_{l} - X'_{l}}{\theta_{l}} \bigg \vert \right) \prod_{i \ne j \ne l} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \text{   for $l \ne j$}
!> \end{cases} \f]
!!$
!>  The entire Hessian can be computed efficiently by storing the following intermediate products as well as the 1-D covariance function values. Define the lower product as:
!>  \f[ L_{j} = \prod_{i<j} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \f],
!>  the upper product as:
!>  \f[ U_{j} = \prod_{i>j} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \f]
!>  and the intermediate product as:
!>  \f[ M_{l,j} = \prod_{l<i<j} k_{i} \left( \bigg \vert \frac{X_{i} - X'_{i}}{\theta_{i}} \bigg \vert \right) \f],
!>  the hessian of the covariance can be computed as:
!>  \f[ \frac{\partial^{2}}{\partial x_{j} \partial x'_{l}} k(\vec{X},\vec{X}'; \theta) = \begin{cases} L_{j} \frac{d^{2}}{dx_{j} dx'_{j}} k_{j} \left(\bigg \vert \frac{X_{j} - X'_{j}}{\theta_{j}} \bigg \vert \right) U_{j} \text{   for $l=j$, } \\\
!> -L_{j} \frac{d}{dx_{j}} k_{j} \left( \bigg \vert \frac{X_{j} - X'_{j}}{\theta_{j}} \bigg \vert \right) M{j,l} \frac{d}{dx_{l}} k_{l} \left( \bigg \vert \frac{X_{l} - X'_{l}}{\theta_{l}} \bigg \vert \right) U_{l} \text{   for $l \ne j$}
!> \end{cases} \f].

!>  The lower, upper and intermediate products can be built up recursively based on the stored 1-D covariance function values. This allows the full Hessian to be calculated in order \f$D^{2}\f$ time. For the intermediate product, the diagonal and sub-diagonal is defined as:
!>  \f[\begin{cases} M_{j,j} = 1 & \text{for $j=1,..,D$, } \\\
!>                   M_{j,j+1} = 1 & \text{for $j=1,..,D-1$, } \end{cases} \f]

!>  Further cost savings are realized by utilizing the fact that the Hessian is symmetric. Hence, only the upper diagonal components of the Hessian are analytically computed and the result of this calculation is copied to the lower diagonal components.

!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
!>         University of Wyoming
!> \date   May 3, 2012
!> \param(in) <b> ndim </b>: The dimension of the problem
!> \param(in) <b> X </b>:  Coordinates of the first point (size=[ndim])
!> \param(in) <b> Xt </b>:  Coordinates of the second point (size=[ndim])
!> \param(in) <b> theta </b>: Length scale in each dimension
!> \param(out) <b> dk </b>: Covariance value

subroutine d2covarfunc(ndim,X,Xt,theta,d2k)
  implicit none
  integer, intent(in) :: ndim
  real(8), intent(in) :: X(ndim),Xt(ndim),theta(ndim)

  real(8), intent(out) :: d2k(ndim,ndim)

  integer i,j

  real(8) Fs(ndim),Gs(ndim),Hs(ndim)

  real(8) L(ndim),U(ndim),M(ndim,ndim)
  
  real(8) dkr,d2kr

  if (ndim==1) then
     d2k=d2covar1d(X(1),Xt(1),theta(1))
  else
     do i=1,ndim
        Fs(i)=covar1d(X(i),Xt(i),theta(i))
     end do
     do i=1,ndim
        Gs(i)=dcovar1d(X(i),Xt(i),theta(i))
     end do
     do i=1,ndim
        Hs(i)=d2covar1d(X(i),Xt(i),theta(i))
     end do

     L(1)=1.D0
     do i=1,ndim-1
        L(i+1)=L(i)*Fs(i)
     end do

     U(ndim)=1.d0
     do i=ndim,2,-1
        U(i-1)=U(i)*Fs(i)
     end do

     do i=1,ndim-1
        M(i,i+1)=1.D0
        do j=i+2,ndim
           M(i,j)=M(i,j-1)*Fs(j-1)
        end do
     end do

     do i=1,ndim
        d2k(i,i)=L(i)*Hs(i)*U(i)
     end do

     do i=1,ndim
        do j=i+1,ndim
           d2k(i,j)=-L(i)*Gs(i)*M(i,j)*Gs(j)*U(j)
           d2k(j,i)=d2k(i,j)
        end do
     end do
  end if

  return
end subroutine d2covarfunc
end module covars
