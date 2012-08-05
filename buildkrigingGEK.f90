!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file buildkrigingGEK.f90
!> \brief Subroutine to Build Kriging Surface based on Function and Gradient Values

!> \detail This subroutine creates a gradient-enhanced Kriging model based on Function and Gradient values supplied at the Training points. The process of creating a gradient-enhanced Kriging model is identical to that of creating a function-only model (described in <a href=buildkriging_8f90.html>buildkriging</a>). The details of creating the Kriging model are not repeated here. To create a gradient-enhanced model, the definitions for the training data, covariance matrix, and collocation matrix must be extended to include the derivative values. Using these modified defintions, the procedure in <a href=buildkriging_8f90.html>buildkriging</a> can be used to create the Kriging model. First, the vector of training data is redefined as a block vector containing function and derivative observations, given as:
!> \f[
!>	\underline{Y} = \left[\begin{array}{c} Y \\\ \delta Y \end{array} \right]
!>    \f]
!> where the vector \f$Y\f$ contains the function observations at the training points and the vector \f$\delta Y\f$ contains the derivative observations for the training points. Based on training data in this form, the covariance matrix must also be extended to include the correlation between derivative values. The block definition of this covariance matrix is given as:
!> \f[
!>    \underline{K} = \left[ \begin{array}{cc}
!>                cov(Y,Y) & cov(Y,\delta Y) \\\
!>                cov(\delta Y,Y) & cov(\delta Y,\delta Y)
!>               \end{array} \right] \f]
!>  where \f$cov(Y,Y)\f$ is the covariance matrix between function values (same as the covariance matrix used in a function-only model), \f$cov(\delta Y, Y) \f$ is a rectangular matrix whose elements represent the covariance between the function values and derivative values at the training points and \f$cov(\delta Y, \delta Y) \f$ represents the covariance between derivative observations. 

!> For a function value at point \f$ \vec{x}' \f$ and a derivative value with respect to \f$x_{k}\f$ evaluated at point \f$\vec{x}\f$, the covariance is found by differentiation of the function describing the covariance between function values. The expression for this covariance between function and derivative values is given as:
!>  \f[
!>      cov(\frac{\partial y}{\partial x_{k}},y')=\frac{\partial}{\partial x_{k}}  k(\vec{x},\vec{x}') 
!>  \f]
!>  For the covariance between a derivative value with respect to \f$x_{k}\f$ evaluated at point \f$\vec{x}\f$ and a derivative value with respect to \f$x_{l}\f$ evaluated at point \f$\vec{x}'\f$, the covariance is found by differentiation of the above expression:
!>   \f[   cov(\frac{\partial y}{\partial x_{k}},\frac{\partial y'}{\partial x'_{l}})=\frac{\partial^2}{\partial x_{k}\partial x'_{l}}  k(\vec{x},\vec{x}') \f]

!>  Hence, to construct a gradient-enhanced Kriging model, the covariance function must be twice differentiable. This limits the choice of covariance function in this work to the Matern function with \f$\nu \geq 3/2\f$. These covariance function and derivatives are computed using the functions in module <a href=covars_8f90.html>covars</a>. 

!>  Finally, the collocation matrix must be extend to include the derivative of the basis functions. The block collocation matrix can be written as:
!>  \f[ \underline{H} = \left[\begin{array}{c} H \\\ G \end{array} \right] \f]
!>  where \f$H\f$ represents the function collocation matrix and \f$G\f$ represents the derivative collocation matrix. The elements of the function collocation matrix are defined as the basis functions for the regression evaluated at the training points. The elements of the derivative collocation matrix are defined as the derivatives of the basis functions for the regression evaluated at the training points. Using these definitions, the mean function values at the training points are defined as:
!> \f[ \bar{Y} = H \beta\f]
!> while the mean derivative values are given as: 
!> \f[ \delta \bar{Y} = G \beta\f]
 
!> For a gradient-enhanced model, the size of the covariance matrix grows as the dimension of the space is increased (size=[(ndim+1)ntot x (ndim+1)ntot] if all components of the gradient are included at every point); hence, the expense of constructing the gradient-enhanced Kriging model grows quickly with dimension. LAPACK has been utilized extensively to help alleviate this cost in large dimension but it is still likely too slow to construct a gradient-enhanced model in more than \f$ \approx 20 \f$ dimensions. 

!> The derivative values themselves are inputted into the subroutine as a type of sparse matrix. Hence, three vectors are supplied to enter the derivative values (\a pts, \a dims, and \a dY). The elements of the vector \a pts give the index of the training point that the derivative value is evaluated at. The elements of the vector \a dims represent the component of the gradient that the derivative value represents and the vector \a dY contains the derivative values themselves. This data structure allows only certain components of the derivative to be enforced at different points. It also allows for only function values to be supplied at certain points. 

!> The final processed data produced by this subroutine is defined as:
!> \f[ \underline{V} = \underline{K}^{-1} \left( \underline{Y} - \underline{H}^{T} \beta) \f]


!> \author Brian Lockwood 
!>         Department of Mechanical Engineering
!>         University of Wyoming
!> \date   May 1, 2012
!> \param(in) <b>ndim </b>: The dimension of the problem
!> \param(in) <b>ntot </b>: The number of training points (also function values)
!> \param(in) <b>X </b>:  The location of the training points (size=[ndimxntot])
!> \param(in) <b>Y </b>:  Funcation values for the training points (size=[ntot])
!> \param(in) <b>gtot</b>: Number of derivative values included in training data (ndim*ntot if all derivatives are included at the training points)
!> \param(in) <b>pts</b>:  List identifying what point the derivative value is enforced at (size=[gtot] with values ranging from 1 to ntot)
!> \param(in) <b>dims</b>:  List identifying the dimension the derivative is taken with respect to (size=[gtot] with values ranging from 1 to ndim)
!> \param(in) <b>dY</b>:  Derivative values (size=[gtot])
!> \param(in) <b>stot </b>: Number of Terms in the regression 
!> \param(in) <b>H</b>: The collocation matrix for the regression including derivative values. (size=[stotxntot+gtot]) <br>
!>                      Columns 1:ntot are the basis evaluated at the training points <br>
!>                      Columns ntot+1:ntot+gtot are the derivative of the basis evaluated at the training points 
!> \param(out) <b> beta</b>: Regression coefficients based on the optimal estimate for the Kriging model (size=[stot]) 
!> \param(out) <b>hyper</b>: Hyperparameters for the Kriging Model (size=[ndim+3]) <br>
!>                   <ul><li>hyper(1:ndim) correspond to the length scale used in each dimension (\f$\theta\f$ in the covariance functions) 
!>                   <li>hyper(ndim+1) is the magnitude of the covariance matrix 
!>                   <li>hyper(ndim+2) is the fitted noise for the function evaluations
!>                   <li>hyper(ndim+3) is the fitted noise for the gradient evaluations </ul>
!> \param(out) <b>V</b>:     Processed Training Data (size=[ntot+gtot]) <br>
!>                   In order to make predictions from the Kriging Model, the inverse of the covariance matrix onto the 
!>                   residual of the training data from the regression is all that is needed. To avoid re-computing and inverting the covariance matrix, the product is stored explicitly
!> \param(in) <b>hyperflagi</b>: Flag to govern how the hyperparameters are selected <br>
!>            <ul><li>hyperflag=0 uses the likelihood formula to determine the length scale only <br>
!>                        The noise in the Function is hard-coded as a small value simply to ensure the covariance matrix is properly conditioned <br>
!>                        The magnitude for the covariance function is determined explicitly based on optimality for the likelihood function (Solve for magnitude when derivative of likelihood is set to zero) 
!>            <li>hyperflag=1 uses the likelihood formula based on all hyperparameters with no parameters determined by explicit relation <br>
!>                        All hyperparameters included noise and magnitude are determined via the optimization. Required when the amount of noise in the function evaluations is to be fitted. However, the dimension of the optimization problem is now higher so this will be slower.</ul>
!>  \param(in) <b>optflagi</b>: Flag to govern the optimization algorithm used to determine the hyperparameters <br>
!>                <ul><li>optflag=0  Simplex Search (Nelder-Mead Method) to determine the hyperparameters <br>
!>                           Local optimization technique using simplex elements to determine search direction and line search to determine step sizes <br>
!>                           Fastest method but inherently sequential and may stall on local optimum
!>                <li>optflag=1  Pattern Search used to determine hyperparameters <br>
!>                           Global optimization method with fixed search directions (forward and backward in each direction) <br>
!>                           Each iteration requires \f$2*ndim\f$ likelihood evaluations; however, these may be performed in parallel using openmp (Each thread performs own likelihood evaluation so OMP_STACK_SIZE must be set large enough) <br> </ul>
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul>
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. <br>
!>              When using gradient values, \f$\nu=1/2\f$ is not differentiable enough so \f$\nu \geq 3/2\f$ must be used<br>
!>              Good rule of thumb is to use the least smooth (\f$\nu=3/2\f$) unless there is a reason to assume more smoothness in the data <br>
!>              For small numbers of training points, higher \f$\nu\f$ can improve accuracy
subroutine buildkrigingGEK(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,beta,V,hyper,hyperflag,optflagi,covarflagi)
  use choleskymod
  use opt,only: optflag
  use covars,only: covarflag
  use hyperparameters_mle_mod
  use hyperparameters_all_mod
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: gtot
  integer, intent(in) :: pts(gtot), dims(gtot)
  real(8), intent(in) :: dY(gtot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot+gtot)
  
  real(8), intent(out) :: V(ntot+gtot),hyper(ndim+3),Beta(stot)

  integer, intent(in) :: optflagi, hyperflag, covarflagi

  real(8) Kmat(ntot+gtot,ntot+gtot)
  real(8) theta(ndim),sigma,sigmaN,sigmaNG

  real(8) B(ntot+gtot),Q(stot)

  real(8) Z(ntot+gtot)
  
  integer i,j,k,l

  real(8) A(stot,stot)
  real(8) C(stot,ntot+gtot)

  real(8) W(ntot+gtot)

  real(8) logpy

  optflag=optflagi
  covarflag=covarflagi

!  Training Data
  do i=1,ntot
     Z(i)=Y(i)
  end do
  do i=ntot+1,ntot+gtot
     Z(i)=dY(i-ntot)
  end do

!  Find Hyper Parameters
  if (hyperflag==0) then
     call hyperparameters_mle(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,sigma,sigmaN,sigmaNG,logpy)
  else
     call hyperparameters_all(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,theta,sigma,sigmaN,sigmaNG,logpy)
  end if

!  Calculate Covariance Matrix
  call covarmatrix(ndim,ntot,gtot,pts,dims,X,theta,Kmat)

!  Add in Noise and magnitude  
  Kmat=sigma**2*Kmat
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+sigmaN**2
  end do
  do i=1,gtot
     Kmat(ntot+i,ntot+i)=Kmat(ntot+i,ntot+i)+sigmaNG**2
  end do

!  Cholesky and Invert
  call invertsymmetric(ntot+gtot,Kmat)

!  Invert Covariance Matrix
  call symmatrixmulttrans(ntot+gtot,Kmat,stot,H,C)
  call matrixmulttrans(stot,ntot+gtot,H,stot,C,A)

  call symmatvec(ntot+gtot,Kmat,Z,W)
  call matvec(stot,ntot+gtot,H,W,Q)
  Beta=Q
  call symmetricsolve(stot,A,Beta)

 !  Compute processed training data
  call matvectrans(stot,ntot+gtot,C,Beta,B)

  do i=1,ntot+gtot
     V(i)=W(i)-B(i)
  end do

  hyper(1:ndim)=theta(:)
  hyper(ndim+1)=sigma
  hyper(ndim+2)=sigmaN
  hyper(ndim+3)=sigmaNG 

  return
end subroutine buildkrigingGEK
