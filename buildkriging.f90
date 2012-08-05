!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file buildkriging.f90
!> \brief Subroutine to Build Kriging Surface based on Function Values

!> \detail This subroutine creates a Kriging surface based on the supplied training data \a X and \a Y, with the parameters H, hyperflag, optflagi, covarflagi using to specify details of the construction process. The creation of the Kriging model is governed by the following steps:
!>  <ol> <li> The fitting of model parameters (hyperparameters) based on maximization of the likelihood equation
!>       <li> Creation of the Covariance Matrix
!>       <li> Determination of the Optimal regression parameters 
!>       <li> Creation of the processed data \a V by inverting the covariance matrix and multiplying the difference between the training data and regression value by the inverse of the covariance matrix </ol>
!>  The Kriging model is built on the assumption that the data \a Y obey a Gaussian process with an assumed form for the mean function and the covariance between data points:
!>  \f[
!>          Y = N(m(\vec{x}), K(\vec{x},\vec{x}))\f]
!>  where \f$m(\vec{x})\f$ is the mean function and \f$K(\vec{x},\vec{x})\f$ represents the covariance between function values. For this work, a regression mean function is assumed. Using this form, the mean function has the following form:
!> \f[ m(x) = h^{T}(\vec{x}) \beta \f]
!> where \f$h^{T}(\vec{x})\f$ represents a column vector containing the basis functions of the basis evaluated at the points \f$\vec{x}\f$. The regression parameters \f$\beta\f$ are treated as part of the Kriging model and are determined while constructioning the model. Using this form of the mean function yields a Universal Kriging model. The case of a \f$p=0\f$ regression (where the vector \f$h(\vec{x})\f$ reduces to unity) is referred to as ordinary Kriging and is also covered by this functional form. The assumption of a vague prior on the regression parameters gives the following closed form for the parameters:
!>\f[
!>   \beta=(H K^{-1} H^{T})^{-1} H^T K^{-1} Y = A^{-1} H^T K^{-1} Y
!>\f]
!> where \f$K\f$ is the covariance matrix between the training data. For a Kriging model, the covariance between function values is assumed to be only a function of the distance between points. The multi-dimension covariance function is constructed using a tensor product of one dimension functions. The multi-dimension covariance is calculated in subroutine <a href=namespacecovars.html#adb9de23578a8355229466a5dcd8e26e6covarfuncs> covarfunc</a>. The elements in ths covariance matrix are given as:
!>  \f[
!>      K_{i,j} = cov(y_{i},y_{j}) = \sigma^{2} k(\vec{X}_{i},\vec{X}_{j}; \theta) + \sigma^{2}_{n} \delta_{i,j}
!>  \f]
!>  The parameters \f$\sigma\f$ and \f$\theta\f$ (and in some cases \f$\sigma_{n}\f$) are denoted as hyperparameters and are determined maximizing the likelihood equation for the Kriging model. This likelihood gives the probability that a Gaussian process with specified hyperparamters describes the training data \a X and \a Y. By picking the hyperparameters that maximize this probability, a Kriging model that best describes the data can be constructed.
!>  The hyperparameters can be determined in two different ways in this code. The first way is based on the Mean square error estimate and is computed in <a href=likelihood__mle_8f90.html>likelihood_mle</a>. For this method, the noise level \f$\sigma_{n}\f$ is prescribed to ensure proper conditioning of the covariance matrix and the optimal covariance magnitude \f$\sigma\f$ is determined in closed form. Hence, only the length scales \f$\theta\f$ are determined through numerical optimization. This optimization is performed using a <a href=simplexsearch_8f90.html>simplex search</a> or <a href=patternsearch_8f90.html>pattern search</a>. 

!> The second way of determining hyperparameters is based on the likelihood equation for a gaussian process with a vague prior on the regression parameters. This likelihood is computed in  <a href=likelihood_8f90.html>likelihood</a>. Using this equation, optimization is used to determine all of the parameters, including the covariance magnitude \f$\sigma\f$ and noise \f$\sigma_{n}\f$. This way of determining hyperparameters should be used when the noise level of the function needs to be fitted. <br>

!> With the regression and covariance parameters determined, the final processed data can be constructed using the inverse of the covariance matrix. To make predictions from the Kriging model, the following vector is required:
!> \f[ V = K^{-1} (Y - H^{T} \beta) \f]
!> where \f$ K \f$ is the covariance matrix, the product \f$H^{T} \beta\f$ represents the mean function evaluated at the training points and \f$Y\f$ represents the function values at the training points. Using this processed data, the regression parameters and covariance parameters, predictions can be made based on the Kriging model using: 
!><ul><li> <a href="krigingfuncpredict_8f90.html">krigingfuncpredict</a> 
!><li><a href="krigingfuncvariance_8f90.html">krigingfuncvariance</a>
!><li><a href="kriginggradpredict_8f90.html">kriginggradpredict</a>
!><li> <a href="kriginggradvariance_8f90.html">kriginggradvariance</a> </ul>

!>
!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 1, 2012
!> \param(in) <b>ndim</b> The dimension of the problem
!> \param(in) <b>ntot </b>: The number of Training points
!> \param(in) <b>X </b>:  The location of the training points (size=[ndimxntot])
!> \param(in) <b>Y </b>:  Funcation values for the training points (size=[ntot])
!> \param(in) <b>stot </b>: Number of Terms in the regression 
!> \param(in) <b>H</b>: The collocation matrix for the regression (size=[stotxntot]) 
!> \param(out) <b> beta</b>: Regression coefficients based on the optimal estimate for the Kriging model (size=[stot]) 
!> \param(out) <b>hyper</b>: Hyperparameters for the Kriging Model (size=[ndim+2]) <br>
!>                   <ul><li>hyper(1:ndim) correspond to the length scale used in each dimension (\f$\theta\f$ in the covariance functions) 
!>                   <li>hyper(ndim+1) is the magnitude of the covariance matrix 
!>                   <li>hyper(ndim+2) is the fitted noise for the function evaluations </ul>
!> \param(out) <b>V</b>:     Processed Training Data (size=[ntot]) <br>
!>                   In order to make predictions from the Kriging Model, the inverse of the covariance matrix onto the 
!>                   residual of the training data from the regression is all that is needed. To avoid re-computing and inverting the covariance matrix, the product is stored explicitly
!> \param(in) <b>hyperflagi</b>: Flag to govern how the hyperparameters are selected <br>
!>            <ul><li>hyperflag=0 uses the likelihood formula to determine the length scale only <br>
!>                        The noise in the Function is hard-coded as a small value simply to ensure the covariance matrix is properly conditioned <br>
!>                        The magnitude for the covariance function is determined explicitly based on optimality for the likelihood function (Solve for magnitude when derivative of likelihood is set to zero) 
!>            <li>hyperflag=1 uses the likelihood formula based on all hyperparameters with no parameters determined by explicit relation <br>
!>                        All hyperparameters including noise and magnitude are determined via the optimization. Required when the amount of noise in the function evaluations is to be fitted. However, the dimension of the optimization problem is now higher so this will be slower.</ul>
!>  \param(in) <b>optflagi</b>: Flag to govern the optimization algorithm used to determine the hyperparameters <br>
!>                <ul><li>optflag=0  Simplex Search (Nelder-Mead Method) to determine the hyperparameters <br>
!>                           Local optimization technique using simplex elements to determine search direction and line search to determine step sizes <br>
!>                           Fastest method but inherently sequential and may stall on local optimum
!>                <li>optflag=1  Pattern Search used to determine hyperparameters <br>
!>                           Global optimization method with fixed search directions (forward and backward in each direction) <br>
!>                           Each iteration requires \f$2*ndim\f$ likelihood evaluations; however, these may be performed in parallel using openmp (Each thread performs own likelihood evaluation so OMP_STACK_SIZE must be set large enough) <br> </ul>
!>   \param(in) <b>covarflagi</b>: Flag to govern which covariance function is used <br>
!>              <ul><li>covarflag==0 Uses Matern function with \f$\nu=1/2\f$
!>              <li>covarflag==1 Uses Matern function with \f$\nu=3/2\f$
!>              <li>covarflag==2 Uses Matern function with \f$\nu=5/2\f$ </ul> <br>
!>              The parameter \f$\nu\f$ governs the smoothness and differentiability of the covariance function. For function only Kriging, all three options are available. <br>
!>              Good rule of thumb is to use the least smooth (\f$\nu=1/2\f$) unless there is a reason to assume more smoothness in the data <br>
!>              For small numbers of training points, higher \f$\nu\f$ can improve accuracy
!
!>  \retval beta Regression coefficients based on the optimal estimate for the Kriging model (size=[stot])
!>  \retval hyper Hyperparameters for the Kriging Model (size=[ndim+2])
!>  \retval V Processed Training Data (size=[ntot])

subroutine buildkriging(ndim,ntot,X,Y,stot,H,beta,V,hyper,hyperflag,optflagi,covarflagi)
  use choleskymod
  use opt,only: optflag
  use covars,only: covarflag
  use hyperparameters_mle_mod
  use hyperparameters_all_mod
  use covarmatrix_mod
  implicit none
  integer, intent(in) :: ndim, ntot
  real(8), intent(in) :: X(ndim,ntot),Y(ntot)
  integer, intent(in) :: stot
  real(8), intent(in) :: H(stot,ntot)
  
  real(8), intent(out) :: V(ntot),hyper(ndim+2),Beta(stot)

  integer, intent(in) :: optflagi, hyperflag, covarflagi

  real(8) Kmat(ntot,ntot)
  real(8) theta(ndim),sigma,sigmaN
 
  integer i,j,k,l

  real(8) A(stot,stot)

  real(8) W(ntot),B(ntot),Q(stot)

  real(8) logpy

  real(8) C(stot,ntot)

  optflag=optflagi
  covarflag=covarflagi

! Determine the Hyperparameters 
  if (hyperflag==0) then
     call hyperparameters_mle(ndim,ntot,X,Y,stot,H,theta,sigma,sigmaN,logpy) 
  else
     call hyperparameters_all(ndim,ntot,X,Y,stot,H,theta,sigma,sigmaN,logpy) 
  end if

!  Calculate Covariance Matrix 
  call covarmatrix(ndim,ntot,X,theta,Kmat)

!  Add in Noise and magnitude  
  Kmat=sigma**2*Kmat
  do i=1,ntot
     Kmat(i,i)=Kmat(i,i)+sigmaN**2
  end do

!  Cholesky and Invert
  call invertsymmetric(ntot,Kmat)

!  Construct Regression Matrix
  call symmatrixmulttrans(ntot,Kmat,stot,H,C)
  call matrixmulttrans(stot,ntot,H,stot,C,A)

!  Cholesky Regression Matrix and solve for Beta
  call symmatvec(ntot,Kmat,Y,W)
  call matvec(stot,ntot,H,W,Q)
  Beta=Q
  call symmetricsolve(stot,A,Beta)

  !  Compute processed training data
  call matvectrans(stot,ntot,C,Beta,B)

  do i=1,ntot
     V(i)=W(i)-B(i)
  end do

  hyper(1:ndim)=theta(:)
  hyper(ndim+1)=sigma
  hyper(ndim+2)=sigmaN

  return
end subroutine buildkriging

