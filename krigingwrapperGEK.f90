!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingwrapperGEK.f90
!> \brief Example Code for building a gradient-enhanced Kriging model and making predictions from it.
program krigingwrapper
  implicit none
  integer, parameter :: ndim=1 ! Dimension of the Space
  integer, parameter :: ntot=100 ! Number of Training Points
  integer, parameter :: stot=1 ! Number of Terms in the regression (stot=1 gives constnt mean function)
  integer, parameter :: mtot=1000 ! Number of Test points
  integer, parameter :: gtot=ndim*ntot ! Number of derivative values in the training data set
  integer, parameter :: htot=ndim*mtot ! Number of derivative values predicted in the test set

  ! Training data
  real(8) X(ndim,ntot),Y(ntot),dY(gtot)
  integer pts(gtot),dims(gtot)
  real(8) dYt(ndim)
  
  !  Collocation Matrix 
  real(8) H(stot,ntot+gtot)

  !  Processed Kriging Data
  real(8) V(ntot+gtot),beta(stot),hyper(ndim+3)

  integer i,k,j,l

  integer ierror

  !  Test Data
  real(8) Xm(ndim,mtot),Ym(mtot),Sm(mtot)
  real(8) Hm(stot,mtot),Gm(stot,htot)
  real(8) dYm(htot), dSm(htot)
  integer ptsm(htot), dimsm(htot)
  
  !  Misc Variables
  real(8) Xl(ndim),Xf(ndim)
  real(8) Fmin

  real(8) rms, Ye,rmsg

  real(8) Lb(ndim), Ub(ndim)
  real(8) Ymax
  real(8) EI(mtot)

  real(8) s,chi,phi,Pi,Ymin

  integer covarflag

  real(8) U(mtot)
  real(8) Sigma(mtot,mtot)
  real(8) Yt(mtot)

  ! Covariance Flag: Matern function with nu = 3/2
  covarflag=1
  
  PI=4.D0*ATAN(1.D0)

  WRITE(*,*) 'Generate Sample Locations'
  ! Random points in the interval [-5.12,5.12] are used 
  do i=1,ntot
     do j=1,ndim
        call random_number(X(j,i))
        X(j,i)=-5.12+10.24*X(j,i)
     end do
  end do
  
  !  Evaluate the Training data. Training data includes function and gradient values for "toyfunc"
  WRITE(*,*) 'Get Training Data'
  !  Test function/gradient is evaluated at each sample point and the linked list that describes the location of the gradient are filled. 
  k=0
  do i=1,ntot
     call toyfunc(ndim,X(:,i),Y(i),dYt(:))
     do j=1,ndim
        k=k+1
        pts(k)=i
        dims(k)=j
        dY(k)=dYt(j)
     end do   
  end do
  H=0.d0
  do i=1,ntot
     H(1,i)=1.D0
  end do

  WRITE(*,*) 'Build Kriging'
  !  Build Gradient-enhanced Kriging model to generate beta, V and Hyper
  !  Uses Pattern search for the optimization and the MLE likelihood formula
  call buildkrigingGEK(ndim,ntot,X,Y,gtot,pts,dims,dY,stot,H,beta,V,hyper,0,1,covarflag)

  WRITE(*,*) 'Hyper Parameters:'
  WRITE(*,*) 'Length scales:', hyper(1:ndim)
  WRITE(*,*) 'Covariance Magnitude:', hyper(ndim+1)
  WRITE(*,*) 'Noise Level:', hyper(ndim+2), hyper(ndim+3)



  ! Generate Test point locations
  k=0
  do i=1,mtot
     do j=1,ndim
        call random_number(Xm(j,i))
        Xm(j,i)=-5.12+10.24*Xm(j,i)
        k=k+1
        ptsm(k)=i
        dimsm(k)=j        
     end do
  end do
  
  do i=1,mtot
     Hm(1,i)=1.D0
  end do
  Gm=0.D0
  
  WRITE(*,*) 'Predict Function Value'
  !  Predict Output from Gradient-enhanced Kriging model
  call krigingfuncpredictGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,V,hyper,mtot,Xm,Hm,Ym,covarflag)

  !  Predict Variance associated with function values from Gradient-enhanced Kriging model
  WRITE(*,*) 'Predict Function Variance'
  call krigingfuncvarianceGEK(ndim,ntot,X,gtot,pts,dims,stot,H,hyper,mtot,Xm,Hm,Sm,covarflag)

  !  Predict derivative values associated from GE Kriging model
  WRITE(*,*) 'Predict Gradient Value'
  call kriginggradpredictGEK(ndim,ntot,X,gtot,pts,dims,stot,H,beta,V,hyper,mtot,Xm,htot,ptsm,dimsm,Gm,dYm,covarflag)  

  !  Predict variance of derivative values from GE Kriging model
  WRITE(*,*) 'Predict Gradient Variance'
  call kriginggradvarianceGEK(ndim,ntot,X,gtot,pts,dims,stot,H,hyper,mtot,Xm,htot,ptsm,dimsm,Gm,dSm,covarflag)  

  ! Output Training data
  open(unit=1,file='train.dat',status='replace',iostat=ierror)
  do i=1,ntot
     !          Location, Function value, Function variance, Derivative, Derivative Variance
     WRITE(1,*) X(1,i), Y(i), 0.D0, dY(i), 0.D0
  end do

  ! Output Test data
  open(unit=1,file='test.dat',status='replace',iostat=ierror)
  do i=1,mtot
     !          Location, Function value, Function variance, Derivative, Derivative Variance
     WRITE(1,*) Xm(1,i), Ym(i), Sm(i), dYm(i), dSm(i)
  end do
99 format(99(F25.15))



  stop
end program krigingwrapper
