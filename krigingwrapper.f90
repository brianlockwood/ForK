!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingwrapper.f90
!> \brief Example Code for building a function-only Kriging model and making predictions from it.
program krigingwrapper
  implicit none
  integer, parameter :: ndim=1 ! Dimension of the Space
  integer, parameter :: ntot=100 ! Number of Training Points
  integer, parameter :: stot=1 ! Number of Terms in the regression (stot=1 gives constnt mean function)
  integer, parameter :: mtot=1000 ! Number of Test points
  integer, parameter :: htot=mtot*ndim ! Number of derivatives predicted at test points

  ! Training data
  real(8) X(ndim,ntot),Y(ntot)
  
  !  Collocation Matrix 
  real(8) H(stot,ntot)

  !  Processed Kriging Data
  real(8) V(ntot),beta(stot),hyper(ndim+2)

  integer i,k,j,l

  integer ierror

  !  Test Data
  real(8) Xm(ndim,mtot),Ym(mtot),Sm(mtot)
  real(8) Hm(stot,mtot)
  real(8) Gm(stot,htot)
  real(8) dYm(htot), dSm(htot)
  integer ptsm(htot), dimsm(htot)
  
  !  Misc Variables
  real(8) dYt(ndim)
  
  integer covarflag
  real(8) Pi

  ! Covariance Flag: Matern function with nu = 3/2
  covarflag=1
  
  PI=4.D0*ATAN(1.D0)

  WRITE(*,*) 'Generate Sample Locations'
  do i=1,ntot
     do j=1,ndim
        call random_number(X(j,i))
        X(j,i)=-5.12+10.24*X(j,i)
     end do
  end do
  
  !  Evaluate the Training data. Training data includes function and gradient values for "toyfunc"
  WRITE(*,*) 'Get Training Data'
  k=0
  do i=1,ntot
     call toyfunc(ndim,X(:,i),Y(i),dYt(:))
  end do
  H=0.d0
  do i=1,ntot
     H(1,i)=1.D0
  end do

  WRITE(*,*) 'Build Kriging'
  !  Build Gradient-enhanced Kriging model to generate beta, V and Hyper
  !  Uses Pattern search for the optimization and the MLE likelihood formula
  call buildkriging(ndim,ntot,X,Y,stot,H,beta,V,hyper,0,1,covarflag)

  WRITE(*,*) 'Hyper Parameters:'
  WRITE(*,*) 'Length scales:', hyper(1:ndim)
  WRITE(*,*) 'Covariance Magnitude:', hyper(ndim+1)
  WRITE(*,*) 'Noise Level:', hyper(ndim+2)


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
  
  WRITE(*,*) 'Predict Output'
  !  Predict Output from Gradient-enhanced Kriging model
  call krigingfuncpredict(ndim,ntot,X,stot,H,beta,V,hyper,mtot,Xm,Hm,Ym,covarflag)

WRITE(*,*) 'Predict Output'
  !  Predict Variance associated with function values from Gradient-enhanced Kriging model
  call krigingfuncvariance(ndim,ntot,X,stot,H,hyper,mtot,Xm,Hm,Sm,covarflag)

WRITE(*,*) 'Predict Output'
  !  Predict derivative values associated from GE Kriging model
  call kriginggradpredict(ndim,ntot,X,stot,H,beta,V,hyper,mtot,Xm,htot,ptsm,dimsm,Gm,dYm,covarflag)  

WRITE(*,*) 'Predict Output'
  !  Predict variance of derivative values from GE Kriging model
  call kriginggradvariance(ndim,ntot,X,stot,H,hyper,mtot,Xm,htot,ptsm,dimsm,Gm,dSm,covarflag)  

  ! Output Training data
  open(unit=1,file='train.dat',status='replace',iostat=ierror)
  do i=1,ntot
     !          Location, Function value, Function variance, Derivative, Derivative Variance
     WRITE(1,99) X(1,i), Y(i)
  end do

  ! Output Test data
  open(unit=2,file='test.dat',status='replace',iostat=ierror)
  do i=1,mtot
     !          Location, Function value, Function Std. Dev., Derivative, Derivative Std. Dev.
     WRITE(2,99) Xm(1,i), Ym(i), sqrt(Sm(i)), dYm(i), sqrt(dSm(i))
  end do
99 format(99(F25.15))

  stop
end program krigingwrapper
