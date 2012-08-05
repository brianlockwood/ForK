!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file krigingfuncsample.f90
!> \brief This subroutine samples from the Kriging model's output distribution. It requires <a href=krigingfunccovar_8f90.html>krigingfunccovar</a> and <a href=krigingfuncpredict_8f90.html>krigingfuncpredict</a> to be call first and cholesky must be performed on the output distribution covariance matrix.

!> \detail This subroutine extracts a sample from the Kriging models output distribution. This sample extraction is given as:
!>\f[ Y_{p} = Y_{m} + L u\f]
!> where \f$Y_{m}\f$ is the mean of the Kriging function's output (the result of <a href=krigingfuncpredict_8f90.html>krigingfuncpredict</a>) and \f$u\f$ is a vector of normally distributioned values of length \f$M\f$ (normal distribution with zero mean and standard deviation of 1). The matrix \f$L\f$ is the result of cholesky performed on the covariance matrix of the Kriging model's output distribution (the result of <a href=krigingfunccovar_8f90.html>krigingfunccovar</a>).

!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 18, 2012
!> \param(in) <b> mtot</b>: The number of test points to be predicted
!> \param(in) <b> Ym </b> : The mean of the output distribution (the result of <a href=krigingfuncpredict_8f90.html>krigingfuncpredict</a>) (size=[mtot])
!> \param(in) <b> S </b>:   The cholesky factorization of the covariance matrix from the output distribution (<a href=krigingfunccovar_8f90.html>krigingfunccovar</a>) stored in the lower triangle. (size=[mtotxmtot])
!> \param(in) <b> U </b> :  A vector of normally distributed random numbers (size=[mtot])
!> \param(out) <b> Yp </b>: The output sample from the Kriging model. This subroutine should be called multiple times to acquire a large number of output samples, allowing for global predictions to be calculated multiple times and a output distribution constructed. 

subroutine krigingfuncsample(mtot,Ym,S,U,Yp)
  use choleskymod
  implicit none
  integer, intent(in) :: mtot
  real(8), intent(in) :: Ym(mtot)
  real(8), intent(in) :: S(mtot,mtot)
  real(8), intent(in) :: U(mtot)
  real(8), intent(out) :: Yp(mtot)
  
  integer i,j

!  U times Cholesky of Covariance gives Yp
  do i=1,mtot
      Yp(i)=Ym(i)
     do j=1,i
        Yp(i)=Yp(i)+S(i,j)*U(j)
     end do
  end do
  
  return
end subroutine krigingfuncsample
