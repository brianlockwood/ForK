!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!>\file toyfunc.f90
!>\brief Subroutine for computing simple functions that can be used to test the Kriging model and the optimization subroutines. Only used within the example programs.

subroutine toyfunc(ndim,X,Y,dY)
  implicit none
  integer, intent(in) :: ndim
  real(8), intent(in) :: X(ndim)
  real(8), intent(out) :: Y, dY(ndim)

  integer i

  real(8) PI
  real(8) alpha

PI=4.D0*ATAN(1.D0)

Y=10.D0*ndim
do i=1,ndim
   Y=Y+X(i)**2-10.D0*cos(2*pi*X(i))
   dy(i)=2.D0*X(i)+20.D0*PI*sin(2.D0*pi*X(i))
end do

  return
end subroutine toyfunc
