!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file opt.f90
!> \brief This module contains some supplemental subroutines used for the <a href=patternsearch_8f90.html>patternsearch</a>. These functions are used simply to allow for more generality in the patternsearch. The functions are hard-coded to trival values for now.

module opt
  implicit none
  integer optflag
contains
   subroutine forcingfunc(t,rho)
   implicit none
   real(8), intent(in) :: t
   real(8), intent(out) :: rho

!   rho=5*t**2
   rho=0.D0
   
   return
 end subroutine forcingfunc

 subroutine expansionfunc(t,theta)
   implicit none
   real(8), intent(in) :: t
   real(8), intent(out) :: theta

   theta=2.D0
   
   return
 end subroutine expansionfunc
   
 subroutine contractionfunc(t,theta)
   implicit none
   real(8), intent(in) :: t
   real(8), intent(out) :: theta

   theta=0.5
   
   return
 end subroutine contractionfunc
end module opt
