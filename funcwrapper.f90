!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file funcwrapper.f90
!> \brief This subroutine unwraps the supplemental arguments from the structure "generalarg" and supplies the appropriate part of the structure to the function "func" defined in the module <a href=funcmod_8f90.html>funcmod</a>. This unwrapping is based on the descriptor component of the "generalarg" structure.

subroutine funcwrapper(n,X,f,args)
  use argument
  use funcmod
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: X(n)
  real(8), intent(out) :: f

  type (generalarg), intent(in) :: args

  character(len=60) descriptor

  descriptor = args%descriptor

  if (descriptor=='Likelihood') then
     call func(n,X,f,args%likearg)
  elseif (descriptor=='Likelihood Grad') then
     call func(n,X,f,args%gradlikearg)
  elseif (descriptor=='AllLikelihood') then
     call func(n,X,f,args%alllikearg)
  elseif (descriptor=='AllLikelihood Grad') then
     call func(n,X,f,args%gradalllikearg)
  elseif (descriptor=='Kriging') then
     call func(n,X,f,args%krigingarg)
  elseif (descriptor=='Kriging Grad') then
     call func(n,X,f,args%gradkrigingarg)
  else
     WRITE(*,*) 'SOMETHING IS WRONG'
  end if
  
  return
end subroutine funcwrapper
