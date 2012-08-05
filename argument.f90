!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file argument.f90
!> \brief This modeule contains the definitions for the structures used to pass supplemental arguments to the function calls within the optimization subroutines (<a href=patternsearch_8f90.html>patternsearch</a> and <a href=simplexsearch_8f90.html>simplexsearch</a>). <br>
!> Each function type (such as likelihood, likelihood_mle, and the associated Kriging subroutines) has its own structure that stores the supplemental arguments required for the function call. Theses structures are themselves elements of the structure "generalarg" which is used as a wrapper for inputing these structures into the optimization subroutines.

module argument
  implicit none
  
  TYPE arg_toy
     character(len=60) descriptor
     real(8) alpha
  END TYPE arg_toy
  TYPE arg_likelihood
     character(len=60) descriptor
     integer ndim
     integer ntot
     real(8), pointer :: X(:,:)
     real(8), pointer :: Y(:)
     integer stot
     real(8), pointer :: H(:,:)
     real(8) sigmaN
  END TYPE arg_likelihood
  TYPE arg_alllikelihood
     character(len=60) descriptor
     integer ndim
     integer ntot
     real(8), pointer :: X(:,:)
     real(8), pointer :: Y(:)
     integer stot
     real(8), pointer :: H(:,:)
  END TYPE arg_alllikelihood
  TYPE arg_likelihood_grad
     character(len=60) descriptor
     integer ndim
     integer ntot
     real(8), pointer :: X(:,:)
     real(8), pointer :: Y(:)
     integer gtot
     integer, pointer :: pts(:)
     integer, pointer :: dims(:)
     real(8), pointer :: dY(:)
     integer stot
     real(8), pointer :: H(:,:)
     real(8) sigmaN
     real(8) sigmaNG
  END TYPE arg_likelihood_grad
  TYPE arg_alllikelihood_grad
     character(len=60) descriptor
     integer ndim
     integer ntot
     real(8), pointer :: X(:,:)
     real(8), pointer :: Y(:)
     integer gtot
     integer, pointer :: pts(:)
     integer, pointer :: dims(:)
     real(8), pointer :: dY(:)
     integer stot
     real(8), pointer :: H(:,:)
  END TYPE arg_alllikelihood_grad
  TYPE arg_kriging
     character(len=60) descriptor
     character(len=60) value
     integer ndim
     integer ntot
     real(8), pointer :: X(:,:)
     integer stot
     real(8), pointer :: H(:,:)
     real(8), pointer :: beta(:)
     real(8), pointer :: V(:)
     real(8), pointer :: hyper(:)
     integer covarflag
     real(8) Ymin
  end type arg_kriging
  TYPE arg_kriging_grad
     character(len=60) descriptor
     character(len=60) value
     integer ndim
     integer ntot
     real(8), pointer :: X(:,:)
     integer gtot
     integer, pointer :: pts(:)
     integer, pointer :: dims(:)
     integer stot
     real(8), pointer :: H(:,:)
     real(8), pointer :: beta(:)
     real(8), pointer :: V(:)
     real(8), pointer :: hyper(:)
     integer covarflag
     real(8) Ymin
  end type arg_kriging_grad
  TYPE generalarg
     character(len=60) descriptor
     type (arg_toy) toyarg
     type (arg_likelihood) likearg
     type (arg_likelihood_grad) gradlikearg
     type (arg_alllikelihood) alllikearg
     type (arg_alllikelihood_grad) gradalllikearg
     type (arg_kriging) krigingarg
     type (arg_kriging_grad) gradkrigingarg
  END TYPE generalarg 
end module argument
  
