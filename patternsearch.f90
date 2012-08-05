!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file patternsearch.f90
!> \brief Pattern Search used for global optimization <br>
!>        Used for max likelihood optimization and determining min/max values from Kriging surface

!> \detail The pattern search is a function-only optimization method suitable for global optimization problems. Given an initial function evaluation, the pattern search performs additional function evaluations in a set of fixed search directions about this initial point (with the same step size in each direction). <br><br>
!>         If one of the new function values is less than the initial function value, the center of the stencil is shifted to the new minimum value and the process is repeated about this new minimum point. <br><br>
!>        If none of the new function evaluations are less than the initial value, the center of the stencil is unchanged and additional function evaluations are performed using half the step size. The process is repeated until the step size is sufficiently small or the relative change in the function value is below a specified tolerance.<br><br>
 
!> The directions used in this pattern search are the positive and negative coordinate directions of the design space. This choice leads to \f$2N\f$ total search directions (and consequently \f$2N\f$ function evaluations per iteration). <br><br>
!> For more details, refer to <a href="http://en.wikipedia.org/wiki/Pattern_search_(optimization)">Pattern_search_(optimization)</a>. <br>
!> <br>
!> This subroutine has shared memory parallelization using <a href="http://openmp.org/wp/">openmp</a>. Because each step in the optimization requires \f$2N\f$ independent function evaluations, these evaluations can be performed simultaneously on different threads. Although this parallelization is efficient from a scaling point of view, each thread requires a large memory footprint. In the case of likelihood evaluations (this also applies for some of the Kriging subroutines), each thread must compute and store the covariance matrix, requiring a large amount of stack memory. Hence, the environmental variable "OMP_STACKSIZE" must be set to a suitably large value to avoid a segmentation fault. A good rule of thumb is to base the STACKSIZE on the size of the covariance matrix. This matrix is \f$NxN\f$ for a function only Kriging model and \f$N(D+1)xN(D+1)\f$ for a gradient-enhanced model (where \f$N\f$ is the number of training points and \f$D\f$ is the dimension of the problem). The elements of this matrix are double precision real numbers and thus require \f$8\f$ bytes of storage for each element. Using these rules, a suitable STACKSIZE can be set. In some cases, the desired number of threads must be reduced to allow for the larger STACKSIZE for each thread to fit in main memory. Because the optimization is the most expensive part of the Kriging construction, determining the proper STACKSIZE and using the maximum number of threads is well worth the possible aggrevation. 


!> \author Brian Lockwood <br>
!>         Department of Mechanical Engineering <br>
!>         University of Wyoming
!> \date   May 15, 2012
!> \param(in)  <b> odim </b>: Dimension of the optimization
!> \param(in)  <b> Xl </b>:   Starting point for the optimization
!> \param(out) <b> Xf </b>:   Final optimal location
!> \param(out) <b> Fmin </b>: Final optimal value
!> \param(in)  <b> L </b>:  Lower bound for each variable
!> \param(in)  <b> U </b>:  Upper bound for each variable
!> \param(in)  <b> flag </b>: Determines whether minimum or maximum value is found <br>
!>                            <ul><li> flag = 0 Subroutine finds the minimum value
!>                                <li> flag = 1 Subroutine finds the maximum value</ul>
!> \param(in)  <b> args </b>: Structure containing additional arguments required for the function call <br>
!>                            These structures are in the module <a href=argument_8f90.html> argument </a> <br>
!>                            The structures for each function type are housed within the structure "generalarg"  and the subroutine <a href=funcwrapper_8f90.html>funcwrapper</a> detects the function type and selects the proper element of "generalarg"
 subroutine patternsearch(odim,Xl,Xf,Fmin,L,U,flag,args)
   use opt
   use argument
   implicit none
   !  Pattern search using General
   !  Optimization Variables
   integer, intent(in) :: odim
   real(8), intent(in) :: Xl(odim)
   real(8), intent(in) :: L(odim),U(odim)
   integer, intent(in) :: flag

   real(8), intent(out) :: Xf(Odim)
   real(8), intent(out) :: Fmin

   !  Variables needed for function call
   type (generalarg) args

  integer ndir

  integer i,j
  real(8) step
  real(8) Xc(odim)
  real(8), allocatable :: V(:,:)
  real(8), allocatable :: Xcand(:,:),Fs(:)
  
  real(8) Fnew,Fold
  real(8) Fcand
  
  integer candloc
  integer count
  
  real(8) phi, rho

  integer d,k

  integer time1,rate,time2

  ndir=2*odim
  allocate(V(odim,ndir),Xcand(odim,ndir),Fs(ndir))
  ! Step Directions
   do i=1,odim
      V(:,i)=0.D0
      V(i,i)=1.D0
   end do
   do i=1,odim
      V(:,i+odim)=0.D0
      V(i,i+odim)=-1.D0
   end do                      
!!$     
!!$   ndir=odim+1
!!$   allocate(V(odim,ndir),Xcand(odim,ndir),Fs(ndir))
!!$  ! Step Directions
!!$   do i=1,odim
!!$      V(:,i)=0.D0
!!$      V(i,i)=1.D0
!!$   end do
!!$   V(:,odim+1)=-1.D0
   
   step=1.D0

   ! Map Optimization Variables to Hyper parameters
   call funcwrapper(odim,Xl,Fnew,args)
   if (flag==1) then
      Fnew=-Fnew
   end if
   Fold=1.D3

   count=0
   Xc(:)=Xl(:)
   call system_clock(time1,rate)
   WRITE(*,98) 'Iteration', ('X',j=1,odim), 'F', 'Step size'
98 format(99(A15))
   do 
      if (step<1.D-6) exit
      count=count+1
      do i=1,ndir
         Xcand(:,i)=Xc(:)+step*V(:,i)
         do j=1,odim
            if (Xcand(j,i)<L(j)) then
               Xcand(j,i)=L(j)
            elseif (Xcand(j,i)>U(j)) then
               Xcand(j,i)=U(j)
            end if
         end do
      end do
      !$omp parallel do private(i) schedule(dynamic,1)
      do i=1,ndir
         call funcwrapper(odim,Xcand(:,i),Fs(i),args)
         if (flag==1) then
            Fs(i)=-Fs(i)
         end if
      end do
      !$omp end parallel do
      call forcingfunc(step,rho)
      Fcand=1.D12
      do i=1,ndir
         if (Fs(i)<(Fcand-rho)) then
            Fcand=Fs(i)
            candloc=i
         end if
      end do
      if (Fcand<Fnew) then
         Fold=Fnew
         Fnew=Fcand
         Xc(:)=Xcand(:,candloc)
         call expansionfunc(step,phi)
      else
         call contractionfunc(step,phi)
      end if
      step=step*phi
      WRITE(*,99) count, Xc(:),-Fnew, step
99 format(I15,99(E15.5))
      if (abs((Fnew-Fold)/Fnew)<1.D-6) exit
   end do
   call system_clock(time2,rate)
   
   Xf(:)=Xc(:)
   call funcwrapper(odim,Xf(:),Fmin,args)
!   WRITE(*,*) ndir*count, dble(time2-time1)/(dble(rate)),Fmin

   return
 end subroutine patternsearch
