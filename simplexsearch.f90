!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file simplexsearch.f90
!> \brief Simplex Search used for optimization <br>
!>        Used for max likelihood optimization and determining min/max values from Kriging surface <br>
!>        Suitable for fast local optimization using function values

!> \detail This subroutine implements the <a href=http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method>Nelder-Mead method</a> for optimization applied to the max likelihood problem and optimization performed on the Kriging surface itself. <br>
!> <br>
!>  For this method, a simplex element in \f$N\f$  dimensions (set by the argument \a odim) is constructed and the search direction for new function values is constructed based on the vertices and faces of this element. The procedure for determining the minimum value is outlined below. For maximization, all function values can be multiplied by \f$-1\f$ and the minimization algorithm can be applied.

!>  <ul><li> Starting with \f$N + 1\f$ function evaluations, the maximum function value in the element is determined (\f$F_{max}\f$)
!>      <li> The search direction is given by the vector between the maximum function value of the element and the centroid of the remaining vertices. <br>
!>           \f[\vec{p} = \vec{X}_{face} - \vec{X}_{max}\f] where \f[\vec{X}_{face} =\frac{1}{N} \sum_{i \ne max} \vec{X}_{i}\f]
!>      <li> A line search is next performed in this direction to determine the new candidate point. The following points are tested.
!>          <ol><li> Reflected point: \f$\vec{X}_{r} = \vec{X}_{face} + \vec{p} \f$
!>              <li> Expanded point: \f$\vec{X}_{e} = \vec{X}_{face} + 2 \vec{p} \f$
!>              <li> Outside Contraction: \f$\vec{X}_{oc} = \vec{X}_{face} + \frac{1}{2} \vec{p} \f$
!>              <li> Inside Contraction: \f$\vec{X}_{ic} = \vec{X}_{face} - \frac{1}{2} \vec{p} \f$</ol>
!>          If the reflected point's function value is the new minimum value, the expanded point is tested. <br>
!>          If the reflected point's function value is less than the previous maximum value in the element but would still correspond to the maximum vertex value in the element, outside contraction is performed. <br>
!>          If the reflected point's function value is greater than the previous maximum value, inside contraction is performed. <br>
!>      <li> Based on the results of the line search, the previous maximum vertex value is replaced by the best candidate determined in the line search (either the reflected point, expanded point, outside contraction or inside contraction). 
!>      <li> In the event that none of the candidate locations tested in the line search are able to reduce the maximum function value of the element, a rescaling is performed around the minimum location.
!>              \f[ \vec{X}'_{i} = \vec{X}_{min} + \frac{1}{2} (\vec{X}_{i} - \vec{X}_{min}) \f] </ul> <br>
!> <br>
!> More details of this algorithm can be found in <i>Numerical Optimization</i> by Nocedal and Wright

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
 subroutine simplexsearch(odim,Xl,Xf,Fmin,L,U,flag,args)
   use argument
   implicit none
   !  Simplex search 

   !  Optimization Variables
   integer, intent(in) :: odim
   real(8), intent(in) :: Xl(odim)
   real(8), intent(in) :: L(odim),U(odim)
   integer, intent(in) :: flag

   real(8), intent(out) :: Xf(Odim)
   real(8), intent(out) :: Fmin

   !  Variables needed for function call
   type (generalarg) args

  real(8) Xcell(odim,odim+1),Xface(odim),Fcell(odim+1)

   real(8) V(odim)

   real(8) alpha,gamma,sigma,rho
   
   integer i,j

   real(8) Fnew,Fold,Fconv
   real(8) Fmax,Fsecond
   real(8) Fr,Fe,Fc

   integer maxloc, minloc
   
   real(8) Xr(odim),Xe(odim),Xco(odim)

   integer count
   
   real(8) maxedge, edgelength

   integer success

   alpha=1.D0
   gamma=2.D0
   sigma=0.5D0
   rho=0.5D0
    
   Xcell(:,odim+1)=Xl(:)
   do i=1,odim
      Xcell(:,i)=Xl(:)
      Xcell(i,i)=Xcell(i,i)+1.D0
      do j=1,odim
         if (Xcell(j,i)<L(j)) then
            Xcell(j,i)=L(j)
         elseif (Xcell(j,i)>U(j)) then
            Xcell(j,i)=U(j)
         end if
      end do
   end do
  
   ! Map Optimization Variables to Hyper parameters
   do i=1,odim+1
      call funcwrapper(odim,Xcell(:,i),Fnew,args)
      if (flag==1) then
         Fnew=-Fnew
      end if
      Fcell(i)=Fnew
   end do
   Fconv=minval(Fcell)
   fold=1.D3
   
   count=0
  WRITE(*,98) 'Iteration', ('X',j=1,odim), 'F', 'Step size'
98 format(99(A15))
   do 
     count=count+1
     !  Determine Min,Max,Second Worst Candidates
     Fmax=-1.D12
     do i=1,odim+1
        if (Fcell(i)>=Fmax) then
           Fmax=Fcell(i)
           maxloc=i
        end if
     end do
     Fmin=1.D12
     do i=1,odim+1
        if (Fcell(i)<Fmin) then
           Fmin=Fcell(i)
           minloc=i
        end if
     end do
     Fsecond=-1.D12
     do i=1,odim+1
        if (i/=maxloc) then
           if (Fcell(i)>Fsecond) then
              Fsecond=Fcell(i)
           end if
        end if
     end do
     
     !  Check for Convergence
     if (Fmin<Fconv) then
        Fold=Fconv
        Fconv=Fmin
     end if
     !  Max Edgelength
     maxedge=0.D0
     do i=1,odim+1
        if (i/=minloc) then
           edgelength=0.D0
           do j=1,odim
              edgelength=edgelength+(Xcell(j,i)-Xcell(j,minloc))**2
           end do
           edgelength=sqrt(edgelength)
           maxedge=max(maxedge,edgelength)
        end if
     end do
     if (maxedge<1.D-6) exit
     if (abs((Fconv-Fold)/Fconv)<1.D-6) exit

     !  Determine Centroid and Search Direction
     Xface=0.D0
     do i=1,odim+1
        if (i/=maxloc) then
           Xface(:)=Xface(:)+Xcell(:,i)
        end if
     end do
     Xface=Xface/(1.D0*odim)
     V(:)=(Xface(:)-Xcell(:,maxloc))
     
!  Check Reflected Point
     Xr(:)=Xface(:)+alpha*V(:)
     do j=1,odim
        if (Xr(j)<L(j)) then
           Xr(j)=L(j)
        elseif (Xr(j)>U(j)) then
           Xr(j)=U(j)
        end if
     end do

      call funcwrapper(odim,Xr,Fr,args)
      if (flag==1) then
         Fr=-Fr
      end if

      success=0
      !  Reflected Point is Better than but not the best
      if ((Fr<Fsecond).and.(Fr>=Fmin)) then
         success=1
         Xcell(:,maxloc)=Xr(:)
         Fcell(maxloc)=Fr
      !  Reflected Point is Best so far, expand further
      elseif (Fr<Fmin) then
         success=1
         Xe(:)=Xface(:)+gamma*V(:)
         do j=1,odim
            if (Xe(j)<L(j)) then
               Xe(j)=L(j)
            elseif (Xe(j)>U(j)) then
               Xe(j)=U(j)
            end if
         end do
         call funcwrapper(odim,Xe,Fe,args)
         if (flag==1) then
            Fe=-Fe
         end if
         ! Expanded point is better than reflected, use it
         if (Fe<Fr) then
            Xcell(:,maxloc)=Xe(:)
            Fcell(maxloc)=Fe
         else
            Xcell(:,maxloc)=Xr(:)
            Fcell(maxloc)=Fr
         end if
         !  Reflected Point is not better
      else
         !  Reflected Point is an improvement over worst candidate, try outside contraction
         if ((Fr>=Fsecond).and.(Fr<Fmax)) then
            Xco(:)=Xface(:)+rho*V(:)
            do j=1,odim
               if (Xco(j)<L(j)) then
                  Xco(j)=L(j)
               elseif (Xco(j)>U(j)) then
                  Xco(j)=U(j)
               end if
            end do
            
            call funcwrapper(odim,Xco,Fc,args)
            if (flag==1) then
               Fc=-Fc
            end if
            
            if (Fc<Fr) then
               success=1
               Xcell(:,maxloc)=Xco(:)
               Fcell(maxloc)=Fc
            end if
         !  Reflected Point is not better than worst candidate, try inside contraction
         else
            Xco(:)=Xface(:)-rho*V(:)
            do j=1,odim
               if (Xco(j)<L(j)) then
                  Xco(j)=L(j)
               elseif (Xco(j)>U(j)) then
                  Xco(j)=U(j)
               end if
            end do
            
            call funcwrapper(odim,Xco,Fc,args)
            if (flag==1) then
               Fc=-Fc
            end if
            if (Fc<Fmax) then
               success=1
               Xcell(:,maxloc)=Xco(:)
               Fcell(maxloc)=Fc
            end if
         end if
      end if
      !  No new candidates were found, rescale around best point
      if (success==0) then
         do i=1,odim+1
            if (i/=minloc) then
               Xcell(:,i)=Xcell(:,minloc)+sigma*(Xcell(:,i)-Xcell(:,minloc))
            end if
         end do
      end if

      WRITE(*,99) count, Xcell(:,minloc),-Fconv
      if (abs((Fconv-Fold)/Fconv)<1.D-6) exit
99 format(I15,99(E15.5))
   end do
  
   WRITE(*,*) count

   Xf(:)=Xcell(:,minloc)
   call funcwrapper(odim,Xf,Fmin,args)

   return
 end subroutine simplexsearch
