!Copyright (C) 2012 Brian A. Lockwood
!
!This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License along with
!this program. If not, see <http://www.gnu.org/licenses/>.

!> \file choleskymod.f90
!> \brief Module containing matrix manipulation subroutines such as Cholesky factorization and matrix multiply. For the most part, the functions contained here are mere wrappers for <a href=http://www.netlib.org/lapack/>LAPACK</a>  calls. Commented out functions contain the "dumb" versions of these <a href=http://www.netlib.org/lapack/>LAPACK</a>  calls and are used primarily for debugging purposes.

module choleskymod
  implicit none
contains
  !> \brief Subroutine to perform Cholesky factorization in place and stored as a lower triangular (uses LAPACK).
#ifdef HAVELAPACK
  subroutine cholesky(n,A)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in out) :: A(n,n)

    integer laerror

    call dpotrf('L',n,A,n,laerror)
    
    return
  end subroutine cholesky
#else
  subroutine cholesky(n,A)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in out) :: A(n,n)

    integer i,j,k

    real(8) Ljj, Lij

    do j=1,n
       Ljj=A(j,j)
       do k=1,j-1
          Ljj=Ljj-A(j,k)**2
       end do
       Ljj=sqrt(Ljj)
       do i=j+1,n
          Lij=A(i,j)
          do k=1,j-1
             Lij=Lij-A(i,k)*A(j,k)
          end do
          A(i,j)=Lij/Ljj
       end do
       A(j,j)=Ljj
    end do

    do i=1,n
       do j=i+1,n
          A(i,j)=0.D0
       end do
    end do
    
    return
  end subroutine cholesky
#endif

!> \brief Subroutine to solve Ax=b if A is a cholesky factorized matrix stored as a lower triangular.
subroutine choleskysolve(n,A,b,x)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n),b(n)
  real(8), intent(out) :: x(n)

  integer i,j

  do i=1,n
     X(i)=B(i)
     do j=1,i-1
        X(i)=X(i)-A(i,j)*X(j)
     end do
     X(i)=X(i)/A(i,i)
  end do
  
  do i=n,1,-1
     do j=i+1,n
        X(i)=X(i)-A(j,i)*X(j)
     end do
     X(i)=X(i)/A(i,i)
  end do
  
  return
end subroutine choleskysolve

!> \brief Subroutine to invert a cholesky factorized matrix.
subroutine invertcholesky(n,A,Ai)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n)
  
  real(8), intent(out) :: Ai(n,n)
  
  integer i,j
  real(8) X(n), b(n)
  
  do i=1,n
     B(:)=0.d0
     B(i)=1.D0
     call choleskysolve(n,A,b,x)
     Ai(:,i)=x(:)
  end do
  
  do i=1,n
     do j=i+1,n
        Ai(i,j)=Ai(j,i)
     end do
  end do
  
  return
end subroutine invertcholesky
  
!> \brief Subroutine to invert a positive definite symmetric matrix and store the inverse in place (uses LAPACK). Result stored as a lower triangular matrix.
#ifdef HAVELAPACK
subroutine invertsymmetric(n,A)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in out) :: A(n,n)

  integer laerror
  integer i,j

  call dpotrf('L',n,A,n,laerror)
  call DPOTRI('L',n,A,n,laerror)

  return
end subroutine invertsymmetric
#else
subroutine invertsymmetric(n,A)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in out) :: A(n,n)

  real(8) Ai(n,n)

  call cholesky(n,A)

  call invertcholesky(n,A,Ai)

  A=Ai

  return
end subroutine invertsymmetric
#endif

!> \brief Solves Ax=b for a symmetric positive definite matrix. Initial value of X is the vector b and the final value is the solution. Uses LAPACK.
#ifdef HAVELAPACK
subroutine symmetricsolve(n,A,x)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in out) :: A(n,n)
  real(8), intent(in out) :: X(n)

  integer laerror

  call DPOSV('L',n,1,A,n,X,n,laerror)

  return
end subroutine symmetricsolve
#else
subroutine symmetricsolve(n,A,x)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in out) :: A(n,n)
  real(8), intent(in out) :: X(n)

  real(8) B(n)

  B=X

  call cholesky(n,A)
  call choleskysolve(n,A,B,x)
  
  return
end subroutine symmetricsolve
#endif

!> \brief Subroutine that performs matrix multiplication \f$A*B=C\f$. Uses LAPACK.
#ifdef HAVELAPACK
subroutine matrixmult(n,m,A,p,B,C)
  implicit none
  integer, intent(in) :: n,m,p
  real(8), intent(in) :: A(n,m), B(m,p)
  real(8), intent(out) :: C(n,p)

  call dgemm('N','N',n,p,m,1.D0,A,n,B,m,0.0,C,n)  

  return
end subroutine matrixmult
#else
subroutine matrixmult(n,m,A,p,B,C)
  implicit none
  integer, intent(in) :: n,m,p
  real(8), intent(in) :: A(n,m), B(m,p)
  real(8), intent(out) :: C(n,p)

  integer i,j,k

  do i=1,n
     do j=1,p
        C(i,j)=0.D0
        do k=1,m
           C(i,j)=C(i,j)+A(i,k)*B(k,j)
        end do
     end do
  end do

  return
end subroutine matrixmult
#endif


!> \brief Subroutine that performs matrix multiplication \f$A*B=C\f$ when A is symmetric. Uses LAPACK.
#ifdef HAVELAPACK
subroutine symmatrixmult(n,A,p,B,C)
  implicit none
  integer, intent(in) :: n,p
  real(8), intent(in) :: A(n,n), B(n,p)
  real(8), intent(out) :: C(n,p)

  integer laerror

  call dsymm('L','L',n,p,1.D0,A,n,B,n,0.D0,C,n)

  return
end subroutine symmatrixmult
#else
subroutine symmatrixmult(n,A,p,B,C)
  implicit none
  integer, intent(in) :: n,p
  real(8), intent(in) :: A(n,n), B(n,p)
  real(8), intent(out) :: C(n,p)

  integer i,j,k

  do i=1,n
     do j=1,p
        C(i,j)=0.D0
        do k=1,i
           C(i,j)=C(i,j)+A(i,k)*B(k,j)
        end do
        do k=i+1,n
           C(i,j)=C(i,j)+A(k,i)*B(k,j)
        end do
     end do
  end do

  return
end subroutine symmatrixmult
#endif


!> \brief Subroutine that performs matrix multiplication as \f$A*B^{T}=C\f$. Uses LAPACK.
#ifdef HAVELAPACK
subroutine matrixmulttrans(n,m,A,p,B,C)
  implicit none
  integer, intent(in) :: n,m,p
  real(8), intent(in) :: A(n,m), B(p,m)
  real(8), intent(out) :: C(p,n)

  call dgemm('N','T',p,n,m,1.D0,B,p,A,n,0.0,C,p)  

  return
end subroutine matrixmulttrans
#else
subroutine matrixmulttrans(n,m,A,p,B,C)
  implicit none
  integer, intent(in) :: n,m,p
  real(8), intent(in) :: A(n,m), B(p,m)
  real(8), intent(out) :: C(p,n)

  integer i,j,k

  do j=1,n
     do i=1,p
        C(i,j)=0.D0
        do k=1,m
           C(j,i)=C(j,i)+A(i,k)*B(j,k)
        end do
     end do
  end do

  return
end subroutine matrixmulttrans
#endif


!> \brief Subroutine that performs matrix multiplication as \f$A*B^{T}=C\f$ when A is symmetric. Uses LAPACK.
#ifdef HAVELAPACK
subroutine symmatrixmulttrans(n,A,p,B,C)
  implicit none
  integer, intent(in) :: n,p
  real(8), intent(in) :: A(n,n), B(p,n)
  real(8), intent(out) :: C(p,n)

  integer laerror

  call dsymm('R','L',p,n,1.D0,A,n,B,p,0.D0,C,p)

  return
end subroutine symmatrixmulttrans
#else
subroutine symmatrixmulttrans(n,A,p,B,C)
  implicit none
  integer, intent(in) :: n,p
  real(8), intent(in) :: A(n,n), B(p,n)
  real(8), intent(out) :: C(p,n)

  integer i,j,k

  do i=1,n
     do j=1,p
        C(j,i)=0.D0
        do k=1,i
           C(j,i)=C(j,i)+A(i,k)*B(j,k)
        end do
        do k=i+1,n
           C(j,i)=C(j,i)+A(k,i)*B(j,k)
        end do
     end do
  end do

  return
end subroutine symmatrixmulttrans
#endif

!> \brief Subroutine that performs matrix-vector multiplication as \f$ A x = b \f$. Uses LAPACK.
#ifdef HAVELAPACK
subroutine matvec(n,m,A,x,B)
  implicit none
  integer, intent(in) :: n,m
  real(8), intent(in) :: A(n,m),x(m)
  real(8), intent(out) :: B(n)

  call dgemv('N',n,m,1.D0,A,n,X,1,0.D0,B,1)

  return
end subroutine matvec
#else
subroutine matvec(n,m,A,x,B)
  implicit none
  integer, intent(in) :: n,m
  real(8), intent(in) :: A(n,m),x(m)
  real(8), intent(out) :: B(n)

  integer i,j

  do i=1,n
     B(i)=0.D0
     do j=1,m
        B(i)=B(i)+A(i,j)*X(j)
     end do
  end do

  return
end subroutine matvec
#endif


!> \brief Subroutine that performs matrix-vector multiplication as \f$ A x = b \f$ when A is symmetric. Uses LAPACK.
#ifdef HAVELAPACK
subroutine symmatvec(n,A,x,B)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n),x(n)
  real(8), intent(out) :: B(n)

  integer laerror

  call dsymv('L',n,1.D0,A,n,X,1,0.D0,B,1)

  return
end subroutine symmatvec
#else
subroutine symmatvec(n,A,x,B)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n),x(n)
  real(8), intent(out) :: B(n)

  integer i,j

  do i=1,n
     B(i)=0.D0
     do j=1,i
        B(i)=B(i)+A(i,j)*X(j)
     end do
     do j=i+1,n
        B(i)=B(i)+A(j,i)*X(j)
     end do
  end do

  return
end subroutine symmatvec
#endif

!> \brief Subroutine that performs matrix-vector multiplication as \f$ A^{T} x = b \f$. Uses LAPACK.
#ifdef HAVELAPACK
subroutine matvectrans(n,m,A,x,B)
  implicit none
  integer, intent(in) :: n,m
  real(8), intent(in) :: A(n,m),x(n)
  real(8), intent(out) :: B(m)

  call dgemv('T',n,m,1.D0,A,n,X,1,0.D0,B,1)

  return
end subroutine matvectrans
#else
subroutine matvectrans(n,m,A,x,B)
  implicit none
  integer, intent(in) :: n,m
  real(8), intent(in) :: A(n,m),x(n)
  real(8), intent(out) :: B(m)

  integer i,j

  do i=1,m
     B(i)=0.D0
     do j=1,n
        B(i)=B(i)+A(j,i)*X(j)
     end do
  end do

  return
end subroutine matvectrans
#endif

end module choleskymod
