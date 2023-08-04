!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_eigen_mod 
!>@author Asdrubal Lozada
!>        Laboratory of Theoretical Chemistry - LQT
!>        Federal University of São Carlos
!>        <http://www.lqt.dq.ufscar.br>
!>@email  aslozada@gmail.com
!-----------------------------------------------------------------------------------
!   Copyright 2018 Asdrubal Lozada
!
!   This program is free software: you can redistribute it and/or modify 
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------
!Written by Asdrubal Lozada
!Laboratory of Theoretical Chemistry, LQT -- UFSCar
!e-mail: aslozada@gmail.com
!-----------------------------------------------------------------------------------
module kanon_eigen_mod
  use kanon_kinds_mod, only:&
    & DP
  use kanon_error_mod, only:&
    & error,&
    & ERROR_ALLOCATION,&
    & ERROR_DEALLOCATION

  implicit none
!@  private

  !>@class eigen  
  type, public :: eigen_t
    real(DP), allocatable, dimension(:)   :: values
    real(DP), allocatable, dimension(:,:) :: vectors
  end type eigen_t

  public :: Get_eigenvalues

  !>@note use jacobi method
contains
  type(eigen_t) function Get_eigenvalues ( A, n ) result(eigen)
  real(DP), intent(inout) :: A(:,:)
  integer, intent(in)     :: n
  real(DP), allocatable ::&
    & v(:,:), d(:), b(:), z(:)
  integer, allocatable :: r(:)
  real(DP) :: sm, thresh, g, h, t, theta, c, s, tau, rot
  integer  :: max_iter, p, q, i, j, k, l, aux, ierr=0

  allocate(v(n,n), d(n), b(n), z(n), r(n), stat=ierr)
  if ( ierr /=0 ) call error%MsgError( 1, ERROR_ALLOCATION )
                
  allocate(eigen%values(n), eigen%vectors(n,n), stat=ierr) 
  if ( ierr /=0 ) call error%MsgError(1, ERROR_ALLOCATION)

  if ( n < 1 ) call error%MsgError( 1, message_error='Invalid lenght matrix.') 

  b = 0.0_DP; d = 0.0_DP; z = 0.0_DP
  rot = 0.0_DP; sm = 0.0_DP; thresh = 0.0_DP
  g = 0.0_DP;  h = 0.0_DP;  t = 0.0_DP
  theta = 0.0_DP; c = 0.0_DP; s = 0.0_DP
  tau = 0.0_DP;  aux = 0

  v = 0.0_DP
                
  do p = 1, n
    do q = 1, n
      if ( p == q ) v(p,q) = 1.0_DP
    end do
  end do
    
  do p = 1, n
     d(p) = a(p,p)
  end do
    
  do p = 1, n
     b(p) = d(p)
  end do

  max_iter = 50

  do i = 1, max_iter

    do p = 1, n
      do q = 1, p-1
        thresh = thresh + a(p,q)**2
      end do
    end do

    thresh = dsqrt(thresh)/real(4*n,DP)

    if ( thresh == 0.0_DP ) exit

    sm = 0.0_DP
    do p = 1, n - 1
       do q = p + 1, n
        sm = sm + dabs(a(p,q))
       end do
    end do


    if ( sm == 0.0_DP ) exit
    if ( i < 4 ) then
       thresh = 0.2_DP * sm / (n * n)
    else
       thresh = 0.0_DP
    end if

    do p = 1, n - 1
       do q = p + 1, n      
         g = 100.0_DP * dabs(a(p,q))


         if ( (i > 4) .and. ( (dabs(d(p))+g) == dabs(d(p))  )&
            & .and. ( (dabs(d(q))+g)  == dabs(d(q)) )   ) then

             a(p,q) = 0.0_DP

         elseif(dabs(a(p,q)) > thresh) then
                h = d(q) - d(p)

                if ( (dabs(h)+g) == dabs(h) ) then
                   t = a(p,q)/h
                else
                   theta = 0.5_DP * h / a(p,q)
                   if ( theta > 0.0_DP ) then
                      t = 1.0_DP / (theta + dsqrt(1.0_DP + theta**2)) 
                   else
                      t = -1.0_DP / (-theta + dsqrt(1.0_DP + theta**2)) 
                   end if
                end if  
            

                c = 1.0_DP / dsqrt(1.0_DP + t**2)
                s = t * c
                tau = s / (1.0_DP + c)
                h = t * a(p,q)

                z(p) = z(p) - h
                z(q) = z(q) + h
                d(p) = d(p) - h
                d(q) = d(q) + h
                a(p,q) = 0.0_DP
          
                do j = 1, p - 1
                   g = a(j,p)
                   h = a(j,q)
                   a(j,p) = g - s * (h + g * tau)
                   a(j,q) = h + s * (g - h * tau)
                end do

                do j = p + 1, q - 1
                   g = a(p,j)
                   h = a(j,q)
                   a(p,j) = g - s * (h + g * tau)
                   a(j,q) = h + s * (g - h * tau)
                end do

                do j = q + 1, n
                   g = a(p,j)
                   h = a(q,j)
                   a(p,j) = g - s * (h + g * tau)
                   a(q,j) = h + s * (g - h * tau)
                end do

                do j = 1, n
                   g = v(j,p)
                   h = v(j,q)
                   v(j,p) = g - s * (h + g * tau)
                   v(j,q) = h + s * (g - h * tau)
                end do

                rot = rot + 1
         end if 
        end do
    end do

    do  p = 1, n
        b(p) = b(p) + z(p)
        d(p) = b(p)
        z(p) = 0.0_DP
    end do

    end do 

    r = 0
    do k = 1, n
       r(k) = k
    end do

    do k = 1, n - 1
       do l = k + 1, n

          if ( d(k) < d(l) ) then
             aux = r(k)
             r(k) = r(l)
             r(l) = aux
          end if
          aux = 0

       end do
    end do  
    
    eigen%values   = 0.0_DP
    eigen%vectors = 0.0_DP

    do k = 1, n
       eigen%values(k) = d(r(k))
    end do

    do j = 1, n
       do k = 1, n
          eigen%vectors(j,k) = v(j,r(k))
       end do  
    end do 
                
    do k = 1, n
       do j = 1, n
          eigen%vectors(j,k) = v(j,r(k))
       end do  
    end do

    deallocate(v, d, r, b, z, stat=ierr)
    if ( ierr /=0 ) call error%MsgError( 1, ERROR_DEALLOCATION )
  end function Get_eigenvalues


!=======================================================
  

! Adapted from https://people.sc.fsu.edu/~jburkardt/

subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

  implicit none

  integer :: n

  real(dp) ::  a(n,n), bw(n), c, d(n), g, gapq, h
  integer :: i, it_max, it_num, j, k, l, m, p, q, rot_num
  real(dp) ::  s, t, tau, term, termp, termq, theta, thresh
  real(dp) ::  v(n,n), w(n), zw(n)

  do j = 1, n
    do i = 1, n
      v(i,j) = 0.0D+00
    end do
    v(j,j) = 1.0D+00
  end do

  do i = 1, n
    d(i) = a(i,i)
  end do

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0D+00
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
      do i = 1, j - 1
        thresh = thresh + a(i,j) ** 2
      end do
    end do

    thresh = sqrt ( thresh ) / real ( 4 * n, kind = dp )

    if ( thresh == 0.0D+00 ) then
      exit 
    end if

    do p = 1, n
      do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then 
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h                  
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
    end do

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

  end do
!
!  Restore upper triangle of input matrix.
!
  do j = 1, n
    do i = 1, j - 1
      a(i,j) = a(j,i)
    end do
  end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
  do k = 1, n - 1

    m = k

    do l = k + 1, n
      if ( d(l) < d(m) ) then
        m = l
      end if
    end do

    if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

    end if

  end do

  return
end subroutine jacobi_eigenvalue

!=======================================================





 end module kanon_eigen_mod
