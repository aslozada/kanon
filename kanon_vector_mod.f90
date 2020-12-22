!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_vector_mod
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
module kanon_vector_mod
  use kanon_kinds_mod, only:&
    & DP
  use kanon_error_mod, only:&
    & error,&
    & ERROR_ALLOCATION
    
  implicit none
  private

  !>@class vector
  type, public :: vector_t
     real(DP)  :: xyz(3)
  contains
     procedure :: norm   => Get_norm_vector  
  end type vector_t

  !>@params[in] type vector
  public :: Scalar_product
  public :: Cross_product 
  public :: Scalar_triple_product
  public :: Define_Rotation
  public :: Rotate_by_quaternion
  public :: Rotate_matrix_quaternion
  public :: Prevent_NAN

contains
  subroutine Rotate_matrix_quaternion( s, vector_V, Matrix, n_cols, qMatrix )
    type(vector_t), intent(in) :: vector_V
    real(DP), intent(in)    :: s
    real(DP), intent(in)    :: Matrix(:,:)
    real(DP), intent(inout) :: qMatrix(:,:)
    integer,  intent(in)    :: n_cols
    real(DP) :: q(3,3)
    integer  :: ierr=0

    real(DP), allocatable, dimension(:,:) :: Matrix_T, qMatrix_T

    if ( n_cols < 1 ) call error%MsgError( 1, message_error='Invalid length matrix.')

    allocate(Matrix_T(3,n_cols), qMatrix_T(3,n_cols), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )
   
    q = 0.0_DP

    q(1,1) = 1.0_DP - (2.0_DP* vector_V%xyz(2)**2) - (2.0_DP * vector_V%xyz(3)**2)
    q(2,2) = 1.0_DP - (2.0_DP* vector_V%xyz(1)**2) - (2.0_DP * vector_V%xyz(3)**2)
    q(3,3) = 1.0_DP - (2.0_DP* vector_V%xyz(1)**2) - (2.0_DP * vector_V%xyz(2)**2)

    q(1,2) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(2)) - 2.0_DP * (s*vector_V%xyz(3))
    q(1,3) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(3)) + 2.0_DP * (s*vector_V%xyz(2))

    q(2,1) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(2)) + 2.0_DP * (s*vector_V%xyz(3))
    q(2,3) = 2.0_DP * (vector_V%xyz(2)*vector_V%xyz(3)) - 2.0_DP * (s*vector_V%xyz(1))

    q(3,1) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(3)) - 2.0_DP * (s*vector_V%xyz(2))  
    q(3,2) = 2.0_DP * (vector_V%xyz(2)*vector_V%xyz(3)) + 2.0_DP * (s*vector_V%xyz(1))


    Matrix_T   = transpose(Matrix)  
    qMatrix_T  = matmul(q,Matrix_T) 
    qMatrix    = transpose(qMatrix_T)
  end subroutine Rotate_matrix_quaternion

  type(vector_t) function Rotate_by_quaternion( s, vector_V, vector_P) result(rotate_P)
    type(vector_t), intent(in) :: vector_V, vector_P
    real(DP), intent(in) :: s
    real(DP) :: q(3,3), r(3)

    q = 0.0_DP

    q(1,1) = 1.0_DP - (2.0_DP* vector_V%xyz(2)**2) - (2.0_DP * vector_V%xyz(3)**2)
    q(2,2) = 1.0_DP - (2.0_DP* vector_V%xyz(1)**2) - (2.0_DP * vector_V%xyz(3)**2)
    q(3,3) = 1.0_DP - (2.0_DP* vector_V%xyz(1)**2) - (2.0_DP * vector_V%xyz(2)**2)

    q(1,2) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(2)) - 2.0_DP * (s*vector_V%xyz(3))
    q(1,3) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(3)) + 2.0_DP * (s*vector_V%xyz(2))

    q(2,1) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(2)) + 2.0_DP * (s*vector_V%xyz(3))
    q(2,3) = 2.0_DP * (vector_V%xyz(2)*vector_V%xyz(3)) - 2.0_DP * (s*vector_V%xyz(1))

    q(3,1) = 2.0_DP * (vector_V%xyz(1)*vector_V%xyz(3)) - 2.0_DP * (s*vector_V%xyz(2))  
    q(3,2) = 2.0_DP * (vector_V%xyz(2)*vector_V%xyz(3)) + 2.0_DP * (s*vector_V%xyz(1))

    r(1) = vector_P%xyz(1); r(2) = vector_P%xyz(2); r(3) = vector_P%xyz(3)

    r = matmul(r,q)
    rotate_P%xyz(:) = r(:)
  end function Rotate_by_quaternion

  subroutine Define_rotation( vector_U, vector_V, vector_W, s ) 
    type(vector_t), intent(inout) :: vector_U, vector_V, vector_W
    real(DP), intent(inout)     :: s
    real(DP) :: nU, nV, nW, UV, theta

    nU = vector_U%norm()
    nV = vector_V%norm()
    UV = Scalar_product( vector_U, vector_V )

    theta = dacos(UV/(nU*nV))
    vector_W = Cross_product( vector_U, vector_V)
    nW = vector_W%norm()

    s = dcos(theta/2.0_DP)

    if ( nW > 0.0_DP ) then
       vector_W%xyz(:) = vector_W%xyz(:)/nW

       vector_W%xyz(1) = dsin(theta/2.0_DP)*vector_W%xyz(1)
       vector_W%xyz(2) = dsin(theta/2.0_DP)*vector_W%xyz(2)
       vector_W%xyz(3) = dsin(theta/2.0_DP)*vector_W%xyz(3)
    else    
       vector_W%xyz(1) = dsin(theta/2.0_DP)*vector_W%xyz(1)
       vector_W%xyz(2) = dsin(theta/2.0_DP)*vector_W%xyz(2)
       vector_W%xyz(3) = dsin(theta/2.0_DP)*vector_W%xyz(3)
    end if    
  end subroutine Define_rotation

  real(DP) function Get_norm_vector( this ) result(norm)
    class(vector_t), intent(in) :: this

    norm = 0.0_DP
    norm = dsqrt(sum(this%xyz(:) * this%xyz(:)))
  end function Get_norm_vector  

  real(DP) function Scalar_product( vector_U, vector_V ) result(dot)
    type(vector_t), intent(in) :: vector_U, vector_V

    dot = 0.0_DP
    dot = sum(vector_U%xyz(:) * vector_V%xyz(:))
  end function Scalar_product

  type(vector_t) function Cross_product( vector_U, vector_V ) result(vector_W)
    type(vector_t), intent(in) :: vector_U, vector_V

     vector_W%xyz(1) =   (vector_U%xyz(2)*vector_V%xyz(3)) - (vector_U%xyz(3)*vector_V%xyz(2))
     vector_W%xyz(2) =   (vector_U%xyz(3)*vector_V%xyz(1)) - (vector_U%xyz(1)*vector_V%xyz(3))
     vector_W%xyz(3) =   (vector_U%xyz(1)*vector_V%xyz(2)) - (vector_U%xyz(2)*vector_V%xyz(1))
  end function Cross_product

  real(DP) function Scalar_triple_product ( vector_U, vector_V, vector_W ) result(dot)  
    type(vector_t), intent(in) :: vector_U, vector_V, vector_W

    dot = 0.0_DP
    dot = Scalar_product(vector_W,Cross_product(vector_U,vector_V)) 
  end function Scalar_triple_product  

  subroutine Prevent_NAN( dot ) 
    real(kind=DP), intent(inout) :: dot
     
    if ( dot < -1.000000000000000_DP ) then                                 
         dot = -1.000000000000000_DP                                         
    elseif( dot > 1.000000000000000_DP) then                                
          dot = 1.000000000000000_DP                                         
    end if             
  end subroutine Prevent_NAN
end module kanon_vector_mod
