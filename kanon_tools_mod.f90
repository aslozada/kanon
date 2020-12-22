!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_tools_mod
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
module kanon_tools_mod
  use kanon_error_mod, only:&
    & error,&
    & ERROR_NOT_FILE,&
    & ERROR_READING
  use kanon_lexical_mod, only:&
    & buffer
  use kanon_kinds_mod
    
    implicit none
    private

	 !>@class file
    type, public :: file_t
      integer :: unit
      character(len=:), allocatable :: name
    contains
      procedure :: Inquire => Request_file
      procedure :: Output  => Set_output_file
    end type file_t

    integer, public  ::&
     & ios  = 0,&
     & ierr = 0

    type(file_t), public ::&
     & file_pattern,&
     & file_image,&
     & file_fragment,&
     & file_single,&
     & file_index1,&
     & file_index2,&
     & file_triangle,&
     & file_trajectory

     public ::&
      & Get_time,&
      & Parallel_matmul,&
      & Get_timing

  !>@params[in/out] this
  !>@params[in] number_unit
  !>@params[in] name_file
contains 
  subroutine Set_output_file( this, number_unit, name_file )
    class(file_t), intent(inout) :: this
    integer, intent(in) :: number_unit
    character(*), intent(in) :: name_file

    this%unit = number_unit
    this%name = trim(name_file)

    close(this%unit)
    open(unit=this%unit, file=this%name, status='unknown')
    
  end subroutine Set_output_file

  subroutine Request_file( this )
    class(file_t), intent(inout) :: this
    logical  ::&
     & unit_open = .false. ,&
     & file_exists = .false.,&
     & file_open = .false.

    close(this%unit)
    inquire ( unit = this%unit, opened = unit_open, iostat = ios )
    if ( ios == 0 ) then
       if ( .not. unit_open ) then
          inquire ( file = trim(this%name), exist = file_exists, iostat = ios,&
           & opened = file_open )
          if ( ios == 0 ) then
             if ( file_exists ) then
                if ( .not. file_open ) then 
                   open ( unit = this%unit, file = trim(this%name), status = 'old',&
                        &    form = 'formatted' )
                        rewind( unit = this%unit)
                else
                   call error%MsgError( 1, message_error='File: '//trim(this%name)//' is in use.&
                     & Check the name file.')         
               end if        
             else  
               ios = 9999
               call error%MsgError( 1, message_error='no such file: '//trim(this%name) )
             end if
          else
            ios = 9999
            call error%MsgError( 1, ERROR_READING )
          end if
       else
         call error%MsgError( 1, message_error='this unit is in use.' )
       end if
    else
      ios = 9999
      call error%MsgError( 1, ERROR_READING )
    end if      
  end subroutine Request_file


  subroutine Get_time
    character(8)  :: date
    character(10) :: time
    
    call date_and_time(date,time)
    write(*,'(a,2x,a,2x,a)') date, time
  end subroutine Get_time

  !>@brief Define parallel matrix multiplication 
  !>       openmp
  subroutine Parallel_matmul(A, B, C, mm, nn, ll)
    real(kind=DP), intent(in) :: A(:,:), B(:,:)
    real(kind=DP), intent(inout) :: C(:,:)
    integer, intent(in), optional :: mm, nn, ll
    integer :: k, i, j
    real(kind=DP) :: tmp_mat

    !$omp parallel do private(tmp_mat,i,j,k)
    do i = 1, mm   
       do j = 1, nn
          tmp_mat = 0.0_DP
          do k = 1, ll
             tmp_mat = tmp_mat + (A(i,k) * B(k,j))
          end do
          C(i,j) = tmp_mat
       end do
    end do
    !$omp end parallel do

  end subroutine Parallel_matmul

  subroutine Get_timing(start, finish)
    real(kind=DP), intent(in) :: start, finish
    write(*,'(a,f12.8,a)') 'Time: ',(finish-start),' sec.'
  end subroutine Get_timing 

end module kanon_tools_mod
