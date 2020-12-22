!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_error_mod
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
module kanon_error_mod
    
  implicit none
  private
  
  !>@brief Error handling  
  integer, parameter, public ::    &
    & ERROR_OK                 = 0,&
    & ERROR_INVALID_OPTION     = 1,&
    & ERROR_NOT_ARGUMENT       = 2,&
    & ERROR_INVALID_CHARACTER  = 3,&
    & ERROR_NOT_FILE           = 4,&
    & ERROR_NOT_UNIT           = 5,&
    & ERROR_ALLOCATION         = 6,&
    & ERROR_DEALLOCATION       = 7,&
    & ERROR_READING            = 8,&
    & ERROR_INVALID_ARGUMENT   = 9

  !>@class error
  type error_t
    character(len=:), allocatable :: message
  contains
    procedure :: MsgError => Print_error_message
  end type error_t

  type(error_t), public :: error

contains
  !>@param[in/out] this	
  !>@param[in] type_error
  !>@param[in] code_error
  !>@param[in] message_error
  subroutine Print_error_message ( this, type_error, code_error, message_error)
    class(error_t), intent(inout)      :: this
    integer, intent(in)                :: type_error
    integer, intent(in), optional      :: code_error
    character(*), intent(in), optional :: message_error    

    this%message = char(0)

    if ( present(code_error) ) then
       select case(code_error)
       case(0); this%message = 'Successful termination.'
       case(1); this%message = 'Invalid option.'
       case(2); this%message  = 'Missing arguments.'
       case(3); this%message = 'Invalid character.'
       case(4); this%message = 'No such file.'
       case(5); this%message = 'Unit in use.'
       case(6); this%message = 'Memory allocation failed.'
       case(7); this%message = 'Memory deallocation failed.'
       case(8); this%message = 'An error ocurred while reading.'
       case(9); this%message = 'Invalid argument in command line.'  
       end select
    end if

    if ( present(message_error) .and. .not. present(code_error) ) then
         this%message = trim(message_error)
    elseif( present(code_error) .and. present(message_error) ) then   
         this%message = this%message//' '//trim(message_error)
    end if

    select case(type_error)
    case(0)
      this%message = this%message
    case(1)
      this%message = 'Error: '//this%message
    case(2)
      this%message = 'Warning: '//this%message
    end select

    write(*,'(a)') trim(this%message)
    if ( type_error == 1 ) stop 
  end subroutine Print_error_message 

end module kanon_error_mod
