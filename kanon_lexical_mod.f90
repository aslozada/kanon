!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_lexical_mod
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
module kanon_lexical_mod

  implicit none
  private

  character(len=256), public ::   &
     & buffer
  character(len=:), allocatable :: foo 

  !>@class string
  type, public :: string_t
    character(len=:), allocatable :: string
  contains
	 !>@params[in/out] this
	 !>@params[in/out] pos Position of current character 
    procedure :: GetChar       => Get_character
    procedure :: IsAlpha       => Character_is_alpha
    procedure :: IsDigit       => Character_is_digit
    procedure :: IsASCII       => Character_is_ascii
    procedure :: IsSpace       => Character_is_space
    procedure :: IsOperator    => Character_is_operator
    procedure :: IsDelimiter   => Character_is_delimiter
    procedure :: IsPrompt      => Character_is_prompt
    procedure :: IsRprompt     => Character_is_restricted_prompt
    procedure :: IsBlank       => Character_is_blank
    procedure :: RmBlank       => Remove_repeated_blank
  end type string_t

  !>@note arg Global variable
  type(string_t), public :: arg

  !>@brief Inquire by regular expresions:
  !>       [a-z|A-Z] [0-9]{./-+*}

contains
  character function Get_character( this, pos ) result(ch)
    class(string_t), intent(inout) :: this
    integer, intent(inout)         :: pos

    this%string = trim(this%string)
    ch = this%string(pos:pos)
    pos = pos + 1 
  end function Get_character  

  logical function Character_is_alpha( this, pos ) result(is)
    class(string_t), intent(inout) :: this
    integer, intent(inout)         :: pos
    character(len=1)               :: ch

    is = .false. 
 
    this%string = trim(this%string)
    ch = this%string(pos:pos)
    is =((lge(ch,'a').and.lle(ch,'z')) .or. (lge(ch,'A').and.lle(ch,'Z')))
    pos = pos + 1
  end function Character_is_alpha  
 
  logical function Character_is_digit( this, pos ) result(is)
    class(string_t), intent(inout) :: this
    integer, intent(inout)         :: pos
    character(len=1)               :: ch

    is = .false.
   
    this%string = trim(this%string)
    ch = this%string(pos:pos)
    is =(lge(ch,'0').and.lle(ch,'9').or.(ch==achar(46))) 
    pos = pos + 1
  end function Character_is_digit  
 
  logical function Character_is_ascii( this, pos ) result(is)
    implicit none
    class(string_t), intent(inout) :: this
    integer, intent(inout)         :: pos
    character(len=1)               :: ch

    is = .false.

    this%string = trim(this%string)
    ch = this%string(pos:pos)
    is =(lge(ch,achar(0)).and.lle(ch,achar(127)))
    pos = pos + 1 
  end function Character_is_ascii 

  logical function Character_is_space( this, pos ) result(is)
    implicit none
    class(string_t), intent(inout) :: this
    integer, intent(inout)         :: pos
    character(len=1)               :: ch

    is = .false.

    this%string = trim(this%string)
    ch = this%string(pos:pos)
    is =(ch == achar(9)  .or. ch == achar(10) .or.&
       &ch == achar(11) .or. ch == achar(11) .or.&
       &ch == achar(13) .or. ch == achar(32))
    pos = pos + 1        
  end function Character_is_space

  logical function Character_is_operator( this, pos ) result(is)
    implicit none
    class(string_t), intent(inout)  :: this
    integer, intent(inout)          :: pos
    character(len=1)                :: ch

    is = .false.

    this%string = trim(this%string)
    ch = this%string(pos:pos)
    is =(ch == achar(42) .or. ch == achar(43) .or.&
       &ch == achar(45) .or. ch == achar(60) .or.&
       &ch == achar(61) .or. ch == achar(62))
    pos = pos + 1  
  end function Character_is_operator

  logical function Character_is_delimiter( this, pos ) result(is)
    implicit none
    class(string_t), intent(inout) :: this
    integer, intent(inout)         :: pos
    character(len=1)               :: ch

    is = .false.

    this%string = trim(this%string)
    ch = this%string(pos:pos)
    is =(ch == achar(44) .or. ch == achar(58) .or.&
        &ch == achar(59) .or. ch == achar(63))
    
    pos = pos + 1
  end function Character_is_delimiter

  logical function Character_is_prompt( this, pos ) result(is)
     class(string_t), intent(inout) :: this
     integer, intent(inout)         :: pos
     character(len=1)               :: ch

     is = .false.

     this%string = trim(this%string)
     ch = this%string(pos:pos)
     is =((lge(ch,'a').and.lle(ch,'z')).or.&
          &(lge(ch,'A').and.lle(ch,'Z')).or.&
          &(lge(ch,'0').and.lle(ch,'9')).or.& 
          &(ch == achar(45)).or.(ch == achar(46)).or.&
          &(ch == achar(95)).or.(ch == achar(63)).or.&
          &(ch == achar(61)).or.(ch == achar(47)))
     pos = pos + 1
  end function Character_is_prompt    

  logical function Character_is_restricted_prompt( this, pos ) result(is)
     class(string_t), intent(inout) :: this
     integer, intent(inout)         :: pos
     character(len=1)               :: ch

     is = .false.

     this%string = trim(this%string)
     ch = this%string(pos:pos)
     is =( (ch == achar(45)).or.(ch == achar(46)).or.&
          &(ch == achar(95)).or.(ch == achar(63)).or.&
          &(ch == achar(61)).or.(ch == achar(47)))
     pos = pos + 1
  end function Character_is_restricted_prompt   

  logical function Character_is_blank( this, pos ) result(is)
    class(string_t), intent(inout) :: this
    integer, intent(inout)         :: pos
    character(len=1)               :: ch

    is = .false.
    ch = this%string(pos:pos)
    is = (ch == achar(32))
    pos = pos + 1
  end function Character_is_blank

  subroutine Remove_repeated_blank( this )
    class(string_t), intent(inout)  :: this
    integer :: k, j
    character(len=1) :: c, cb

    foo = ''

    foo = trim(this%string)
    k = scan(foo, ' ')
    this%string = foo(1:k)

    do j = k+1, len_trim(foo)
       c  = foo(j:j)
       cb = foo(j-1:j-1)

       if ( c /= ' ' .or. (c == ' ' .and. cb /= ' ')) then
          this%string = this%string//c
       end if
    end do
  end subroutine Remove_repeated_blank

end module kanon_lexical_mod
