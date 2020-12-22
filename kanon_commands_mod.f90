!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_commands_mod 
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
module kanon_commands_mod
  use kanon_error_mod, only: &
    & error,&
    & ERROR_INVALID_OPTION,&
    & ERROR_INVALID_CHARACTER
  use kanon_lexical_mod, only: &
    & string_t,&
    & buffer
  
  implicit none
  private

  integer, public ::&
    & argc = 1,&
    & pos  = 2
  character(len=1), public ::&
    & copt

  type, public :: option_t
    character(len=40) :: string
    logical           :: has
    character(len=1)  :: value
  end type option_t

  !>@class command
  !>@brief Get options from command line
  type, public :: command_t
    type(option_t), allocatable, dimension(:) :: options
  contains
    procedure :: Getopt => Get_option 
    procedure :: Longs  => Check_long_option 
  end type command_t

  type(command_t), public :: type_command

  type, public :: keyword_t
     type(option_t), allocatable, dimension(:) :: options
  end type keyword_t

contains

  !>@params[in/out] this
  !>@params[in/out] string type string
  !>@params[in] short
  character function Get_option( this, string, short ) result(opt)
    class(command_t), intent(inout) :: this
    type(string_t), intent(inout)   :: string
    character(*), intent(in)        :: short
    logical                         :: is
    integer                         :: n, j, m

    call get_command_argument( argc, buffer )
    string%string = trim(buffer)

    n = 1
    do
      if( n > len_trim(string%string) ) exit
      is = string%IsPrompt(n) 
      if (.not. is) call error%MsgError(&
       & 1, message_error='Invalid character in command line: '//string%string(n-1:n-1) )
    end do

    if ( argc == 1 ) then
       if ( string%string(1:1) /= '-' ) call error%MsgError( 1, &
         & message_error='Invalid option in command line: '//string%string)
    end if

    n = 2; j = 3; m = 1
    if ( len_trim(string%string) >= 2 ) then

       if( string%string(1:2) == '--' .and. string%IsAlpha(j) ) then
         opt = this%Longs(string)
       elseif( string%string(1:1) == '-' .and. string%IsAlpha(n) ) then

         if ( len_trim(string%string) > 2 ) call error%MsgError( 1, &
           & message_error='Invalid option in commman line: '//trim(string%string))

         opt = Check_short_option( string, short )
       else
         call error%MsgError( 1, message_error='Invalid option in command line: '//trim(string%string))  
       end if

    elseif(len_trim(string%string)==1 .and. string%IsRprompt(m)) then
      call error%Msgerror( 1, message_error='Reserved character: '//string%string(1:1))
    elseif(len_trim(string%string)==1) then
      call error%MsgError( 1, message_error='Invalid option in command line: '//trim(string%string))  
    end if
  end function Get_option

  character function Check_long_option( this, string ) result(opt)
    class(command_t), intent(inout) :: this
    type(string_t),   intent(inout) :: string
    integer :: n, j, m

    opt = '?'

    do n = 1, size(this%options)

       if ( string%string(3:len_trim(string%string)) == this%options(n)%string ) then
          opt = this%options(n)%value

          if ( this%options(n)%has ) then

             if ( command_argument_count() < (argc + 1) ) then
                call error%MsgError( &
                 & 1, message_error='The option '//trim(this%options(n)%string)//' requires an argument.')
             end if 

             argc = argc + 1
             call get_command_argument(argc, buffer)
             string%string = buffer
             
             if ( len_trim(string%string) < 1 ) then
                call error%MsgError( 1, message_error='too few arguments' )
             end if

             m = 1

             if( len_trim(string%string)==1 .and. string%IsRprompt(m)) then
                call error%Msgerror( 1, message_error='Reserved character: '//string%string(1:1))
             end if

             m = 1
             do 
               if ( m > len_trim(string%string) ) exit
                 
               if ( .not. string%IsPrompt(m) ) then
                   call error%MsgError( 1, ERROR_INVALID_CHARACTER )
               end if
             end do

             j = 2

             if ( len_trim(string%string) >=2 ) then
                if ( (string%string(1:1) == '-' .and. string%isalpha(j)) .or.&
                     & string%string(1:2) == '--' ) then
                     call error%MsgError( 1, message_error='too few arguments. '&
                     &//trim(string%string)//':'//' reserved expression.' )
                end if
             end if

             argc = argc + 1
             exit

           else 
             argc = argc + 1
          end if
       end if

       if ( opt /= '?' ) exit
        
    end do

    if ( opt == '?' ) then
       call error%MsgError( 1, ERROR_INVALID_OPTION )
    end if
    
  end function Check_long_option

  character function Check_short_option( string, short ) result(opt)
    type(string_t), intent(inout) :: string
    character(*), intent(in)      :: short
    integer :: k, n

    opt = string%string(2:2)

    k = index(short,opt)

    if ( k == 0 ) then
       argc = argc + 1
       opt = '?'
    else

       if ( k < len_trim(short) ) then
           

          if ( short(k+1:k+1) == ':' ) then
             
              if ( command_argument_count() < (argc + 1) ) then
                 call error%MsgError( 1, message_error='The option '//&
                 &opt//' requires an argument.')
              end if

              argc = argc + 1
             
               call get_command_argument(argc, buffer)
               string%string = trim(buffer)

               if ( len_trim(string%string) < 1 ) then
                  call error%MsgError( 1, message_error='too few arguments.' )
               end if
               
               n = 1

               if( len_trim(string%string)==1 .and. string%IsRprompt(n)) then
                  call error%Msgerror( 1, message_error='Reserved character: '//string%string(1:1))
               end if

               n = 1

               do 
                 if ( n > len_trim(string%string) ) exit
                 
                 if ( .not. string%IsPrompt(n) ) then
                    call error%MsgError( 1, ERROR_INVALID_CHARACTER )
                 end if
               
               end do

               n = 2
               if ( string%string(1:1) == '-' .and. string%IsAlpha(n) ) then
                  call error%MsgError( 1, message_error='too few arguments. '&
                  &//trim(string%string)//':'//' reserved expression.')
               end if

               argc = argc + 1

          else

               argc = argc + 1

               call get_command_argument(argc, buffer)
               string%string = trim(buffer)

               if( len_trim(string%string) >= 1 ) then
                 if (string%string(1:1) /= '-') call error%MsgError( 1,&
                  & message_error='Invalid argument in command line.' )
               end if

          end if
       else
           argc = argc + 1
           if ( k == len_trim(short) ) then   
              call get_command_argument(argc,buffer)
              string%string = trim(buffer)

              if (len_trim(string%string) >= 1 ) then
                 if (string%string(1:1) /= '-') call error%MsgError(1,&
                   & message_error='Invalid argument in command line.')
              end if

           end if    
       end if
    end if

  end function Check_short_option

end module kanon_commands_mod
