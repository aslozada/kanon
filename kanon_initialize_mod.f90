!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_initialize_mod
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
module kanon_initialize_mod
  use kanon_error_mod, only:&
    & error,&
    & ERROR_INVALID_OPTION,&
    & ERROR_ALLOCATION,&
    & ERROR_INVALID_ARGUMENT
  use kanon_tools_mod, only:&
    & file_t,&
    & file_pattern, file_image,&
    & file_fragment, file_single,&
    & file_index1, file_index2,&
    & file_triangle,&
    & file_trajectory
  use kanon_commands_mod
  use kanon_lexical_mod

  implicit none
  private
  
  public :: Get_command_options 

  !>@group Initialize commands
  logical, public ::&
   & is_fragment    = .false.,&
   & is_image       = .false.,&
   & is_pattern     = .false.,&
   & is_single      = .false.,&
   & is_torsion     = .false.,&
   & is_trajectory  = .false.,&
   & is_monolayer   = .false.,&
   & is_bilayer     = .false.,&
   & is_delaunay    = .false.,&
   & is_bohr        = .false.,&
   & is_verbose     = .false.,&
   & is_patraj      = .false.,&
   & is_onlytraj    = .false.,&
   & is_performance = .false.

   integer, public ::&
     & type_torsion = 0
	 

   character(len=5), parameter, public :: PROG_NAME = 'kanon'
   character(len=6), parameter, public :: VERSION   = '1.2.0'
   character(len=3), public    :: func  = ''

contains
  subroutine Get_command_options()
    character(len=1) :: opt
    character(len=:), allocatable :: foo
    integer :: ios = 0

    call Build_long_options( type_command )

    if ( argc > command_argument_count() ) then
        call error%MsgError( 1, message_error='no options: '//'try '//PROG_NAME//' --help' )
    end if

    do
       if ( argc > command_argument_count() ) exit
          opt = type_command%getopt(arg,'hvf:i:p:s:t:d:mBT:bVP')

          select case(opt)
          case('f')  
            is_fragment  = .true. 
            file_fragment%name = trim(arg%string)
            file_fragment%unit = 40
          case('h') 
            call Show_help()
          case('i') 
            is_image = .true. 
            file_image%name = trim(arg%string)
            file_image%unit = 20
            call Check_name_file(file_image%name)
          case('p') 
            is_pattern = .true.
            file_pattern%name = trim(arg%string)
            file_pattern%unit = 10
            call Check_name_file(file_pattern%name)
          case('s')
            is_single = .true.
            file_single%name = trim(arg%string) 
            file_single%unit = 30
            call Check_name_file(file_single%name)
          case('T') 
            is_torsion = .true. 
            foo = trim(arg%string)
            read(foo,*, iostat=ios) type_torsion
            if ( ios /= 0 ) call error%MsgError( 1, ERROR_INVALID_ARGUMENT )
          case('m') 
            is_monolayer = .true.
          case('b')  
            is_bilayer = .true.
          case('d')
            is_delaunay = .true.
            file_triangle%name = trim(arg%string)
            file_triangle%unit = 12
          case('t') 
            is_trajectory = .true.
            file_trajectory%name = trim(arg%string) 
            file_trajectory%unit = 93
            call Check_name_file(file_trajectory%name)
          case('B')  
            is_bohr = .true.
          case('V')
            is_verbose = .true. 
          case('P')
            is_performance = .true.  
          case('v') 
            call Show_license
          case('?')
            call error%MsgError( 1, ERROR_INVALID_OPTION ) 
          case default   
            call error%MsgError( 1, ERROR_INVALID_OPTION ) 
          end select

    end do

    call Check_arguments()
    call Get_header()
  end subroutine Get_command_options

  subroutine Check_arguments()

    if ( is_single .and. is_pattern )&
      & call error%MsgError( 1, message_error='Ambiguous pattern molecule.')

    if ( is_trajectory .and. is_image )&
      & call error%MsgError( 1, message_error='Duplicated image. Use trajectory xor image')

    if ( is_trajectory .and. is_single )&
      & call error%MsgError( 1, message_error='Ambiguous image definition.')

    if ( is_trajectory .and. is_pattern )&
      & is_patraj =.true.

    if ( is_trajectory .and. .not. is_pattern )&
      & is_onlytraj = .true.
    
    if ( is_bohr .and. (.not. is_pattern))&
     & call error%MsgError( 1, message_error='missing arguments in command line.'&
            &//' try: '//PROG_NAME//' --help' )
    
    if ( is_image .and. (.not. is_pattern))&
     & call error%MsgError( 1, message_error='missing arguments in command line.'&
            &//' try: '//PROG_NAME//' --help' )
    
    if ( is_delaunay .and. (.not. is_pattern) )&
     & call error%MsgError( 1, message_error='missing arguments in command line.'&

            &//' try: '//PROG_NAME//' --help' )
    
    if ( is_torsion .and. (.not. is_pattern .and. .not. is_image)) then
       call error%MsgError( 1, message_error='missing arguments in command line.'&
            &//' try: '//PROG_NAME//' --help' )
    end if

    if ( is_monolayer .and. (.not. is_pattern .and. .not. is_image )) then
       call error%MsgError( 1, message_error='missing arguments in command line.'&
            &//' try: '//PROG_NAME//' --help' )
    end if

    if ( is_bilayer .and. (.not. is_pattern .and. .not. is_image)) then
       call error%MsgError( 1, message_error='missing arguments in command line.'&
            &//' try: '//PROG_NAME//' --help' )
    end if
    
    if ( is_torsion .and. is_pattern .and. is_image) then
       file_index1%unit = 75
       file_index1%name = 'index1.dat'
       call file_index1%inquire()
 
       if ( type_torsion == 2 ) then
          file_index2%unit = 76
          file_index2%name = 'index2.dat'
          call file_index2%inquire()
       end if
       
    end if

  end subroutine Check_arguments

  subroutine Check_name_file(string)
    character(*), intent(inout) :: string

    string = trim(string)
    select case(string)
    case('Pattern.xyz')
        call error%MsgError(1, message_error='Reserved name file: '//trim(string))
    case('Rotate_image.xyz')
        call error%MsgError(1, message_error='Reserved name file: '//trim(string))
    end select

  end subroutine Check_name_file

  subroutine Build_long_options( this )
    class(command_t), intent(inout) :: this
    integer :: ierr = 0

    allocate(this%options(14), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    this%options(1)   = option_t('fragment',   .true., 'f')
    this%options(2)   = option_t('help',       .false.,'h')
    this%options(3)   = option_t('image',      .true., 'i')
    this%options(4)   = option_t('pattern',    .true., 'p')
    this%options(5)   = option_t('single',     .true., 's')
    this%options(6)   = option_t('version',    .false.,'v')
    this%options(7)   = option_t('torsion',    .true., 'T')
    this%options(8)   = option_t('trajectory', .true.,'t')
    this%options(9)   = option_t('monolayer',  .false.,'m')
    this%options(10)  = option_t('bilayer',    .false.,'b')
    this%options(11)  = option_t('delaunay',   .true., 'd')
    this%options(12)  = option_t('bohr',       .false.,'B')
    this%options(13)  = option_t('verbose',    .false.,'V')
    this%options(14)  = option_t('performance',.false.,'P')
  end subroutine Build_long_options  

  subroutine Show_help()

    write(*,'(a)')'Usage: '//PROG_NAME//' [options]'//' <arguments>'
    write(*,'(a)')'Example: '//PROG_NAME//' -p'//' pattern.xyz'//' -i'//' image.xyz'
    write(*,'(a)')
    write(*,'(a)')'Mandatory arguments to long options are mandatory to short options too.'
    write(*,'(a)')'Options:'
    write(*,'(a)') '  --help,-h                    '//'Display this help.'
    write(*,'(a)') '  --version,-v                 '//'Display the version of '//PROG_NAME//'.'
    write(*,'(a)') '  --single,-s <file>           '//'Calculate the chirality index from <file>.'
    write(*,'(a)') '  --pattern, -p <file>         '//'Input the pattern molecule xyz file.'
    write(*,'(a)') '  --image, -i <file>           '//'Input an image molecule xyz file.'

	 write(*,'(a)') 'Advanced use'
    write(*,'(a)') '  --torsion, -T <type>         '//'Active torsion surface calculation  [INT]'
    write(*,'(a)') '  --monolayer, -m              '//'Active monolayer option for torsion calculation.'             
    write(*,'(a)') '  --bilayer, -b                '//'Active bilayer option for torsion calculation' 
    write(*,'(a)') '  --delaunay -d                '//'Build a set triangles on surface (2D).'
    write(*,'(a)') '  --trajectory, -t             '//'Read a cartesian coordiantes trajectory.'
    write(*,'(a)') '  --verbose, -V                '//'Display extra information.'

    if ( is_verbose ) then
        write(*,'(a)')'----------------------------------------------------------------------------------------'
        write(*,'(a)')'KANON this a program to Hausdorff distance-derived calculation an chirality Index'
        write(*,'(a)')'the Hausdorff chirality measure is defined as'
        write(*,'(a)')'                                H(Q) = hmin(QQ*)/d(Q)'
        write(*,'(a)')'where hmin is the Hausdorff to the optimal overlap and d(Q) represent'
        write(*,'(a)')'the diameter of molecule. See: Buda A. B and Mislow K. J. Am. Chem. Soc. 1992, 114.'
        write(*,'(a)')'---------------------------------------------------------------------------------------'
        write(*,'(a)')'KANON too calculates some symmetry properties and performs inertia moments calculation,'
        write(*,'(a)')'build inertia tensor, and performs alignment between O and O*.'
        write(*,'(a)')'KANON build handedness measure that verify torsion/deflection'
        write(*,'(a)')'of molecular surfaces.'
        write(*,'(a)')'---------------------------------------------------------------------------------------'

    end if

    if ( .not. is_verbose ) then
       write(*,'(a)')
       write(*,'(a)')'For further informations and examples use: '//PROG_NAME//' [--verbose|-V] [--help|-h]'
    end if   
    write(*,'(a)')
    write(*,'(a)')'Report bugs to: aslozada@gmail.com'
    write(*,'(a)')'kanon: <http://www.lqt.dq.ufscar.br>'

    stop
  end subroutine Show_help

  subroutine Show_license()

    write(*,'(a)') PROG_NAME//' '//VERSION
    write(*,'(a)') 'Copyright 2018 Asdrubal Lozada'
    write(*,'(a)') 'License GPLv3+: GNU GPL version 3 or later &
                     &<http://gnu.org/license/gpl.html>.'
    write(*,'(a)') 'This a free software: you are free to change&
                     & and redistribute it.'
    write(*,'(a)') 'There is NO WARRANTY, to the extent permited by law.'
    write(*,*)
    write(*,'(a)') 'Written by Asdrubal Lozada'
    write(*,'(a)') 'Laboratory of Theoretical Chemistry, LQT -- UFSCar'
	 write(*,'(a)') '<http://www.lqt.dq.ufscar.br>'
    write(*,'(a)') 'e-mail: aslozada@gmail.com'

    stop

  end subroutine Show_license

  subroutine Get_header()
 
   write(*,'(a)') '---------------------------------------------------------------------------------'
   write(*,'(a)') '                                       KANON'
   write(*,'(a)')''
   write(*,'(a)') '                                by Asdrubal Lozada'
   write(*,'(a)')'                         Laboratory of Theoretical Chemistry'
   write(*,'(a)')'                           Federal University of São Carlos'
   write(*,'(a)')'                                       Brazil'
   write(*,'(a)')''
   write(*,'(a)')'                             Program Version 1.2.0 RELEASE'
   write(*,'(a)') '---------------------------------------------------------------------------------'
   write(*,'(a)')
   write(*,'(a)')'This program calculates an Hausdorff`s distance-based chirality Index'
   write(*,'(a)')'(see Mol. Phys. 107, 281)*.'
   write(*,'(a)')'The exclude intersection algorithm is used to calculate the Hausdorff`s distance'
   write(*,'(a)')'(see Trans. Patt. A. Mach. Int. 37, 2153)'
   write(*,'(a)')'This program uses J. Burkardt`s source code: table_delaunay.f90'
   write(*,'(a)')'see: https://people.sc.fsu.edu/~jburkardt/f_src/table_delaunay/table_delaunay.html'
   write(*,'(a)')
   write(*,'(a)') '---------------------------------------------------------------------------------'
  end subroutine Get_header
end module kanon_initialize_mod
