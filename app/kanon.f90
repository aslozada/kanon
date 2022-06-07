!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   Main program
!>@author Asdrubal Lozada
!>        Laboratory of Theoretical Chemistry - LQT
!>        Federal University of São Carlos
!>        <http://www.lqt.dq.ufscar.br>
!>@email  aslozada@gmail.com
!-----------------------------------------------------------------------------------
!   Copyright 2018-2022 Asdrubal Lozada
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
program kanon
  use kanon_kinds_mod
  use kanon_initialize_mod
  use kanon_tools_mod 
  use kanon_molecule_mod
  use kanon_delaunay_mod

  implicit none

  type(molecule_t) :: pattern
  type(molecule_t) :: image
  type(molecule_t) :: rotate_image
  type(molecule_t), allocatable :: trajectory(:)

  type(file_t) :: aligned_pattern
  type(file_t) :: aligned_image

  call Get_command_options()
  call Get_info_molecule()
  call Use_single_option()

  call aligned_pattern%output(33,"Aligned_pattern.xyz")
  call aligned_image%output(34,"Aligned_image.xyz")
  call Use_pattern_image()
  
  call Get_pattern_from_trajectory()
  call Get_image_from_trajectory()
  
  call Get_triangles()
  call Calc_torsion()

contains
  subroutine Get_info_molecule()
    real(kind=DP) :: diameter

    if ( is_pattern .and. (.not. is_image&
       & .and. .not. is_delaunay )) then
       call pattern%readmol( file_pattern, 1)
       diameter = Calc_Diameter( pattern )
       write(*,'(a)')'-----------------------------------------------------------------------------'
       write(*,'(a,f10.5,a)')'Diameter of molecule: ', diameter, ' ang.'

       call pattern%Point(1)
    end if
  end subroutine Get_info_molecule

  subroutine Use_single_option()
    real(kind=DP) :: diameter
    
    if ( is_single ) then
      call pattern%readmol( file_single, 1)
      diameter = Calc_diameter( pattern )
      
      write(*,'(a)')'-----------------------------------------------------------------------------'
      write(*,'(a,f10.5,a)')'Diameter of molecule: ', diameter, ' ang.'
      write(*,'(a)')'-----------------------------------------------------------------------------'
      
		call pattern%axis(1)

      image = Build_image( pattern )
      rotate_image = Optimal_rotation( pattern, image )

    end if
  end subroutine Use_single_option

  subroutine Use_pattern_image()
    real(kind=DP) :: diameter
    integer :: k

    if ( is_pattern .and. is_image ) then
       call pattern%readmol( file_pattern, 1)
       diameter = Calc_Diameter( pattern )
       write(*,'(a)')'-----------------------------------------------------------------------------'
       write(*,'(a,f10.5,a)')'Diameter of molecule: ', diameter, ' ang.'
       write(*,'(a)')'-----------------------------------------------------------------------------'
       call pattern%axis(1)

       call image%readmol( file_image, 1)
       diameter = Calc_Diameter( image )
       write(*,'(a)')'-----------------------------------------------------------------------------'
       write(*,'(a,f10.5,a)')'Diameter of image molecule: ', diameter, ' ang.'
       write(*,'(a)')'-----------------------------------------------------------------------------'
       call image%axis(1)

       rotate_image = Optimal_rotation( pattern, image )

       write(*,'(a)')'-----------------------------------------------------------------------------'
       write(*,'(a)')'Aligned molecules: Aligned_pattern.xyz / Aligned_image.xyz'
       write(*,'(a)')'-----------------------------------------------------------------------------'

       write(aligned_pattern%unit,'(i6)') pattern%natm
       write(aligned_pattern%unit,'(a)')'Aligned molecule: pattern. Builded by KANON.'
       write(aligned_image%unit,'(i6)') rotate_image%natm
       write(aligned_image%unit,'(a)')'Aligned molecule: image. Builded by KANON.'

       do k = 1, pattern%natm
          write(aligned_pattern%unit,'(a,3f10.5)') pattern%atoms(k)%sym,&
            & pattern%atoms(k)%xyz(:)
       end do

       do k = 1, rotate_image%natm
          write(aligned_image%unit,'(a,3f10.5)') rotate_image%atoms(k)%sym,&
            & rotate_image%atoms(k)%xyz(:)
       end do

    end if
    
  end subroutine Use_pattern_image

  subroutine Get_pattern_from_trajectory()
    use kanon_error_mod
    real(kind=DP) :: diameter
    integer :: frames, k


    if ( is_onlytraj ) then

        close(26)
        open(26, file='CHI_traj.dat', status='unknown')
        write(26,'(a)') '#Result CHI values along the trajectory.'
        write(26,'(a)')'#------------------------------------------------------------'
        write(26,'(a)')'# Frame     CHI'
        write(26,'(a)')'#'
        write(*,'(a)')'---------------------------------------------------------------------------------'
       
        write(*,'(a)')'Chirality index calculation along the trajetory: '//trim(file_trajectory%name)
        call Read_trajectory_xyz(trajectory, file_trajectory, frames) 
        write(*,'(a,i6)')'Total number of frames: ', frames
        write(*,'(a)')
        write(*,'(a)')'Calculating...'

          do k = 1, frames
             diameter = Calc_Diameter(trajectory(k))

             if ( is_verbose  ) then
                 write(*,'(a,i5,f10.5,a)')'Diameter of image: ',k, diameter, ' ang.'
             end if
             
             call trajectory(k)%axis(1)
             
             image = Build_image( trajectory(k) )

                rotate_image = Optimal_rotation(&
                & trajectory(k), image, frame=k, is_traj=1, utraj=26, measure=0)

                deallocate(rotate_image%atoms, stat=ierr)
                if ( ierr /= 0 ) call error%MsgError( 1, ERROR_DEALLOCATION )

          end do
    write(*,'(a)')'Check the output file: CHI_traj.dat'

    call Get_time()
    call error%MsgError(0, ERROR_OK)
    end if
    
  end subroutine Get_pattern_from_trajectory
  
  subroutine Get_image_from_trajectory()
    integer :: k, ierr=0
    integer :: frames
    real(kind=DP) :: diameter
    
	 if ( is_trajectory ) then
       if ( is_pattern) then
          call pattern%readmol( file_pattern, 1)
          call pattern%axis(1)
          
          call Read_trajectory_xyz(trajectory, file_trajectory, frames)

          close(26)
          open(26, file='RMSD_traj.dat', status='unknown')
          write(26,'(a)') '#Result RMSD values along the trajectory.'
          write(26,'(a)')'#----------------------------------------------------------------------------'
          write(26,'(a)')'# Frame     RMSD'
          write(26,'(a)')'#'
          
          write(*,'(a)')'-------------------------------------------------------'
          do k = 1, frames
             diameter = Calc_Diameter(trajectory(k))
             write(*,'(a,i5,f10.5,a)')'Diameter of image: ',k, diameter, ' ang.'
             call trajectory(k)%axis(1)

             rotate_image = Optimal_rotation(&
              & pattern, trajectory(k), frame=k, is_traj=1, utraj=26, measure=1)
             deallocate(rotate_image%atoms, stat=ierr)
             if ( ierr /= 0 ) call error%MsgError( 1, ERROR_DEALLOCATION )
          end do

       end if
    end if
  end subroutine Get_image_from_trajectory

  subroutine Get_triangles()
    logical :: wrong
    if (is_delaunay) then
       call pattern%readmol(file_pattern,0)
       call Get_delaunay_triangles(pattern,wrong)

       if (.not. wrong ) then
          write(*,'(a)')'Copy '//trim(file_triangle%name)//' to index1 or to index2&
            & for pattern or image, respectively.' 
       end if
       call error%MsgError( 0, message_error='Sucessful termination.')
    end if
  end subroutine Get_triangles


  subroutine Calc_torsion ( )
    
    if ( is_torsion) then
      select case(type_torsion)
      case(1)
        write(*,'(a)')'Torsion on surface / Monolayer mode'
        call Calc_surface_torsion( pattern, image, file_index1)
      case(2)
        write(*,'(a)')'Torsion on surface / Bilayer mode'
        call Calc_surface_torsion( pattern, image, file_index1)
        call Calc_surface_torsion( pattern, image, file_index2)
      end select
    end if

  end subroutine Calc_torsion

end program kanon
