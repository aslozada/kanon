!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_molecule_mod 
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
module kanon_molecule_mod
  use kanon_kinds_mod, only:&
    & DP
  use kanon_constants_mod
  use kanon_error_mod, only:&
    & error,&
    & ERROR_ALLOCATION,&
    & ERROR_DEALLOCATION,&
    & ERROR_READING,&
    & ERROR_INVALID_CHARACTER
  use kanon_lexical_mod, only:&
    & buffer
  use kanon_tools_mod
  use kanon_eigen_mod
  use kanon_vector_mod
  use kanon_initialize_mod!, only:&
 !   & is_bohr

  implicit none
  private

  type atom_t
    real(DP)         :: xyz(3)
    character(len=2) :: sym
    real(DP)         :: mass
  end type atom_t

  !>@class molecule
  type, public  :: molecule_t
    type(atom_t), allocatable, dimension(:) :: atoms
    integer   :: natm
    real(DP)  :: com(3) 
  !>@params[in/out] type molecule
  !>@brief molecule is enter as XYZ file format
  contains
    procedure :: Malloc   => Allocate_molecule
    procedure :: Readmol  => Read_molecule_xyz 
    procedure :: Getcom   => Get_center_of_masses  
    procedure :: Inertia  => Get_inertia_moments
    procedure :: Axis     => Rotate_molecule_to_axis
    procedure :: Point    => Get_point_group
  end type molecule_t

  public ::&
   & Build_image,&
   & Optimal_rotation,&
   & Calc_RMSD,&
   & Calc_Hausdorff,& 
   & Calc_Diameter,&
   & Calc_surface_torsion,&
   & Read_trajectory_xyz

  integer, public ::&
  & natm_on,&
  & frames_on

contains

  subroutine Calc_surface_torsion( pattern, image, index_list )
    type(molecule_t), intent(in) :: pattern, image
    type(file_t), intent(inout) :: index_list
    
    type(vector_t) :: U, V, W, VU, WU, ONP 
    type(vector_t) :: U2, V2, W2, VU2, WU2, ONP2 
    integer :: i, j, m, n, ios=0
    integer, allocatable, dimension(:,:) :: ver
    real(DP) :: nVU, nWU, nONP, nVU2, nWU2, nONP2, dot

    type(vector_t) ::  me1, ame1, bme1, cfme1, VR, VRR
    real(DP) :: name1, nbme1, ncfme1, dot1, nVR, nVRR, dot2

    integer :: s1, s2, s3, ss, mn, itera, ierr=0
    integer, allocatable, dimension(:) :: signs, idlm2
    real(kind=DP), allocatable, dimension(:) :: angle, draw
    real(kind=DP) :: sum, counter, diff, eps0
    character(len=4) :: foo
    integer :: ref(3)

    n = pattern%natm; m = image%natm 

    rewind(index_list%unit)
    read(index_list%unit,*,iostat=ios) foo, mn
    if ( ios /= 0 ) call error%MsgError( 1, message_error='Error reading > Check index file: '//index_list%name )
    read(index_list%unit,*,iostat=ios) foo, ref(:)
    if ( ios /= 0 ) call error%MsgError( 1, message_error='Error reading > Check index file: '//index_list%name )

    write(*,'(a)')'---------------------------------------------------------------------------------'
    write(*,'(a,i5)')'Number of triangles: ', mn
    write(*,'(a,3i3)')'Reference vertices: ', ref(:)

    if ( ref(1) > n .or. ref(2) > n .or. ref(3) > n) &
      &call error%MsgError( 1, message_error='Invalid reference vertice > Check index file: '//index_list%name)
    
    if ( ref(1) == 0 .or. ref(2) == 0 .or. ref(3) == 0 ) &
      &call error%MsgError( 1, message_error='Invalid reference vertice > Check index file: '//index_list%name)

    if ( (ref(1)==ref(2)) .or. (ref(1)==ref(3)) .or. (ref(2)==ref(3)) ) &
      &call error%MsgError( 1, message_error='Collinear system > Check index file: '//index_list%name)

    s1 = 0; s2 = 0; s3 = 0

    allocate(idlm2(n), stat=ierr)
    if ( ierr /=0 ) call error%MsgError( 1, ERROR_ALLOCATION )
    idlm2 = 0

    allocate(signs(n), stat=ierr)
    if ( ierr /=0 ) call error%MsgError( 1, ERROR_ALLOCATION )
    signs = 0

    allocate(draw(n), stat=ierr)
    if ( ierr /=0 ) call error%MsgError( 1, ERROR_ALLOCATION )
    draw = 0.0_DP

    me1%xyz(:)  = (pattern%atoms(ref(1))%xyz(:)+pattern%atoms(ref(2))%xyz(:))/2.0_DP
    me1%xyz(:) = (me1%xyz(:) + pattern%atoms(ref(3))%xyz(:))/2.0_DP

    ame1%xyz(:) = (pattern%atoms(ref(1))%xyz(:) - me1%xyz(:))
    bme1%xyz(:) = (pattern%atoms(ref(3))%xyz(:) -  me1%xyz(:))

    name1 = ame1%norm()
    nbme1 = bme1%norm()

    ame1%xyz(:) = ame1%xyz(:)/name1
    bme1%xyz(:) = bme1%xyz(:)/nbme1

    cfme1 = cross_product( ame1, bme1 )
    ncfme1 = cfme1%norm()

    cfme1%xyz(:) = cfme1%xyz(:)/ncfme1

    eps0 = 0.00005_DP
    itera = 0
        
    do i = 1, n
       VRR%xyz(:) = pattern%atoms(i)%xyz(:) - me1%xyz(:)
       nVRR = VRR%norm()

       dot2 = Scalar_product( cfme1, VRR )

       if ( dabs(dot2/(ncfme1*nVRR)) < eps0) then
          itera = itera + 1
          idlm2(itera) = i
       end if

    end do

    do i = 1, m
       VR%xyz(:) = (image%atoms(i)%xyz(:) - me1%xyz(:))
       nVR = VR%norm()

       VRR%xyz(:) = (pattern%atoms(i)%xyz(:) - me1%xyz(:))
       nVRR =  VRR%norm()

       dot2 = Scalar_product( cfme1,VRR )
       dot2 = dot2/( ncfme1*nVRR )

       call Prevent_NAN(dot2)

       dot1 = Scalar_product(cfme1,VR)
       dot1 = dot1/(ncfme1*nVR)

       call Prevent_NAN(dot1)

       dot1 = dacos(dot1)
       dot2 = dacos(dot2)

       diff = dot1-dot2

       if (dabs(diff) < eps0) then
          signs(i) = 0
       elseif( diff < eps0) then
             signs(i) = -1
       elseif( diff > eps0) then
             signs(i) = 1
       end if

    end do

    allocate(ver(mn,3), angle(mn), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )
    ver = 0

    do i = 1, mn
       read(index_list%unit,*,iostat=ios) ver(i,:)
       if ( ios /= 0 ) call error%MsgError( 1, message_error='Error reading > Check index file: '//index_list%name )
       
       if ( ver(i,1) > n .or. ver(i,2) > n .or. ver(i,3) > n) &
          &call error%MsgError( 1, message_error='Invalid index vertice > Check index file: '//index_list%name)
    end do 

    do i = 1, mn
      do j = 1, n

       if ( j == ver(i,1) ) then
          U%xyz(:) = pattern%atoms(j)%xyz(:)
          U2%xyz(:) = image%atoms(j)%xyz(:)
          s1 = signs(j)
       end if  
              
       if ( j == ver(i,2) ) then
          V%xyz(:) = pattern%atoms(j)%xyz(:)
          V2%xyz(:) = image%atoms(j)%xyz(:)
          s2 = signs(j)
       end if
              
       if ( j == ver(i,3) ) then
          W%xyz(:) = pattern%atoms(j)%xyz(:)
          W2%xyz(:) = image%atoms(j)%xyz(:)
          s3 = signs(j)
       end if
            
      end do 

      VU%xyz(:)  = V%xyz(:)  - U%xyz(:) 
      VU2%xyz(:) = V2%xyz(:) - U2%xyz(:) 
      WU%xyz(:)  = W%xyz(:)  - U%xyz(:)
      WU2%xyz(:) = W2%xyz(:) - U2%xyz(:)

      nVU  = VU%norm() 
      nVU2  = VU2%norm() 
      nWU  = WU%norm() 
      nWU2  = WU2%norm() 

      VU%xyz(:)  = VU%xyz(:)/nVU
      VU2%xyz(:) = VU2%xyz(:)/nVU2
      WU%xyz(:)  = WU%xyz(:)/nWU
      WU2%xyz(:) = WU2%xyz(:)/nWU2

      ONP = cross_product(VU,WU)
      nONP = ONP%norm() 
      
      ONP%xyz(:) = ONP%xyz(:)/nONP

      ONP2 = cross_product(VU2,WU2)
      nONP2 = ONP2%norm() 

      ONP2%xyz(:) = ONP2%xyz(:)/nONP2

      dot = Scalar_product(ONP,ONP2)

      call Prevent_NAN(dot)
 
      ss = s1 + s2 + s3

      if (ss >= 0) then
          ss = 1
      else
          ss = -1
      end if

      angle(i) = ((dacos(dot)*180.0_DP)/PI)

      write(69,*)'Triangle ', i, ((dacos(dot)*180.0_DP)/PI)*ss
           
    end do 

    write(*,'(a)')'Angles save in Angles_values.dat'
    write(*,'(a)')'--------------------------------------------------------------------------------'
    open(19, file='Angle_values.dat', status='unknown')

    do j = 1, n
       sum = 0.0_DP
       counter = 0
       
       do i = 1, mn
          if ( ver(i,1) == j .or. ver(i,2) == j .or. ver(i,3) == j ) then
             sum = sum + angle(i)
             counter = counter + 1
           end if

       end do

       draw(j) = (sum/counter) * signs(j)

       do i = 1, n
          if ( j == idlm2(i) ) then
             write(19,'(i2,2x,f10.5)') j, draw(j)
          end if
       end do

    end do

    deallocate(angle, ver, signs, draw)
    close(index_list%unit)
  end subroutine Calc_surface_torsion
    
  real(DP) function Calc_Hausdorff( pattern, image, diameter ) result(chi)
    type(molecule_t), intent(in) :: pattern
    type(molecule_t), intent(in) :: image
    real(DP), intent(in)       :: diameter
    integer :: n, m, k, i, j, ierr=0, ik
    real(DP), allocatable, dimension(:,:) :: A, B, E
    real(DP), allocatable, dimension(:,:) :: Br, Er
    real(DP) :: dx, dy, dz, p2(3), tmp(3), cmax, cmin, d2, xi
    real(DP), parameter :: thresh = 0.0_DP

    n = pattern%natm
    m = image%natm
    
    if ( n /= m ) call error%MsgError( 1, message_error='Size pattern/image dont match.' )

    allocate(A(3,n), B(3,n), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    allocate(E(3,n), Br(3,n), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )
 
    do k = 1, n
       A(1,k) = pattern%atoms(k)%xyz(1)
       A(2,k) = pattern%atoms(k)%xyz(2)
       A(3,k) = pattern%atoms(k)%xyz(3)
    end do

    do k = 1, n
       B(1,k) = image%atoms(k)%xyz(1)
       B(2,k) = image%atoms(k)%xyz(2)
       B(3,k) = image%atoms(k)%xyz(3)
    end do

!@====================================================    

    o1: do i = 1, n
        cmin = 999999.9
        i1: do j = 1, n
               d2 = (A(1,i)-B(1,j))**2 + &
                  &(A(2,i)-B(2,j))**2 + &
                  &(A(3,i)-B(3,j))**2

               if(d2 < cmin) then
                 cmin = d2     
               end if
             end do i1

             if(cmin > cmax) then
               cmax = cmin
             end if
    end do o1    

   chi = sqrt(cmax) / diameter



!@====================================================

!    k = 0

!    c1: do i = 1, n
!        d1: do j = 1, n
!               dx = dabs(A(1,i) - B(1,j))    
!               dy = dabs(A(2,i) - B(2,j))    
!               dz = dabs(A(3,i) - B(3,j))  
!               if ( dx < thresh .and. dy < thresh .and. dz < thresh ) then 
!                  cycle d1
!               end if 
!        end do d1
!        k = k + 1
!        E(:,k) = A(:,i)
!    end do c1
    
!    allocate(Er(3,k), stat=ierr)
!    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

!    do i = 1, k
!       Er(:,i) = E(:,i)
!    end do

!    do i = 1, k
!       xi = Get_Random()
!       ik = int(xi*k)+1

!       tmp(:) = Er(:,i)
!       p2(:)  = Er(:,ik)
!
!       Er(:,i) = p2(:)
!       Er(:,ik) = tmp(:)
!    end do

!    Br = B

!    do i = 1, n
!       xi = Get_Random()
!       ik = int( xi*n) + 1

!       tmp(:) = Br(:,i)
!       p2(:)  = Br(:,ik)

!       Br(:,i)  = p2(:)
!       Br(:,ik) = tmp(:)
!    end do

!    cmax = 0.0_DP

!    o1: do i = 1, k
!        cmin = 9999.0_DP

!        i1: do j = 1, n
!            d2 = (Er(1,i)-Br(1,j))*(Er(1,i)-Br(1,j))+&
!                 &(Er(2,i)-Br(2,j))*(Er(2,i)-Br(2,j))+&
!                 &(Er(3,i)-Br(3,j))*(Er(3,i)-Br(3,j))
!            d2 = dsqrt(d2)   
            
!            if ( d2 < cmin ) then
!               cmin = d2
!            end if

!       end do i1

!       if ( cmin > cmax ) then
!          cmax = cmin
!       end if
!   end do o1

!   chi = cmax/diameter
  end function Calc_Hausdorff

  real(DP) function Calc_Diameter( molecule ) result(phi)
    type(molecule_t), intent(in) :: molecule
    real(DP), allocatable, dimension(:,:) ::  r2
    real(DP) :: rx, ry, rz
    integer :: n, ierr=0, i, j

     n = molecule%natm

     allocate(r2(n,n), stat=ierr)
     if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

     r2 = 0.0_DP 

     do i = 1, n-1
        do j = i+1, n
           rx = molecule%atoms(i)%xyz(1) - molecule%atoms(j)%xyz(1)
           ry = molecule%atoms(i)%xyz(2) - molecule%atoms(j)%xyz(2)
           rz = molecule%atoms(i)%xyz(3) - molecule%atoms(j)%xyz(3)

           rx = rx*rx
           ry = ry*ry
           rz = rz*rz

           r2(i,j) = rx + ry + rz
        end do
     end do

     phi = maxval(r2)
     phi = dsqrt(phi)

     deallocate(r2, stat=ierr)
     if( ierr /= 0 ) call error%MsgError( 1, ERROR_DEALLOCATION )
  end function Calc_Diameter

  real(DP) function Get_Random() result(xi)
    integer :: n, ierr=0
    integer, allocatable, dimension(:) :: seed

    call random_seed(size=n)
    allocate(seed(n), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    call random_seed(get=seed)
    call random_number(xi)
  end function Get_Random

  real(DP) function Calc_RMSD( pattern, image ) result(rmsd)
    type(molecule_t), intent(in) :: pattern, image
    real(DP) :: x, y, z
    integer :: n, m, k

    n = pattern%natm
    m = image%natm

    if ( n /= m ) call error%MsgError( 1, message_error='Size pattern/image dont match.' )

    rmsd = 0.0_DP

    do k = 1, n
       x = pattern%atoms(k)%xyz(1) - image%atoms(k)%xyz(1)
       y = pattern%atoms(k)%xyz(2) - image%atoms(k)%xyz(2)
       z = pattern%atoms(k)%xyz(3) - image%atoms(k)%xyz(3)

       x = x**2
       y = y**2
       z = z**2

      rmsd = rmsd + ( x + y + z )
    end do

    rmsd = rmsd / n
    rmsd = dsqrt(rmsd)
  end function Calc_RMSD

  type(molecule_t) function Optimal_rotation( &
      &pattern, image, frame, is_traj, utraj, measure) result(rimage)
    type(molecule_t), intent(in) :: pattern
    type(molecule_t), intent(inout) :: image 
    integer, intent(in), optional :: is_traj, utraj, frame, measure
    real(DP), allocatable, dimension(:,:) ::&
     & X0, X1, X0T, X1T
    integer :: n, m, ierr=0, k, j, kk
    real(DP), dimension(3,3) ::&
     & V10, V01, V10V01, I, A, VV, V00, V11, V, VT, V10P, SumV
    real(DP), dimension(3)   :: c   
    real(DP) :: T, T0, T1, Tr, s, theta, nW
    real(DP), dimension(4,4) :: B

    real(DP) :: D02, D2

    real(DP) :: va(3,3),vb(4,4), da(3), db(4)

    type(eigen_t)  :: Aeig, Beig
    type(vector_t) :: W
    type(molecule_t)  :: backup

    real(DP) :: rmsd(4,181), mrmsd(4)
    real(DP) :: chi(4,181), mchi(4)
    real(DP) :: chia(0:180,4)

    real(DP) :: chi1, rmsd1
    real(DP) :: diameter, mrms
    integer  :: mlrmsd(4), mlchi(4), id, mid, iter

    real(DP), allocatable :: Ma(:,:), MaR(:,:)

    real(DP) :: start, finish

    real(DP) :: mineig, maxeig

    integer :: it_max, it_num, rot_num

    integer :: inL1, inL4

    call cpu_time(start)

    n = pattern%natm
    m = image%natm

    allocate(Ma(m,3), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )
    allocate(MaR(m,3), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1,ERROR_ALLOCATION )

    backup = image
    rmsd = 999.9_DP; chi = 999.9_DP

    close(66); close(77); close(88)
    open(66,file='Pattern.xyz',status='unknown')
    open(77,file='Rotate_image.xyz',status='unknown')
    open(88,file='RMSD_CHI.dat',status='unknown')

     if ( n /= m ) call error%MsgError( 1, message_error='Size pattern/image dont match.' )

     allocate(X0(n,3), X1(m,3), stat=ierr) 
     if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )
     allocate(X0T(3,n), X1T(3,m), stat=ierr)
     if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

     allocate(rimage%atoms(m), stat=ierr)
     if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

     X0  = 0.0_DP
     X1  = 0.0_DP
     X0T = 0.0_DP
     X1T = 0.0_DP

     do k = 1, n
        X0(k,1) = pattern%atoms(k)%xyz(1)
        X0(k,2) = pattern%atoms(k)%xyz(2)
        X0(k,3) = pattern%atoms(k)%xyz(3)
     end do

     do k = 1, m
        X1(k,1) = image%atoms(k)%xyz(1)
        X1(k,2) = image%atoms(k)%xyz(2)
        X1(k,3) = image%atoms(k)%xyz(3)
     end do

     X0T = transpose(X0)
     X1T = transpose(X1)

     V   = 0.0_DP; VT  = 0.0_DP; VV     = 0.0_DP
     V00 = 0.0_DP; V00 = 0.0_DP; V11    = 0.0_DP
     V10 = 0.0_DP; V01 = 0.0_DP; V10V01 = 0.0_DP

     if ( m < 20 ) then
        V00 = matmul(X0T,X0)
        V11 = matmul(X1T,X1)
     else   
        call Parallel_matmul(X0T,X0,V00, mm=3, nn=3, ll=m)
     end if   

     T0 = V00(1,1) + V00(2,2) + V00(3,3)
     T1 = V11(1,1) + V11(2,2) + V11(3,3)
     T  = (T0 + T1)/2.0_DP

     if ( m < 20 ) then
        V10 = matmul(X1T,X0)
        V01 = matmul(X0T,X1)

        V10P = matmul(X1T,X0)
     else
        call Parallel_matmul(X1T,X0,V10, mm=3, nn=3, ll=m)
        call Parallel_matmul(X0T,X1,V01, mm=3, nn=3, ll=m)
        call Parallel_matmul(X1T,X0,V10P, mm=3, nn=3, ll=m)
     end if   

     c(1) = V10(2,3) - V10(3,2)
     c(2) = V10(3,1) - V10(1,3)
     c(3) = V10(1,2) - V10(2,1)

     SumV = V10 + V01
     Tr = SumV(1,1)+SumV(2,2)+SumV(3,3)

     I = 0.0_DP
     I(1,1) = 1.0_DP; I(2,2) = 1.0_DP; I(3,3) = 1.0_DP

     A = SumV - Tr * I

     B(1,1) = 0.0_DP; B(1,2) = c(1); B(1,3) = c(2); B(1,4) = c(3)
     B(2,1) = c(1);   B(2,2) = A(1,1); B(2,3) = A(1,2); B(2,4) = A(1,3)
     B(3,1) = c(2);   B(3,2) = A(2,1); B(3,3) = A(2,2); B(3,4) = A(2,3)
     B(4,1) = c(3);   B(4,2) = A(3,1); B(4,3) = A(3,2); B(4,4) = A(3,3)
    
     if ( .not. is_trajectory .or. ( is_trajectory .and. is_verbose ) ) then
        write(*,'(a)')'---------------------------------------------------------------------------------'
        write(*,'(a)') 'Symmetrical matrix in Diamond`s method to optimal rotation: B'
        write(*,'(a)')'---------------------------------------------------------------------------------'
        do j = 1, 4
           write(*,'(4f10.3)') B(j,:)
        end do
     end if   

     !========================================
       it_max = 100

      call jacobi_eigenvalue ( 3, A, it_max, va, da, it_num, rot_num )
      call jacobi_eigenvalue ( 4, B, it_max, vb, db, it_num, rot_num )
     !========================================

       allocate(Beig%values(4),Beig%vectors(4,4))
       allocate(Aeig%values(3),Aeig%vectors(3,3))

       Beig%values(:) = db(:)
       Beig%vectors(:,:) =  vb(:,:)

     if ( .not. is_trajectory .or. ( is_trajectory .and. is_verbose ) ) then
        write(*,'(a)')'---------------------------------------------------------------------------------'
        write(*,'(a)')' Eigenvalues   -   Eigenvectors'
        write(*,'(a)')'---------------------------------------------------------------------------------'
        write(*,'(5f10.5)') Beig%values(1),Beig%vectors(:,1)
        write(*,'(5f10.5)') Beig%values(2),Beig%vectors(:,2)
        write(*,'(5f10.5)') Beig%values(3),Beig%vectors(:,3)
        write(*,'(5f10.5)') Beig%values(4),Beig%vectors(:,4)
        write(*,'(a)')'--------------------------------------------------------------------------------'
      end if

      mineig = Beig%values(1)
      maxeig = Beig%values(1)

      do j = 2, 4
         if(Beig%values(j) < mineig) then
            mineig = Beig%values(j)
            inL4 = j
         end if
         if(Beig%values(j) > maxeig) then
           maxeig = Beig%values(j)
           inL1 = j
         end if
      end do


      write(*,'(a)')'---------------------------------------------------------------------------------'
      
      ! use molecule inertia
  !    diameter = 4.0_DP * T / 3.0_DP
      diameter = Calc_Diameter( pattern )

     do j = 1, m
        Ma(j,:) = image%atoms(j)%xyz(:)
     end do

   !@  if ( .not. is_trajectory .or. ( is_trajectory .and. is_verbose ) ) then
   !@     write(*,'(a)')'Check files: Pattern.xyz / Rotate_image.xyz'
   !@  end if

     write(66,*) m
     write(66,'(a)')'Pattern: Build by KANON.' 
     do j = 1, m
        write(66,'(a,3f10.5)') pattern%atoms(j)%sym, pattern%atoms(j)%xyz(:)
     end do 

     write(88,'(a)') '#Trajectory RMSD  CHI' 

     iter = 0

     !$omp parallel do
     do kk = 1,4 ! loop over quaternion space

      do k = 0, 180
        theta = real(k,kind=DP) *PI / 180.0_DP
        
        W%xyz(1) =  Beig%vectors(2,kk)
        W%xyz(2) =  Beig%vectors(3,kk)
        W%xyz(3) =  Beig%vectors(4,kk)

        nW = W%norm()

        s = dcos(theta/2.0_DP)

        W%xyz(:) = W%xyz(:)/nW
        W%xyz(1) =  dsin(theta/2.0_DP) * W%xyz(1)
        W%xyz(2) =  dsin(theta/2.0_DP) * W%xyz(2)
        W%xyz(3) =  dsin(theta/2.0_DP) * W%xyz(3)

        call Rotate_matrix_quaternion( s, W, Ma, m, MaR )

        do j = 1, m
           image%atoms(j)%xyz(1) = MaR(j,1)
           image%atoms(j)%xyz(2) = MaR(j,2)
           image%atoms(j)%xyz(3) = MaR(j,3)
        end do
        
        chi1  = Calc_Hausdorff( pattern, image, diameter )
        iter = iter + 1

        chia(k,kk) = chi1

     end do
   end do

   !$omp end parallel do

    !==========================================

      
     if ( .not. is_trajectory .or. ( is_trajectory .and. is_verbose ) ) then
         write(*,'(a)')'---------------------------------------------------------------------------------'
         write(*,'(a)')'Hausdorff-derived chirality Index'
         write(*,'(a,f12.8,2x,f12.8)')'CHI:  ', minval(chia)
         write(*,'(a)')'--------------------------------------------------------------------------------'
     end if    

     ! TODO
    !@ if ( present(is_traj) .and. present(frame) .and. present(utraj) ) then
      !@  if ( is_traj == 1 ) then
       !@    if ( present(measure) .and. measure==0) then
       !@        write(utraj,'(i6,2f10.5)') frame, minval(mchi)
       !@    elseif( present(measure) .and. measure==1) then
       !@        write(utraj,'(i6,2f10.5)') frame,  minval(mrmsd)
       !@    end if    
       !@ end if
   !@  end if

   !@  mid = mlrmsd(1)
   !@  mrms = mrmsd(1)
   !@  id = 1
   !@  do k = 2, 4
   !@     if ( mrmsd(k) < mrms ) then
   !@        mrms = mrmsd(k)
   !@        mid = mlrmsd(k)
   !@        id = k
   !@     end if  
   !@  end do
     
   !@  image = backup
     
   !@  theta  = (real(mid-1,DP)*PI) / 180.0_DP
   !@  W%xyz(1) =  Beig%vectors(2,id)
   !@  W%xyz(2) =  Beig%vectors(3,id)
   !@  W%xyz(3) =  Beig%vectors(4,id)

   !@  nW = W%norm()

    !@ s = dcos(theta/2.0_DP)
   !@  if ( nw > 0.0_DP ) then
   !@     W%xyz(:) = W%xyz(:)/nW
   !@     W%xyz(1) =  dsin(theta/2.0_DP) * W%xyz(1)
   !@     W%xyz(2) =  dsin(theta/2.0_DP) * W%xyz(2)
   !@     W%xyz(3) =  dsin(theta/2.0_DP) * W%xyz(3)
       
    !@    call Rotate_matrix_quaternion( s, W, Ma, m, MaR )

     !@   rimage%natm = m

     !@   do j = 1, m
     !@      rimage%atoms(j)%xyz(1) = MaR(j,1)
     !@      rimage%atoms(j)%xyz(2) = MaR(j,2)
     !@      rimage%atoms(j)%xyz(3) = MaR(j,3)

     !!@      rimage%atoms(j)%sym = image%atoms(j)%sym
     !@   end do

   !@  else    
   !@     rimage = image
     !@end if

     call cpu_time(finish)

     if ( is_performance) call Get_timing(start, finish)
  end function Optimal_rotation   

  type(molecule_t) function Build_image( pattern ) result(image)
    type(molecule_t), intent(in) :: pattern
    real(DP), allocatable, dimension(:,:) ::&
     & I, P, X0, X0T, X1, X1T
    integer :: k, n, ierr=0
                    
    n = pattern%natm

    allocate(I(n,n), P(n,n), stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    I = 0.0_DP; P = 0.0_DP

    do k = 1, n
      I(k,k) = 1.0_DP
    end do  

    P = I

    allocate(X0(n,3), X0T(3,n), X1(n,3), X1T(3,n), stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    do k = 1, n
       X0(k,1) = pattern%atoms(k)%xyz(1)
       X0(k,2) = pattern%atoms(k)%xyz(2)
       X0(k,3) = pattern%atoms(k)%xyz(3)
    end do

    if ( n < 20 ) then
       X1 = matmul(-P,X0)
    else   
       call Parallel_matmul(-P,X0,X1, mm=n, nn=3, ll=n)
    end if   

    image%natm = n
    allocate(image%atoms(n), stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    do k = 1, n
       image%atoms(k)%sym = pattern%atoms(k)%sym
       image%atoms(k)%xyz(1) = X1(k,1)
       image%atoms(k)%xyz(2) = X1(k,2)
       image%atoms(k)%xyz(3) = X1(k,3)
    end do

    deallocate(X0, X0T, X1, X1T, I, P, stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_DEALLOCATION )
  end function Build_image  

  subroutine Rotate_molecule_to_axis( this, out )
    class(molecule_t), intent(inout) :: this
    type(vector_t) :: Z, Y, V, W, rP, F
    type(vector_t) :: e1, e2, e3, e11, e22, e33
    integer, intent(in) :: out
    real(DP) :: s
    integer  :: k
    type(eigen_t) :: eig

    call this%Getcom()
    call this%Inertia(eig)
    
    if ( out == 1 ) then
     if ( .not. is_trajectory .or. ( is_trajectory .and. is_verbose ) ) then
         write(*,'(a,3f10.5)') 'Center of masses: ', this%com(:)
         write(*,'(a)')'------------------------------------------------------------------------------'
         write(*,'(a,3(3x,f10.5))')'Moments of inertia: ', eig%values(:)
         write(*,'(a)')'Eigenvectors'
         do k = 1, 3
            write(*,'(3f10.3)') eig%vectors(k,:)
         end do
     end if
   end if


    Z%xyz(1) = 0.0_DP; Z%xyz(2) = 0.0_DP; Z%xyz(3) = 1.0_DP
    Y%xyz(1) = 0.0_DP; Y%xyz(2) = 1.0_DP; Y%xyz(3) = 0.0_DP

    V%xyz(1) = eig%vectors(1,3)
    V%xyz(2) = eig%vectors(2,3)
    V%xyz(3) = eig%vectors(3,3)

    call Define_rotation( Z, V, W, s )

    e3 = Rotate_by_quaternion( s, W, V )

    V%xyz(1) = eig%vectors(1,2)
    V%xyz(2) = eig%vectors(2,2)
    V%xyz(3) = eig%vectors(3,2)
    
    e2 = Rotate_by_quaternion( s, W, V )
        
    V%xyz(1) = eig%vectors(1,1)
    V%xyz(2) = eig%vectors(2,1)
    V%xyz(3) = eig%vectors(3,1)
        
    e1 = Rotate_by_quaternion( s, W, V )

    do k = 1, this%natm
       F%xyz(1) = this%atoms(k)%xyz(1)
       F%xyz(2) = this%atoms(k)%xyz(2)
       F%xyz(3) = this%atoms(k)%xyz(3)

       rP = Rotate_by_quaternion( s, W, F )
       this%atoms(k)%xyz(1) = rP%xyz(1)
       this%atoms(k)%xyz(2) = rP%xyz(2)
       this%atoms(k)%xyz(3) = rP%xyz(3)
    end do

    call Define_rotation( Y, e2, W, s )
    e22 = Rotate_by_quaternion( s, W, e2 )
    e33 = Rotate_by_quaternion( s, W, e3 )
    e11 = Rotate_by_quaternion( s, W, e1 )

    do k = 1, this%natm
       F%xyz(1) = this%atoms(k)%xyz(1)
       F%xyz(2) = this%atoms(k)%xyz(2)
       F%xyz(3) = this%atoms(k)%xyz(3)

       rP = Rotate_by_quaternion(s,W,F)
       this%atoms(k)%xyz(1) = rP%xyz(1)
       this%atoms(k)%xyz(2) = rP%xyz(2)
       this%atoms(k)%xyz(3) = rP%xyz(3)
    end do
  end subroutine Rotate_molecule_to_axis

  subroutine Get_inertia_moments( this, eaxis )
    class(molecule_t), intent(in) :: this 
    real(DP) :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz
    real(DP) :: tensor(3,3)
    type(eigen_t), intent(inout) :: eaxis

    Ixx =  sum(this%atoms(:)%mass*(this%atoms(:)%xyz(2)**2+this%atoms(:)%xyz(3)**2))
    Iyy =  sum(this%atoms(:)%mass*(this%atoms(:)%xyz(1)**2+this%atoms(:)%xyz(3)**2))
    Izz =  sum(this%atoms(:)%mass*(this%atoms(:)%xyz(1)**2+this%atoms(:)%xyz(2)**2))
    Ixy = -sum(this%atoms(:)%mass*(this%atoms(:)%xyz(1)*this%atoms(:)%xyz(2)))
    Ixz = -sum(this%atoms(:)%mass*(this%atoms(:)%xyz(1)*this%atoms(:)%xyz(3)))
    Iyz = -sum(this%atoms(:)%mass*(this%atoms(:)%xyz(2)*this%atoms(:)%xyz(3)))

    tensor(1,1) = Ixx; tensor(1,2) = Ixy; tensor(1,3) = Ixz
    tensor(2,1) = Ixy; tensor(2,2) = Iyy; tensor(2,3) = Iyz
    tensor(3,1) = Ixz; tensor(3,2) = Iyz; tensor(3,3) = Izz

    if (is_verbose) then
       write(*,'(a)')'-----------------------------------------------------------------------------'
       write(*,'(a)')'Inertia tensor'
       write(*,'(3f10.5)') tensor(1,:)
       write(*,'(3f10.5)') tensor(2,:)
       write(*,'(3f10.5)') tensor(3,:)
       write(*,'(a)')'-----------------------------------------------------------------------------'
    end if
    
    eaxis = Get_eigenvalues(tensor,3)
  end subroutine Get_inertia_moments

  subroutine Get_center_of_masses( this )
    class(molecule_t), intent(inout) :: this
    real(DP) :: mass_mol

    mass_mol = sum(this%atoms(:)%mass)

    this%com(1) = sum(this%atoms(:)%xyz(1) * this%atoms(:)%mass)
    this%com(2) = sum(this%atoms(:)%xyz(2) * this%atoms(:)%mass)
    this%com(3) = sum(this%atoms(:)%xyz(3) * this%atoms(:)%mass)

    this%com(:) = this%com(:) / mass_mol

    this%atoms(:)%xyz(1) = this%atoms(:)%xyz(1) - this%com(1)
    this%atoms(:)%xyz(2) = this%atoms(:)%xyz(2) - this%com(2)
    this%atoms(:)%xyz(3) = this%atoms(:)%xyz(3) - this%com(3)
  end subroutine Get_center_of_masses

  subroutine Read_molecule_xyz( this, xyz, out )
    class(molecule_t), intent(inout) :: this
    type(file_t), intent(inout)      :: xyz
    integer, intent(in)              :: out
    integer                          :: line, k

    line = 0

    call xyz%inquire()
    rewind(xyz%unit)
    
    do
      read(xyz%unit,'(a)',iostat=ios) buffer
      if ( ios /= 0 ) exit
      line = line + 1
    end do

    if (line < 1) call error%MsgError(  1, message_error='Missing data in:'//xyz%name )
       
    rewind(xyz%unit)

    read(xyz%unit,*,iostat=ios) this%natm
    if( ios/=0 ) call error%MsgError( 1, ERROR_READING )
    if (line < this%natm+2) call error%MsgError( 1, message_error='Missing data in:'//xyz%name )
    
    read(xyz%unit,'(a)')
    call this%Malloc( this%natm )

    do k = 1, this%natm
       read(xyz%unit,*,iostat=ios) this%atoms(k)%sym,this%atoms(k)%xyz(:)
       if( ios/=0 ) call error%MsgError( 1, ERROR_READING )
          this%atoms(k)%mass = Add_mass_atom( this%atoms(k)%sym )  
    end do
      
    if ( out == 1 ) then

       if( is_bohr ) then
         write(*,'(a)')'---Converting bohr to angstroms---'
         write(*,'(a)')''
         do k = 1, this%natm
            this%atoms(k)%xyz(:) = BOHR2ANG * this%atoms(k)%xyz(:)
         end do
       end if

       write(*,'(a)') 'Cartesian coordinates input file: '//trim(xyz%name)
       write(*,'(a,i5)')'Number atoms: ', this%natm 
       do k = 1, this%natm
       write(*,'(a3,3f10.5)') this%atoms(k)%sym, this%atoms(k)%xyz(:)
       end do
       write(*,'(a)')''
    end if

  end subroutine Read_molecule_xyz

  subroutine Read_trajectory_xyz( trajectory, file_traj, frames_on )
    type(file_t), intent(inout) :: file_traj
    type(molecule_t), allocatable, intent(inout) :: trajectory(:)
    integer, intent(out) :: frames_on
    integer :: j, k, line
 
    call file_traj%inquire()
    rewind(file_traj%unit)

    line = 0
    do 
       read(file_traj%unit,*,iostat=ios) buffer
       if ( ios /= 0 ) exit
       line = line+1
    end do

    if ( line < 1 ) call error%MsgError( 1, &
      &message_error='Missing data in '//trim(file_trajectory%name)//'. File empty.')

    rewind(file_traj%unit)

    read(file_traj%unit,*,iostat=ios) natm_on
    
    if ( ios /= 0 ) call error%MsgError( 1,&
      & ERROR_READING)

    write(*,*) line, natm_on, natm_on+2

    if ( mod(line,natm_on+2) /= 0 ) call error%MsgError( 1, &
      &message_error='Missing data in '//trim(file_trajectory%name))

    frames_on = (line/(natm_on+2))

    if ( .not. is_trajectory .or. ( is_trajectory .and. is_verbose ) ) then
        write(*,'(a,i6)') 'Number of frames on trajectory ',frames_on
    end if    

    allocate(trajectory(frames_on), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1,&
      & ERROR_ALLOCATION) 
   
    do k = 1, frames_on
       allocate(trajectory(k)%atoms(natm_on), stat=ierr)
       if ( ierr /= 0 ) call error%MsgError( 1,&
          & ERROR_ALLOCATION) 
    end do
 
    rewind(file_traj%unit)

    do j = 1, frames_on
       read(file_traj%unit,*) 
       read(file_traj%unit,*) 

       trajectory(j)%natm = natm_on 

       do k = 1, natm_on
          read(file_traj%unit,*) trajectory(j)%atoms(k)%sym, trajectory(j)%atoms(k)%xyz(:)
          trajectory(j)%atoms(k)%mass = Add_mass_atom( trajectory(j)%atoms(k)%sym )  
       end do
    end do
  end subroutine Read_trajectory_xyz

  subroutine Allocate_molecule( this, n ) 
    class(molecule_t), intent(inout) :: this  
    integer, intent(in) :: n
    integer :: ierr=0

    if ( allocated(this%atoms) ) deallocate(this%atoms)

    allocate(this%atoms(n), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION)

  end subroutine Allocate_molecule

  real(DP) function Add_mass_atom( symbol ) result(mass)
     character(*), intent(in) :: symbol

     select case(symbol)
     case('H');  mass = MASS_HYDROGEN
     case('He'); mass = MASS_HELIUM
     case('Li'); mass = MASS_LITHIUM
     case('Be'); mass = MASS_BERYLLIUM
     case('B');  mass = MASS_BORON
     case('C');  mass = MASS_CARBON
     case('N');  mass = MASS_NITROGEN
     case('O');  mass = MASS_OXYGEN
     case('F');  mass = MASS_FLUORINE
     case('Ne'); mass = MASS_NEON
     case('Na'); mass = MASS_SODIUM
     case('Mg'); mass = MASS_MAGNESIUM
     case('Al'); mass = MASS_ALUMINIUM
     case('Si'); mass = MASS_SILICON
     case('P');  mass = MASS_PHOSPHORUS
     case('S');  mass = MASS_SULFUR
     case('Cl'); mass = MASS_CHLORINE
     case('Ar'); mass = MASS_ARGON
     case('K');  mass = MASS_POTASSIUM
     case('Ca'); mass = MASS_CALCIUM
     case('Sc'); mass = MASS_SCANDIUM
     case('Ti'); mass = MASS_TITANIUM
     case('V');  mass = MASS_VANADIUM
     case('Cr'); mass = MASS_CHROMIUM
     case('Mn'); mass = MASS_MANGANESE
     case('Fe'); mass = MASS_IRON
     case('Co'); mass = MASS_COBALT
     case('Ni'); mass = MASS_NICKEL
     case('Cu'); mass = MASS_COPPER
     case('Zn'); mass = MASS_ZINC
     case('Ga'); mass = MASS_GALLIUM
     case('Ge'); mass = MASS_GERMANIUM
     case('As'); mass = MASS_ARSENIC
     case('Se'); mass = MASS_SELENIUM
     case('Br'); mass = MASS_BROMINE
     case('Kr'); mass = MASS_KRYPTON
     case('Rb'); mass = MASS_RUBIDIUM
     case('Sr'); mass = MASS_STRONTIUM
     case('Y');  mass = MASS_YTTRIUM
     case('Zr'); mass = MASS_ZIRCONIUM
     case('Nb'); mass = MASS_NIOBIUM
     case('Mo'); mass = MASS_MOLYBDENUM
     case('Ru'); mass = MASS_RUTHENIUM
     case('Rh'); mass = MASS_RHODIUM
     case('Pd'); mass = MASS_PALLADIUM
     case('Ag'); mass = MASS_SILVER
     case('Cd'); mass = MASS_CADMIUM
     case('In'); mass = MASS_INDIUM
     case('Sn'); mass = MASS_TIN
     case('Sb'); mass = MASS_ANTIMONY
     case('Te'); mass = MASS_TELLURIUM
     case('I');  mass = MASS_IODINE
     case('Xe'); mass = MASS_XEON
     case('Cs'); mass = MASS_CESIUM
     case('Ba'); mass = MASS_BARIUM
     case('La'); mass = MASS_LANTHANUM
     case('Hf'); mass = MASS_HAFNIUM
     case('Ta'); mass = MASS_TANTALUM
     case('W');  mass = MASS_TUNGSTEN
     case('Re'); mass = MASS_RHENIUM
     case('Os'); mass = MASS_OSMIUM
     case('Ir'); mass = MASS_IRIDIUM
     case('Pt'); mass = MASS_PLATINUM
     case('Au'); mass = MASS_GOLD
     case('Hg'); mass = MASS_MERCURY
     case('Tl'); mass = MASS_THALLIUM
     case('Pb'); mass = MASS_LEAD
     case('Bi'); mass = MASS_BISMUTH
     case('X');  mass = MASS_DUMMY
     case default
       call error%MsgError( 1, message_error='Unknown symbol. Mass will be set to 1.0' )  
    end select
  end function Add_mass_atom  

  subroutine Get_point_group( molecule, show )
    class(molecule_t), intent(inout) :: molecule
	 integer, intent(in) :: show
    type( eigen_t) :: eigen
    real(kind=DP)  :: a, b, c, threshold

    character(len=2), parameter ::&
      & IDENTITY  = 'E',&
      & PLANE     = 's',&
      & INVERSION = 'i',&
      & PROPER    = 'Cn',&
      & IMPROPER  = 'Sn'

    logical ::&
     & is_spherical =.false.,&
     & is_oblate    = .false.,&
     & is_prolate   = .false.,&
     & is_assymetric= .false.,&
     & is_linear    = .false.

    threshold = 0.05_DP

    call molecule%Getcom()  
    write(*,'(a)')'--------------------------------------------------------------------------------'
    write(*,'(a)')'Center of masses' 
    write(*,'(3f10.5)') molecule%com(:)
    write(*,'(a)')'--------------------------------------------------------------------------------'

    call molecule%Inertia( eigen )

	 if ( show == 1 ) then
       write(*,'(a)')'--------------------------------------------------------------------------------'
       write(*,'(a)')'Principal moments of inertia'
       write(*,'(a)')'     Ia           Ib           Ic'
       write(*,'(f10.5,3x,f10.5,3x,f10.5)') eigen%values(3), eigen%values(2), eigen%values(1)
       write(*,'(a)')'--------------------------------------------------------------------------------'
       write(*,'(a)')'Inertia axis'
       write(*,'(3f10.5)') eigen%vectors(1,:) 
       write(*,'(3f10.5)') eigen%vectors(2,:) 
       write(*,'(3f10.5)') eigen%vectors(3,:) 
	 end if	 


    write(*,'(a)')'--------------------------------------------------------------------------------'
    write(*,'(a,f10.5)')'Symmetry elements. Tolerance: ',threshold
    write(*,'(a)')'--------------------------------------------------------------------------------'
    
    c = eigen%values(1)
    b = eigen%values(2)
    a = eigen%values(3)

    if ( dabs(b-c) < threshold .and. dabs(a) < threshold ) then
       write(*,'(a)') 'Linear molecule.'
       is_linear = .true.
    elseif ( dabs(a-b) < threshold .and. dabs(b-c) < threshold ) then
       write(*,'(a)') 'Spherical molecule.'
       is_spherical = .true.
    elseif( dabs(b-c) < threshold .and.&
     & ( dabs(a-b) > threshold .and. dabs(a-c) > threshold ) ) then   
       write(*,'(a)') 'Prolate symmetric top' 
       is_prolate = .true.
    elseif( dabs(a-b) < threshold .and.&
     & ( dabs(a-c) > threshold .and. dabs(b-c) > threshold ) ) then   
       write(*,'(a)') 'Oblate symmetric top'
      is_oblate = .true.
    elseif( dabs(a-b) > threshold .and.&
      &  dabs(a-c) > threshold .and. dabs(b-c) > threshold ) then  
       write(*,'(a)') 'Assymetric top'
       is_assymetric = .true.
    end if   
    write(*,'(a)')'--------------------------------------------------------------------------------'
  end subroutine Get_point_group

  type(molecule_t) function Check_symmetry_operation( pattern, symm ) result(image)
    type(molecule_t), intent(in) :: pattern
    character(*), intent(in), optional  :: symm
    integer :: n, k, ierr=0
    real(kind=DP), allocatable :: I(:,:), P(:,:)
    real(kind=DP), allocatable :: X0(:,:), X1(:,:)

    n = pattern%natm

    allocate(I(n,n), stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    allocate(P(n,n), stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    allocate(X0(n,3), stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    allocate(X1(n,3), stat=ierr)
    if( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    I = 0.0_DP; P = 0.0_DP

    do k = 1, n
      I(k,k) = 1.0_DP
    end do 

    do k = 1, n
       X0(k,1) = pattern%atoms(k)%xyz(1)
       X0(k,2) = pattern%atoms(k)%xyz(2)
       X0(k,3) = pattern%atoms(k)%xyz(3)
    end do

    P = I

    select case(symm)
    case('E')
      X1 = X0
    case('s')
    case('i')
      X1 = matmul(-P,X0)
    case('Cn')
    case('Sn')
    end select

    allocate(image%atoms(n), stat=ierr)
    if ( ierr /= 0 ) call error%MsgError( 1, ERROR_ALLOCATION )

    image%natm = n


    do k = 1, n
       image%atoms(k)%sym = pattern%atoms(k)%sym
       image%atoms(k)%xyz(1) = X1(k,1)
       image%atoms(k)%xyz(2) = X1(k,2)
       image%atoms(k)%xyz(3) = X1(k,3)
    end do

    deallocate(X0, X1, P, I, stat=ierr )
    if (ierr/=0) call error%MsgError( 1, ERROR_DEALLOCATION )

  end function Check_symmetry_operation

end module kanon_molecule_mod
