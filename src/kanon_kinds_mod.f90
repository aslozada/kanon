!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_kinds_mod
!>@brief  Defines portable type variables:
!>        SP Single Precision
!>        DP Double Precision
!>@author Asdrubal Lozada
!>        Laboratory of Theoretical Chemistry - LQT
!>        <http://www.lqt.dq.ufscar.br>
!>        Federal University of São Carlos
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
module kanon_kinds_mod
    
  implicit none
  private

  integer, parameter, public ::       &
    & SP = selected_real_kind(6,37),  &
    & DP = selected_real_kind(15,300),&
    & I8 = selected_int_kind(15)
end module kanon_kinds_mod
