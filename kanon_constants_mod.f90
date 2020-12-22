!-----------------------------------------------------------------------------------
!KANON: A program for chirality calculation and other elements of molecular symmetry
!-----------------------------------------------------------------------------------
!>@file   kanon_constants_mod
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
module kanon_constants_mod
  use kanon_kinds_mod, only:&
      & DP
  
  implicit none
  private

  !>@note Standard atomic weights <https://iupac.org>  
  real(kind=DP), public, parameter ::  &
    &MASS_HYDROGEN   = 1.00794_DP,   &
    &MASS_HELIUM     = 4.00260_DP,   &
    &MASS_LITHIUM    = 6.941_DP,     & 
    &MASS_BERYLLIUM  = 9.012182_DP,  &
    &MASS_BORON      = 10.811_DP,    &
    &MASS_CARBON     = 12.0107_DP,   &
    &MASS_NITROGEN   = 14.0067_DP,   &
    &MASS_OXYGEN     = 15.9994_DP,   &
    &MASS_FLUORINE   = 18.9984032_DP,& 
    &MASS_NEON       = 20.1797_DP,   &
    &MASS_SODIUM     = 22.989770_DP, & 
    &MASS_MAGNESIUM  = 24.305_DP,    &
    &MASS_ALUMINIUM  = 26.981538_DP, &
    &MASS_SILICON    = 28.0855_DP,   &
    &MASS_PHOSPHORUS = 30.973761_DP, &
    &MASS_SULFUR     = 32.065_DP,    &
    &MASS_CHLORINE   = 35.453_DP,    &
    &MASS_ARGON      = 39.948_DP,    &
    &MASS_POTASSIUM  = 39.0983_DP,   &
    &MASS_CALCIUM    = 40.078_DP,    &
    &MASS_SCANDIUM   = 44.955910_DP, &
    &MASS_TITANIUM   = 47.867_DP,    &
    &MASS_VANADIUM   = 50.9415_DP,   &
    &MASS_CHROMIUM   = 51.9961_DP,   &
    &MASS_MANGANESE  = 54.938044_DP, &
    &MASS_IRON       = 55.845_DP,    &
    &MASS_COBALT     = 58.933194_DP, &
    &MASS_NICKEL     = 58.6934_DP,   &
    &MASS_COPPER     = 63.546_DP,    &
    &MASS_ZINC       = 65.38_DP,     &
    &MASS_GALLIUM    = 69.723_DP,    &
    &MASS_GERMANIUM  = 72.630_DP,    &
    &MASS_ARSENIC    = 74.921595_DP, &
    &MASS_SELENIUM   = 78.971_DP,    &
    &MASS_BROMINE    = 79.901_DP,    &
    &MASS_KRYPTON    = 83.798_DP,    &
    &MASS_RUBIDIUM   = 85.4678_DP,   &
    &MASS_STRONTIUM  = 87.62_DP,     &
    &MASS_YTTRIUM    = 88.90584_DP,  &
    &MASS_ZIRCONIUM  = 91.224_DP,    &
    &MASS_NIOBIUM    = 92.90637_DP,  &
    &MASS_MOLYBDENUM = 95.95_DP,     &
    &MASS_RUTHENIUM  = 101.07_DP,    &
    &MASS_RHODIUM    = 102.90550_DP, &
    &MASS_PALLADIUM  = 106.42_DP,    &
    &MASS_SILVER     = 107.8682_DP,  &
    &MASS_CADMIUM    = 112.414_DP,   &
    &MASS_INDIUM     = 114.818_DP,   &
    &MASS_TIN        = 118.710_DP,   &
    &MASS_ANTIMONY   = 121.760_DP,   &
    &MASS_TELLURIUM  = 127.60_DP,    &
    &MASS_IODINE     = 126.90447_DP, & 
    &MASS_XEON       = 131.293_DP,   &
    &MASS_CESIUM     = 132.90545196_DP,&
    &MASS_BARIUM     = 137.327_DP,   &
    &MASS_LANTHANUM  = 138.90547_DP, &
    &MASS_CERIUM     = 140.116_DP,   &
    &MASS_PRASEODYMIUM = 140.90766_DP, &
    &MASS_NEODYMIUM  = 144.242_DP,   &
    &MASS_SAMARIUM   = 150.36_DP,    &
    &MASS_EUROPIUM   = 151.964_DP,   &
    &MASS_GADOLIMIUM = 157.25_DP,    &
    &MASS_TERBIUM    = 158.92535_DP, &
    &MASS_DYSPROSIUM = 162.500_DP,   &
    &MASS_HOLMIUM    = 164.93033_DP, &
    &MASS_ERBIUM     = 167.259_DP,   &
    &MASS_THULIUM    = 168.93422_DP, &
    &MASS_YTTERBIUM  = 173.054_DP,   &
    &MASS_LUTETIUM   = 174.9668_DP,  &
    &MASS_HAFNIUM    = 178.49_DP,    &
    &MASS_TANTALUM   = 180.94788_DP, &
    &MASS_TUNGSTEN   = 183.84_DP,    &
    &MASS_RHENIUM    = 186.207_DP,   &
    &MASS_OSMIUM     = 190.23_DP,    &
    &MASS_IRIDIUM    = 192.217_DP,   &
    &MASS_PLATINUM   = 195.084_DP,   &
    &MASS_GOLD       = 196.96655_DP, &
    &MASS_MERCURY    = 200.592_DP,   &
    &MASS_THALLIUM   = 204.382_DP,   &
    &MASS_LEAD       = 207.2_DP,     &
    &MASS_BISMUTH    = 208.98040_DP, &
    &MASS_DUMMY      =   1.00000_DP

  real(DP), parameter, public ::&
    & PI       = 3.1415926535897932_DP,&
    & BOHR2ANG = 0.52917721067_DP
end module kanon_constants_mod
