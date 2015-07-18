! ----------------------------------------------------------------
! QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  
!
! QUantum Mechanical VIbratioal Anlalysis or QUMVIA, for short, is
! an original implementation of Vibrational Self-consistent Field
! (VSCF) and Vibrational Configuration Interaction (VCI) using 
! distributed gaussian basis set with analytical integrals and a
! configuration selection algorithm.
!
! Implementation by Diego Javier Alonso de Armino
! Copyright(C) 2015 Diego J. Alonso de Armiño. All rights reserved
!
! This file is part of QUMVIA.
!
! QUMVIA is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! QUMVIA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! A copy of the licence can be found in the root directory of
! QUMVIA.  If not, see <http://www.gnu.org/licenses/>.
!
! QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  
! ----------------------------------------------------------------


!     PARAMETERS
!     -----------------------------------------------------------------------------------
!     The following array contains all atomic masses ordered by atomic number, 
!     so that atomic_masses(<atomic_number>) -where <atomic_number> may be 1 to 18-
!     returns the atomic mass corresponding to <atomic_number>. 
!     All masses correspond to the most abundant isotope and are in Daltons. Atomic 
!     masses for H, He, and O are obtanied from  Mohr, Taylor, Newell, Rev. Mod. Phys. 
!     80 (2008) 633-730. The rest are extracted from NIST: 
!     http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl. 
!     -----------------------------------------------------------------------------------
      real*8, PARAMETER :: EMASS_CODATA08 = 0.0005485799111D0 ! Electron mass in Daltons
      real*8, PARAMETER :: AMU_TO_AU = 1.0d0/0.0005485799111D0    ! Conversion factor Da to AU
      real*8, PARAMETER :: atomic_masses(1:18) =&
          (/  1.00866491574D0, 4.002603254153D0,  7.016004D0,  9.012182D0, 11.009305D0, &
             12.0D0,          14.0030740048D0,   15.99491461956D0, 18.998403D0,         &
             19.992440D0,     22.989770D0,       23.985042D0,      26.981538D0,         &
             27.976927D0,     30.973762D0,       31.972071D0,      34.968853D0,         &
             39.962383D0 /)
      real*8, PARAMETER :: atomic_masses_au(1:18) = atomic_masses * AMU_TO_AU
      real*8, PARAMETER :: invsqrt_atomic_masses_au(1:18) = 1.0d0/SQRT(atomic_masses_au)
      real*8, PARAMETER :: sqrt_atomic_masses_au(1:18) = SQRT(atomic_masses_au)

