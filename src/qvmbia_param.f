
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
      real*8, PARAMETER :: D1=1.0D00, D2=2.0D00

!     Gauss-Hermite quadrature points for N=16.
      real*8,parameter,dimension(16) ::  ghx16=(/-4.688738939305818d+0,&
                                    -3.869447904860123d+0,&
                                    -3.176999161979960d+0,&
                                    -2.546202157847480d+0,&
                                    -1.951787990916254d+0,&
                                    -1.380258539198882d+0,&
                                   -0.8229514491446558d+0,&
                                   -0.2734810461381524d+0,&
                                    0.2734810461381525d+0,&
                                    0.8229514491446560d+0,&
                                     1.380258539198880d+0,&
                                     1.951787990916254d+0,&
                                     2.546202157847482d+0,&
                                     3.176999161979959d+0,&
                                     3.869447904860125d+0,&
                                     4.688738939305819d+0 /)

!     Factor for converting frequencies from atomic units to cm-1
      real*8  :: e, a0, c, e0, xl, pi, fact1, wave
!     -------------------------------------------------------------------------
!
!     Set up some constants
!     First, WAVE: factor for converting frequencies from AU to cm-1
!
      E   =   1.6021892D00         ! unit charge
      A0  =   0.5291771D00         ! bohr radius
      C   =   2.9979245D00         ! speed of light
      E0  =   8.85418782D00        ! permittivity of vacuum
      XL  =   6.022045D00          ! Avagadro's number
      PI  =   D2*ACOS(0.0D00)      ! PI
!
      FACT1 = E*E*1.0D02/(D2*D2*PI*E0*A0)  ! HARTREE to ERG
      FACT1 = FACT1*XL/(A0*A0)             ! hartree/bohr**2/amu to sec**-2
      FACT1 = SQRT(FACT1)                  ! to sec**-1
      FACT1 = FACT1/(D2*PI*C)              ! to wavenumbers
      WAVE  = FACT1*1.0D04                 ! include all powers of ten
