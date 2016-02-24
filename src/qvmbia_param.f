
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

 
!     Isotopic atomic masses.
      
!     Atomic masses (u.m.a.) of most common isotopes
      real*8,parameter,dimension(216) :: isomass = (/ &
!      H-1             H-2             H-3
      &  1.007825037D0,  2.014101787D0,  3.016049286D0,  0d0, &
!      He-4            He-3
      &  4.00260325D0,   3.016029297D0,  0d0,            0d0, &
!      Li-7            Li-6
      &  7.0160045D0,    6.0151232D0,    0d0,            0d0, &
!      Be-9
      &  9.0121825D0,    0d0,            0d0,            0d0, &
!       B-11            B-10
      & 11.0093053D0,   10.0129380d0,    0d0,            0d0, &
!       C-12            C-13
      & 12.000000000d0, 13.003354839D0,  0d0,            0d0, &
!       N-14            N-15
      & 14.003074008D0, 15.000108978D0,  0d0,            0d0, &
!       O-16            O-18            O-17
      & 15.99491464D0,  17.99915939D0,  16.9991306D0,    0d0, &
!       F-19
      & 18.99840325D0,   0d0,            0d0,            0d0, &
!      Ne-20           Ne-22           Ne-21
      & 19.9924391D0,   21.9913837D0,   20.9938453D0,    0d0, &
!      Na-23
      & 22.9897697D0,    0d0,            0d0,            0d0, &
!      Mg-24           Mg-26           Mg-25
      & 23.9850450d0,   25.9825954D0,   24.9858392D0,    0d0, &
!      AL-27
      & 26.9815413D0,    0d0,            0d0,            0d0, &
!      Si-28           Si-29           Si-30
      & 27.9769284D0,   28.9764964D0,   29.9737717D0,    0d0, &
!       P-31
      & 30.9737634D0,    0d0,            0d0,            0d0, &
!       S-32            S-34            S-33            S-36
      & 31.9720718D0,   33.96786774D0,  32.9714591D0,   35.9670790d0, &
!      Cl-35           Cl-37
      & 34.968852729D0, 36.965902624D0,  0d0,            0d0, &
!      Ar-40           Ar-36           Ar-38
      & 39.9623831D0,   35.967545605D0, 37.9627322D0,    0d0, &
!       K-39            K-41            K-40
      & 38.9637079D0,   40.9618254D0,   39.9639988D0,    0d0, &
!      Ca-40           Ca-44           Ca-42           Ca-48
      & 39.9625907D0,   43.9554848D0,   41.9586218D0,   47.952532D0, &
!      Sc-45
      & 44.9559136D0,    0d0,            0d0,            0d0, &
!      Ti-48           Ti-46           Ti-47           Ti-49
      & 47.9479467D0,   45.9526327D0,   46.9517649D0,   48.9478705D0, &
!       V-51            V-50
      & 50.9439625D0,   49.9471613D0,    0d0,            0d0, &
!      Cr-52           Cr-53           Cr-50           Cr-54
      & 51.9405097D0,   52.9406510d0,   49.9460463D0,   53.9388822D0, &
!      Mn-55
      & 54.9380463D0,    0d0,            0d0,            0d0, &
!      Fe-56           Fe-54           Fe-57           Fe-58
      & 55.9349393D0,   53.9396121D0,   56.9353957D0,   57.9332778D0, &
!      Co-59
      & 58.9331978D0,    0d0,            0d0,            0d0, &
!      Ni-58           Ni-60           Ni-62           Ni-61
      & 57.9353471D0,   59.9307890d0,   61.9283464D0,   60.9310586D0, &
!      Cu-63           Cu-65
      & 62.9295992D0,   64.9277924D0,    0d0,            0d0, &
!      Zn-64           Zn-66           Zn-68           Zn-67
      & 63.9291454D0,   65.9260352D0,   67.9248458D0,   66.9271289D0, &
!      Ga-69           Ga-71
      & 68.9255809D0,   70.9247006D0,    0d0,            0d0, &
!      Ge-74           Ge-72           Ge-70           Ge-73
      & 73.9211788D0,   71.9220800d0,   69.9242498D0,   72.9234639D0, &
!      As-75
      & 74.9215955D0,    0d0,            0d0,            0d0, &
!      Se-80           Se-78           Se-82           Se-76
      & 79.9165205D0,   77.9173040d0,   81.916709D0,    75.9192066D0, &
!      Br-79           Br-81
      & 78.9183361D0,   80.916290d0,     0d0,            0d0, &
!      Kr-84           Kr-86           Kr-82           Kr-83
      & 83.9115064D0,   85.910614D0,    81.913483D0,    82.914134D0, &
!      Rb-85
      & 84.9117D0,      0d0,             0d0,             0d0, &
!      Sr-88           Sr-84           Sr-86           Sr-87
      & 87.9056D0,      83.9134d0,      85.9094d0,      86.9089d0, &
!      Y-89
      & 88.9054D0,      0d0,             0d0,             0d0, &
!      Zr-90           Zr-91           Zr-92           Zr-94
      & 89.9043D0,      90.9053D0,      91.9046D0,      93.9061D0, &
!      Nb-93
      & 92.9060d0,      0d0,             0d0,             0d0, &
!      Mo-98           Mo-92           Mo-95           Mo-96
      & 97.9055D0,      91.9063D0,      94.90584D0,     95.9046D0, &
!      Tc
      & 98.0d0,         0d0,             0d0,             0d0, &
!      Ru-102          Ru-99           Ru-100          Ru-104
      & 101.9037D0,     98.9061D0,      99.9030d0,      103.9055D0, &
!      Rh-103
      & 102.9048D0,     0d0,             0d0,             0d0, &
!      Pd-106          Pd-104           Pd-105         Pd-108
      & 105.9032D0,     103.9036D0,      104.9046D0,    107.90389D0, &
!      Ag-107          Ag-109
      & 106.90509d0,    108.9047D0,      0d0,             0d0, &
!      Cd-114          Cd-110           Cd-111         Cd-112
      & 113.9036D0,     109.9030d0,      110.9042D0,    111.9028D0, &
!      In-115          In-113
      & 114.9041D0,     112.9043D0,      0d0,             0d0, &
!      Sn-118          Sn-116           Sn-117         Sn-119
      & 117.9018D0,     115.9021D0,      116.9031D0,    118.9034D0, &
!      Sb-121          Sb-123
      & 120.9038D0,     122.9041D0,      0d0,             0d0, &
!      Te-130          Te-125           Te-126         Te-128
      & 129.9067D0,     124.9044D0,      125.9032D0,    127.9047D0, &
!      I-127
      & 126.9004D0,     0d0,             0d0,             0d0, &
!      Xe-132          Xe-129           Xe-131         Xe-134
      & 131.9042D0,     128.9048D0,      130.9051D0,    133.9054D0/)


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

      
