program qumvia_main

! ----------------------------------------------------------------
! QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA
!
!                        WELCOME TO
!
!  .d8888b. 0000 0000  .00..00. .00. 0000  000 00000   d888b
! d88P^ Y88b 888  888 `88bP"b8bP"88b  888  888  888   d88 88b
! 888    888 888  888  888  888  888  888  888  888   888 888
! Y88b..d88P 888  888  888  888  888  Y88bd88Y  888   888o888
!  Y888888P  'V8V"V8P o888 o888 o888    Y88Y   00000 o888 888o
!    `Yp.
!
!                A VSCF/VCI IMPLEMENTATION
!                            by
!               Diego Javier Alonso de Armino
!                   diegoarmino@gmail.com
!
!
! QUantum Mechanical VIbratioal Anlalysis or QUMVIA, for short, is
! an original implementation of Vibrational Self-consistent Field
! (VSCF) and Vibrational Configuration Interaction (VCI) using
! distributed gaussian basis set with analytical integrals and a
! configuration selection algorithm.
!
!
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

#ifdef qvacpu
   use qvamod_cpu
#else
   use qvamod_lio
   use qvamod_lioexcl
   use garcha_mod, only: natom, nsol, Iz, basis_set, fitting_set, &
                         int_basis, omit_bas, verbose, writeforces,&
                         r, rqm
#endif
   use qvamod_common
   use M_kracken

   implicit none

   integer   :: nqmatoms,ierr
   type(qva_cli_type), save  :: qva_cli
   type(qva_nml_type), save  :: qva_nml
#ifdef qvalio
   integer   :: nclatoms, qmcharge
   integer, allocatable :: at_numbers(:)
   real*8, allocatable :: qvageom(:,:)
   type(lio_nml_type), save  :: lio_nml
   double precision :: a0  =   0.5291771D00         ! bohr radius
   integer :: i,j
#endif

!  PARSING COMMAND LINE ARGUMENTS
!  The parsing (and only the parsing) is performed using M_kracken
!  module from John S. Urban.
!  Copyright(C) 1989,1996,2013 John S. Urban   all rights reserved
!  http://home.comcast.net/~urbanjost/LIBRARY/libCLI/arguments/krackenhelp.html

   call kracken('qva','-i input.iqv    &
                &      -o output.oqv   &
                &      -c geom.cqv     &
                &      -n nmodes.nqv   &
                &      -h hess.hqv     &
                &      -q qff.qqv     &
                &      -s stencil.sqv')
   qva_cli%inp= sget('qva_i')
   qva_cli%out= sget('qva_o')
   qva_cli%geo= sget('qva_c')
   qva_cli%nmo= sget('qva_n')
   qva_cli%hes= sget('qva_h')
   qva_cli%qff= sget('qva_q')
   qva_cli%ste= sget('qva_s')

!  OPEN OUTPUT FILE
   open(UNIT=77,FILE=qva_cli%out,ACTION='WRITE',iostat=ierr)
   if (ierr /= 0) then
      write(*,'(A)') 'ERROR OPENNING OUTPUT FILE'
      STOP
   end if

write(77,'(A)') ' QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  '
write(77,'(A)') '        '
write(77,'(A)') '                        WELCOME TO '
write(77,'(A)') '        '
write(77,'(A)') '  .d8888b. 0000 0000  .00..00. .00. 0000 0000 00000   d888b'
write(77,'(A)') ' d88P^ Y88b 888  888  888P"88bP"88b  888  888  888   d88 88b'
write(77,'(A)') ' 888    888 888  888  888  888  888  888  888  888   888 888'
write(77,'(A)') ' Y88b..d88P 888  888  888  888  888  Y88bd88Y  888   888o888'
write(77,'(A)') '  Y888888P  `V8V"V8P o888 o888 o888    Y88Y   00000 o888 888o'
write(77,'(A)') '    `Yp.  '
write(77,'(A)') ''
write(77,'(A)') ''
write(77,'(A)') ' Copyright(C) 2015 Diego J. Alonso de Armiño. All rights reserved'
write(77,'(A)') ''
write(77,'(A)') ' QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  '

!  READ NAMELIST AND FIRST LINE OF GEOMETRY FILE'
   call get_qva_nml(qva_cli%inp,qva_nml)
   call print_qva_nml(qva_nml)
   call readnqmatoms(qva_cli%geo,nqmatoms)
#ifdef qvalio
!  Read geometry file.
   allocate(at_numbers(nqmatoms),qvageom(3,nqmatoms))
   call readgeom(qva_cli,nqmatoms,qvageom,at_numbers)

   call lio_defaults()
   call read_options(qva_cli%inp,qmcharge)
   call read_coords_lio(qva_cli%geo)

   call get_lio_nml(qva_cli%inp,lio_nml)
   call print_lio_nml(lio_nml)
   call init_lio_common(nqmatoms,at_numbers,nsol,qmcharge,0)

   deallocate(at_numbers,qvageom)
#endif
!   close(unit=10)

   IF (qva_nml%nhess .eq. 1) THEN

!     STATE-SPECIFIC VSCF / CONFIGURATION SELECTION VCI
      call run_vscfvci(qva_cli,qva_nml,nqmatoms)
#ifdef qvalio
      call lio_finalize()
#endif
      STOP

   ELSE IF (qva_nml%nhess .eq. 2) THEN

!     DISTORTED GEMETRIES GENERATION FOR HESSIAN CALCULATIONS USING
!     GAUSSIAN OR GAMESS, FOR QUARTIC FORCE FIELD CALCULATION.
      call geoms4qff(qva_cli,qva_nml,nqmatoms)
      STOP

   END IF

end program qumvia_main
