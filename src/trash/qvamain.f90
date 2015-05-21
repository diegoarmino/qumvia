program qumvia_main

! ----------------------------------------------------------------
! QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  
!
! Implementation by Diego Javier Alonso de Armino
!
! Date Febraury 2015
!
! An original implementation of Vibrational Self-consistent Field
! (VSCF) and Vibrational Configuration Interaction (VCI) using 
! distributed gaussian basis set with analytical integrals and a
! configuration selection algorithm.
!
! QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  QUMVIA  
! ----------------------------------------------------------------

#ifdef qvacpu
#warning 'USING CPU MODULES IN MAIN PROGRAM'
   use qvamod_cpu
#elif qvalio
#warning 'USING LIO MODULES IN MAIN PROGRAM'
   use qvamod_lio
#endif
   use qvamod_common

   implicit none

!   type qva_nml_type
!     integer :: nhess
!     real*8  :: vscf_gauswidth
!     integer :: vci_qmax1
!     integer :: vci_qmax2
!     integer :: vci_qmax3
!     integer :: vci_qmax4
!     integer :: qva_naddref
!     integer :: qumvia_qff
!     integer :: qumvia_nmc
!     integer :: qva_extprog
!     real*8  :: ethresh
!     real*8  :: resthresh
!     real*8  :: selcut1
!     real*8  :: selcut2
!     real*8  :: qva_dstep
!   end type qva_nml_type
!
   type(qva_nml_type), save  :: qva_nml
   integer   :: nqmatoms

!  READ NAMELIST AND FIRST LINE OF GEOMETRY FILE
   call get_qva_nml(qva_nml)
   call readnqmatoms(nqmatoms)


   IF (qva_nml%nhess .eq. 1) THEN

!     STATE-SPECIFIC VSCF / CONFIGURATION SELECTION VCI           
      call run_vscfvci(qva_nml,nqmatoms)
      STOP

   ELSE IF (qva_nml%nhess .eq. 2) THEN

!     DISTORTED GEMETRIES GENERATION FOR HESSIAN CALCULATIONS USING
!     GAUSSIAN OR GAMESS, FOR QUARTIC FORCE FIELD CALCULATION.
!     call geoms4qff(qva_extprog,nmodes,eig,dy,nqmatoms,qvageom,at_numbers,ndf,nvdf)
      call geoms4qff(qva_nml,nqmatoms)
      STOP

   END IF

end program qumvia_main

