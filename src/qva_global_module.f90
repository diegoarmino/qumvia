module qva_global_module
   implicit none
!  -----------------------------------------------------------
   type qva_nml_type
     integer :: nhess
     real*8  :: vscf_gauswidth
     integer :: vci_qmax1
     integer :: vci_qmax2
     integer :: vci_qmax3
     integer :: vci_qmax4
     integer :: qva_naddref
     integer :: qumvia_qff
     integer :: qumvia_nmc
     integer :: qva_extprog
     integer :: doconfsel
     integer :: csdepth
     integer :: nmorse
     integer :: nsinh
     integer :: hess_norder
     integer :: natharm
     real*8  :: csiterfactor
     real*8  :: ethresh
     real*8  :: resthresh
     real*8  :: selcut1
     real*8  :: selcut2
     real*8  :: qva_dstep
     real*8  :: hess_h
     character(99) :: harmask
   end type qva_nml_type
  
   type qva_cli_type
      character(99) :: inp
      character(99) :: out
      character(99) :: geo
      character(99) :: nmo
      character(99) :: hes
      character(99) :: qff
      character(99) :: ste
   end type qva_cli_type

   type lio_nml_type
      character(len=20) :: basis
      character(len=20) :: output
      character(len=20) :: fcoord
      character(len=20) :: fmulliken
      character(len=20) :: frestart
      character(len=20) :: frestartin
      logical :: verbose
      logical :: OPEN
      integer :: NMAX
      integer :: NUNP
      logical :: VCINP
      real*8  :: GOLD
      real*8  :: told
      real*8  :: rmax
      real*8  :: rmaxs
      logical :: predcoef
      integer :: idip
      logical :: writexyz
      logical :: intsoldouble
      logical :: DIIS
      integer :: ndiis
      real*8  :: dgtrig
      integer :: Iexch
      logical :: integ
      logical :: DENS
      integer :: IGRID
      integer :: IGRID2
      integer :: timedep
      real*8  :: tdstep
      integer  :: ntdstep
      logical :: field
      logical :: exter
      real*8  :: a0
      real*8  :: epsilon
      real*8  :: Fx
      real*8  :: Fy
      real*8  :: Fz
      integer  :: NBCH
      integer :: propagator
      logical :: writedens
      logical :: tdrestart
      integer :: qmcharge
   end type lio_nml_type

!   type(qva_cli_type)  :: qva_cli             ! File I/O info.
!   type(qva_nml_type)  :: qva_nml             ! Qumvia namelist
!   integer             :: nqmatoms,nclatoms   ! Number of qm an mm atoms respect.
!   real*8,allocatable  :: qvageom(:,:)        ! QM coords
!   real*8,allocatable  :: clgeom(:,:)         ! MM coords and charges.
!   integer,allocatable :: atnum_qm(:)         ! Atomic numbers for QM atoms.
!   integer,allocatable :: atnum(:)            ! At numbers for all atoms.
!  -----------------------------------------------------------
end module
