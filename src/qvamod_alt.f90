#ifdef qvacpu
module qvamod_cpu
#else
module qvamod_lio
#endif

   use qvamod_common

#ifdef qvalio
   use qvamod_lioexcl
#endif


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
! Copyright(C) 2015 Diego J. Alonso de Armiño. All rights reserved.
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

   implicit none
   private
   public :: run_vscfvci,selectqff,ssvscf_csvci2


contains

      subroutine run_vscfvci(qva_cli,qva_nml,nqmatoms)
         implicit none
!        -----------------------------------------------------------
         type(qva_cli_type), intent(in) :: qva_cli
         type(qva_nml_type), intent(in) :: qva_nml
         integer,intent(in) :: nqmatoms

!        LOCAL VARIABLES
         integer :: i
         integer :: at_numbers(nqmatoms)
         integer :: nmcoup
         integer :: ndf
         integer :: nvdf
         integer :: ngaus
         integer :: nconf
         integer :: qmaxx1
         integer :: qmaxx2
         integer :: qmaxx3
         integer :: qmaxx4
         integer :: bdim
         integer :: nrst
         integer :: qumvia_qff
         integer :: qumvia_nmc
         integer :: qva_extprog
         integer :: naddref
         integer,parameter :: nclatoms = 0

         real*8 :: dy
         real*8 :: ethresh
         real*8 :: resthresh
         real*8 :: selcut1
         real*8 :: selcut2
         real*8 :: gwidth
         real*8 :: escf
         real*8 :: clcoords(4,nclatoms)
         real*8 :: Evscf
         real*8 :: qvageom(3,nqmatoms)
         real*8 :: Evirt1(3*nqmatoms-6)
         real*8 :: eig(3*nqmatoms)
         real*8 :: nmodes(3*nqmatoms,3*nqmatoms)
         real*8 :: atmass(nqmatoms)
         real*8 :: L(3*nqmatoms,3*nqmatoms-6)
         real*8 :: freq(3*nqmatoms-6)
   
!        LOCAL ALLOCATABLE VARIABLES
         real*8,allocatable :: Emod(:,:)
         real*8,allocatable :: P(:,:,:)
         real*8,allocatable :: Po(:,:,:)
         real*8,allocatable :: Scho(:,:,:)
         real*8,allocatable :: Q1(:,:,:)
         real*8,allocatable :: Q2(:,:,:)
         real*8,allocatable :: Q3(:,:,:)
         real*8,allocatable :: Hcore(:,:,:)
         real*8,allocatable :: GDmtrx(:,:,:)
         real*8,allocatable :: GTmtrx(:,:,:)
         real*8,allocatable :: hii(:)
         real*8,allocatable :: tiii(:)
         real*8,allocatable :: uiiii(:)
         real*8,allocatable :: tiij(:,:)
         real*8,allocatable :: tjji(:,:)
         real*8,allocatable :: uiiij(:,:)
         real*8,allocatable :: ujjji(:,:)
         real*8,allocatable :: uiijj(:,:)
         real*8,allocatable :: tijk(:,:,:)
         real*8,allocatable :: uiijk(:,:,:)
         real*8,allocatable :: uijjk(:,:,:)
         real*8,allocatable :: uijkk(:,:,:)
         real*8,allocatable :: ECI(:)
         real*8,allocatable :: Cci(:,:)
!        -----------------------------------------------------------

!        Read geometry file.
         call readgeom(qva_cli,nqmatoms,qvageom,at_numbers)

!        Some aliases
         qumvia_qff = qva_nml%qumvia_qff
         qumvia_nmc = qva_nml%qumvia_nmc
         qva_extprog = qva_nml%qva_extprog
         dy=qva_nml%qva_dstep
         gwidth=qva_nml%vscf_gauswidth

         ngaus=16
         ndf=3*nqmatoms
         nvdf=ndf-6

         allocate ( &
            P(ngaus,ngaus,nvdf),     &
            Po(ngaus,ngaus,nvdf),    &
            Scho(ngaus,ngaus,nvdf),  &
            Q1(ngaus,ngaus,nvdf),    &
            Q2(ngaus,ngaus,nvdf),    &
            Q3(ngaus,ngaus,nvdf),    &
            Hcore(ngaus,ngaus,nvdf), &
            hii(nvdf),               &
            tiii(nvdf),              &
            uiiii(nvdf),             &
            tiij(nvdf,nvdf),         &
            tjji(nvdf,nvdf),         &
            uiiij(nvdf,nvdf),        &
            ujjji(nvdf,nvdf),        &
            tijk(nvdf,nvdf,nvdf),    &
            uiijk(nvdf,nvdf,nvdf),   &
            uijjk(nvdf,nvdf,nvdf),   &
            uijkk(nvdf,nvdf,nvdf),   &
            uiijj(nvdf,nvdf)         &
         )
   
         P=0.0d0
         Po=0.0d0
         Evscf=0.0d0
         Evirt1=0.0d0
         Q1=0.0d0
         Q2=0.0d0
         Q3=0.0d0
         Hcore=0.0d0
         hii=0.0d0
         tiii=0.0d0
         uiiii=0.0d0
         tiij=0.0d0
         tjji=0.0d0
         uiiij=0.0d0
         ujjji=0.0d0
         uiijj=0.0d0
         tijk=0.0d0
         uiijk=0.0d0
         uijjk=0.0d0
         uijkk=0.0d0
         ECI=0.0d0
         Cci=0.0d0
   
!        Decide whether if computing hessian from scratch or reading it.
         if (qumvia_qff < 3) then
#ifdef qvalio
            call hessian(qvageom,nqmatoms,at_numbers,nmodes,eig)
#else
            write(77,'(A)') 'FATAL ERROR: qumvia_qff<3 is only valid '
            write(77,'(A)') 'for QUMVIA_LIO. '
            STOP
#endif
   
         else if (qumvia_qff == 5) then
            call readgaunmodes(qva_cli,nqmatoms,ndf,nvdf,L,freq,hii,atmass)
            nmodes=0d0
            nmodes(:,7:ndf)=L
            eig=0d0
            eig(7:ndf)=hii
         end if
   
!        ---------------------------------------------------------------
!        STATE-SPECIFIC VSCF/VCI           
!        ---------------------------------------------------------------
!        Number of additional reference configurations
         naddref=qva_nml%qva_naddref
!        Max excitation in singles
         qmaxx1=qva_nml%vci_qmax1
!        Max excitqvaon in doubles
         qmaxx2=qva_nml%vci_qmax2
!        Max excitation in triples
         qmaxx3=qva_nml%vci_qmax3
!        Max excitation in quadruples
         qmaxx4=qva_nml%vci_qmax4
!        Overall Energy Threshold
         ethresh=qva_nml%ethresh
!        Resonance cutoff
         resthresh=qva_nml%resthresh
!        Selection cutoff for 1st order resonance
         selcut1=qva_nml%selcut1
!        Selection cutoff for 2nd order resonance
         selcut2=qva_nml%selcut2
!        Total number of configurations
         nconf = 1 + &
               & nvdf*qmaxx1 + &
               & nvdf*(nvdf-1)*qmaxx2**2/2 + &
               & nvdf*(nvdf-1)*(nvdf-2)*qmaxx3**3/6 + &
               & nvdf*(nvdf-1)*(nvdf-2)*(nvdf-3)*qmaxx4**4/24
         write(77,'(A,I15)') 'NUMBER OF INITIAL CONFIGURATIONS: ',nconf
         if (nconf > 500000) nconf=500000
!        Number of reference states to compute
         write(*,'(A)') '***BEGINNING STATE-SPECIFIC VSCF/csVCI CALCULATION***'
         write(77,'(A)') '***BEGINNING STATE-SPECIFIC VSCF/csVCI CALCULATION***'
         nrst=1+nvdf+naddref
!        DANGER version 2 of the following subroutine
         call ssvscf_csvci2(qva_cli,qumvia_nmc,qumvia_qff,ethresh,resthresh,selcut1,&
                         &selcut2,ndf,nvdf,ngaus,nmcoup,nqmatoms,nclatoms,&
                         &qmaxx1,qmaxx2,qmaxx3,qmaxx4,nconf,at_numbers,dy,gwidth,&
                         &eig,nmodes,qvageom,clcoords,naddref,nrst)
   
         deallocate(P,Po,Q1,Q2,Q3,Hcore,hii,tiii,uiiii,tiij,tjji,&
                   &uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk)
   
      end subroutine




!######################################################################
!     FORCE FIELD SELECTION
!######################################################################
      
      subroutine selectqff(qva_cli,qumvia_qff,nmodes,eig,dy,nqmatoms,&
                     &nclatoms,qmcoords,clcoords,at_numbers,ndf,nvdf,&
                     &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj,&
                     &tijk,uiijk,uijjk,uijkk)
      implicit none

      type(qva_cli_type), intent(in) :: qva_cli
      integer, intent(in) :: qumvia_qff   
      integer, intent(in) :: ndf                  ! Number of deg of freedm
      integer, intent(in) :: nvdf                 ! Number of vib degrees of
      integer, intent(in) :: nqmatoms             ! Number of QM atoms
      integer, intent(in) :: nclatoms             ! Number of classical atoms
      integer, intent(in) :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
      real*8, intent(in)  :: nmodes(ndf,ndf)      ! mw Hessian eigenvectors.
      real*8, intent(in)  :: eig(ndf)             ! mw Hessian eigenvalues.
      real*8, intent(in)  :: dy                   ! Step size factor dy.
      real*8, intent(in)  :: qmcoords(3,nqmatoms) ! QM atom coordinates
      real*8, intent(in)  :: clcoords(4,nclatoms) ! MM atom coordinates and
      real*8, intent(out) :: hii(nvdf)           ! Diagonal cubic coupling
      real*8, intent(out) :: tiii(nvdf)           ! Diagonal cubic coupling
      real*8, intent(out) :: uiiii(nvdf)          ! Diagonal quartic coupling
      real*8, intent(out) :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
      real*8, intent(out) :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
      real*8, intent(out) :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling
      real*8, intent(out) :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling
      real*8, intent(out) :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling
      real*8, intent(out) :: tijk(nvdf,nvdf,nvdf)      ! 3-mode cubic coupling
      real*8, intent(out) :: uiijk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling
      real*8, intent(out) :: uijjk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling
      real*8, intent(out) :: uijkk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling


      write(*,*) 'BEGINNING QUARTIC FORCE FIELD COMPUTATION'

      select case (qumvia_qff)
         case (1)
#ifdef qvalio
         call fullnumqff(nmodes,eig,dy,nqmatoms,nclatoms,&
                     &qmcoords,clcoords,at_numbers,ndf,nvdf,&
                     &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj,&
                     &tijk,uiijk,uijjk,uijkk)
#else 
            write(77,'(A)') 'FATAL ERROR: qumvia_qff<3 is only valid '
            write(77,'(A)') 'for QUMVIA_LIO. '
            STOP
#endif
         case (2)
#ifdef qvalio
         call seminumqff(nmodes,eig,dy,nqmatoms,nclatoms,&
                     &qmcoords,clcoords,at_numbers,ndf,nvdf,&
                     &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj, &
                     &tijk,uiijk,uijjk,uijkk)
#else
            write(77,'(A)') 'FATAL ERROR: qumvia_qff<3 is only valid '
            write(77,'(A)') 'for QUMVIA_LIO. '
            STOP
#endif
         case (3)
         call readqff(qva_cli,nvdf,hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
                 & tijk,uiijk,uijjk,uijkk)
         case (4)
         call readgamessqff(qva_cli,nvdf,hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
                 & tijk,uiijk,uijjk,uijkk)
         case (5)
         call hseminumqff(qva_cli,dy,nqmatoms,nclatoms,qmcoords,&
                         &clcoords,at_numbers,ndf,nvdf,&
                         &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,&
                         &uiijj,tijk,uiijk,uijjk,uijkk)
         case default
         call hseminumqff(qva_cli,dy,nqmatoms,nclatoms,qmcoords,&
                         &clcoords,at_numbers,ndf,nvdf,&
                         &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,&
                         &uiijj,tijk,uiijk,uijjk,uijkk)
      end select

      end subroutine



       subroutine ssvscf_csvci2(qva_cli,qumvia_nmc,qumvia_qff,ethresh,resthresh,selcut1,&
                          &selcut2,ndf,nvdf,ngaus,nmcoup,nqmatoms,nclatoms,&
                          &qmaxx1,qmaxx2,qmaxx3,qmaxx4,nconf,at_numbers,dy,gwidth,&
                          &eig,nmodes,qmcoords,clcoords,naddref,nrst)
 !     -----------------------------------------------------------------
 !     THIS SUBROUTINE PERFORMS VIBRATIONAL SELF-CONSISTENT FIELD 
 !     FOLLOWED BY CONFIGURATION INTERACTION CALCULATION IN A 
 !     STATE-SPECIFIC FASHION.
 !     A separate VSCF and VCI computation is performed for each state
 !     of interest, by default the ground state and the first excited 
 !     states of each normal mode, using that state as the reference.
 !     The reference state is the vscf configuration over which the
 !     effective potential is computed.
 !     -----------------------------------------------------------------
       implicit none
 
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in)  :: nrst
       integer,intent(in)  :: qumvia_nmc           ! # of coupled mode in Hci
       integer,intent(in)  :: qumvia_qff           ! # keyword for qff read/calc
       integer,intent(in)  :: naddref              ! # of additional ref confs.
       integer,intent(in)  :: ndf                  ! Number of classical atoms
       integer,intent(in)  :: nvdf                 ! Number of classical atoms
       integer,intent(in)  :: ngaus                ! Number of classical atoms
       integer,intent(in)  :: nmcoup               ! Number of normal modes to couple in QFF.
       integer,intent(in)  :: nqmatoms             ! Number of QM atoms
       integer,intent(in)  :: nclatoms             ! Number of classical atoms
       integer,intent(in)  :: qmaxx1               ! Max excitation for singles.
       integer,intent(in)  :: qmaxx2               ! Max excitation for doules.
       integer,intent(in)  :: qmaxx3               ! Max excitation for triples
       integer,intent(in)  :: qmaxx4               ! Max excitation for quadruples
       integer,intent(in)  :: nconf                ! Dimension of CI basis set.
       integer,intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
       real*8,intent(in)   :: dy                   ! Step size factor for num derivatives.
       real*8,intent(in)   :: resthresh            ! Step size factor for num derivatives.
       real*8,intent(in)   :: ethresh              ! Step size factor for num derivatives.
       real*8,intent(in)   :: selcut1              ! Step size factor for num derivatives.
       real*8,intent(in)   :: selcut2              ! Step size factor for num derivatives.
       real*8,intent(in)   :: gwidth               ! Width factor for gaussian primitives.
       real*8,intent(in)   :: eig(ndf)             ! Eigenvalues of the hessian matrix.
       real*8,intent(in)   :: nmodes(ndf,ndf)      ! Eigenvectors of the hessian matrix.
       real*8,intent(in)   :: qmcoords(3,nqmatoms) ! QM atom coordinates
       real*8,intent(in)   :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
 
       integer            :: i
       integer            :: ref(8)
       integer            :: cnf
       integer            :: bdim
       integer            :: trash
       integer            :: vscfstat
       real*8             :: P(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
       real*8             :: Po(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
       real*8             :: Scho(ngaus,ngaus,nvdf)  ! Cholesky factor.
       real*8             :: Evscf                ! VSCF energy.
       real*8             :: Q1(ngaus,ngaus,nvdf) !  <gi|Qa^1|gj> matrix
       real*8             :: Q2(ngaus,ngaus,nvdf) !  <gi|Qa^2|gj> matrix
       real*8             :: Q3(ngaus,ngaus,nvdf) !  <gi|Qa^3|gj> matrix
       real*8             :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian (non-orthogonal basis)
       real*8             :: Essvhf(nrst)
       real*8             :: Pss(ngaus,nrst)
       real*8             :: Poss(ngaus,nrst)
       real*8             :: GDmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8             :: GTmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8             :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8             :: hii(nvdf)             ! Diagonal Hessian eigenvalues.
       real*8             :: tiii(nvdf)           ! Diagonal cubic coupling terms
       real*8             :: uiiii(nvdf)          ! Diagonal quartic coupling terms
       real*8             :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
       real*8             :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
       real*8             :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8             :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8             :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8             :: tijk(nvdf,nvdf,nvdf)      ! 3-mode cubic coupling terms
       real*8             :: uiijk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
       real*8             :: uijjk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
       real*8             :: uijkk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
       real*8             :: Evci(nrst)
       real*8             :: Eref
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
 
 !     Allocatable variables.
       integer,allocatable,dimension(:,:) :: addrefs
       integer,allocatable,dimension(:,:) :: fundrefs
       integer,allocatable,dimension(:,:) :: refstat
 
 !     -----------------------------------------------------------------
 !     GENERATING REFERENCE FONFIGURATIONS
      
       write(77,'(A,I3)') 'NUMBER OF ADDITIONAL REFERENCE CONFIGS:',naddref
       if (naddref > 0) then 
          allocate ( addrefs(8,naddref),fundrefs(8,nvdf+1) )
          addrefs=0
          fundrefs=0
          call readaddref(naddref,addrefs)
          call genConf3(0,0,0,0,hii,ethresh,nvdf,nvdf+1,trash,fundrefs)
          allocate ( refstat(8,nrst) )
          refstat=0
          refstat(:,1:nvdf+1) = fundrefs
          refstat(:,nvdf+2:nrst) = addrefs
          deallocate (addrefs,fundrefs)
       else
          allocate ( refstat(8,nrst) )
          refstat=0
          call genConf3(0,0,0,0,hii,ethresh,nvdf,nvdf+1,trash,refstat)
       end if
       write(77,'(A)') 'REFERENCE CONFIGURATIONS'
       do i=1,nrst+naddref
          write(77,'(8I3)') refstat(:,i)
       end do
 
 
 !     QUARTIC FORCE FIELD
 !     ---------------------------------------------------------------
       call selectqff(qva_cli,qumvia_qff,nmodes,eig,dy,nqmatoms,nclatoms,&
                      &qmcoords,clcoords,at_numbers,ndf,nvdf,&
                      &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj,&
                      &tijk,uiijk,uijjk,uijkk)
 
 !     COMPUTING VSCF/VCI OVER THESE CONFIGURATIONS      
       do cnf=1,nrst
 !        Computing VSCF
          ref(:) = refstat(:,cnf)
          call ssvscf2(ref,nmodes,eig,dy,nmcoup,gwidth,nqmatoms,nclatoms,&
                    &ngaus,qmcoords,clcoords,at_numbers,ndf,nvdf,&
                    &P,Po,Scho,Evscf,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Emod,&
                    &hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk,vscfstat)
          if (vscfstat /=0) CYCLE
          Essvhf(cnf)=Evscf
          Pss(:,cnf)=P(:,cnf,cnf)
          Poss(:,cnf)=Po(:,cnf,cnf)
 
 !        Computing VCI
          write(77,'(A,8I3)') 'COMPUTING VCI FOR REFERENCE STATE ',ref
          bdim=qmaxx1+1
          call csVCI2(ref,qumvia_nmc,ethresh,resthresh,selcut1,selcut2,&
                    &Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,Emod,&
                    &hii,tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nconf,ngaus,nvdf,qmaxx1,qmaxx2,qmaxx3,qmaxx4,bdim,Eref)
          Evci(cnf)=Eref*h2cm
       end do
 
       write(77,'(A)') 'ssVSCF AND csVCI ENERGIES AND TRANSITIONS (CM-1)'
       write(77,'(A)') 'No.       E(ssVSCF)      dE(ssVSCF)        E(csVCI)        dE(csVCI)'
       do cnf=1,nrst
          write(77,'(I3,4F16.2)') cnf-1, Essvhf(cnf), Essvhf(cnf)-Essvhf(1),&
                                 & Evci(cnf),Evci(cnf)-Evci(1)
       end do
 
 
       end subroutine

#ifdef qvacpu
end module qvamod_cpu
#else
end module qvamod_lio
#endif

 
