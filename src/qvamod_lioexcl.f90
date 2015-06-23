module qvamod_lioexcl

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

   use qvamod_common

   implicit none
   private
   public :: fullnumqff,seminumqff,hessian,lio_nml_type,get_lio_nml


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


contains

!######################################################################
!    LIO INITIALIZATION SECTION
!######################################################################


     subroutine get_lio_nml(qvain,lio_nml)
        implicit none
     
        character(99),intent(in)        :: qvain
        type(lio_nml_type), intent(out) :: lio_nml

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
        logical :: field
        logical :: exter
        real*8  :: tdstep
        integer  :: ntdstep
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

        namelist /lio/ basis,output,fcoord,fmulliken,frestart,frestartin &
        ,verbose,OPEN,NMAX,Nunp,VCINP, &
        GOLD,told,rmax,rmaxs,predcoef, &
        idip,writexyz,intsoldouble,DIIS,ndiis,dgtrig, &
        Iexch,integ,dens,igrid,igrid2,timedep, tdstep, ntdstep, &
        propagator,NBCH, writedens, tdrestart, &
        field,exter,a0,epsilon,Fx,Fy,Fz,qmcharge
   
        integer :: ifind, ierr
        !defaults
   
        basis='basis'  ! name of the base file
        output='output'
        fcoord='qm.xyz'
        fmulliken='mulliken'
        frestart='restart.out'
        frestartin='restart.in'
        verbose=.false.
        OPEN=.false.
        NMAX= 100
        NUNP= 0
        VCINP= .false.
        GOLD= 10.
        told=1.0D-6
        rmax=16
        rmaxs=5
        predcoef=.false.
        idip=1
        writexyz=.true.
        intsoldouble=.true.
        DIIS=.true.
        ndiis=30
        dgtrig=100.
        Iexch=9
        integ=.true.
        DENS = .true.
        IGRID = 2
        IGRID2 = 2
        timedep = 0
        tdstep = 2.D-3
        ntdstep= 1
        field=.false.
        exter=.false.
        a0=1000.0
        epsilon=1.D0
        Fx=0.05
        Fy=0.05
        Fz=0.05
        NBCH=10
        propagator=1
        writedens=.true.
        tdrestart=.false.
        qmcharge=0
   
        ! Read namelist
!        open(UNIT=10,FILE=qvain,action='READ',iostat=ierr)
!        if (ierr /= 0) then
!           write(*,*) 'ERROR OPENING LIO NAMELIST'
!           STOP
!        end if

        rewind 10
        read(10,nml=lio,iostat=ierr)
   
        if ( ierr > 0 ) then
           write(*,'(A)') 'LIO NAMELIST READ ERROR'
        else if ( ierr < 0 ) then
           write(*,'(A)') '&lio namelist read success'
        end if

        close(unit=10,iostat=ierr)
        if (ierr /= 0) then
           write(*,*) 'ERROR CLOSING LIO NAMELIST'
           STOP
        end if
   
   
        lio_nml%basis=basis
        lio_nml%output=output
        lio_nml%fcoord=fcoord
        lio_nml%fmulliken=fmulliken
        lio_nml%frestart=frestart
        lio_nml%frestartin=frestartin
        lio_nml%verbose=verbose
        lio_nml%OPEN=OPEN
        lio_nml%NMAX=NMAX
        lio_nml%NUNP=NUNP
        lio_nml%VCINP=VCINP
        lio_nml%GOLD=GOLD
        lio_nml%told=told
        lio_nml%rmax=rmax
        lio_nml%rmaxs=rmaxs
        lio_nml%predcoef=predcoef
        lio_nml%idip=idip
        lio_nml%writexyz=writexyz
        lio_nml%intsoldouble=intsoldouble
        lio_nml%DIIS=DIIS
        lio_nml%ndiis=ndiis
        lio_nml%dgtrig=dgtrig
        lio_nml%Iexch=Iexch
        lio_nml%integ=integ
        lio_nml%DENS =DENS
        lio_nml%IGRID =IGRID
        lio_nml%IGRID2 =IGRID2
        lio_nml%timedep =timedep
        lio_nml%tdstep = tdstep
        lio_nml%ntdstep = ntdstep
        lio_nml%field=field
        lio_nml%exter=exter
        lio_nml%a0=a0
        lio_nml%epsilon=epsilon
        lio_nml%Fx=Fx
        lio_nml%Fy=Fy
        lio_nml%Fz=Fz
        lio_nml%NBCH=NBCH
        lio_nml%propagator=propagator
        lio_nml%writedens=writedens
        lio_nml%tdrestart=tdrestart
        lio_nml%qmcharge=qmcharge
     end subroutine get_lio_nml
 
!######################################################################
!    HESSIAN AND HARMONIC OSCILLATOR SECTION
!######################################################################

     subroutine projHessian(hess,qmcoords,at_numbers,ndf,nqmatoms,hessp)
     implicit none
!    ------------------------------------------------------------------
     integer,intent(in)  :: ndf
     integer,intent(in)  :: nqmatoms
     integer,intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
     real*8,intent(in)   :: hess(ndf,ndf)
     real*8,intent(in)   :: qmcoords(3,nqmatoms)
     real*8,intent(out)  :: hessp(ndf,ndf)
!    Internal variables
     integer             :: i,j,k,l,n1,n2,nz,cnt
     integer             :: ia,ib,ig,ja,jb,jg,iat,jat
     integer             :: II,JJ
     integer             :: zeros(2),nonz(2)
     integer             :: ipiv(3)
     real*8              :: levic(3,3,3)
     real*8              :: Mass(ndf)
     real*8              :: X0(3,nqmatoms)
     real*8              :: Itsr(3,3)
     real*8              :: innp(3,3)
     real*8              :: outp(3,3)
     real*8              :: ai(3)
     real*8              :: atmp
     real*8              :: cutoff
     real*8              :: det
     real*8              :: trp
     real*8              :: summ
     real*8              :: totM
     real*8              :: P(ndf,ndf)
     real*8              :: tmp(ndf,ndf)
     real*8              :: cmass(3)

     real*8,external  :: ddot

     real*8,  dimension(:), allocatable :: WORK
     integer                            :: LWORK,INFO
     include "qvmbia_param.f"
!    ------------------------------------------------------------------

     cutoff = 1.0d-08
     levic = reshape((/  0.0D+00,  0.0D+00,  0.0D+00,&
                 0.0D+00,  0.0D+00, -1.0D+00,&
                 0.0D+00,  1.0D+00,  0.0D+00,&
                 0.0D+00,  0.0D+00,  1.0D+00,&
                 0.0D+00,  0.0D+00,  0.0D+00,&
                -1.0D+00,  0.0D+00,  0.0D+00,&
                 0.0D+00, -1.0D+00,  0.0D+00,&
                 1.0D+00,  0.0D+00,  0.0D+00,&
                 0.0D+00,  0.0D+00,  0.0D+00  /),  &
                 (/3,3,3/))
!    Converting initial coordinates to AU
     X0=qmcoords/a0
!     write(77,'(A)') 'CONVERSION FACTOR a0 -> Angs'
!     write(77,'(D12.3)') a0
      
!    Computing total mass and atomic mass array.
     totM=0.0d0
     do i=1,nqmatoms
         do j=1,3
             k=at_numbers(i)
             Mass(3*(i-1)+j) = sqrt_atomic_masses_au(k)
         end do
         totM=totM+atomic_masses_au(k)
     end do

!    Computing center of mass
     do i=1,nqmatoms
        do j=1,3
           cmass(j)=cmass(j)+Mass(3*(i-1)+j)*X0(j,i)
        end do
     end do
     cmass = cmass/totM

!    Translating to center of mass and mass-weighting.
     do i=1,nqmatoms
        do j=1,3
           X0(j,i)=(X0(j,i)-cmass(j))*Mass(3*(i-1)+j)
        end do
     end do


!    Moment of inertia tensor.
!    Itsr = Sum_i [ai**T.ai - ai.ai**T]
     Itsr=0.0d0
     do i=1,nqmatoms
        ai=X0(:,i)
        innp=0.0d0
        outp=0.0d0
        innp(1,1) = ddot(3,ai,1,ai,1)
        innp(2,2) = innp(1,1)
        innp(3,3) = innp(1,1)
        call dger(3,3,1d0,ai,1,ai,1,outp,3)
!        call dsyr('U',3,1d0,ai,1,outp,3)
        Itsr=Itsr+innp-outp
     end do

!    Symmetrizing inertia tensor.
     do i=1,2
        do j=i+1,3
          Itsr(j,i)=Itsr(i,j)
        end do
     end do

!     write(77,'(A)') 'INERTIA TENSOR'
!     do i=1,3
!        write(77,'(3D11.2)') Itsr(i,:)
!     end do

!    Invert moment of inertia matrix.
!    We first handle cases with zeroes 
     atmp = Itsr(1,1)*Itsr(2,2)*Itsr(3,3)
     if (abs(atmp) < cutoff) then
        cnt=0
        nz=0
        do i=1,3
           if (Abs(Itsr(i,i)) < cutoff) then
              cnt=cnt+1
              zeros(cnt)=i
           else
              nz=nz+1
              nonz(nz)=i
           end if
        end do
        if (cnt == 2) then
           n1 = nonz(1)
           Itsr(n1,n1)=1d0/Itsr(n1,n1)
        else if (cnt == 1) then
           n1 = nonz(1)
           n2 = nonz(2)
           det = Itsr(n1,n1)*Itsr(n2,n2)-Itsr(n1,n2)*Itsr(n2,n1)
           trp = Itsr(n1,n1)
           Itsr(n1,n1)=Itsr(n2,n2)/det
           Itsr(n2,n2)=trp/det
           Itsr(n1,n2)=-Itsr(n1,n2)/det
           Itsr(n2,n1)=-Itsr(n2,n1)/det
        end if
     else
        if (allocated(WORK)) deallocate(WORK)
        allocate(WORK(100))

        INFO=0
        LWORK=-1
        call dsytrf('U',3,Itsr,3,ipiv,WORK,LWORK,INFO)
        if (INFO /= 0) STOP('ERROR IN INERTIA FACTORIZATION 1')

        INFO=0
        LWORK=WORK(1)
        deallocate(WORK)
        allocate(WORK(LWORK))
        call dsytrf('U',3,Itsr,3,ipiv,WORK,LWORK,INFO)
        if (INFO /= 0) STOP('ERROR IN INERTIA FACTORIZATION 2')
        deallocate(WORK)

        INFO=0
        allocate(WORK(3))
        call dsytri('U',3,Itsr,3,ipiv,WORK,INFO)
        deallocate(WORK)
        if (INFO /= 0) STOP('ERROR IN INERTIA TENSOR INVERSION')
     end if
     do i=1,2
        do j=i+1,3
           Itsr(j,i)=Itsr(i,j)
        end do
     end do
!     write(77,'(A)') 'INVERSE INERTIA TENSOR'
!     do i=1,3
!        write(77,'(3D11.2)') Itsr(i,:)
!     end do

!    WE NOW COMPUTE THE PROJECTOR MATRIX P
!    THE FORMULA CAN BE FOUND IN Miller, Handy and Adams, JCP 72:99,(1980) eq 4.11
!    We then project the hessian matrix using equation 1.5a

!                  PHess = (1-P).Hess.(1-P)

!    BUILDING PROJECTOR P
     do iat=1,nqmatoms
        do ig=1,3
           do jat=1,nqmatoms
              do jg=1,3
                 II = 3*(iat-1)+ig
                 JJ = 3*(jat-1)+jg

                 summ=0d0
                 do ia=1,3
                    do ib=1,3
!                       if (levic(ia,ib,ig)==0) CYCLE
                       do ja=1,3
                          do jb=1,3
!                          if (levic(ja,jb,jg)==0) CYCLE 
                             summ=summ+levic(ia,ib,ig)*levic(ja,jb,jg) &
                                  & *Itsr(ia,ja)*X0(ib,iat)*X0(jb,jat)
                          end do
                       end do
                    end do
                 end do
                 P(II,JJ) = summ

                 if (ig == jg) then
                    P(II,JJ)=P(II,JJ)+Mass(II)*Mass(JJ)/totM
                 end if

              end do
           end do
        end do
     end do
!     write(77,'(A)') 'PROJECTOR P'
!     do i=1,ndf
!        write(77,'(999D9.2)') P(i,:)
!     end do


!    COMPUTING (1-P)
     P=-P
     do i=1,ndf
        P(i,i) = 1.0d0 + P(i,i)
     end do
     do i=1,ndf-1
        do j=i+1,ndf
           if (abs(P(i,j)) < cutoff) P(i,j)=0d0
           P(j,i)=P(i,j)
        end do
     end do
!     write(77,'(A)') 'PROJECTOR 1-P'
!     do i=1,ndf
!        write(77,'(999D9.2)') P(i,:)
!     end do


!    PROJECTION
!    HESSP = (1-P).HESS.(1-P)
     hessp=0d0
     call dsymm('L','U',ndf,ndf,1d0,hess,ndf,P,ndf,0d0,tmp,ndf)
     call dgemm('N','N',ndf,ndf,ndf,1d0,P,ndf,tmp,ndf,0d0,hessp,ndf)
     do i=1,ndf
        do j=1,ndf
           if (abs(hessp(i,j)) < 1d-10) hessp(i,j)=0d0
        end do
     end do
!     write(77,'(A)') 'PROJECTED HESSIAN'
!     do i=1,ndf
!        write(77,'(999D10.2)') hessp(i,:)
!     end do
    
     return

     end subroutine

     subroutine fixphase(hess,qmcoords,ndf,nqmatoms)
!    ------------------------------------------------------------------------ 
!    FIXES THE PHASE OF THE EIGENVECTORS OF THE HESSIAN BY COMPUTING THE 
!    DOT PRODUCT BETWEEN THESE AND THE INITIAL GEOMETRY, WHICH PROVIDE 
!    A FIXED REFERENCE. IF THE DOT PRODUCT IS NEGATIVE, THE SIGN OF THE 
!    EIGENVECTOR IS CHANGED.
!    ------------------------------------------------------------------------ 

     implicit none
!    ------------------------------------------------------------------------ 
     integer,intent(in)   :: ndf
     integer,intent(in)   :: nqmatoms
     real*8,intent(in)    :: qmcoords(3,ndf)
     real*8,intent(inout) :: hess(ndf,ndf)
!    Internal variables
     integer       :: a
     integer       :: at
     integer       :: x
     real*8        :: vec(ndf)
     real*8        :: prod
     real*8        :: X0(ndf)
!    External functions
     real*8,external :: ddot
!    ------------------------------------------------------------------------ 

     do at=1,nqmatoms
        do x=1,3
           X0(3*(at-1)+x) = qmcoords(x,at)
        end do 
     end do
     
     do a=1,ndf
        vec = 0d0
        prod = 0d0
        vec = hess(:,a)
        prod = ddot(ndf,X0,1,vec,1)
        if (prod < 0d0) then
           hess(:,a) = -vec
        end if
     end do

     end subroutine


     subroutine hessian(qmcoords,nqmatoms,at_numbers,nmodes,eig)

      implicit none

!     -------------------------------------------------------------------------
      real*8,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      integer, intent(in) :: nqmatoms             ! Number of QM atoms
      integer, intent(in) :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
      real*8,  intent(out):: nmodes(3*nqmatoms,3*nqmatoms)  ! Hessian eigenvectors
      real*8,  intent(out):: eig(3*nqmatoms)      ! Hessian eigenvalues
!     LOCAL VARIABLES
      integer,parameter   :: nclatoms=0
      integer             :: ncoords              ! 3*nqmatoms
      integer             :: i,j,k,r,s
      real*8              :: dxyzcl(3,nclatoms)   ! SCF MM force
      real*8              :: dxyzqm(3,nqmatoms)   ! SCF QM force
      real*8              :: h                    ! Step for num differentiation.
      real*8              :: hessf                ! 1/(2*h)

      real*8              :: escf                 ! SCF energy
      real*8              :: dipxyz(3)            ! Dipole moment
      real*8              :: grad(3*nqmatoms+1,3*nqmatoms)
      real*8              :: hess(3*nqmatoms,3*nqmatoms)
      real*8              :: hessp(3*nqmatoms,3*nqmatoms)
      real*8              :: at_masses(3*nqmatoms)
      real*8              :: sign_eig(3*nqmatoms)
      real*8              :: freq(3*nqmatoms)
      real*8              :: qmxyz(3,nqmatoms)
      real*8              :: clcoords(4,nqmatoms)

      type(qva_nml_type), save     :: qva_nml

!     Workspace for diagonalization. 
      character(len=117) :: format1
      character(len=117) :: format2
      real*8,  dimension(:), allocatable :: WORK, WORK2
      integer, dimension(:), allocatable :: IWORK, IWORK2
      integer :: LWORK, LIWORK, INFO

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
      real*8, PARAMETER :: EMASS_CODATA08 = 0.0005485799111D0 ! Electron mass in Daltons (me * Na). The mass of 1 mol of electrons.
      real*8, PARAMETER :: AMU_TO_AU = 1.0d0/EMASS_CODATA08   ! Conversion factor Da to AU
      real*8, PARAMETER :: atomic_masses(18) =&
          (/  1.00866491574D0, 4.002603254153D0,  7.016004D0,  9.012182D0, 11.009305D0, &
             12.0D0,          14.0030740048D0,   15.99491461956D0, 18.998403D0,         &
             19.992440D0,     22.989770D0,       23.985042D0,      26.981538D0,         &
             27.976927D0,     30.973762D0,       31.972071D0,      34.968853D0,         &
             39.962383D0 /)
      real*8, PARAMETER :: atomic_masses_au(18) = atomic_masses * AMU_TO_AU
      real*8, PARAMETER :: sqrt_atomic_masses_au(18) = SQRT(atomic_masses_au)
      real*8, PARAMETER :: invsqrt_atomic_masses_au(18) = 1.0d0/SQRT(atomic_masses_au)
      real*8, PARAMETER :: D1=1.0D00, D2=2.0D00
      real*8, PARAMETER :: PI  =   D2*ACOS(0.0D00)      ! PI
      real*8, PARAMETER :: h2cm=219475.64d0  ! Convert Hartree to cm-1

!     -------------------------------------------------------------------------
!
!     Set up some constants
!     First, WAVE: factor for converting frequencies from AU to cm-1
!
!     -------------------------------------------------------------------------
!     DEBUG
!     -------------------------------------------------------------------------
!     write(77,'(A)') 'AMU_TO_AU'
!     write(77,'(D14.4)') AMU_TO_AU
!     write(77,'(A)') 'atomic masses in AU'
!     write(77,'(18D14.4)') atomic_masses_au
!     write(77,'(A)') 'sqrt_atomic_masses_au'
!     write(77,'(18D14.4)') sqrt_atomic_masses_au
!     write(77,'(A)') 'invsqrt_atomic_masses_au'
!     write(77,'(18D14.4)') invsqrt_atomic_masses_au
!     -------------------------------------------------------------------------
!
      write(*,'(A)') 'BEGINING HESSIAN CALCULATION'
      h=1D-3                           ! Step for numerical diferentiation.
      hessf=0.5D0/h                   ! 1/(2*h)
      ncoords=3*nqmatoms

      dxyzqm=0.0D0
      dxyzcl=0.0D0
      grad=0.0D0
      hess=0.0D0
      eig=0.0D0
      freq=0.0D0
      sign_eig=1.0D00
      at_masses=0.0D0

      write(*,'(A)') 'BEFORE CALLING LIO'
!     Calculating energy and gradient at the initial geometry.
      call SCF_in(escf,qmcoords,clcoords,nclatoms,dipxyz)
      write(*,'(A)') 'AFTER FIRST CALL TO LIO'
      call dft_get_qm_forces(dxyzqm)
      write(*,'(A)') 'AFTER SECOND CALL TO LIO'
!      call dft_get_mm_forces(dxyzcl,dxyzqm)

      do i=1,nqmatoms
          do j=1,3
              grad(1,3*(i-1)+j)=dxyzqm(j,i)
          end do
      end do
!      write(77,'(99D12.3)') grad(1,:)

!     Calculating gradients in displaced positions along the 
!     QM coordinates.
      do i=1,nqmatoms
          do j=1,3
              k=3*(i-1)+j+1
              qmxyz=qmcoords   ! In Angstroms. Transfromed in SCF_in
              qmxyz(j,i)=qmxyz(j,i)+h
              call SCF_in(escf,qmxyz,clcoords,nclatoms,dipxyz)
              call dft_get_qm_forces(dxyzqm)
!              call dft_get_mm_forces(dxyzcl,dxyzqm)
              do r=1,nqmatoms
               do s=1,3
                 grad(k,3*(r-1)+s)=dxyzqm(s,r)
               end do
              end do
          end do
      end do
      write(*,'(A)') 'JUST COMPUTED ALL GRADIENTS'


!     Create a vector of inverse square-root masses corresponding to each 
!     cartesian coordinate.
      do i=1,nqmatoms
          do j=1,3
              k=at_numbers(i)
              at_masses(3*(i-1)+j) = invsqrt_atomic_masses_au(k)
          end do
      end do

!     Calculating non diagonal elements of the Hessian    
      do i=1,ncoords
          do j=i,ncoords
              hess(i,j)=grad(i+1,j)-grad(1,j)+grad(j+1,i)-grad(1,i)
              hess(i,j)=at_masses(i)*at_masses(j)*hessf*hess(i,j)
          end do
      end do

!     Projecting out translational and rotaional modes form hessian
      call projHessian(hess,qmcoords,at_numbers,3*nqmatoms,nqmatoms,hessp)
      hess = hessp

!     Write hessian to file
!      write(77,*), '     WRITING HESSIAN MATRIX     '
!      write(77,*), '     ----------------------     '
!      do i=1,ncoords
!          write(77,'(I3,999D11.3)'), i, ( hess(i,j), j=1,ncoords )
!      end do

!     Scaling hessian to avoid numerical problems.
      hess=hess*1.0D03


!     DIAGONALIZATION OF THE HESSIAN
      allocate ( work(1000), IWORK(1000) )
      LWORK=-1
      call dsyevd('V','U',ncoords,hess,ncoords,eig,WORK,LWORK,IWORK,LWORK,INFO)
      LWORK=WORK(1)
      LIWORK=IWORK(1)
      if(allocated(WORK2)) deallocate (WORK2,IWORK2)
      allocate (WORK2(LWORK),IWORK2(LIWORK))
      call dsyevd('V','U',ncoords,hess,ncoords,eig,WORK2,LWORK,IWORK2,LIWORK,INFO)
      deallocate( WORK, IWORK, WORK2, IWORK2 )


      eig=eig/AMU_TO_AU
!     Convert frequecies to wavenumbers
      do i = 1,ncoords
         freq(i) = sign(sign_eig(i),eig(i))*h2cm*sqrt(abs(eig(i)))
      ENDDO

!     Fix the phase of the eigenvectors
      call fixphase(hess,qmcoords,ncoords,nqmatoms)


      format1="(9999F9.2)"
      format2="(9999D11.3)"
!     Write eigenvector and eigenvalues to file
      write(77,*) '     CONVERSION FACTOR        '
      write(77,*) '     --------------------     '
      write(77,format2) h2cm
      write(77,*) '     EIGENVALUES AU           '
      write(77,*) '     --------------------     '
      write(77,format2) eig
      write(77,*) '     FREQUENCIES CM-1         '
      write(77,*) '     --------------------     '
      write(77,format1) freq
      write(77,*) '     WRITING EIGENVECTORS     '
      write(77,*) '     --------------------     '
      write(77,*)
      do k=1,nqmatoms
         r=3*(k-1)+1
         s=r+2
         write(77,format1)  (hess(r:s,j),j=7,ncoords)
      end do

      nmodes=hess

      return
      end subroutine
 






!######################################################################
!     FORCE FIELD CALCULATION
!######################################################################
      


      subroutine fullnumqff(nmodes,eig,dy,nqmatoms,nclatoms,qmcoords,clcoords,&
                     &at_numbers,ndf,nvdf,&
                     &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj,&
                     &tijk,uiijk,uijjk,uijkk)
      implicit none

      integer, intent(in) :: ndf                  ! Number of deg of freedm (3*nqmatoms)
      integer, intent(in) :: nvdf                 ! Number of vib degrees of freedom (ndf-6)
      integer, intent(in) :: nqmatoms             ! Number of QM atoms
      integer, intent(in) :: nclatoms             ! Number of classical atoms
      integer, intent(in) :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
      real*8, intent(in)  :: nmodes(ndf,ndf)      ! mw Hessian eigenvectors.
      real*8, intent(in)  :: eig(ndf)             ! mw Hessian eigenvalues.
      real*8, intent(in)  :: dy                   ! Step size factor dy. Default=0.5d0
      real*8, intent(in)  :: qmcoords(3,nqmatoms) ! QM atom coordinates
      real*8, intent(in)  :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
      real*8, intent(out) :: hii(nvdf)           ! Diagonal cubic coupling terms
      real*8, intent(out) :: tiii(nvdf)           ! Diagonal cubic coupling terms
      real*8, intent(out) :: uiiii(nvdf)          ! Diagonal quartic coupling terms
      real*8, intent(out) :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
      real*8, intent(out) :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
      real*8, intent(out) :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling terms
      real*8, intent(out) :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling terms
      real*8, intent(out) :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling terms
      real*8, intent(out) :: tijk(nvdf,nvdf,nvdf)      ! 3-mode cubic coupling terms
      real*8, intent(out) :: uiijk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
      real*8, intent(out) :: uijjk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
      real*8, intent(out) :: uijkk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms

      real*8              :: L(ndf,nvdf)          ! Array of VIBRATIONAL nmodes.
      real*8              :: omega(nvdf)          ! Harmonic frequencies omega(nm)=Sqrt(hii(nm))
      real*8              :: X0(3,nqmatoms)       ! Initial geometry in cartesian and Angstroms.
      real*8              :: Xm(3,nqmatoms)       ! Displaced geometry in cartesian and Angstroms.
      real*8              :: dX(3,nqmatoms)       ! Displacement vector (cart and Angs).
      real*8              :: dX1(3,nqmatoms)      ! Displacement vector (cart and Angs).
      real*8              :: dX2(3,nqmatoms)      ! Displacement vector (cart and Angs).
      real*8              :: dX3(3,nqmatoms)      ! Displacement vector (cart and Angs).
      real*8              :: dQ(nvdf)             ! dQ(mode) = dy/Sqrt(omega_i) in Atomic units.
      real*8              :: Mass(ndf)            ! Sqrt mass weights matrix.
      real*8              :: Minv(ndf)            ! Inverse sqrt mass weights matrix.
      real*8              :: escf                 ! Energy computed by SCF_in()
      real*8              :: dipxyz(3)            ! Dipole moment computed by SCF_in()
      real*8              :: V0                   ! Energy at initial geometry.
      real*8              :: V1(6,nvdf)           ! Energies at 1-mode grid points
      real*8              :: V2(12,nvdf,nvdf)     ! Energies at 2-mode coupling gridpoints
      real*8              :: V3(-1:1,-1:1,-1:1,nvdf,nvdf,nvdf) ! Energies at 3-mode coupling gridpoints
      integer             :: i, j, k
      integer             :: m, nm, p
      integer             :: nm1, nm2, nm3, gp
      integer             :: step_count
      integer             :: dd(6)
      integer             :: grdi(12)
      integer             :: grdj(12)
      integer             :: grti(8)
      integer             :: grtj(8)
      integer             :: grtk(8)
    
      include "qvmbia_param.f"
!     Eliminating rotational and translational degrees of freedom.
      do nm=7,ndf
         L(:,nm-6)=nmodes(:,nm)
         hii(nm-6)=eig(nm)
      end do

!     We define the step size acording to J.Chem.Phys 121:1383
!     Using a dimensionless reduced coordinate y=sqrt(omega_i/hbar)Qi
!     where omega_i is the freq of mode i and Qi is normal coordinate 
!     of mode i. Then the step size dQi is given by
!
!                    dQi = dy/sqrt(omega_i)
!
!     Where omega_i = Sqrt(Ki/mi) in atomic units (hbar=1). 
!     dy is arbitrary and set to 0.5 by default.
!     To obtain the displacement vector in cartesian coordinates dXi
!                  
!                   dXi = M^(-1) * Li * dQi
!
!     Where M^(-1) is the inverse mass weights matrix diag(1/Sqrt(mi))
!     Li is the normal mode i (a column eigenvector of the hessian).

      omega = Sqrt(hii)
      dQ = dy/Sqrt(omega)

!     Energy at initial geometry
      call SCF_in(escf,qmcoords,clcoords,nclatoms,dipxyz)
      V0=escf
      write(77,'(D14.6)') V0

!     Converting units from Angstroms to Bohrs. We work in Bohrs.
!     Convert back to Angstroms before calling SCF_in
      X0=qmcoords/a0 

!     Building mass weights matrix Minv = diag(1/Sqrt(m_i))
      do i=1,nqmatoms
          do j=1,3
              k=at_numbers(i)
              Mass(3*(i-1)+j) = sqrt_atomic_masses_au(k)
              Minv(3*(i-1)+j) = invsqrt_atomic_masses_au(k)
          end do
      end do
!-----------------------------------------------------------------------
!     DIAGONAL QFF TERMS
!-----------------------------------------------------------------------
!     Calculating displaced coordinates and energy at 6 grid points 
!     along each normal mode.
      dd=(/-3,-2,-1,1,2,3/)
      tiii=0.0d0
      uiiii=0.0d0
      step_count=1
      write(*,*) 'ENERGY STEP No  ', step_count
      do nm=1,nvdf
         do j=1,6
            do k=1,nqmatoms
               do m=1,3
                  p = 3*(k-1)+m
                  dX(m,k) = Minv(p)*L(p,nm)*dQ(nm)
                  Xm(m,k) = X0(m,k) + dd(j)*dX(m,k)         ! Moving along nmode
               end do
            end do
            Xm = Xm*a0  ! Convert to Angstroms.
            call SCF_in(escf,Xm,clcoords,nclatoms,dipxyz)   ! Calc energy
            V1(j,nm)=escf
            step_count=step_count+1
            write(*,*) 'ENERGY STEP No  ', step_count
         end do
!        Calculating QFF diagonal terms tiii and uiiii
         tiii(nm)=V1(6,nm)-3.0D0*V1(4,nm)+3.0D0*V1(3,nm)-V1(1,nm)
         tiii(nm)=tiii(nm)*0.125D0/dQ(nm)**3

         uiiii(nm)=V1(5,nm)-4.0D0*V1(4,nm)+6.0D0*V0-4.0D0*V1(3,nm)+V1(2,nm)
         uiiii(nm)=uiiii(nm)/dQ(nm)**4
      end do
      
!-----------------------------------------------------------------------
!     2-MODE COUPLING QFF TERMS
!-----------------------------------------------------------------------
!            1  2  3  4  5  6  7  8  9 10 11 12
      grdi=(/1,-1, 1,-1, 3, 3, 1,-1,-3,-3,-1, 1/)
      grdj=(/1, 1,-1,-1, 1,-1,-3,-3,-1, 1, 3, 3/)         
      V2=0.0D0
!     Calculating energy at stencil 2-d points {grdi(gp),grdj(gp)}
      do nm1=1,nvdf-1
         do nm2=nm1+1,nvdf
            do gp=1,12

               do k=1,nqmatoms
                  do m=1,3
                     p = 3*(k-1)+m
                     dX1(m,k)=Minv(p)*L(p,nm1)*dQ(nm1)
                     dX2(m,k)=Minv(p)*L(p,nm2)*dQ(nm2)
                     Xm(m,k)=X0(m,k) + grdi(gp)*dX1(m,k) + grdj(gp)*dX2(m,k)
                  end do
               end do 

               Xm = Xm*a0 ! Convert to Angstroms.
               call SCF_in(escf,Xm,clcoords,nclatoms,dipxyz)
               V2(gp,nm1,nm2)=escf
               step_count=step_count+1
               write(*,*) 'ENERGY STEP No  ', step_count
            end do 
         end do
      end do
     
!     Calculating QFF 2-mode coupling terms Tiij Tjji Uiiij Ujjji Uiijj 
      tiij=0.0d0
      tjji=0.0d0
      uiiij=0.0d0
      ujjji=0.0d0
      uiijj=0.0d0
      do nm1=1,nvdf-1
         do nm2=nm1+1,nvdf
            tiij(nm1,nm2)=V2(1,nm1,nm2)+V2(2,nm1,nm2)-V2(3,nm1,nm2)-V2(4,nm1,nm2)&
                     &-2.0D0*V1(4,nm2)+2.0D0*V1(3,nm2)
            tiij(nm1,nm2)=tiij(nm1,nm2)*0.5D0/dQ(nm1)**2/dQ(nm2)

            tjji(nm1,nm2)=V2(1,nm1,nm2)+V2(3,nm1,nm2)-V2(2,nm1,nm2)-V2(4,nm1,nm2)&
                     &-2.0D0*V1(4,nm1)+2.0D0*V1(3,nm1)
            tjji(nm1,nm2)=tjji(nm1,nm2)*0.5D0/dQ(nm2)**2/dQ(nm1)

            uiiij(nm1,nm2)= V2(5,nm1,nm2)-3.0D0*V2(1,nm1,nm2)+3.0D0*V2(2,nm1,nm2)-V2(10,nm1,nm2)&
                      &-V2(6,nm1,nm2)+3.0D0*V2(3,nm1,nm2)-3.0D0*V2(4,nm1,nm2)+V2(9,nm1,nm2)
            uiiij(nm1,nm2)=uiiij(nm1,nm2)*0.0625D0/dQ(nm1)**3/dQ(nm2)

            ujjji(nm1,nm2)= V2(12,nm1,nm2)-3.0D0*V2(1,nm1,nm2)+3.0D0*V2(3,nm1,nm2)-V2(7,nm1,nm2) &
                      &-V2(11,nm1,nm2)+3.0D0*V2(2,nm1,nm2)-3.0D0*V2(4,nm1,nm2)+V2(8,nm1,nm2)
            ujjji(nm1,nm2)=ujjji(nm1,nm2)*0.0625D0/dQ(nm2)**3/dQ(nm1)

            uiijj(nm1,nm2)= V2(1,nm1,nm2)+V2(2,nm1,nm2)+V2(3,nm1,nm2)+V2(4,nm1,nm2)+4.0D0*V0 &
                      &-2.0D0*( V1(3,nm1)+V1(4,nm1)+V1(3,nm2)+V1(4,nm2) )
            uiijj(nm1,nm2)=uiijj(nm1,nm2)/dQ(nm1)**2/dQ(nm2)**2
         end do 
      end do
!-----------------------------------------------------------------------
!     3-MODE COUPLING QFF TERMS
!-----------------------------------------------------------------------
!            1  2  3  4  5  6  7  8 
      grti=(/1, 1,-1,-1, 1, 1,-1,-1/)
      grtj=(/1,-1, 1,-1, 1,-1, 1,-1/)
      grtk=(/1,-1,-1, 1,-1, 1, 1,-1/)
!     Calculating energy at stencil 2-d points {grdi(gp),grdj(gp)}
      do nm1=1,nvdf-2
         do nm2=nm1+1,nvdf-1
            do nm3=nm2+1,nvdf
               do gp=1,12
   
                  do k=1,nqmatoms
                     do m=1,3
                        p = 3*(k-1)+m
                        dX1(m,k)=Minv(p)*L(p,nm1)*dQ(nm1)
                        dX2(m,k)=Minv(p)*L(p,nm2)*dQ(nm2)
                        dX3(m,k)=Minv(p)*L(p,nm2)*dQ(nm3)
                        Xm(m,k)=X0(m,k) + grti(gp)*dX1(m,k) + &
                        &grtj(gp)*dX2(m,k) + grtj(gp)*dX3(m,k)
                     end do
                  end do 
   
                  Xm = Xm*a0 ! Convert to Angstroms.
                  call SCF_in(escf,Xm,clcoords,nclatoms,dipxyz)
                  V3(grti(gp),grtj(gp),grtk(gp),nm1,nm2,nm3)=escf
                  step_count=step_count+1
                  write(*,*) 'ENERGY STEP No  ', step_count
               end do
            end do 
         end do
      end do
!     Calculating QFF 2-mode coupling terms Tiij Tjji Uiiij Ujjji Uiijj 
      tijk=0.0d0
      uiijk=0.0d0
      uijjk=0.0d0
      uijkk=0.0d0
      do nm1=1,nvdf-2
         do nm2=nm1+1,nvdf-1
            do nm3=nm2+1,nvdf
               tijk(nm1,nm2,nm3)= &
                 &  V3( 1, 1, 1,nm1,nm2,nm3)+V3( 1,-1,-1,nm1,nm2,nm3) &
                 & +V3(-1, 1,-1,nm1,nm2,nm3)+V3(-1,-1, 1,nm1,nm2,nm3) &
                 & -V3( 1, 1,-1,nm1,nm2,nm3)-V3( 1,-1, 1,nm1,nm2,nm3) &
                 & -V3(-1, 1, 1,nm1,nm2,nm3)-V3(-1,-1,-1,nm1,nm2,nm3)
               tijk(nm1,nm2,nm3)=tijk(nm1,nm2,nm3)/(8d0*dQ(nm1)*dQ(nm2)*dQ(nm3))
   
               uiijk(nm1,nm2,nm3) = &
                &   V3( 1, 1, 1,nm1,nm2,nm3) + V3( 1,-1,-1,nm1,nm2,nm3) &
                & + V3(-1, 1, 1,nm1,nm2,nm3) + V3(-1,-1,-1,nm1,nm2,nm3) &
                & - V3( 1, 1,-1,nm1,nm2,nm3) - V3( 1,-1, 1,nm1,nm2,nm3) &
                & - V3(-1, 1,-1,nm1,nm2,nm3) - V3(-1,-1, 1,nm1,nm2,nm3) &
                & + 2d0*V2(2,nm2,nm3) + 2d0*V2(3,nm2,nm3) &
                & - 2d0*V2(1,nm2,nm3) - 2d0*V2(4,nm2,nm3)
               uiijk(nm1,nm2,nm3)=uiijk(nm1,nm2,nm3)/(4d0*dQ(nm2)*dQ(nm3)*dQ(nm1)**2)
   
               uijjk(nm1,nm2,nm3) = &
                &   V3( 1, 1, 1,nm1,nm2,nm3) - V3( 1,-1,-1,nm1,nm2,nm3) &
                & - V3(-1, 1, 1,nm1,nm2,nm3) + V3(-1,-1,-1,nm1,nm2,nm3) &
                & - V3( 1, 1,-1,nm1,nm2,nm3) + V3( 1,-1, 1,nm1,nm2,nm3) &
                & + V3(-1, 1,-1,nm1,nm2,nm3) - V3(-1,-1, 1,nm1,nm2,nm3) &
                & + 2d0*V2(2,nm1,nm3) + 2d0*V2(3,nm1,nm3) &
                & - 2d0*V2(1,nm1,nm3) - 2d0*V2(4,nm1,nm3)
               uijjk(nm1,nm2,nm3)=uijjk(nm1,nm2,nm3)/(4d0*dQ(nm1)*dQ(nm3)*dQ(nm2)**2)
   
               uijkk(nm1,nm2,nm3) = &
                &   V3( 1, 1, 1,nm1,nm2,nm3) - V3( 1,-1,-1,nm1,nm2,nm3) &
                & - V3(-1, 1, 1,nm1,nm2,nm3) + V3(-1,-1,-1,nm1,nm2,nm3) &
                & + V3( 1, 1,-1,nm1,nm2,nm3) - V3( 1,-1, 1,nm1,nm2,nm3) &
                & - V3(-1, 1,-1,nm1,nm2,nm3) + V3(-1,-1, 1,nm1,nm2,nm3) &
                & + 2d0*V2(2,nm1,nm2) + 2d0*V2(3,nm1,nm2) &
                & - 2d0*V2(1,nm1,nm2) - 2d0*V2(4,nm1,nm2)
               uijkk(nm1,nm2,nm3)=uijjk(nm1,nm2,nm3)/(4d0*dQ(nm1)*dQ(nm2)*dQ(nm3)**2)
            end do
         end do 
      end do

      write(*,*) '----------------------------------------------'
      write(*,*) '   FINISHED WITH QUARTIC FORCE FIELD COMP     '
      write(*,*) '----------------------------------------------'
      write(77,'(A)') '             QFF PARAMETERS'
      write(77,'(A)') '----------------------------------------------'
      write(77,'(A)') '                DIAGONAL'
      write(77,'(A)') '----------------------------------------------'
      do nm1=1,nvdf
         write(77,'(I6,3D16.5)') nm1,hii(nm1),tiii(nm1),uiiii(nm1)
      end do
      write(77,'(A)') '----------------------------------------------'
      write(77,'(A)') '           OFF-DIAGONAL DOUBLES               '
      write(77,'(A)') '----------------------------------------------'
      do nm1=1,nvdf-1
         do nm2=nm1+1,nvdf
            write(77,'(I4,I4,D16.5,D16.5,D16.5,D16.5,D16.5)') nm1,nm2,tiij(nm1,nm2), &
                     &tjji(nm1,nm2), uiiij(nm1,nm2), ujjji(nm1,nm2), uiijj(nm1,nm2)
         end do
      end do
      write(77,'(A)') '----------------------------------------------'
      write(77,'(A)') '           OFF-DIAGONAL TRIPLES               '
      write(77,'(A)') '----------------------------------------------'
      do nm1=1,nvdf-2
         do nm2=nm1+1,nvdf-1
            do nm3=nm2+1,nvdf
               write(77,'(3I4,99D16.5)') nm1,nm2,nm3,tijk(nm1,nm2,nm3), &
                     &uiijk(nm1,nm2,nm3), uijjk(nm1,nm2,nm3), uijkk(nm1,nm2,nm3)
            end do
         end do
      end do

      write(77,'(A)') 'HESSIAN EIGENVALUES'
      write(77,'(9D14.6)') eig
      end subroutine

      subroutine seminumqff(nmodes,eig,dy,nqmatoms,nclatoms,qmcoords,clcoords,&
                     &at_numbers,ndf,nvdf,&
                     &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj, &
                     &tijk,uiijk,uijjk,uijkk)
!     ------------------------------------------------------------------
!     COMPUTES A SEMINUMERICAL QUARTIC FORCE FIELD USING ANALYTICAL 
!     GRADIENT CALLS UP TO 3-MODES COUPLINGS.
!     ------------------------------------------------------------------
      implicit none

      integer, intent(in) :: ndf                  ! Number of deg of freedm (3*nqmatoms)
      integer, intent(in) :: nvdf                 ! Number of vib degrees of freedom (ndf-6)
      integer, intent(in) :: nqmatoms             ! Number of QM atoms
      integer, intent(in) :: nclatoms             ! Number of classical atoms
      integer, intent(in) :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
      real*8, intent(in)  :: nmodes(ndf,ndf)      ! mw Hessian eigenvectors.
      real*8, intent(in)  :: eig(ndf)             ! mw Hessian eigenvalues.
      real*8, intent(in)  :: dy                   ! Step size factor dy. Default=0.5d0
      real*8, intent(in)  :: qmcoords(3,nqmatoms) ! QM atom coordinates
      real*8, intent(in)  :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
      real*8, intent(out) :: hii(nvdf)           ! Diagonal cubic coupling terms
      real*8, intent(out) :: tiii(nvdf)           ! Diagonal cubic coupling terms
      real*8, intent(out) :: uiiii(nvdf)          ! Diagonal quartic coupling terms
      real*8, intent(out) :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
      real*8, intent(out) :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
      real*8, intent(out) :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling terms
      real*8, intent(out) :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling terms
      real*8, intent(out) :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling terms
      real*8, intent(out) :: tijk(nvdf,nvdf,nvdf)      ! 3-mode cubic coupling terms
      real*8, intent(out) :: uiijk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
      real*8, intent(out) :: uijjk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
      real*8, intent(out) :: uijkk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms

      real*8              :: L(ndf,nvdf)          ! Array of VIBRATIONAL nmodes.
      real*8              :: omega(nvdf)          ! Harmonic frequencies omega(nm)=Sqrt(hii(nm))
      real*8              :: dxyzcl(3,nclatoms)   ! SCF MM force
      real*8              :: dxyzqm(3,nqmatoms)   ! SCF QM force
      real*8              :: X0(3,nqmatoms)       ! Initial geometry in cartesian and Angstroms.
      real*8              :: Xm(3,nqmatoms)       ! Displaced geometry in cartesian and Angstroms.
      real*8              :: dX(3,nqmatoms)       ! Displacement vector (cart and Angs).
      real*8              :: dX1(3,nqmatoms)      ! Displacement vector (cart and Angs).
      real*8              :: dX2(3,nqmatoms)      ! Displacement vector (cart and Angs).
      real*8              :: dQ(nvdf)             ! dQ(mode) = dy/Sqrt(omega_i) in Atomic units.
      real*8              :: Mass(ndf)            ! Sqrt mass weights matrix.
      real*8              :: Minv(ndf)            ! Inverse sqrt mass weights matrix.
      real*8              :: gtmp(ndf)
      real*8              :: gnc(nvdf)
      real*8              :: grad0(nvdf)
      real*8              :: grad1(nvdf,4,nvdf)
      real*8              :: grad2(nvdf,4,nvdf,nvdf)
      real*8              :: escf                 ! Energy computed by SCF_in()
      real*8              :: dipxyz(3)            ! Dipole moment computed by SCF_in()
      real*8              :: tmp1
      real*8              :: tmp2
      real*8              :: tmp3
      integer             :: i, j, k
      integer             :: m, nm, p
      integer             :: nm1, nm2, nm3, gp
      integer             :: step_count
      integer             :: dd(4)
      integer             :: grdi(4),grdj(4)
    
      include "qvmbia_param.f"
!     Eliminating rotational and translational degrees of freedom.
      do nm=7,ndf
         L(:,nm-6)=nmodes(:,nm)
         hii(nm-6)=eig(nm)
      end do

!     We define the step size acording to J.Chem.Phys 121:1383
!     Using a dimensionless reduced coordinate y=sqrt(omega_i/hbar)Qi
!     where omega_i is the freq of mode i and Qi is normal coordinate 
!     of mode i. Then the step size dQi is given by
!
!                    dQi = dy/sqrt(omega_i)
!
!     Where omega_i = Sqrt(Ki/mi) in atomic units (hbar=1). 
!     dy is arbitrary and set to 0.5 by default.
!     To obtain the displacement vector in cartesian coordinates dXi
!                  
!                   dXi = M^(-1) * Li * dQi
!
!     Where M^(-1) is the inverse mass weights matrix diag(1/Sqrt(mi))
!     Li is the normal mode i (a column eigenvector of the hessian).
      write(77,'(A)')'------------------------------------------------------------------------' 
      write(77,'(A)')'         STARTING SEMINUMERICAL QFF COMPUTATION USING GRADIENTS         ' 
      write(77,'(A)')'------------------------------------------------------------------------' 

      omega = Sqrt(hii)
      dQ = dy/Sqrt(omega)


!     Building mass weights matrix Minv = diag(1/Sqrt(m_i))
      do i=1,nqmatoms
          do j=1,3
              k=at_numbers(i)
              Mass(3*(i-1)+j) = sqrt_atomic_masses_au(k)
              Minv(3*(i-1)+j) = invsqrt_atomic_masses_au(k)
          end do
      end do


!     Energy at initial geometry
      call SCF_in(escf,qmcoords,clcoords,nclatoms,dipxyz)
      call dft_get_qm_forces(dxyzqm)
      call dft_get_mm_forces(dxyzcl,dxyzqm)

!     Mass weight gradient
      do i=1,nqmatoms
          do j=1,3
              gtmp(3*(i-1)+j)=Minv(3*(i-1)+j)*dxyzqm(j,i)
          end do
      end do

!     Convert to normal coordinates.
      call dgemv('T',ndf,nvdf,1d0,L,ndf,gtmp,1,0d0,gnc,1)
      grad0=gnc



!-----------------------------------------------------------------------
!     DIAGONAL QFF TERMS
!-----------------------------------------------------------------------
!     Converting units from Angstroms to Bohrs. We work in Bohrs.
!     Convert back to Angstroms before calling SCF_in
      X0=qmcoords/a0 

!     Calculating displaced coordinates and energy at 6 grid points 
!     along each normal mode.
      dd=(/-2,-1,1,2/)
      write(77,'(A)') 'LABELS FOR 1d STENCIL POINTS'
      write(77,'(4I3)') (dd(i),i=1,4)
      tiii=0.0d0
      uiiii=0.0d0
      grad1=0d0
      step_count=1
      write(*,*) 'ENERGY STEP No  ', step_count
      do nm=1,nvdf
         do j=1,4
            do k=1,nqmatoms
               do m=1,3
                  p = 3*(k-1)+m
                  dX(m,k) = Minv(p)*L(p,nm)*dQ(nm)   ! Scaling/Mass unweighting
                  Xm(m,k) = X0(m,k) + dd(j)*dX(m,k)         ! Moving along nmode
               end do
            end do
            Xm = Xm*a0  ! Convert to Angstroms.
            call SCF_in(escf,Xm,clcoords,nclatoms,dipxyz)   ! Calc energy
            call dft_get_qm_forces(dxyzqm)
            call dft_get_mm_forces(dxyzcl,dxyzqm)

      !     Mass weight gradient
            do k=1,nqmatoms
                do m=1,3
                    gtmp(3*(k-1)+m)=Minv(3*(k-1)+m)*dxyzqm(m,k)
                end do
            end do
      
      !     Convert to normal coordinates.
            call dgemv('T',ndf,nvdf,1d0,L,ndf,gtmp,1,0d0,gnc,1)
            grad1(:,j,nm)=gnc

            step_count=step_count+1
            write(77,'(A,I6,A,I6)') 'GRAD STEP No ', step_count,'/',1+nvdf*4+2*nvdf*(nvdf-1)
         end do
!        Calculating QFF diagonal terms tiii and uiiii
         tiii(nm)=grad1(nm,4,nm)-2d0*grad0(nm)+grad1(nm,1,nm)
         tiii(nm)=tiii(nm)*0.25D0/dQ(nm)**2

         uiiii(nm)=grad1(nm,4,nm)+2D0*(grad1(nm,2,nm)-grad1(nm,3,nm))-grad1(nm,1,nm)
         uiiii(nm)=uiiii(nm)*0.5d0/dQ(nm)**3
      end do
      
!-----------------------------------------------------------------------
!     2-MODE COUPLING QFF TERMS
!-----------------------------------------------------------------------
!            1  2  3  4
      grdi=(/1,-1, 1,-1/)
      grdj=(/1, 1,-1,-1/)         
      grad2=0.0D0
!     Calculating energy at stencil 2-d points {grdi(gp),grdj(gp)}
      do nm1=1,nvdf-1
         do nm2=nm1+1,nvdf
            do gp=1,4

               do k=1,nqmatoms
                  do m=1,3
                     p = 3*(k-1)+m
                     dX1(m,k)=Minv(p)*L(p,nm1)*dQ(nm1)
                     dX2(m,k)=Minv(p)*L(p,nm2)*dQ(nm2)
                     Xm(m,k)=X0(m,k) + grdi(gp)*dX1(m,k) + grdj(gp)*dX2(m,k)
                  end do
               end do 

               Xm = Xm*a0 ! Convert to Angstroms.
               call SCF_in(escf,Xm,clcoords,nclatoms,dipxyz)
               call dft_get_qm_forces(dxyzqm)
               call dft_get_mm_forces(dxyzcl,dxyzqm)
   
               gtmp=0d0
               gnc=0d0

         !     Mass weight gradient
               do k=1,nqmatoms
                   do m=1,3
                       gtmp(3*(k-1)+m)=Minv(3*(k-1)+m)*dxyzqm(m,k)
                   end do
               end do
         
         !     Convert to normal coordinates.
               call dgemv('T',ndf,nvdf,1d0,L,ndf,gtmp,1,0d0,gnc,1)
               grad2(:,gp,nm1,nm2)=gnc
   
               step_count=step_count+1
               write(77,'(A,I6,A,I6)') 'GRAD STEP No ', step_count,'/',1+nvdf*4+2*nvdf*(nvdf-1)
            end do 
         end do
      end do
     
!     Calculating QFF 2-mode coupling terms Tiij Tjji Uiiij Ujjji Uiijj 
      tiij=0.0d0
      tjji=0.0d0
      uiiij=0.0d0
      ujjji=0.0d0
      uiijj=0.0d0
      do nm1=1,nvdf-1
         do nm2=nm1+1,nvdf
!           KEY
!           dd=(/-2,-1,1,2/)
!           grdi=(/1,-1, 1,-1/)
!           grdj=(/1, 1,-1,-1/)         

            tmp1=0d0
            tmp2=0d0
            tmp1=grad1(nm2,4,nm1)-2d0*grad0(nm2)+grad1(nm2,1,nm1)
            tmp1=tmp1*0.25d0/dQ(nm1)**2
            tmp2=grad2(nm1,1,nm1,nm2)-grad2(nm1,2,nm1,nm2)-grad2(nm1,3,nm1,nm2)+grad2(nm1,4,nm1,nm2)
            tmp2=tmp2/(4d0*dQ(nm1)*dQ(nm2))
            tiij(nm1,nm2)=(tmp1+tmp2)*0.5d0

            tmp1=0d0
            tmp2=0d0
            tmp1=grad1(nm1,4,nm2)-2d0*grad0(nm1)+grad1(nm1,1,nm2)
            tmp1=tmp1*0.25d0/dQ(nm2)**2
            tmp2=grad2(nm2,1,nm1,nm2)-grad2(nm2,2,nm1,nm2)-grad2(nm2,3,nm1,nm2)+grad2(nm2,4,nm1,nm2)
            tmp2=tmp2/(4d0*dQ(nm1)*dQ(nm2))
            tjji(nm1,nm2)=(tmp1+tmp2)*0.5d0

            tmp1=0d0
            tmp2=0d0
            tmp1=grad1(nm2,4,nm1)-2d0*grad1(nm2,3,nm1)+2d0*grad1(nm2,2,nm1)-grad1(nm2,1,nm1)
            tmp1=tmp1*0.5d0/dQ(nm1)**3
            tmp2=grad2(nm1,1,nm1,nm2)-grad2(nm1,3,nm1,nm2)+grad2(nm1,2,nm1,nm2)-grad2(nm1,4,nm1,nm2)
            tmp2=(tmp2-2d0*grad1(nm1,3,nm2)+2d0*grad1(nm1,2,nm2))/(2d0*dQ(nm2)*dQ(nm1)**2)
            uiiij(nm1,nm2)=(tmp1+tmp2)*0.5d0

            tmp1=0d0
            tmp2=0d0
            tmp1=grad1(nm1,4,nm2)-2d0*grad1(nm1,3,nm2)+2d0*grad1(nm1,2,nm2)-grad1(nm1,1,nm2)
            tmp1=tmp1*0.5d0/dQ(nm2)**3
            tmp2=grad2(nm2,1,nm1,nm2)+grad2(nm2,3,nm1,nm2)-grad2(nm2,2,nm1,nm2)-grad2(nm2,4,nm1,nm2)
            tmp2=(tmp2-2d0*grad1(nm2,3,nm1)+2d0*grad1(nm2,2,nm1))/(2d0*dQ(nm1)*dQ(nm2)**2)
            ujjji(nm1,nm2)=(tmp1+tmp2)*0.5d0

            tmp1=0d0
            tmp2=0d0
            tmp1=grad2(nm2,1,nm1,nm2)-grad2(nm2,3,nm1,nm2)+grad2(nm2,2,nm1,nm2)-grad2(nm2,4,nm1,nm2)
            tmp1=(tmp1-2d0*grad1(nm2,3,nm2)+2d0*grad1(nm2,2,nm2))/(2d0*dQ(nm2)*dQ(nm1)**2)
            tmp2=grad2(nm1,1,nm1,nm2)+grad2(nm1,3,nm1,nm2)-grad2(nm1,2,nm1,nm2)-grad2(nm1,4,nm1,nm2)
            tmp2=(tmp2-2d0*grad1(nm1,3,nm1)+2d0*grad1(nm1,2,nm1))/(2d0*dQ(nm1)*dQ(nm2)**2)
            uiijj(nm1,nm2)= (tmp1+tmp2)*0.5d0
         end do 
      end do

!     BEGINING COMPUTATION OF 3-MODE COUPLING TERMS.
      tijk=0d0
      uiijk=0d0
      uijjk=0d0
      uijkk=0d0
      do nm1=1,nvdf-2
         do nm2=nm1+1,nvdf-1
            do nm3=nm2+1,nvdf
!              KEY
!              dd=(/-2,-1,1,2/)
!              grdi=(/1,-1, 1,-1/)
!              grdj=(/1, 1,-1,-1/)         
               tmp1=0d0
               tmp2=0d0
               tmp3=0d0

               tmp1=grad2(nm1,1,nm2,nm3)-grad2(nm1,3,nm2,nm3)
               tmp1=tmp1-grad2(nm1,2,nm2,nm3)+grad2(nm1,4,nm2,nm3)
               tmp1=tmp1/(2d0*dQ(nm2)*dQ(nm3))

               tmp2=grad2(nm2,1,nm1,nm3)-grad2(nm2,3,nm1,nm3)
               tmp2=tmp2-grad2(nm2,2,nm1,nm3)+grad2(nm2,4,nm1,nm3)
               tmp2=tmp2/(2d0*dQ(nm1)*dQ(nm3))

               tmp3=grad2(nm3,1,nm1,nm2)-grad2(nm3,3,nm1,nm2)
               tmp3=tmp3-grad2(nm3,2,nm1,nm2)+grad2(nm3,4,nm1,nm2)
               tmp3=tmp3/(2d0*dQ(nm1)*dQ(nm2))

               tijk(nm1,nm2,nm3)=(tmp1+tmp2+tmp3)/3d0

             
               tmp1=0d0
               tmp2=0d0

               tmp1=grad2(nm3,1,nm1,nm2)-grad2(nm3,3,nm1,nm2)
               tmp1=tmp1+grad2(nm3,2,nm1,nm2)-grad2(nm3,4,nm1,nm2)
               tmp1=tmp1-2d0*grad1(nm3,3,nm2)+2d0*grad1(nm3,2,nm2)
               tmp1=tmp1/(2d0*dQ(nm2)*dQ(nm1)**2)

               tmp2=grad2(nm2,1,nm1,nm3)-grad2(nm2,3,nm1,nm3)
               tmp2=tmp2+grad2(nm2,2,nm1,nm3)-grad2(nm2,4,nm1,nm3)
               tmp2=tmp2-2d0*grad1(nm2,3,nm3)+2d0*grad1(nm2,2,nm3)
               tmp2=tmp2/(2d0*dQ(nm3)*dQ(nm1)**2)

               uiijk(nm1,nm2,nm3)=(tmp1+tmp2)*0.5d0


               tmp1=0d0
               tmp2=0d0

               tmp1=grad2(nm1,1,nm2,nm3)-grad2(nm1,3,nm2,nm3)
               tmp1=tmp1+grad2(nm1,2,nm2,nm3)-grad2(nm1,4,nm2,nm3)
               tmp1=tmp1-2d0*grad1(nm1,3,nm3)+2d0*grad1(nm1,2,nm3)
               tmp1=tmp1/(2d0*dQ(nm3)*dQ(nm2)**2)

               tmp2=grad2(nm3,1,nm1,nm2)+grad2(nm3,3,nm1,nm2)
               tmp2=tmp2-grad2(nm3,2,nm1,nm2)-grad2(nm3,4,nm1,nm2)
               tmp2=tmp2-2d0*grad1(nm3,3,nm1)+2d0*grad1(nm3,2,nm1)
               tmp2=tmp2/(2d0*dQ(nm1)*dQ(nm2)**2)

               uijjk(nm1,nm2,nm3)=(tmp1+tmp2)*0.5d0


               tmp1=0d0
               tmp2=0d0

               tmp1=grad2(nm2,1,nm1,nm3)+grad2(nm2,3,nm1,nm3)
               tmp1=tmp1-grad2(nm2,2,nm1,nm3)-grad2(nm2,4,nm1,nm3)
               tmp1=tmp1-2d0*grad1(nm2,3,nm1)+2d0*grad1(nm2,2,nm1)
               tmp1=tmp1/(2d0*dQ(nm1)*dQ(nm3)**2)

               tmp2=grad2(nm1,1,nm2,nm3)+grad2(nm1,3,nm2,nm3)
               tmp2=tmp2-grad2(nm1,2,nm2,nm3)-grad2(nm1,4,nm2,nm3)
               tmp2=tmp2-2d0*grad1(nm1,3,nm2)+2d0*grad1(nm1,2,nm2)
               tmp2=tmp2/(2d0*dQ(nm2)*dQ(nm3)**2)

               uijkk(nm1,nm2,nm3)=(tmp1+tmp2)*0.5d0

            end do
         end do
      end do   
      write(*,*) '----------------------------------------------'
      write(*,*) '   FINISHED WITH QUARTIC FORCE FIELD COMP     '
      write(*,*) '----------------------------------------------'
      write(77,'(A)') '             QFF PARAMETERS'
      write(77,'(A)') '----------------------------------------------'
      write(77,'(A)') '                DIAGONAL'
      write(77,'(A)') '----------------------------------------------'
      do nm1=1,nvdf
         write(77,'(I6,3D16.5)') nm1,hii(nm1),tiii(nm1),uiiii(nm1)
      end do
      write(77,'(A)') '----------------------------------------------'
      write(77,'(A)') '           OFF-DIAGONAL DOUBLES               '
      write(77,'(A)') '----------------------------------------------'
      do nm1=1,nvdf-1
         do nm2=nm1+1,nvdf
            write(77,'(I4,I4,D16.5,D16.5,D16.5,D16.5,D16.5)') nm1,nm2,tiij(nm1,nm2), &
                     &tjji(nm1,nm2), uiiij(nm1,nm2), ujjji(nm1,nm2), uiijj(nm1,nm2)
         end do
      end do
      write(77,'(A)') '----------------------------------------------'
      write(77,'(A)') '           OFF-DIAGONAL TRIPLES               '
      write(77,'(A)') '----------------------------------------------'
      do nm1=1,nvdf-2
         do nm2=nm1+1,nvdf-1
            do nm3=nm2+1,nvdf
               write(77,'(3I4,99D16.5)') nm1,nm2,nm3,tijk(nm1,nm2,nm3), &
                     &uiijk(nm1,nm2,nm3), uijjk(nm1,nm2,nm3), uijkk(nm1,nm2,nm3)
            end do
         end do
      end do
      write(77,'(A)') 'HESSIAN EIGENVALUES'
      write(77,'(9D14.6)') eig
      end subroutine


end module qvamod_lioexcl

 