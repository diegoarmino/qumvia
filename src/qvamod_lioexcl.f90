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
   public :: fullnumqff,seminumqff,hessian,lio_nml_type,get_lio_nml,&
           & print_lio_nml, read_coords_lio


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

     subroutine print_lio_nml(lio_nml)
        implicit none
        type(lio_nml_type), intent(in) :: lio_nml

        write(77,'(A)') ' LIO NAMELIST PARAMETERS '
        write(77,'(A)') ' ----------------------- '
        write(77,*)      '  OPEN ', lio_nml%OPEN
        write(77,*)     '  NMAX ', lio_nml%NMAX
        write(77,*)     '  NUNP ', lio_nml%NUNP
        write(77,*)     '  VCINP ', lio_nml%VCINP
!       write(77,*)     '  GOLD ', GOLD
        write(77,*)     '  told ', lio_nml%told
!       write(77,*)     '  rmax ', rmax
!       write(77,*)     '  rmaxs ', rmaxs
!       write(77,*)     '  predcoef ', predcoef
!       write(77,*)     '  idip ', idip
        write(77,*)     '  writexyz ', lio_nml%writexyz
        write(77,*)     '  DIIS ', lio_nml%DIIS
        write(77,*)     '  ndiis ', lio_nml%ndiis
        write(77,*)     '  Iexch ', lio_nml%Iexch
!       write(77,*)     '  integ ', integ
!       write(77,*)     '  DENS ' ,  DENS
        write(77,*)     '  IGRID ', lio_nml%IGRID
        write(77,*)     '  IGRID2 ', lio_nml%IGRID2
        write(77,*)     '  timedep ', lio_nml%timedep
        if(lio_nml%timedep.gt.0) then
           write(77,*)     'Runing TDDFT calculation with Lio'
           write(77,*)     '  tdstep ', lio_nml%tdstep
           write(77,*)     '  ntdstep ', lio_nml%ntdstep
           write(77,*)     '  field ', lio_nml%field
           write(77,*)     '  exter ', lio_nml%exter
           write(77,*)     '  a0 ', lio_nml%a0
           write(77,*)     '  epsilon ', lio_nml%epsilon
           write(77,*)     '  Fx ', lio_nml%Fx
           write(77,*)     '  Fy ', lio_nml%Fy
           write(77,*)     '  Fz ', lio_nml%Fz
           write(77,*)     '  NBCH ', lio_nml%NBCH
           write(77,*)     '  propagator ', lio_nml%propagator
           write(77,*)     '  writedens ', lio_nml%writedens
           write(77,*)     '  writedens ', lio_nml%tdrestart
        endif

     end subroutine print_lio_nml


     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
     !%% READ_COORDS  Reads atoms' coordinates from an input file for LIO           !
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
     subroutine read_coords_lio(inputCoord)

         use garcha_mod, only : natom, ntatom, nsol, iz, r, rqm, pc

         implicit none

         character(len=20), intent(in) :: inputCoord

         integer :: ios,i
         logical :: fileExists

         inquire(file=inputCoord,exist=fileExists)
         if(fileExists) then
             open(unit=101,file=inputCoord,iostat=ios)
         else
             write(*,*) 'Input coordinates file ',adjustl(inputCoord),' not found.'
             stop
         endif

         ! Reads coordinates file.
         ntatom = natom + nsol
         allocate (iz(natom), r(ntatom,3), rqm(natom,3), pc(ntatom))
         read(101,*)
         do i=1,natom
             read(101,*) iz(i), r(i,1:3)
             rqm(i,1:3) = r(i,1:3)
         enddo
         do i=natom+1,ntatom
             read(101,*) pc(i), r(i,1:3)
         enddo
         r  = r   / 0.529177D0
         rqm= rqm / 0.529177D0

         return
     end subroutine read_coords_lio

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

     subroutine eckart(hess,qmcoords,at_numbers,ndf,nqmatoms,hessp)
     implicit none
!    ------------------------------------------------------------------
     integer,intent(in)  :: ndf
     integer,intent(in)  :: nqmatoms
     integer,intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
     real*8,intent(in)   :: hess(ndf,ndf)
     real*8,intent(in)   :: qmcoords(3,nqmatoms)
     real*8,intent(out)  :: hessp(ndf,ndf)
!    Internal variables
     integer             :: i,j,k,iat
     integer             :: zeros(2),nonz(2)
     integer             :: ipiv(3)
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
     real*8              :: TRmat(ndf,6)
     real*8              :: unitv(3,3)
     real*8              :: SV(6)

     real*8,external  :: ddot

     real*8,  dimension(:,:), allocatable :: VT,U

     real*8,  dimension(:), allocatable :: WORK,IWORK
     integer                            :: LWORK,INFO
     include "qvmbia_param.f"
!    ------------------------------------------------------------------

     unitv=reshape((/ 1d0, 0d0, 0d0, &
                  &   0d0, 1d0, 0d0, &
                  &   0d0, 0d0, 1d0 /),(/3,3/))

!    Converting initial coordinates to AU
     X0=qmcoords/a0
!     write(77,'(A)') 'CONVERSION FACTOR a0 -> Angs'
!     write(77,'(D12.3)') a0

!    Computing total mass and atomic mass array and center of mass.
     totM=0.0d0
     do i=1,nqmatoms
         do j=1,3
             k=at_numbers(i)
             Mass(3*(i-1)+j) = SQRT(isomass((k-1)*4+1))
             cmass(j)=cmass(j)+isomass((k-1)*4+1)*X0(j,i)
!             version with atomic masses in AU
!             Mass(3*(i-1)+j) = sqrt_atomic_masses_au(k)
!             cmass(j)=cmass(j)+atomic_masses_au(k)*X0(j,i)
         end do
         totM=totM+atomic_masses_au(k)
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


!    WE NOW COMPUTE THE PROJECTOR MATRIX P
!    THE FORMULA CAN BE FOUND IN
!    Szalay, JCP 140:234107, (2014), eq 9, and eq. 17
!    Note that although eq. 12 is only an orthonormalized version of
!    eq. 9. Also, we only implemented Eckart conditions. Sayvetz
!    conditions are left for another time.

!    Also, the following paper may be of interest.
!    Miller, Handy and Adams, JCP 72:99,(1980) eq 4.11

!    The imposition of Eckart conditions is made by projection.

!                  PHess = P.Hess.P

!    Initializing TRmat
     Trmat=0d0

!    Computing translation and rotation vectors.
     do i=1,3
        do iat=1,nqmatoms
           TRmat(3*(iat-1)+i,i)=Mass(3*(iat-1)+i)/Sqrt(TotM)  ! Translation
           TRmat(3*(iat-1)+i,4:6) = crossp(X0(:,iat),unitv(i,:)) ! Rotation
        end do
     end do

!    Orthonormalization of translation and rotation vectors by singular value decomposition method.
     LWORK=-1
     allocate(VT(6,6),U(0,0),WORK(1000),IWORK(8*6))
     call dgesdd('O',ndf,6,TRmat,ndf,SV,U,ndf,VT,6,WORK,LWORK,IWORK,INFO)
     LWORK=WORK(1)
     deallocate(WORK)
     allocate(WORK(LWORK))
     call dgesdd('O',ndf,6,TRmat,ndf,SV,U,ndf,VT,6,WORK,LWORK,IWORK,INFO)
     deallocate(WORK,IWORK,U,VT)
     if (INFO /= 0) STOP ('Error during singular value decomposition in eckart subroutine')

!    Building projection matrix by external product of TRmat.
     call dgemm('n','t',ndf,ndf,6,1d0,TRmat,ndf,TRmat,ndf,0d0,P,ndf)

     P=-P
     do i=1,ndf
        P(i,i) = 1.0d0 + P(i,i)
     end do

!     cutoff = 1.0d-08
!     do i=1,ndf-1
!        do j=i+1,ndf
!           if (abs(P(i,j)) < cutoff) P(i,j)=0d0
!           P(j,i)=P(i,j)
!        end do
!     end do
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

     function crossp(v1,v2)
!       ---------------------------------------------------------------------
!       Computes cross product v1 x v2 = vr, with all vectors dimension 3
!       ---------------------------------------------------------------------
        real*8  :: v1(3), v2(3)
        integer :: i

        real*8,dimension(3) :: crossp

        crossp(1) = v1(2)*v2(3)-v1(3)*v2(2)
        crossp(2) = v1(3)*v2(1)-v1(1)*v2(3)
        crossp(3) = v1(1)*v2(2)-v1(2)*v2(1)

        return
     end function crossp


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


     subroutine hessian(qmcoords,nqmatoms,h,hess_norder,at_numbers,nmodes,eig)

      implicit none

!     -------------------------------------------------------------------------
      real*8,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      real*8,  intent(in) :: h                    ! Finite differences displacement factor.
      integer, intent(in) :: hess_norder          ! Finite differences displacement factor.
      integer, intent(in) :: nqmatoms             ! Number of QM atoms
      integer, intent(in) :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
      real*8,  intent(out):: nmodes(3*nqmatoms,3*nqmatoms)  ! Hessian eigenvectors
      real*8,  intent(out):: eig(3*nqmatoms)      ! Hessian eigenvalues
!     LOCAL VARIABLES
      integer,parameter   :: nclatoms=0
      integer             :: ncoords              ! 3*nqmatoms
      integer             :: i,j,k,r,s,m
      real*8              :: dxyzcl(3,nclatoms)   ! SCF MM force
      real*8              :: dxyzqm(3,nqmatoms)   ! SCF QM force
      real*8              :: hau                  ! h in bohrs

      real*8              :: escf                 ! SCF energy
      real*8              :: tmp1,tmp2            ! Auxiliary variables for hessian calc.
      real*8              :: hhi,hhj
      real*8              :: dipxyz(3)            ! Dipole moment
      real*8              :: grad0(3*nqmatoms)
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
      real*8,  dimension(:), allocatable    :: WORK, WORK2
      real*8,  dimension(:,:,:),allocatable :: grad
      integer, dimension(:), allocatable    :: IWORK, IWORK2
      integer :: LWORK, LIWORK, INFO

!     PARAMETERS
      real*8, PARAMETER :: h2cm=219475.64d0  ! Convert hartee/bohr^2/emu to cm-1
      real*8, PARAMETER :: freq2cm=5141.1182009496200901D0 ! Convert freq to cm-1

      include "qvmbia_param.f"
!     -------------------------------------------------------------------------
!
      hau=h/a0                         ! step for num differentiation in bohrs
      ncoords=3*nqmatoms

      dxyzqm=0.0D0
      dxyzcl=0.0D0
      hess=0.0D0
      eig=0.0D0
      freq=0.0D0
      sign_eig=1.0D00
      at_masses=0.0D0

!     Create a vector of inverse square-root masses corresponding to each
!     cartesian coordinate.
      do i=1,nqmatoms
         do j=1,3
            k=at_numbers(i)
            at_masses(3*(i-1)+j) = isomass((k-1)*4+1)
         end do
      end do

      if (hess_norder==1) then
         allocate (grad(3*nqmatoms,-1:1,3*nqmatoms))
      else
         allocate (grad(3*nqmatoms,-2:2,3*nqmatoms))
      end if
      grad=0.0D0

!     Calculating energy and gradient at the initial geometry.
 !     PRINT THE GEOMETRY JUST READ
       write(77,'(A)') 'INPUT GEOMETRY'
       write(77,'(A)') '-----------------------------------------------------------------'
       do i=1,nqmatoms
          write(77,'(I5,3D20.11)') at_numbers(i),qmcoords(:,i)
       end do
       write(77,'(A)') '-----------------------------------------------------------------'
       write(77,*)

      call SCF_in(escf,qmcoords,clcoords,nclatoms,dipxyz)
      call dft_get_qm_forces(dxyzqm)
      call dft_get_mm_forces(dxyzcl,dxyzqm)

      do i=1,nqmatoms
         do j=1,3
            grad0(3*(i-1)+j)=dxyzqm(j,i)
         end do
      end do

!     Calculating gradients in displaced positions along the
!     QM coordinates.

      do i=1,nqmatoms
         do j=1,3

            k=3*(i-1)+j

!           NOTE:
!           qmxyz must be in Angstroms. It will be later transfromed into AU in SCF_in.
!           The displacement factor h has units of Angs*amu^1/2, so in order to use it
!           for simple cartesians in Angs we must divide by sqrt(mi), with mi the mass
!           of the atom being displaced. In this way heavy atoms are less displaced than
!           light atoms, and the numerical differentiation should be better balanced.

!           Forward displacement
!           --------------------
            qmxyz=qmcoords
            qmxyz(j,i)=qmcoords(j,i)+h/sqrt(at_masses(k))
            call SCF_in(escf,qmxyz,clcoords,nclatoms,dipxyz)
            call dft_get_qm_forces(dxyzqm)
            call dft_get_mm_forces(dxyzcl,dxyzqm)


            do r=1,nqmatoms
             do s=1,3
               grad(k,1,3*(r-1)+s)=dxyzqm(s,r)
             end do
            end do

!           Backward displacement
!           ---------------------
            qmxyz=qmcoords
            qmxyz(j,i)=qmcoords(j,i)-h/sqrt(at_masses(k))
            call SCF_in(escf,qmxyz,clcoords,nclatoms,dipxyz)
            call dft_get_qm_forces(dxyzqm)
            call dft_get_mm_forces(dxyzcl,dxyzqm)

            do r=1,nqmatoms
             do s=1,3
               grad(k,-1,3*(r-1)+s)=dxyzqm(s,r)
             end do
            end do

            if (hess_norder==2) then

!              Forward 2x displacement
!              -----------------------
               qmxyz=qmcoords
               qmxyz(j,i)=qmcoords(j,i)+2d0*h/sqrt(at_masses(k))
               call SCF_in(escf,qmxyz,clcoords,nclatoms,dipxyz)
               call dft_get_qm_forces(dxyzqm)
               call dft_get_mm_forces(dxyzcl,dxyzqm)
               do r=1,nqmatoms
                do s=1,3
                  grad(k,2,3*(r-1)+s)=dxyzqm(s,r)
                end do
               end do

!              Backward 2x displacement
!              ------------------------
               qmxyz=qmcoords
               qmxyz(j,i)=qmcoords(j,i)-2d0*h/sqrt(at_masses(k))
               call SCF_in(escf,qmxyz,clcoords,nclatoms,dipxyz)
               call dft_get_qm_forces(dxyzqm)
               call dft_get_mm_forces(dxyzcl,dxyzqm)
               do r=1,nqmatoms
                do s=1,3
                  grad(k,-2,3*(r-1)+s)=dxyzqm(s,r)
                end do
               end do

            end if

         end do
      end do
      write(77,'(A)') 'FINISHED ALL GRADIENTS COMPUTATIONS, NOW COMPUTING HESSIAN'


!     Calculating elements of the Hessian and mass weighting.
      IF (hess_norder==1) THEN

         do i=1,ncoords
             do j=i,ncoords
                tmp1=0d0
                tmp2=0d0
                hhi=hau/sqrt(at_masses(i))
                hhj=hau/sqrt(at_masses(j))
                if (i==j) then
                   tmp1=(grad(i,1,j)-grad(i,-1,j))*0.5d0/hhi    ! Finite difference.
                   hess(i,j)=tmp1/sqrt(at_masses(i)*at_masses(j))          ! Mass weighting hessian
                else
                   tmp1=(grad(i,1,j)-grad(i,-1,j))*0.5d0/hhi
                   tmp2=(grad(j,1,i)-grad(j,-1,i))*0.5d0/hhj
                   hess(i,j)=(tmp1+tmp2)*0.5d0/sqrt(at_masses(i)*at_masses(j))
                   hess(j,i)=hess(i,j)
                end if
             end do
         end do

      ELSE

         do i=1,ncoords
             do j=i,ncoords
                tmp1=0d0
                tmp2=0d0
                hhi=hau/sqrt(at_masses(i))
                hhj=hau/sqrt(at_masses(j))
                if (i==j) then
                   tmp1=(grad(i,-2,j)+8d0*grad(i,1,j)-8d0*grad(i,-1,j)-grad(i,2,j))/(12d0*hhi)    ! Finite difference.
                   hess(i,j)=tmp1/sqrt(at_masses(i)*at_masses(j))          ! Mass weighting hessian
                else
                   tmp1=(grad(i,-2,j)+8d0*grad(i,1,j)-8d0*grad(i,-1,j)-grad(i,2,j))/(12d0*hhi)    ! Finite difference.
                   tmp2=(grad(j,-2,i)+8d0*grad(j,1,i)-8d0*grad(j,-1,i)-grad(j,2,i))/(12d0*hhj)    ! Finite difference.
                   hess(i,j)=(tmp1+tmp2)*0.5d0/sqrt(at_masses(i)*at_masses(j))
                   hess(j,i)=hess(i,j)
                end if
             end do
         end do

      END IF

!     Projecting out translational and rotaional modes from hessian (Eckart conditions).
      call eckart(hess,qmcoords,at_numbers,3*nqmatoms,nqmatoms,hessp)

!     Write hessian to file
!      write(77,*), '     WRITING HESSIAN MATRIX     '
!      write(77,*), '     ----------------------     '
!      do i=1,ncoords
!          write(77,'(I3,999D11.3)'), i, ( hess(i,j), j=1,ncoords )
!      end do

!     Scaling hessian to avoid numerical problems.
      hess=hess*1000d0


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


!     Descaling eigenvalues.
      eig=eig/1000d0

!     Convert frequecies to wavenumbers
      do i = 1,ncoords
         freq(i) = sign(sign_eig(i),eig(i))*freq2cm*sqrt(abs(eig(i)))
      ENDDO

!     Fix the phase of the eigenvectors
      call fixphase(hess,qmcoords,ncoords,nqmatoms)


      format1="(9999F9.2)"
      format2="(9999D11.3)"
!     Write eigenvector and eigenvalues to file
      write(77,*)
      write(77,*) '     VIBRATIONAL FREQUENCIES WITHOUT ECKART CONDITIONS'
      write(77,*) '     -------------------------------------------------'
      write(77,*)
      write(77,format1) freq
      write(77,*)
      write(77,*)

      hess = hessp

!     Scaling hessian to avoid numerical problems.
      hess=hess*1000D0

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

!     Descaling eigenvalues.
      eig=eig/1000d0

!     Convert frequecies to wavenumbers
      do i = 1,ncoords
         freq(i) = sign(sign_eig(i),eig(i))*freq2cm*sqrt(abs(eig(i)))
      ENDDO

!     Fix the phase of the eigenvectors
      call fixphase(hess,qmcoords,ncoords,nqmatoms)


      format1="(9999F9.2)"
      format2="(9999D11.3)"
!     Write eigenvector and eigenvalues to file
      write(77,*) '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*) '     VIBRATIONAL RESULTS WITH ECKART CONDITIONS'
      write(77,*) '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*)
      write(77,*) ' LIST OF FREQUENCIES (CM-1)'
      write(77,format1) freq
      write(77,*)
      write(77,*)             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*)             '           WRITING NORMAL MODES            '
      write(77,*)             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*)
      do j=7,ncoords
         write(77,'(A,i3)')   'MODE =',j-6
         write(77,'(A,F9.2)') 'FREQUENCY (CM-1) =',freq(j)
         write(77,'(A)')      '-------------------------------------------'
         write(77,'(A)')      ' AT       X            Y            Z      '
         write(77,'(A)')      '-------------------------------------------'
         do k=1,nqmatoms
            r=3*(k-1)+1
            s=r+2
            write(77,'(i3,3F13.8)')  k, hess(r:s,j)
         end do
         write(77,*)
         write(77,*)
      end do

      nmodes=hess
!     Converting eigenvectors form Hartree/bohr^2/amu to Hartree/bohr^2/emu
      eig=eig/AMU_TO_AU

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
