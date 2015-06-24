 module qvamod_common

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
    public :: geoms4qff,get_qva_nml,readnqmatoms,qva_nml_type, &
           &  readgeom, readgaunmodes, readqff, readgamessqff, &
           &  hseminumqff, readaddref, genconf3, &
           &  ssvscf2, csVCI2,qva_cli_type
 
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
      real*8  :: ethresh
      real*8  :: resthresh
      real*8  :: selcut1
      real*8  :: selcut2
      real*8  :: qva_dstep
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

 
 contains
 
       subroutine get_qva_nml(qvain,qva_nml)
          implicit none
       
          character(99),intent(in)        :: qvain
          type(qva_nml_type), intent(out) :: qva_nml
       
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
          real*8  :: ethresh
          real*8  :: resthresh
          real*8  :: selcut1
          real*8  :: selcut2
          real*8  :: qva_dstep
       
          namelist /qva/ nhess,vscf_gauswidth, &
          vci_qmax1,vci_qmax2,qumvia_qff,qumvia_nmc,vci_qmax3,ethresh,&
          resthresh,selcut1,selcut2,vci_qmax4,qva_naddref,qva_dstep,qva_extprog
       
          integer :: ifind, ierr
       
       !  DEFAULT VALUES FOR QUMVIA
          nhess=0
          vscf_gauswidth=0.5
          vci_qmax1=5
          vci_qmax2=3
          vci_qmax3=3
          vci_qmax4=3
          qva_naddref=0
          qumvia_qff=2
          qumvia_nmc=3
          ethresh=15000
          resthresh=3000
          selcut1=1d-6
          selcut2=1d-6
          qva_dstep=0.5d0
          qva_extprog=2
       
       !  READ NAMELIST
       !  DANGER: Is it necessary to open file? I guess so.
          write(*,*) 'OPENING INPUT FILE'
          open(UNIT=10,FILE=qvain,action='READ',iostat=ierr)
          if (ierr /= 0) then
             write(*,*) 'ERROR OPENNING INPUT FILE'
             STOP
          end if
          rewind 10
          read(10,nml=qva,iostat=ierr)
       
          if ( ierr > 0 ) then
             STOP('ERROR READING INPUT FILE')
          end if
       
          qva_nml%nhess=nhess
          qva_nml%vscf_gauswidth=vscf_gauswidth
          qva_nml%vci_qmax1=vci_qmax1
          qva_nml%vci_qmax2=vci_qmax2
          qva_nml%qumvia_qff=qumvia_qff
          qva_nml%qumvia_nmc=qumvia_nmc
          qva_nml%qva_extprog=qva_extprog
          qva_nml%vci_qmax3=vci_qmax3
          qva_nml%vci_qmax4=vci_qmax4
          qva_nml%qva_naddref=qva_naddref
          qva_nml%ethresh=ethresh
          qva_nml%resthresh=resthresh
          qva_nml%selcut1=selcut1
          qva_nml%selcut2=selcut2
          qva_nml%qva_dstep=qva_dstep
       
       end subroutine get_qva_nml
 
 
       subroutine print_namelist(qva_nml)
           implicit none
           type(qva_nml_type), intent(in) :: qva_nml
 
       end subroutine print_namelist
 
 !######################################################################
 !     READ NQMATOMS FROM GEOMETRY FILE
 !######################################################################
 
       subroutine readnqmatoms(qvageom,nqmatoms)
 
       implicit none
 !     -----------------------------------------------------------
 !     READS THE NQMATOMS VARIABLE FROM GEOMETRY INPUT FILE 
 !     -----------------------------------------------------------
       character(99),intent(in)   :: qvageom
       integer,intent(out)        :: nqmatoms
 !     Local variables
       integer   :: openstatus
       integer   :: closestatus
 !     -----------------------------------------------------------
 
       open(UNIT=15, FILE=qvageom, ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A,A)') 'COULD NOT OPEN GEOMETRY FILE ',qvageom
          STOP
       end if
 
 !     READING GEOMETRY FROM FILE
       READ(15,*) nqmatoms
       write(77,'(A,I3)') 'nqmatoms=',nqmatoms
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A,A)') 'ERROR: COULD NOT CLOSE FILE',qvageom
          STOP
       end if
 
       end subroutine
 
 !######################################################################
 !    READ GEOMETRY FILE
 !######################################################################
 
       subroutine readgeom(qva_cli,nqmatoms,qvageom,at)
 
       implicit none
 !     -----------------------------------------------------------
 !     READS GEOMETRY INPUT FILE IN QUMVIA FORMAT (SIMILAR TO XYZ)
 
 !     FORMAT:
 !     ~~~~~~~~~~~~~~
 !     Nat
 !     AtN(1) x  y  z
 !     AtN(2) x  y  z
 !     ......
 !     AtN(N) x  y  z
 !     ~~~~~~~~~~~~~~
 
 !     WHERE AtN(i) IS THE ATOMIC NUMBER OF ATOM i, AND N IS THE
 !     NUMBER OF ATOMS IN THE MOLECULE. FOR NOW THE FORMAT IS 
 !     VERY THIGHTLY FIXED AS I3 FOR AtN(i) AND F15.10 FOR X, Y 
 !     AND Z. THIS WILL PROBABLY CHANGE IN THE FUTURE.
 !     -----------------------------------------------------------
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in)             :: nqmatoms
       real*8,intent(out)             :: qvageom(3,nqmatoms)
       integer,intent(out)            :: at(nqmatoms)
 !     Local variables
       integer   :: n
       integer   :: i
       integer   :: nat
       integer   :: openstatus
       integer   :: closestatus
       real*8    :: x,y,z
 !     -----------------------------------------------------------
 
       open(UNIT=15, FILE=qva_cli%geo, ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN GEOMETRY FILE' 
          STOP
       end if
 
 !     READING GEOMETRY FROM FILE
       READ(15,*) nat
       write(77,'(A,I3)') 'NAT=',nat
       at=0
       qvageom=0d0
 
       DO n=1,nqmatoms
          READ(15,'(I3,3F15.10)') at(n),x,y,z
          qvageom(1,n)=x
          qvageom(2,n)=y
          qvageom(3,n)=z
          write(77,'(A,I3)') 'READING LINE ',n
       END DO
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE'
          STOP
       end if
 
 !     PRINT THE GEOMETRY JUST READ
       write(77,'(A)') 'INPUT GEOMETRY'
       do i=1,nqmatoms
          write(77,'(I6,3D18.11)') at(i),qvageom(:,i)
       end do
 
       end subroutine
 
       subroutine geoms4qff(qva_cli,qva_nml,nqmatoms)
 !     ------------------------------------------------------------------
 !     GENERATES GEOMETRIES FOR SEMINUMERICAL QFF USING HESSIANS
 !     TO BE COMPUTED USING AN EXTERNAL ELECTRONIC STRUCTURE SOFTWARE.
 !     ------------------------------------------------------------------
       implicit none
 
       type(qva_nml_type), intent(in) :: qva_nml
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in) :: nqmatoms
 
       integer             :: ngaus
       integer             :: ndf                  ! Number of deg of freedm (3*nqmatoms)
       integer             :: nvdf                 ! Number of vib degrees of freedom (ndf-6)
       integer             :: qva_extprog          ! External program gam=1,gaus=2
       integer             :: qumvia_nmc
       integer             :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
       real*8              :: dy                   ! Step size factor dy. Default=0.5d0
       real*8              :: qvageom(3,nqmatoms) ! QM atom coordinates
 
       real*8              :: nmodes(3*nqmatoms,3*nqmatoms)      ! mw Hessian eigenvectors.
       real*8              :: freq(3*nqmatoms)            ! Harmonic frequencies
       real*8              :: eig(3*nqmatoms)             ! Harmonic force constants.
       real*8              :: atmass(nqmatoms)     ! Atomic masses
       real*8              :: L(3*nqmatoms,3*nqmatoms-6)          ! Array of VIBRATIONAL nmodes.
       real*8              :: hii(3*nqmatoms-6)            ! Harmonic force constants.
       real*8              :: omega(3*nqmatoms-6)          ! Harmonic frequencies omega(nm)=Sqrt(hii(nm))
       real*8              :: X0(3,nqmatoms)       ! Initial geometry in cartesian and Angstroms.
       real*8              :: Xm(3,nqmatoms)       ! Displaced geometry in cartesian and Angstroms.
       real*8              :: dX(3,nqmatoms)       ! Displacement vector (cart and Angs).
       real*8              :: dX1(3,nqmatoms)      ! Displacement vector (cart and Angs).
       real*8              :: dQ(3*nqmatoms-6)             ! dQ(mode) = dy/Sqrt(omega_i) in Atomic units.
       real*8              :: Mass(3*nqmatoms)            ! Sqrt mass weights matrix.
       real*8              :: Minv(3*nqmatoms)            ! Inverse sqrt mass weights matrix.
       integer             :: nat
       integer             :: an
       integer             :: i, j, k
       integer             :: m, nm, p
       integer             :: step
       integer             :: dd(2)
       integer             :: openstatus,closestatus
       character(len=2),dimension(36),parameter :: atsym=(/"H ","He",&
       & "Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg","Al","Si",&
       & "P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn","Fe",&
       & "Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"/)
     
       include "qvmbia_param.f"
 
 !     Read geometry file.
       call readgeom(qva_cli,nat,qvageom,at_numbers)
 
 !     Some aliases
       qumvia_nmc = qva_nml%qumvia_nmc
       qva_extprog = qva_nml%qva_extprog
       dy=qva_nml%qva_dstep
 
       ngaus=16
       ndf=3*nqmatoms
       nvdf=ndf-6
 
 !     Read normal modes, frequencies and atomic masses from input file.
       if (qva_extprog == 1) then
          call readgamnmodes(qva_cli,nqmatoms,ndf,nmodes,freq,eig,atmass)
 !        Eliminating rotational and translational degrees of freedom.
          do nm=7,ndf
             L(:,nm-6)=nmodes(:,nm)
             hii(nm-6)=eig(nm)
          end do
       else if (qva_extprog == 2) then
          call readgaunmodes(qva_cli,nqmatoms,ndf,nvdf,L,freq,hii,atmass)
       end if
 
 
 !     -----------------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------------
       write(77,'(A)') 'NORMAL MODES: Output from readgamnmodes'
       do i=1,ndf
          write(77,'(99F15.10)') L(i,:)
       end do
       write(77,'(A)') 'FORCE CONSTANTS: Output from readgamnmodes'
       do i=1,nvdf
          write(77,'(F15.10)') hii(i)
       end do
 !     -----------------------------------------------------------------
 
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
       write(77,'(A)')'-----------------------------------------------'
       write(77,'(A)')' GENERATING GEOMETRIES FOR HESSIAN COMPUTATION '
       write(77,'(A)')'-----------------------------------------------'
 
       omega = Sqrt(hii)
       dQ = dy/Sqrt(omega)
 
       write(77,'(A)')'SCALE FACTOR FOR DISPLACEMENTS'
       write(77,'(99F15.10)') dQ
 
 
 !     Building mass weights matrix Minv = diag(1/Sqrt(m_i))
       do i=1,nqmatoms
          do j=1,3
             Mass(3*(i-1)+j) = sqrt(atmass(i))
             Minv(3*(i-1)+j) = 1d0/sqrt(atmass(i))
          end do
       end do
 
       write(77,'(A)')'MASS MATRICES'
       write(77,'(99F15.10)') Mass
       write(77,*) 
       write(77,'(99F15.10)') Minv
 !-----------------------------------------------------------------------
 !     GENERATING STENCIL GEOMETRIES
 !-----------------------------------------------------------------------
 !     Converting units from Angstroms to Bohrs. 
       X0=qvageom/a0
 
 !     Calculating displaced coordinates and energy at 6 grid points 
 !     along each normal mode.
       dd=(/-1,1/)
       write(77,'(A)') 'LABELS FOR 1d STENCIL POINTS'
       write(77,'(4I3)') (dd(i),i=1,2)
 
       open(UNIT=16, FILE=qva_cli%ste, ACTION='WRITE', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN FILE geoms.qva'
          STOP
       end if
 !     ------------------------------------------------------------
 !     DEBUG
 !     ------------------------------------------------------------
       write(77,'(A)') 'EQUILIBRIUM GEOMETRY (X0) (ANGS)'
       do i=1,nqmatoms
          write(77,'(3F15.10)') X0(:,i)*a0
       end do
 !     ------------------------------------------------------------
 
       write(77,'(A)') 'DISPLACEMENT VECTORS (dX) (ANGS)'
       step=1
       do nm=1,nvdf
          do j=1,2
             do k=1,nqmatoms
                do m=1,3
                   p = 3*(k-1)+m
                   dX(m,k) = Minv(p)*L(p,nm)*dQ(nm)   ! Scaling/Mass unweighting
                end do
             end do
             Xm = X0 + dd(j)*dX  ! Displacing along nmode
             Xm=Xm*a0            ! BACK TO ANGSTROMS
 !           ------------------------------------------------------------
 !           DEBUG
 !           ------------------------------------------------------------
             write(77,'(A)') 'dX'
             do i=1,nqmatoms
                write(77,'(3F15.10)') dd(j)*dX(:,i)*a0
             end do
             write(77,'(A)') 'DISPLACED GEOM'
             do i=1,nqmatoms
                write(77,'(3F15.10)') Xm(:,i)
             end do
 !           ------------------------------------------------------------
 !           WRITING DISPLACED COORDINATES TO OUTPUT FILE
             write(16,'(2I3,I5)') nm,dd(j),step
             do k=1,nqmatoms
                an=at_numbers(k)
                if (qva_extprog==1) then
                   write(16,'(A2,I5,3F15.10)') atsym(an),an,Xm(:,k)
                elseif (qva_extprog==2) then
                   write(16,'(A2,3F15.10)') atsym(an),Xm(:,k)
                end if
             end do
             step = step + 1
          end do
       end do
 
       close(16,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE'
          STOP
       end if
 
       end subroutine
 
 
 
  
 
 
 
 
 
 
 !######################################################################
 !     FORCE FIELD SECTION
 !######################################################################
       
 
       subroutine hseminumqff(qva_cli,dy,nqmatoms,nclatoms,qmcoords,&
                             &clcoords,at_numbers,ndf,nvdf,&
                             &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,&
                             &uiijj,tijk,uiijk,uijjk,uijkk)
 !     ------------------------------------------------------------------
 !     COMPUTES A SEMINUMERICAL QUARTIC FORCE FIELD USING ANALYTICAL 
 !     HESSIANS COMPUTED USING AN EXTERNAL PROGRAM (GAUSSIAN).
 !     ------------------------------------------------------------------
       implicit none
 
       type(qva_cli_type), intent(in) :: qva_cli
       integer, intent(in) :: ndf                  ! Number of deg of freedm (3*nqmatoms)
       integer, intent(in) :: nvdf                 ! Number of vib degrees of freedom (ndf-6)
       integer, intent(in) :: nqmatoms             ! Number of QM atoms
       integer, intent(in) :: nclatoms             ! Number of classical atoms
       integer, intent(in) :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
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
 
       real*8              :: nmodes(ndf,ndf)      ! mw Hessian eigenvectors.
       real*8              :: freq(ndf)             ! mw Hessian eigenvalues.
       real*8              :: eig(ndf)             ! mw Hessian eigenvalues.
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
       real*8              :: atmass(nqmatoms)     ! Atomic masses
       real*8              :: Mass(ndf)            ! Sqrt mass weights matrix.
       real*8              :: Minv(ndf)            ! Inverse sqrt mass weights matrix.
       real*8              :: gtmp(ndf)
       real*8              :: gnc(nvdf)
       real*8              :: hess0(nvdf,nvdf)
       real*8              :: hess(nvdf,-1:1,nvdf,nvdf)
       real*8              :: tmp1
       real*8              :: tmp2
       real*8              :: tmp3
       integer             :: i, j, k
       integer             :: m, nm, p
       integer             :: nm1, nm2, nm3, gp
     
       include "qvmbia_param.f"
       write(77,'(A)')'------------------------------------------------------------------------' 
       write(77,'(A)')'   STARTING SEMINUMERICAL QUARTIC FORCE FIELD USING EXTERNAL HESSIANS   '
       write(77,'(A)')'        COMPUTED USING  EXTERNAL ELECTRONIC STRUCTURE SOFTWARE          ' 
       write(77,'(A)')'------------------------------------------------------------------------' 
 
 
 !     Read normal modes, frequencies and atomic masses from input file.
       call readgaunmodes(qva_cli,nqmatoms,ndf,nvdf,L,freq,hii,atmass)
 
 
 !     Building mass weights 
       do i=1,nqmatoms
           do j=1,3
               Mass(3*(i-1)+j) = Sqrt(atmass(i))
               Minv(3*(i-1)+j) = 1d0/Sqrt(atmass(i))
           end do
       end do
 
 !     Eliminating rotational and translational degrees of freedom.
 !      do nm=7,ndf
 !         L(:,nm-6)=nmodes(:,nm)
 !         hii(nm-6)=eig(nm)
 !      end do
 
       omega = Sqrt(hii)
       dQ = dy/Sqrt(omega)
 
       write(77,'(A,99D14.6)') 'dQ(nm) = ',dQ
 
 !     READ, MASS-WEIGHT AND CONVERT TO NORMAL COORDINATES  HESSIANS 
 !     AT DISTORTED GEOMETRIES 
 !      call readgamhess(ndf,nvdf,Minv,L,hess,hess0)
       call readgauhess(qva_cli,ndf,nvdf,Minv,L,hess,hess0)
 
 !-----------------------------------------------------------------------
 !     DIAGONAL QFF TERMS
 !-----------------------------------------------------------------------
       tiii=0.0d0
       uiiii=0.0d0
       do nm=1,nvdf
          tiii(nm)=(hess(nm,1,nm,nm)-hess(nm,-1,nm,nm))*0.5d0/dQ(nm)
          uiiii(nm)=(hess(nm,1,nm,nm)-2d0*hess0(nm,nm)+hess(nm,-1,nm,nm))/dQ(nm)**2
       end do
       
 !-----------------------------------------------------------------------
 !     2-MODE COUPLING QFF TERMS
 !-----------------------------------------------------------------------
 !     Calculating QFF 2-mode coupling terms Tiij Tjji Uiiij Ujjji Uiijj 
       tiij=0.0d0
       tjji=0.0d0
       uiiij=0.0d0
       ujjji=0.0d0
       uiijj=0.0d0
       do nm1=1,nvdf-1
          do nm2=nm1+1,nvdf
 
 !           TIIJ
 !           -----------------------------------------------------------
             tmp1=0d0
             tmp1=hess(nm1,1,nm1,nm2)-hess(nm1,-1,nm1,nm2)
             tmp1=tmp1*0.5d0/dQ(nm1)
             tiij(nm1,nm2)=tmp1
 
 !           TJJI
 !           -----------------------------------------------------------
             tmp1=0d0
             tmp1=hess(nm2,1,nm1,nm2)-hess(nm2,-1,nm1,nm2)
             tmp1=tmp1*0.5d0/dQ(nm2)
             tjji(nm1,nm2)=tmp1
 
 !           UIIIJ
 !           -----------------------------------------------------------
             tmp1=0d0
             tmp1=hess(nm1,1,nm1,nm2)-2d0*hess0(nm1,nm2)+hess(nm1,-1,nm1,nm2)
             tmp1=tmp1/dQ(nm1)**2
             uiiij(nm1,nm2)=tmp1
               
 !           UJJJI
 !           -----------------------------------------------------------
             tmp1=0d0
             tmp1=hess(nm2,1,nm1,nm2)-2d0*hess0(nm1,nm2)+hess(nm2,-1,nm1,nm2)
             tmp1=tmp1/dQ(nm2)**2
             ujjji(nm1,nm2)=tmp1
 
 !           UIIJJ
 !           -----------------------------------------------------------
             tmp1=0d0
             tmp2=0d0
             tmp1=hess(nm1,1,nm2,nm2)-2d0*hess0(nm2,nm2)+hess(nm1,-1,nm2,nm2)
             tmp1=tmp1/dQ(nm1)**2
             tmp2=hess(nm2,1,nm1,nm1)-2d0*hess0(nm1,nm1)+hess(nm2,-1,nm1,nm1)
             tmp2=tmp2/dQ(nm2)**2
             uiijj(nm1,nm2)=0.5d0*(tmp1+tmp2)
 
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
 
 !              TIJK
 !              ---------------------------------------------------------
                tmp1=0d0
                tmp2=0d0
                tmp3=0d0
                tmp1=hess(nm1,1,nm2,nm3)-hess(nm1,-1,nm2,nm3)
                tmp1=tmp1*0.5d0/dQ(nm1)
                tmp2=hess(nm2,1,nm1,nm3)-hess(nm2,-1,nm1,nm3)
                tmp2=tmp2*0.5/dQ(nm2)
                tmp3=hess(nm3,1,nm1,nm2)-hess(nm3,-1,nm1,nm2)
                tmp3=tmp3*0.5/dQ(nm3)
                tijk(nm1,nm2,nm3)=(tmp1+tmp2+tmp3)/3d0
 
              
 !              UIIJK
 !              ---------------------------------------------------------
                tmp1=0d0
                tmp1=hess(nm1,1,nm2,nm3)-2d0*hess0(nm2,nm3)+hess(nm1,-1,nm2,nm3)
                tmp1=tmp1/dQ(nm1)**2
                uiijk(nm1,nm2,nm3)=tmp1
 
 
 !              UIJJK
 !              ---------------------------------------------------------
                tmp1=0d0
                tmp1=hess(nm2,1,nm1,nm3)-2d0*hess0(nm1,nm3)+hess(nm2,-1,nm1,nm3)
                tmp1=tmp1/dQ(nm2)**2
                uijjk(nm1,nm2,nm3)=tmp1
 
 !              UIJKK
 !              ---------------------------------------------------------
                tmp1=0d0
                tmp1=hess(nm3,1,nm1,nm2)-2d0*hess0(nm1,nm2)+hess(nm3,-1,nm1,nm2)
                tmp1=tmp1/dQ(nm3)**2
                uijkk(nm1,nm2,nm3)=tmp1
 
             end do
          end do
       end do   
 
       write(77,'(A)') '##############################################'
       write(77,'(A)') '   FINISHED WITH QUARTIC FORCE FIELD COMP     '
       write(77,'(A)') '##############################################'
       write(77,*)
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
       write(77,'(9D14.6)') hii
       end subroutine
 
 
       subroutine readgamhess(ndf,nvdf,Minv,L,hess,hess0)
 
       implicit none
 !     -----------------------------------------------------------
 !     READS HESSIAN IN GAMESS FORMAT ($HESS GROUP IN PUNCHFILE)
 !     -----------------------------------------------------------
       integer,intent(in)   :: ndf
       integer,intent(in)   :: nvdf
       real*8,intent(in)    :: Minv(ndf)
       real*8,intent(in)    :: L(ndf,nvdf)
       real*8,intent(out)   :: hess(nvdf,-1:1,nvdf,nvdf)
       real*8,intent(out)   :: hess0(nvdf,nvdf)
 !     Local variables
       integer   :: a,b
       integer   :: i,j
       integer   :: n,d
       integer   :: minc,maxc
       integer   :: openstatus
       integer   :: closestatus
       integer   :: nm,dsp,nh
       real*8    :: thess(ndf,ndf)
       real*8    :: tmp(ndf,nvdf)
       real*8    :: nhess(nvdf,nvdf)
       character(len=5)   :: txt1
       character(len=*),parameter :: FMT1 = "(I2,I3,5D15.8)"
       character(len=*),parameter :: FMT2 = "(A5,2I3,I5)"
 !     -----------------------------------------------------------
 
       open(UNIT=15, FILE='hess.gam', ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN FILE hess.gam'
          STOP
       end if
 
       hess=0d0
       DO n=0,nvdf
          DO d=-1,1,2
 
 !           READING HESSIAN FROM FILE
 !            READ(15,*)
             READ(15,FMT2) txt1,nm,dsp,nh
             READ(15,*)
             if (n==0 .AND. nm/=0) STOP('ERROR: THE FIRST HESSIAN READ DOES NOT&
                                       &CORRESPOND TO THE EQUILIBRIUM GEOMETRY')
             DO i=1,ndf
                DO minc = 1,ndf,5
                   maxc=minc+4
                   IF(maxc.GT.ndf) maxc=ndf
                   READ(15,FMT1) a,b,(thess(i,j),j=minc,maxc)
                END DO
             END DO
 
             write(77,'(2I2)') nm,dsp
             do i=1,ndf
                write(77,'(99D15.8)') (thess(i,j),j=1,ndf)
             end do
 
 !           MASS-WEIGHTING HESSIAN MATRIX
             do i=1,ndf
                do j=i,ndf
                   thess(i,j)=thess(i,j)*Minv(i)*Minv(j)
                   if (i/=j) thess(j,i)=thess(i,j)
                end do
             end do
 
             write(77,'(A)') 'MASS WEIGHTED HESSIAN'
             write(77,'(2I2)') nm,dsp
             do i=1,ndf
                write(77,'(99D15.8)') (thess(i,j),j=1,ndf)
             end do
 
 !           CONVERTINV TO NORMAL COORDINATES
 !           This procedure reduces the dimension of the hessian matrix
 !           from 3Nx3N to 3N-6x3N-6.
             tmp=0d0
             nhess=0d0
             call dsymm('L','U',ndf,nvdf,1d0,thess,ndf,L,ndf,0d0,tmp,ndf)
             call dgemm('T','N',nvdf,nvdf,ndf,1d0,L,ndf,tmp,ndf,0d0,nhess,nvdf)
 
             write(77,'(A)') 'MASS WEIGHTED HESSIAN IN NORMAL COORDINATES'
             write(77,'(2I2)') nm,dsp
             do i=1,nvdf
                write(77,'(99D15.8)') (nhess(i,j),j=1,nvdf)
             end do
 
             if (n==0) then
                hess0=nhess
                EXIT
             else
                hess(nm,dsp,:,:)=nhess
             end if
 
          END DO
       END DO
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE'
          STOP
       end if
 
 !     ------------------------------------------------------------------
 !     DEBUG
 !     ------------------------------------------------------------------
       write(77,'(A)') '--------------------------------------------'
       write(77,'(A)') 'MASS WEIGHTED HESSIANS IN NORMAL COORDINATES'
       write(77,'(A)') '--------------------------------------------'
       write(77,'(A)') 'EQUILIBRIUM GEOMETRY'
       DO i=1,nvdf
          write(77,'(99D14.6)') (hess0(i,j),j=1,nvdf)
       END DO
       DO nm=1,nvdf
          DO dsp=-1,1,2
             write(77,*)
             write(77,'(2I3)') nm,dsp
             DO i=1,nvdf
                write(77,'(99D14.6)') (hess(nm,dsp,i,j),j=1,nvdf)
             END DO
          END DO
       END DO
 !     ------------------------------------------------------------------
 
       end subroutine
 
       subroutine readgauhess(qva_cli,ndf,nvdf,Minv,L,hess,hess0)
 
       implicit none
 !     -----------------------------------------------------------
 !     READS HESSIAN IN GAMESS FORMAT ($HESS GROUP IN PUNCHFILE)
 !     -----------------------------------------------------------
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in)   :: ndf
       integer,intent(in)   :: nvdf
       real*8,intent(in)    :: Minv(ndf)
       real*8,intent(in)    :: L(ndf,nvdf)
       real*8,intent(out)   :: hess(nvdf,-1:1,nvdf,nvdf)
       real*8,intent(out)   :: hess0(nvdf,nvdf)
 !     Local variables
       integer   :: ne
       integer   :: a,b
       integer   :: i,j
       integer   :: n,d
       integer   :: minc,maxc
       integer   :: openstatus
       integer   :: closestatus
       integer   :: nm,dsp,nh
       real*8    :: phess(ndf+ndf*(ndf-1)/2)
       real*8    :: thess(ndf,ndf)
       real*8    :: tmp(ndf,nvdf)
       real*8    :: nhess(nvdf,nvdf)
       character(len=5)   :: txt1
       character(len=*),parameter :: FMT1 = "(5D16.8)"
       character(len=*),parameter :: FMT2 = "(A5,2I3,I5)"
 !     -----------------------------------------------------------
 
       open(UNIT=15, FILE=qva_cli%hes, ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN FILE hess.gau'
          STOP
       end if
 
       ne=ndf+ndf*(ndf-1)/2
       hess=0d0
       DO n=0,nvdf
          DO d=-1,1,2
 
 !           READING HESSIAN FROM FILE
             READ(15,FMT2) txt1,nm,dsp,nh
             if (n==0 .AND. nm/=0) then
                STOP('ERROR: THE FIRST HESSIAN READ DOES NOT&
                     &CORRESPOND TO THE EQUILIBRIUM GEOMETRY')
             end if
             READ(15,FMT1) (phess(j),j=1,ne)
 
             do i=1,ndf
                do j=i,ndf
                   thess(i,j) = phess(i+j*(j-1)/2)
                   if (i/=j) thess(j,i) = thess(i,j)
                end do
             end do
 
             write(77,'(2I2)') nm,dsp
             do i=1,ndf
                write(77,'(99D15.8)') (thess(i,j),j=1,ndf)
             end do
 
 !           MASS-WEIGHTING HESSIAN MATRIX
             do i=1,ndf
                do j=i,ndf
                   thess(i,j)=thess(i,j)*Minv(i)*Minv(j)
                   if (i/=j) thess(j,i)=thess(i,j)
                end do
             end do
 
             write(77,'(A)') 'MASS WEIGHTED HESSIAN'
             write(77,'(2I2)') nm,dsp
             do i=1,ndf
                write(77,'(99D15.8)') (thess(i,j),j=1,ndf)
             end do
 
 !           CONVERTING TO NORMAL COORDINATES
 !           This procedure reduces the dimension of the hessian matrix
 !           from 3Nx3N to 3N-6x3N-6.
             tmp=0d0
             nhess=0d0
             call dsymm('L','U',ndf,nvdf,1d0,thess,ndf,L,ndf,0d0,tmp,ndf)
             call dgemm('T','N',nvdf,nvdf,ndf,1d0,L,ndf,tmp,ndf,0d0,nhess,nvdf)
 
             write(77,'(A)') 'MASS WEIGHTED HESSIAN IN NORMAL COORDINATES'
             write(77,'(2I2)') nm,dsp
             do i=1,nvdf
                write(77,'(99D15.8)') (nhess(i,j),j=1,nvdf)
             end do
 
             if (n==0) then
                hess0=nhess
                EXIT
             else
                hess(nm,dsp,:,:)=nhess
             end if
 
          END DO
       END DO
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE'
          STOP
       end if
 
 !     ------------------------------------------------------------------
 !     DEBUG
 !     ------------------------------------------------------------------
       write(77,'(A)') '--------------------------------------------'
       write(77,'(A)') 'MASS WEIGHTED HESSIANS IN NORMAL COORDINATES'
       write(77,'(A)') '--------------------------------------------'
       write(77,'(A)') 'EQUILIBRIUM GEOMETRY'
       DO i=1,nvdf
          write(77,'(99D14.6)') (hess0(i,j),j=1,nvdf)
       END DO
       DO nm=1,nvdf
          DO dsp=-1,1,2
             write(77,*)
             write(77,'(2I3)') nm,dsp
             DO i=1,nvdf
                write(77,'(99D14.6)') (hess(nm,dsp,i,j),j=1,nvdf)
             END DO
          END DO
       END DO
 !     ------------------------------------------------------------------
 
       end subroutine
 
 
 
 
       subroutine readnmodes(ndf,nmdes)
 !     ------------------------------------------------------------------
 !     THIS SUBROUTINE READS NORMAL MODES FROM FILE NAMED nmodes.qba
 !     ------------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: ndf
       real*8,intent(out)   :: nmdes(ndf,ndf)
 
       integer   :: openstatus
       integer   :: closestatus
       integer   :: i,j
 !     ------------------------------------------------------------------
 
       open(UNIT=15, FILE='nmodes.qba', ACTION='READ',iostat=openstatus)
 
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN FILE nmodes.qba'
          STOP
       end if
 
       do i=1,ndf
          read(15,*) (nmdes(i,j),j=1,ndf)
       end do
 
       close(15,iostat=closestatus)
 
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE nmodes.qba'
          STOP
       end if 
 
       end subroutine
 
       subroutine readqff(qva_cli,nvdf,hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
                  & tijk,uiijk,uijjk,uijkk)
 !     ------------------------------------------------------------------
 !     THIS SUBROUTINE READS QUARTIC FORCE FIELD FROM FILE NAMED qff.qba
 !     ------------------------------------------------------------------
       implicit none
 
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in) :: nvdf
       real*8,intent(out) :: hii(nvdf)
       real*8,intent(out) :: tiii(nvdf)
       real*8,intent(out) :: uiiii(nvdf)
       real*8,intent(out) :: tiij(nvdf,nvdf)
       real*8,intent(out) :: tjji(nvdf,nvdf)
       real*8,intent(out) :: uiiij(nvdf,nvdf)
       real*8,intent(out) :: ujjji(nvdf,nvdf)
       real*8,intent(out) :: uiijj(nvdf,nvdf)
       real*8,intent(out) :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: uijkk(nvdf,nvdf,nvdf)
 
       integer   :: i,j,k
       integer   :: nm, nm1, nm2, nm3
 !     ------------------------------------------------------------------
 
       open(UNIT=15, FILE=qva_cli%qff, ACTION='READ')
 
       do i=1,nvdf
          read(15,*) nm, hii(nm), tiii(nm), uiiii(nm)
          if (abs(tiii(nm)) < 1d-9) tiii(nm)=0d0
          if (abs(uiiii(nm)) < 1d-9) uiiii(nm)=0d0
       end do
 
       do i=1,nvdf-1
          do j=i+1,nvdf
             read(15,*) nm1, nm2, tiij(nm1,nm2), tjji(nm1,nm2),&
                   & uiiij(nm1,nm2), ujjji(nm1,nm2), uiijj(nm1,nm2)
             if (abs(tiij(nm1,nm2)) < 1d-9) tiij(nm1,nm2)=0d0
             if (abs(tjji(nm1,nm2)) < 1d-9) tjji(nm1,nm2)=0d0
             if (abs(uiiij(nm1,nm2)) < 1d-9) uiiij(nm1,nm2)=0d0
             if (abs(ujjji(nm1,nm2)) < 1d-9) ujjji(nm1,nm2)=0d0
             if (abs(uiijj(nm1,nm2)) < 1d-9) uiijj(nm1,nm2)=0d0
          end do
       end do
 
       do i=1,nvdf-2
          do j=i+1,nvdf-1
             do k=j+1,nvdf
                read(15,*) nm1, nm2, nm3, tijk(nm1,nm2,nm3), uiijk(nm1,nm2,nm3),&
                         & uijjk(nm1,nm2,nm3),uijkk(nm1,nm2,nm3)
                if (abs(tijk(nm1,nm2,nm3)) < 1d-9) tijk(nm1,nm2,nm3)=0d0
                if (abs(uiijk(nm1,nm2,nm3)) < 1d-9) uiijk(nm1,nm2,nm3)=0d0
                if (abs(uijjk(nm1,nm2,nm3)) < 1d-9) uijjk(nm1,nm2,nm3)=0d0
                if (abs(uijkk(nm1,nm2,nm3)) < 1d-9) uijkk(nm1,nm2,nm3)=0d0
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
 
 
       end subroutine
 
       subroutine readgamessqff(qva_cli,nvdf,hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
                  & tijk,uiijk,uijjk,uijkk)
 !     ------------------------------------------------------------------
 !     THIS SUBROUTINE READS QUARTIC FORCE FIELD FROM FILE NAMED qff.qba
 !     ------------------------------------------------------------------
       implicit none
 
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in) :: nvdf
       real*8,intent(out) :: hii(nvdf)
       real*8,intent(out) :: tiii(nvdf)
       real*8,intent(out) :: uiiii(nvdf)
       real*8,intent(out) :: tiij(nvdf,nvdf)
       real*8,intent(out) :: tjji(nvdf,nvdf)
       real*8,intent(out) :: uiiij(nvdf,nvdf)
       real*8,intent(out) :: ujjji(nvdf,nvdf)
       real*8,intent(out) :: uiijj(nvdf,nvdf)
       real*8,intent(out) :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: uijkk(nvdf,nvdf,nvdf)
 
       integer   :: i,j,k
       integer   :: nm, nm1, nm2
       integer   :: openstatus
       integer   :: closestatus
       character*6 :: chaff
       real*8,parameter  :: amu2au = 1d0/0.0005485799111d0
 !     ------------------------------------------------------------------
 
       open(UNIT=15, FILE=qva_cli%qff, ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN QFF FILE '
          STOP
       end if
 
       do i=1,nvdf
          read(15,*) nm
          nm=nm-6
          read(15,*) chaff, hii(nm)
          read(15,*) chaff, tiii(nm)
          read(15,*) chaff, uiiii(nm)
          hii(nm)=hii(nm)/amu2au
          tiii(nm)=tiii(nm)/Sqrt(amu2au**3)
          uiiii(nm)=uiiii(nm)/amu2au**2
       end do
 
       do i=1,nvdf-1
          do j=i+1,nvdf
             read(15,*) nm2, nm1
             nm1=nm1-6
             nm2=nm2-6
             read(15,*) chaff, tjji(nm1,nm2)
             read(15,*) chaff, tiij(nm1,nm2)
             read(15,*) chaff, ujjji(nm1,nm2)
             read(15,*) chaff, uiiij(nm1,nm2)
             read(15,*) chaff, uiijj(nm1,nm2)
             tjji(nm1,nm2)=tjji(nm1,nm2)/Sqrt(amu2au**3)
             tiij(nm1,nm2)=tiij(nm1,nm2)/Sqrt(amu2au**3)
             ujjji(nm1,nm2)=ujjji(nm1,nm2)/amu2au**2
             uiiij(nm1,nm2)=uiiij(nm1,nm2)/amu2au**2
             uiijj(nm1,nm2)=uiijj(nm1,nm2)/amu2au**2
          end do
       end do
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE QFF FILE '
          STOP
       end if 
 
       tijk=0d0
       uiijk=0d0
       uijjk=0d0
       uijkk=0d0
 
       write(*,*) '----------------------------------------------'
       write(*,*) ' FINISHED READING GAMESS QUARTIC FORCE FIELD  '
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
 
       end subroutine
 
  
 
       subroutine readgamnmodes(qva_cli,nqmatoms,ndf,nmodes,freq,eig,atmass)
       implicit none
 
 !     -----------------------------------------------------------
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in)             :: ndf
       integer,intent(in)             :: nqmatoms
       real*8,intent(out)             :: nmodes(ndf,ndf)
       real*8,intent(out)             :: freq(ndf)
       real*8,intent(out)             :: eig(ndf)
       real*8,intent(out)             :: atmass(nqmatoms)
 
       integer   :: igr
       integer   :: nat
       integer   :: i,j,k,l
       integer   :: m,at,nm
       integer   :: openstatus
       integer   :: closestatus
       character(len=*),parameter :: FMT1 = "(3D17.9)"
       character(len=*),parameter :: FMT2 = "(A,I5,A,F10.5,A)"
       character(len=4)  :: txt1
       character(len=13) :: txt2
       character(len=9)  :: txt3
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
       real*8,parameter  :: amu2au = 1d0/0.0005485799111d0
 !     -----------------------------------------------------------
       nat=nqmatoms
 
       open(UNIT=15, FILE=qva_cli%nmo, ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN FILE'
          STOP
       end if
 
 !     READING THE FILE
       read(15,*)
       read(15,'(3F12.5)') atmass
       do m=1,3*nat
          read(15,FMT2) txt1,igr,txt2,freq(m),txt3
          do i=1,nat
             read(15,FMT1) (nmodes(j+3*(i-1),m),j=1,3)
          end do
       end do
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE'
          STOP
       end if
 
 !     MASS WEIGHT NORMAL MODES
       do nm=1,3*nat
          do at=1,nat
             do j=1,3
                k=j+3*(at-1)
                nmodes(k,nm) = nmodes(k,nm)*Sqrt(atmass(at))
             end do
          end do
       end do
 
 !     CONVERT FREQUENCIES TO AU
       eig = (freq/h2cm)**2
 
 !     CONVERT MASSES TO AU
       atmass = atmass*amu2au
 
 !     -----------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------
       write(77,'(A)') 'NORMAL MODES'
       DO i=1,nat
          DO j=1,3
             k=j+3*(i-1)
             write(77,'(A,9F12.8,A)') '{',(nmodes(k,l),l=1,3*nat),'}'
          END DO
       END DO
 
       write(77,'(A)') 'FREQUENCIES (cm-1)'
       write(77,'(99F12.5)') (freq(j),j=1,3*nat)
 
       write(77,'(A)') 'FREQUENCIES (cm-1)'
       write(77,'(99D12.5)') (eig(j),j=1,3*nat)
 
       write(77,'(A)') 'ATOMIC MASSES'
       write(77,'(99F14.5)') (atmass(j),j=1,nat)
 !     -----------------------------------------------------------
 
       end subroutine
 
 
       subroutine readgaunmodes(qva_cli,nqmatoms,ndf,nvdf,nmodes,freq,fc,atmass)
       implicit none
 
 !     -----------------------------------------------------------
       type(qva_cli_type), intent(in) :: qva_cli
       integer,intent(in)             :: ndf
       integer,intent(in)             :: nvdf
       integer,intent(in)             :: nqmatoms
       real*8,intent(out)             :: nmodes(ndf,nvdf)
       real*8,intent(out)             :: freq(nvdf)
       real*8,intent(out)             :: fc(nvdf)
       real*8,intent(out)             :: atmass(nqmatoms)
 
       integer   :: ne
       integer   :: igr
       integer   :: nat
       integer   :: i,j,k,l
       integer   :: m,at,nm
       integer   :: openstatus
       integer   :: closestatus
       real*8    :: nmt(ndf)
       real*8    :: orth(nvdf,nvdf)
       real*8    :: nmd(nvdf*ndf)
       real*8    :: norm
       character(len=49) :: txt1
       character(len=*),parameter :: FMT1 = "(5D16.8)"
       character(len=*),parameter :: FMT2 = "(A49,I12)"
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
       real*8,parameter  :: amu2au = 1d0/0.0005485799111d0
       real*8,external :: ddot
       real*8,external :: dnrm2
 !     -----------------------------------------------------------
       nat=nqmatoms
 
       open(UNIT=15, FILE=qva_cli%nmo, ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN FILE'
          STOP
       end if
 
 !     READING THE FILE
       READ(15,*)
       READ(15,'(5D16.8)') (freq(i),i=1,nvdf)
       READ(15,*)
       READ(15,'(5D16.8)') (atmass(i),i=1,nat)
       READ(15,FMT2) txt1,ne
       READ(15,FMT1) (nmd(j),j=1,ne)
 
       k=1
       do m=1,nvdf
          do i=1,3*nat
             nmodes(i,m) = nmd(k)
             k=k+1
          end do
       end do
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE'
          STOP
       end if
 
 !     MASS WEIGHT NORMAL MODES
       do nm=1,nvdf
          do at=1,nat
             do j=1,3
                k=j+3*(at-1)
                nmodes(k,nm) = nmodes(k,nm)*Sqrt(atmass(at))
             end do
          end do
       end do
 
 !     NORMALIZE MASS-WEIGHTED NORMAL MODES
       do nm=1,nvdf
          nmt=nmodes(:,nm)
          write(77,'(99D16.8)') nmodes(:,nm)
          write(77,'(99D16.8)') nmt
          write(77,'(F16.8)') dnrm2(ndf,nmt,1)
          nmt=nmt/dnrm2(ndf,nmt,1)
          write(77,'(99D16.8)') nmt
          write(77,'(F16.8)') dnrm2(ndf,nmt,1)
          nmodes(:,nm)=nmt
       end do
 
       call dgemm('T','N',nvdf,nvdf,ndf,1d0,nmodes,ndf,nmodes,ndf,0d0,orth,nvdf)
       do i=1,nvdf
          write(77,'(99D15.6)') orth(i,:)
       end do
 
        write(77,*) freq
 !     CONVERT FREQUENCIES TO AU
       fc = 0d0
       fc = freq
       fc = (fc/h2cm)**2
 
 !     CONVERT MASSES TO AU
       atmass = atmass*amu2au
 
 !     -----------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------
       write(77,'(A)') 'MASS-WEIGHTED NORMAL MODES'
       DO i=1,nat
          DO j=1,3
             k=j+3*(i-1)
             write(77,'(99F12.8)') (nmodes(k,l),l=1,nvdf)
          END DO
       END DO
 
       write(77,'(A)') 'FREQUENCIES (cm-1)'
       write(77,'(99F12.5)') (freq(j),j=1,nvdf)
 
       write(77,'(A)') 'FREQUENCIES (cm-1)'
       write(77,'(99D12.5)') (fc(j),j=1,nvdf)
 
       write(77,'(A)') 'ATOMIC MASSES'
       write(77,'(99F14.5)') (atmass(j),j=1,nat)
 !     -----------------------------------------------------------
 
       end subroutine
 
 
 !#######################################################################
 !     VIBRATIONAL SELF-CONSISTENT FIELD SECTION
 !#######################################################################
 
 
 
       subroutine buildoperators(hii,tiii,uiiii,gwidth,ngaus,nvdf,&
                     &S,Q1,Q2,Q3,Hcore)
       implicit none
 
 !     ------------------------------------------------------------------
       integer, intent(in) :: ngaus  ! Number of primitive gaussian basis functions.
       integer, intent(in) :: nvdf   ! Number of vibrational degrees of freedom
       real*8, intent(in)  :: hii(nvdf)   ! Hessian eigenvalues in AU.(d^2 V/dQi^2)
       real*8, intent(in)  :: tiii(nvdf)  ! d^3 V/dQi^3
       real*8, intent(in)  :: uiiii(nvdf) ! d^4 V/dQi^4
       real*8, intent(in)  :: gwidth    ! Gaussian primitives width parameter. Default 0.5.
       real*8, intent(out) :: S(ngaus,ngaus,nvdf)     ! Overlap Matrix
       real*8, intent(out) :: Q1(ngaus,ngaus,nvdf)    ! Q^1 Operator
       real*8, intent(out) :: Q2(ngaus,ngaus,nvdf)    ! Q^2 Operator
       real*8, intent(out) :: Q3(ngaus,ngaus,nvdf)    ! Q^3 Operator
       real*8, intent(out) :: Hcore(ngaus,ngaus,nvdf) ! Core hamiltonian.
 !     PARAMETERS     
       real*8,parameter    ::  a0 = 0.5291771D00      ! bohr radius
       real*8,parameter    ::  pi = 2.0d0*acos(0.0d0) ! PI
 !     Gauss-Hermite quadrature points for N=16.
       real*8,parameter    ::  ghx16(16)=(/-4.688738939305818d+0,&
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
 !     LOCAL VARIABLES
       character(len=117) :: format1,format2
       real*8   :: A              ! Auxiliary variable for integrals computation.
       real*8   :: B              ! Auxiliary variable for integrals computation.
       real*8   :: C              ! Auxiliary variable for integrals computation.
       real*8   :: D              ! Auxiliary variable for integrals computation.
       real*8   :: Ai(ngaus,nvdf) ! Auxiliary variable for integrals computation.
       real*8   :: omega(nvdf)    ! Vibrational frequencies in AU.
       real*8   :: T(ngaus,ngaus) ! Kinetic energy matrix.
       real*8   :: ghQ(ngaus,nvdf)! Scaled and translated Gauss-Hermite points.
       real*8   :: Q4(ngaus,ngaus,nvdf) ! Auxiliary variable for integrals computation.
       integer  :: i,j,nm         ! Indexes
 !     ------------------------------------------------------------------
 
       format1='(16D11.3)'
       format2='(16F11.5)'
 
 !     DISTRIBUTED GAUSSIAN BASIS SET
 !     Gaussians are placed in Gauss-Hermite quadrature points, and 
 !     scaled so their centers coincide with the nodes of an harmonic
 !     oscillator wavefunction of vibrational quantum number equal to 
 !     the number of gaussians specified by the user (only 16 for now). 
 !     GH quadrature points were calculated using John Burkardt's code
 !     which can be found at
 !     people.sc.fsu.edu/~jburkardt/f_src/hermite_rule/hermite_rule.html
 !     GH points ghx are scaled and displaced acording to
 !     ghQ = ghx/Sqrt(omega_a)
 !     Distributed Gaussians have the form
 !     g(j,nm) = (2 A(j,nm)/pi)^1/4 exp[-A(j,nm)*(Qi-ghQ(j,nm))^2]
 !     So to define the gaussian basis set we must define a set of
 !     centers (GH points, ghQ) and widths (A(j,nm)).
 !     For details see Hamilton & Light (1986)
 
 !     Scaling and displacing GH quadrature points
       ghQ=0.0D0
       omega = Sqrt(hii)
       do nm=1,nvdf
          ghQ(:,nm)=ghx16/Sqrt(omega(nm))
       end do
 
 !     Building gaussian parameters A(j,nm)
       do nm=1,nvdf
          Ai(1,nm) = gwidth*gwidth/(ghQ(2,nm)-ghQ(1,nm))**2
          do j=2,ngaus-1
             Ai(j,nm) = 4.0D0*gwidth*gwidth/(ghQ(j+1,nm)-ghQ(j-1,nm))**2
          end do 
          Ai(ngaus,nm) = gwidth*gwidth/(ghQ(ngaus,nm)-ghQ(ngaus-1,nm))**2
       end do
 
 !     Computing core hamiltonian and Q^n integrals.
       S = 0.0D0
       Q1 = 0.0D0
       Q2 = 0.0D0
       Q3 = 0.0D0
       Q4 = 0.0D0
 
       do nm=1,nvdf
          T=0.0D0
          do i=1,ngaus
             do j=i,ngaus
                A=Sqrt(Sqrt(4.0D0*Ai(i,nm)*Ai(j,nm)/PI**2))
                B=Sqrt(Ai(i,nm)+Ai(j,nm))
                C=Ai(i,nm)*Ai(j,nm)*(ghQ(i,nm)-ghQ(j,nm))**2/B**2
                D=(Ai(i,nm)*ghQ(i,nm)+Ai(j,nm)*ghQ(j,nm))/B
                if (i == j) C=0d0
 
 !              Overlap matrix
                if (i /= j) then 
                    S(i,j,nm)=Sqrt(PI)*A*Exp(-C)/B
                    S(j,i,nm)=S(i,j,nm)
                else
                    S(i,i,nm) = 1d0
                end if
 
 !              Kinetic Energy matrix
                T(i,j)=2.0d0*S(i,j,nm)*Ai(i,nm)*Ai(j,nm)*(-1.0D0+2.0D0*C)/B**2
                if (i /= j) T(j,i)=T(i,j)
 
 !              Q matrix elements
                Q1(i,j,nm)=S(i,j,nm)*D/B
                Q2(i,j,nm)=S(i,j,nm)*(1.0D0+2.0D0*D**2)*0.5d0/B**2
                Q3(i,j,nm)=S(i,j,nm)*(3.0D0*D + 2.0d0*D**3)*0.5d0/B**3
                Q4(i,j,nm)=S(i,j,nm)*(3.0D0+12.0D0*D**2+4.0D0*D**4)*0.25d0/B**4
 
                if (i /= j) Q1(j,i,nm)=Q1(i,j,nm)
                if (i /= j) Q2(j,i,nm)=Q2(i,j,nm)
                if (i /= j) Q3(j,i,nm)=Q3(i,j,nm)
                if (i /= j) Q4(j,i,nm)=Q4(i,j,nm)
             end do
          end do
 
 !        Building core hamiltonian
          Hcore(:,:,nm) = -0.5d0*T(:,:) + 0.5d0*hii(nm)*Q2(:,:,nm) + &
                       & tiii(nm)*Q3(:,:,nm)/6.0d0 + uiiii(nm)*Q4(:,:,nm)/24.0d0
       end do
 !     -----------------------------------------------------------------
 !     DEBUG
 !     Hcore=Hcore*1.0d01
 !     -----------------------------------------------------------------
 !      write(77,'(A)') 'PRINTING CORE HAMILTONIAN'
 !      write(77,'(A)') 'normal mode 1'
 !      do i=1,ngaus
 !         write(78,format2) Hcore(i,:,1)
 !      end do
 !      write(77,'(A)') 'normal mode 2'
 !      do i=1,ngaus
 !         write(77,format1) Hcore(i,:,2)
 !      end do
 !      write(77,'(A)') 'normal mode 3'
 !      do i=1,ngaus
 !         write(77,format1) Hcore(i,:,2)
 !      end do
 !
 !
 !      write(77,'(A)') 'PRINTING OVERLAP MATRIX'
 !      do nm=1,1
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         do i=1,ngaus
 !            write(99,format2) S(i,:,nm)
 !         end do
 !      end do
 !      write(77,'(A)') 'PRINTING Q1 MATRIX'
 !      do nm=1,nvdf
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         write(77,format1) Q1(:,:,nm)
 !      end do
 !      write(77,'(A)') 'PRINTING Q2 MATRIX'
 !      do nm=1,nvdf
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         write(77,format1) Q2(:,:,nm)
 !      end do
 !      write(77,'(A)') 'PRINTING Q3 MATRIX'
 !      do nm=1,nvdf
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         write(77,format1) Q3(:,:,nm)
 !      end do
 !      write(77,'(A)') 'PRINTING Q4 MATRIX'
 !      do nm=1,nvdf
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         write(77,format1) Q4(:,:,nm)
 !      end do
 !     ------------------------------------------------------------------
 
       end subroutine
 
       subroutine cholesky(S,nvdf,ngaus,Scho,IScho)
 !     ------------------------------------------------------------------
 !     This subroutine computes S**(-1/2) using the following steps
 !     1) Diagonalize S
 !        S U = U D (where D is diagonal and U orthogonal)
 !     2) Compute D**(-1/2)
 !     3) Compute S**(-1/2) = U D**(-1/2) U**T
 !     ------------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: nvdf                  ! Number of classical atoms
       integer,intent(in)  :: ngaus                 ! Number of gaussian primitives
       real*8,intent(in)   :: S(ngaus,ngaus,nvdf)   ! Overlap matrix
       real*8,intent(out)  :: Scho(ngaus,ngaus,nvdf)   ! Overlap matrix
       real*8,intent(out)  :: IScho(ngaus,ngaus,nvdf)   ! Overlap matrix
 
       character(len=117) :: format1
       integer  ::  i,nm
       real*8   ::  Mtrx(ngaus,ngaus)
       integer  ::  INFO
 
       Scho=0.0D0 
       IScho=0.0D0 
       do nm=1,nvdf
 
          Mtrx=S(:,:,nm)
 
          call dpotrf('U',ngaus,Mtrx,ngaus,INFO)
          if (INFO /= 0) STOP('Error in Cholesky decomposition')
          Scho(:,:,nm) = Mtrx
 
          call dtrtri('U','N',ngaus,Mtrx,ngaus,INFO)
          if (INFO /= 0) STOP('Error in Cholesky factor inversion')
          IScho(:,:,nm) = Mtrx
 
       end do
 
 !     ------------------------------------------------------------------
 !     DEBUG
 !     ------------------------------------------------------------------
 !      format1='(16D12.4)'
 !
 !      write(77,'(A)') 'PRINTING CHOLESKY FACTORS'
 !      do nm=1,nvdf
 !         write(77,'(A,I3)') 'normal mode ',nm
 !         do i=1,ngaus
 !            write(77,format1) Scho(i,:,nm)
 !         end do
 !      end do
 !
 !      write(77,'(A)') 'PRINTING CHOLESKY INVERSES'
 !      do nm=1,nvdf
 !         write(77,'(A,I3)') 'normal mode ',nm
 !         do i=1,ngaus
 !            write(77,format1) IScho(i,:,nm)
 !         end do
 !      end do
 !     ------------------------------------------------------------------
 
       end subroutine
 
       subroutine computSsqrt(S,nvdf,ngaus,Ssqrt)
 !     ------------------------------------------------------------------
 !     This subroutine computes S**(-1/2) using the following steps
 !     1) Diagonalize S
 !        S U = U D (where D is diagonal and U orthogonal)
 !     2) Compute D**(-1/2)
 !     3) Compute S**(-1/2) = U D**(-1/2) U**T
 !     ------------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: nvdf                  ! Number of classical atoms
       integer,intent(in)  :: ngaus                 ! Number of gaussian primitives
       real*8,intent(in)   :: S(ngaus,ngaus,nvdf)   ! Overlap matrix
       real*8,intent(out)  :: Ssqrt(ngaus,ngaus,nvdf)   ! Overlap matrix
 
       character(len=117) :: format1
       integer  ::  i,nm
       real*8   ::  U(ngaus,ngaus)
       real*8   ::  US(ngaus,ngaus)
       real*8   ::  USUt(ngaus,ngaus)
       real*8   ::  Seval(ngaus)
 !     Workspace for diagonalization. 
       real*8,  dimension(:), allocatable :: WORK, WORK2
       integer, dimension(:), allocatable :: IWORK, IWORK2
       integer :: LWORK, LIWORK, INFO
 
       Ssqrt=0.0d0
       do nm=1,nvdf
          U=S(:,:,nm)  ! Initialization of 
          US=0.0D0     ! allocatable
          USUt=0.0D0   ! variables
          Seval=0.0D0  !
 
 !        Diagonalizing overlap matrix
          allocate ( WORK(10000), IWORK(10000) )
          LWORK=-1
          call dsyevd('V','U',ngaus,U,ngaus,Seval,WORK,LWORK,      &
             &IWORK, LWORK, INFO)
          LWORK=WORK(1)
          LIWORK=IWORK(1)
          if(allocated(WORK2)) deallocate (WORK2,IWORK2)
          allocate (WORK2(LWORK),IWORK2(LIWORK))
          call dsyevd('V','U',ngaus,U,ngaus,Seval,WORK2,LWORK,     &
             &IWORK2, LIWORK, INFO)
          deallocate( WORK, IWORK, WORK2, IWORK2 )
 
 !        Computing S**(-1/2) (Ssqrt)
          Seval = 1.0D0/Sqrt( Seval )
          do i=1,ngaus
             US(:,i)=U(:,i)*Seval(i)
          end do
          call dgemm('N','T',ngaus,ngaus,ngaus,1.0D0,US,           &
               & ngaus,U,ngaus,0.0D0,USUt,ngaus)
 
          Ssqrt(:,:,nm)=USUt(:,:)
          
       end do
 !     -----------------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------------
 !      format1='(16D11.3)'
 !      write(77,'(A)') 'PRINTING S^(-1/2)'
 !      write(77,'(A)') 'normal mode 1'
 !      write(77,format1) Ssqrt(:,:,1)
 !      write(77,'(A)') 'normal mode 2'
 !      write(77,format1) Ssqrt(:,:,2)
 !      write(77,'(A)') 'normal mode 3'
 !      write(77,format1) Ssqrt(:,:,3)
 !     -----------------------------------------------------------------
       end subroutine
  
       subroutine solveGEV(Hcore,Scho,IScho,nvdf,ngaus,P,Po,Emod)
 !     ------------------------------------------------------------------
 !
 !     ------------------------------------------------------------------
       implicit none
 
       integer, intent(in)    :: nvdf
       integer, intent(in)    :: ngaus
       real*8, intent(inout)  :: Hcore(ngaus,ngaus,nvdf)
       real*8, intent(in)     :: Scho(ngaus,ngaus,nvdf)
       real*8, intent(in)     :: IScho(ngaus,ngaus,nvdf)
       real*8, intent(out)    :: P(ngaus,ngaus,nvdf)
       real*8, intent(out)    :: Po(ngaus,ngaus,nvdf)
       real*8, intent(out)    :: Emod(ngaus,nvdf)
 
       character(len=117) :: format1
       integer :: nm,i,ex
       integer :: halfg
       real*8  :: x
       real*8  :: norm
       real*8  :: csum
       real*8  :: eval(ngaus)
       real*8  :: vec(ngaus)
       real*8  :: Mtrx(ngaus,ngaus)
       real*8  :: Mtmp(ngaus,ngaus)
       real*8  :: tScho(ngaus,ngaus)
       real*8  :: tiScho(ngaus,ngaus)
       real*8  :: Ptmp(ngaus,ngaus)
 !     External functions
       real*8,external  :: dnrm2
 !     Workspace for diagonalization. 
       real*8,  dimension(:), allocatable :: WORK, WORK2
       integer, dimension(:), allocatable :: IWORK, IWORK2
       integer :: LWORK, LIWORK, INFO
 
       format1='(16D11.3)'
       do nm=1,nvdf
          eval=0.0d0
          vec=0.0d0
          norm=0.0d0
          Ptmp=0.0d0
          Mtmp=0.0d0
          Mtrx(:,:)=Hcore(:,:,nm)
          tScho(:,:)=Scho(:,:,nm)
          tiScho(:,:)=IScho(:,:,nm)
 
 !        Computing transformed core hamiltonian.
          call dsygst(1,'U',ngaus,Mtrx,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
 
 !        ------------------------------------------------------------------
 !        DEBUG
 !        ------------------------------------------------------------------
 !         write(77,'(A)') 'PRINTING TRANSFORMED HAMILTONIAN'
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         do i=1,ngaus
 !            write(77,format1) Mtrx(i,:)
 !         end do
 !        ------------------------------------------------------------------
 
 !        Diagonalizing transformed core hamiltonian.
          allocate ( WORK(1000) )
          LWORK=-1
          call dsyev('V','U',ngaus,Mtrx,ngaus,eval,WORK,LWORK,INFO)
          LWORK=WORK(1)
          if(allocated(WORK2)) deallocate (WORK2)
          allocate (WORK2(LWORK))
          call dsyev('V','U',ngaus,Mtrx,ngaus,eval,WORK2,LWORK,INFO)
          deallocate( WORK, WORK2 )
          if (INFO /= 0) STOP ('Error during diagonalization')
 
 !        Normalize eigenvectors
          do i=1,ngaus
             vec=Mtrx(:,i)
             norm=dnrm2(ngaus,vec,1)
             vec=vec/norm
             Mtmp(:,i)=vec
          end do
          
 
 !        Convert coeficients P into the non-orthonormal basis representation.
 !        C=S**(-1/2) . C'
          Ptmp=Mtmp
          call dtrmm('L','U','N','N',ngaus,ngaus,1.0D0,tiScho,ngaus,Ptmp,ngaus)
 
 !        Normalize eigenvectors
          do i=1,ngaus
             vec=Ptmp(:,i)
             norm=dnrm2(ngaus,vec,1)
             vec=vec/norm
             Ptmp(:,i)=vec
          end do
 
 !        Setting up a consistent phase between eigenvectors.
          if (mod(ngaus,2) == 0) then
             halfg=ngaus/2
          else if (mod(ngaus,2) == 1) then
             x=ngaus/2
             halfg=floor(x)+1
          else
             STOP ('ERROR DURING PHASE CONSISTENCY CHECK')
          end if
          do ex=1,ngaus
             csum=0.0d0
             do i=1,halfg
                csum=csum+Mtmp(i,ex)
             end do
             if ((csum < 0.0d0) .AND. (mod(ex,2) == 1)) then
                Mtmp(:,ex)=-Mtmp(:,ex)
                Ptmp(:,ex)=-Ptmp(:,ex)
             else if ((csum > 0.0d0) .AND. (mod(ex,2) == 0)) then
                Mtmp(:,ex)=-Mtmp(:,ex)
                Ptmp(:,ex)=-Ptmp(:,ex)
             end if
          end do
                
 !        ------------------------------------------------------------------
 !        DEBUG
 !        ------------------------------------------------------------------
  !        write(77,'(A)') 'PRINTING TRANSFORMED COEFICIENTS'
  !        write(77,'(A,I2)') 'normal mode ',nm
  !        do i=1,ngaus
  !           write(77,format1) Mtrx(i,:)
  !        end do 
 !
 !         write(77,'(A)') 'PRINTING BACK-TRANSFORMED COEFICIENTS'
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         do i=1,ngaus
 !            write(77,format1) Ptmp(i,:)
 !         end do
 !
 !         write(77,'(A)') 'PRINTING MODAL ENERGIES'
 !         write(77,'(A,I2)') 'normal mode ',nm
 !         write(77,format1) eval
 !        ------------------------------------------------------------------
 
 
          Po(:,:,nm)=Mtmp(:,:)
          P(:,:,nm)=Ptmp(:,:)
          Emod(:,nm)=eval(:)
       end do
       
       end subroutine
 
 
       subroutine calcVeff(Po,Q1,Q2,Q3,Scho,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk,nvdf,ngaus,Gmtrx,GDmtrx,GTmtrx)
 !     --------------------------------------------------------------
 !     COMPUTE VSCF EFFECTIVE POTENTIAL MATRIX "G"
 !     --------------------------------------------------------------
       implicit none
       integer,intent(in) :: nvdf
       integer,intent(in) :: ngaus
       real*8,intent(in)  :: Po(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q1(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q2(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q3(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Scho(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: tiij(nvdf,nvdf)
       real*8,intent(in)  :: tjji(nvdf,nvdf)
       real*8,intent(in)  :: uiiij(nvdf,nvdf)
       real*8,intent(in)  :: ujjji(nvdf,nvdf)
       real*8,intent(in)  :: uiijj(nvdf,nvdf)
       real*8,intent(in)  :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)  :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)  :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)  :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: Gmtrx(ngaus,ngaus,nvdf)
       real*8,intent(out) :: GDmtrx(ngaus,ngaus,nvdf)
       real*8,intent(out) :: GTmtrx(ngaus,ngaus,nvdf)
 
       integer :: nm,nm1,nm2,nm3,n1,n2,n3
       integer :: INFO
       real*8  :: Coef1
       real*8  :: Coef2
       real*8  :: Coef3
       real*8  :: tmp(ngaus)
       real*8  :: Q1avg(nvdf)
       real*8  :: Q2avg(nvdf)
       real*8  :: Q3avg(nvdf)
       real*8  :: tScho(ngaus,ngaus)
       real*8  :: Qtmp(ngaus,ngaus)
       real*8  :: Pt(ngaus)
       real*8,external :: ddot
 !     --------------------------------------------------------------
 
 !     ------------------------------------------------------------------
 !     COMPUTING AVERAGES OF QM OPERATORS OVER MODALS
 !     ------------------------------------------------------------------
 !     Computing average normal coordinate operators over ground state modals
 !                           <phi_a|(Q_a)**n|phi_a>
       do nm=1,nvdf
          Pt=Po(:,1,nm)
          tScho=Scho(:,:,nm)
 
          tmp=0.0d0
          Qtmp=Q1(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Pt,1,0.0D0,tmp,1)
          Q1avg(nm) = ddot(ngaus,Pt,1,tmp,1)
 
          tmp=0.0d0
          Qtmp=Q2(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Pt,1,0.0D0,tmp,1)
          Q2avg(nm) = ddot(ngaus,Pt,1,tmp,1)
 
          tmp=0.0d0
          Qtmp=Q3(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Pt,1,0.0D0,tmp,1)
          Q3avg(nm) = ddot(ngaus,Pt,1,tmp,1)
       end do
 
 !     Computing vscf effective potential matrix G
 !     ------------------------------------------------------------------
 !     BEGINING COMPUTATION OF 2 MODE COUPLING TERMS
 !     ------------------------------------------------------------------
       GDmtrx=0.0d0
       GTmtrx=0.0d0
       do nm1=1,nvdf
          do nm2=1,nvdf
             if (nm1 == nm2) CYCLE
             if (nm1 > nm2) then
                n1=nm2
                n2=nm1
 
                Coef1=0.0d0
                Coef2=0.0d0
                Coef3=0.0d0
 
                Coef1=tiij(n1,n2)*Q2avg(n1)/2d0+uiiij(n1,n2)*Q3avg(n1)/6d0
                Coef2=tjji(n1,n2)*Q1avg(n1)/2d0+uiijj(n1,n2)*Q2avg(n1)/4d0
                Coef3=ujjji(n1,n2)*Q1avg(n1)/6d0
 
                GDmtrx(:,:,nm1)=GDmtrx(:,:,n2)+&
                             &Coef1*Q1(:,:,n2)+&
                             &Coef2*Q2(:,:,n2)+&
                             &Coef3*Q3(:,:,n2)
             else
                n1=nm1
                n2=nm2
 
                Coef1=0.0d0
                Coef2=0.0d0
                Coef3=0.0d0
 
                Coef1=tjji(n1,n2)*Q2avg(n2)/2d0+ujjji(n1,n2)*Q3avg(n2)/6d0
                Coef2=tiij(n1,n2)*Q1avg(n2)/2d0+uiijj(n1,n2)*Q2avg(n2)/4d0
                Coef3=uiiij(n1,n2)*Q1avg(n2)/6d0
 
                GDmtrx(:,:,nm1)=GDmtrx(:,:,n1)+&
                             &Coef1*Q1(:,:,n1)+&
                             &Coef2*Q2(:,:,n1)+&
                             &Coef3*Q3(:,:,n1)
             end if
          end do
       end do
 
 !     ------------------------------------------------------------------
 !     BEGINING COMPUTATION OF 3 MODE COUPLING TERMS
 !     ------------------------------------------------------------------
       do nm1=1,nvdf
          do nm2=1,nvdf-1
             do nm3=nm2+1,nvdf
                if (nm1 == nm2) CYCLE
                if (nm1 == nm3) CYCLE
                if (nm1<nm2) then
                   n1=nm1
                   n2=nm2
                   n3=nm3
 
                   coef1=0d0 
                   coef2=0d0
    
                   coef1=tijk(n1,n2,n3)*Q1avg(n2)*Q1avg(n3)+&
                        &uijjk(n1,n2,n3)*Q2avg(n2)*Q1avg(n3)/2d0+&
                        &uijkk(n1,n2,n3)*Q1avg(n2)*Q2avg(n3)/2d0
    
                   coef2=uiijk(n1,n2,n3)*Q1avg(n2)*Q1avg(n3)/2d0
    
                   GTmtrx(:,:,n1)=GTmtrx(:,:,n1)+&
                               &coef1*Q1(:,:,n1)+&
                               &coef2*Q2(:,:,n1)
 
                else if (nm1>nm2 .AND. nm1<nm3) then
                   n1=nm2
                   n2=nm1
                   n3=nm3
 
                   coef1=0d0 
                   coef2=0d0
    
                   coef1=tijk(n1,n2,n3)*Q1avg(n1)*Q1avg(n3)+&
                        &uiijk(n1,n2,n3)*Q2avg(n1)*Q1avg(n3)/2d0+&
                        &uijkk(n1,n2,n3)*Q1avg(n1)*Q2avg(n3)/2d0
    
                   coef2=uijjk(n1,n2,n3)*Q1avg(n1)*Q1avg(n3)/2d0
    
                   GTmtrx(:,:,n2)=GTmtrx(:,:,n2)+&
                               &coef1*Q1(:,:,n2)+&
                               &coef2*Q2(:,:,n2)
                else if (nm1>nm3) then
                   n1=nm2
                   n2=nm3
                   n3=nm1
 
                   coef1=0d0 
                   coef2=0d0
    
                   coef1=tijk(n1,n2,n3)*Q1avg(n1)*Q1avg(n2)+&
                        &uiijk(n1,n2,n3)*Q2avg(n1)*Q1avg(n2)/2d0+&
                        &uijjk(n1,n2,n3)*Q1avg(n1)*Q2avg(n2)/2d0
    
                   coef2=uijkk(n1,n2,n3)*Q1avg(n1)*Q1avg(n2)/2d0
    
                   GTmtrx(:,:,n3)=GTmtrx(:,:,n3)+&
                               &coef1*Q1(:,:,n3)+&
                               &coef2*Q2(:,:,n3)
 
                end if
                
             end do
          end do
       end do 
 
       Gmtrx=GDmtrx+GTmtrx
                   
       end subroutine
 
       subroutine transformMtrx(Ssqrt,Gmtrx,nvdf,ngaus)
 !     ------------------------------------------------------------------
 !     TRANSFORMS MATRIX M INTO M' = S**(-1/2) . M . S**(-1/2)
 !     ------------------------------------------------------------------
       implicit none
 
       integer,intent(in)    :: nvdf,ngaus
       real*8,intent(in)     :: Ssqrt(ngaus,ngaus,nvdf)
       real*8,intent(inout)  :: Gmtrx(ngaus,ngaus,nvdf)
 
       integer :: nm
 !     ------------------------------------------------------------------
       do nm=1,nvdf
 !        For each normal mode, compute S**-1/2 . Gmtrx . S**-1/2
          call dgemm('N','N',ngaus,ngaus,ngaus,1.0D0,Ssqrt(:,:,nm),       &
                  & ngaus,Gmtrx(:,:,nm),ngaus,0.0D0,Gmtrx(:,:,nm),ngaus)
          call dgemm('N','N',ngaus,ngaus,ngaus,1.0D0,Gmtrx(:,:,nm),       &
                  & ngaus,Ssqrt(:,:,nm),ngaus,0.0D0,Gmtrx(:,:,nm),ngaus)
       end do
 
       end subroutine
 
 
       subroutine calcEvirt2(Emod,Po,Scho,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
                  & tijk,uiijk,uijjk,uijkk,ngaus,nvdf,Evirt1)
 !     --------------------------------------------------------------
 !     COMPUTE VSCF VIRTUAL STATES ENERGIES
 !     Only single quanta excited states for later computation of 
 !     fundamental transition energies.
 !     --------------------------------------------------------------
       implicit none
       integer,intent(in) :: nvdf
       integer,intent(in) :: ngaus
       real*8,intent(in)  :: Emod(ngaus,nvdf)
       real*8,intent(in)  :: Po(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Scho(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q1(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q2(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q3(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: tiij(nvdf,nvdf)    
       real*8,intent(in)  :: tjji(nvdf,nvdf)   
       real*8,intent(in)  :: uiiij(nvdf,nvdf) 
       real*8,intent(in)  :: ujjji(nvdf,nvdf)     
       real*8,intent(in)  :: uiijj(nvdf,nvdf)     
       real*8,intent(in)  :: tijk(nvdf,nvdf,nvdf)  
       real*8,intent(in)  :: uiijk(nvdf,nvdf,nvdf) 
       real*8,intent(in)  :: uijjk(nvdf,nvdf,nvdf) 
       real*8,intent(in)  :: uijkk(nvdf,nvdf,nvdf) 
       real*8,intent(out) :: Evirt1(nvdf)
 
       integer :: ex
       integer :: nm
       real*8  :: corr
       real*8  :: corrD
       real*8  :: corrT
 !     --------------------------------------------------------------
 
       Evirt1=0.0D0
       do ex=1,nvdf
          do nm=1,nvdf
             if (nm == ex) then
                Evirt1(ex)=Evirt1(ex)+Emod(2,nm)
             else
                Evirt1(ex)=Evirt1(ex)+Emod(1,nm)
             end if
          end do
       end do
 
       do ex=1,nvdf
          corr=0.0D0
 !        DANGER
 !        call vscfcorr2(ex,Scho,Po,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
 !                 & tijk,uiijk,uijjk,uijkk,nvdf,ngaus,corr,corrD,corrT)
          call vscfcorr2(0,Scho,Po,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
                   & tijk,uiijk,uijjk,uijkk,nvdf,ngaus,corr,corrD,corrT)
          Evirt1(ex)=Evirt1(ex) - corr
       end do
 
       end subroutine
 
       subroutine calcEvirt(Po,Scho,GDmtrx,GTmtrx,Emod,ngaus,nvdf,Evirt1)
 !     --------------------------------------------------------------
 !     COMPUTE VSCF VIRTUAL STATES ENERGIES
 !     Only single quanta excited states for later computation of 
 !     fundamental transition energies.
 !     --------------------------------------------------------------
       implicit none
       integer,intent(in) :: nvdf
       integer,intent(in) :: ngaus
       real*8,intent(in)  :: Po(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Scho(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: GDmtrx(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: GTmtrx(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Emod(ngaus,nvdf)
       real*8,intent(out) :: Evirt1(nvdf)
 
       integer :: ex
       integer :: nm
       integer :: INFO
       real*8  :: corr
       real*8  :: corrD
       real*8  :: corrT
       real*8  :: tmp(ngaus)
       real*8  :: GDtmp(ngaus,ngaus)
       real*8  :: GTtmp(ngaus,ngaus)
       real*8  :: tScho(ngaus,ngaus)
       real*8  :: Ptmp(ngaus)
       real*8,external  :: ddot
 !     --------------------------------------------------------------
 
       Evirt1=0.0D0
       do ex=1,nvdf
          do nm=1,nvdf
             if (nm == ex) then
                Evirt1(ex)=Evirt1(ex)+Emod(2,nm)
             else
                Evirt1(ex)=Evirt1(ex)+Emod(1,nm)
             end if
          end do
       end do
 
       Ptmp=0.0d0
       GDtmp=0.0d0
       GTtmp=0.0d0
       do ex=1,nvdf
          corr=0.0D0
          corrD=0.0D0
          corrT=0.0D0
          do nm=1,nvdf
             tScho(:,:) = Scho(:,:,nm)
             GDtmp(:,:) = GDmtrx(:,:,nm)
             GTtmp(:,:) = GTmtrx(:,:,nm)
 
             call dsygst(1,'U',ngaus,GDtmp,ngaus,tScho,ngaus,INFO)
             if (INFO /= 0) STOP('Error during Matrix transformation')
             call dsygst(1,'U',ngaus,GTtmp,ngaus,tScho,ngaus,INFO)
             if (INFO /= 0) STOP('Error during Matrix transformation')
 
             if (nm == ex) then
                Ptmp(:) = Po(:,2,nm)
             else
                Ptmp(:) = Po(:,1,nm)
             end if
 
             call dsymv('U',ngaus,1.0D0,GDtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
             corrD = corrD + ddot(ngaus,Ptmp,1,tmp,1)
 
             call dsymv('U',ngaus,1.0D0,GTtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
             corrT = corrT + ddot(ngaus,Ptmp,1,tmp,1)
 
          end do
 
          Evirt1(ex)=Evirt1(ex) - 0.5d0*corrD - (2.0d0*corrT/3.0d0)
       end do
 
       end subroutine
 
       subroutine vscfcorr2(ex,Scho,Po,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
                   & tijk,uiijk,uijjk,uijkk,nvdf,ngaus,corr,corrD,corrT)
 !     ------------------------------------------------------------------------ 
 !     COMPUTES VSCF ENERGY CORRECTION
 !     CORR = <PSI|VcD|PSI> + 2 <PSI|VcT|PSI>
 !     WHERE VcD AND VcT ARE THE DOUBLES AND TRIPLES COUPLING POTENTIAL,
 !     RESPECTIVELY.
 !     ------------------------------------------------------------------------ 
       implicit none
 
       integer,intent(in)     :: ex
       integer,intent(in)     :: nvdf
       integer,intent(in)     :: ngaus
       real*8,intent(in)      :: Scho(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Po(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Q1(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Q2(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Q3(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: tiij(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: tjji(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uiiij(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: ujjji(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uiijj(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: tijk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uiijk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uijjk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uijkk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(out)     :: corr
       real*8,intent(out)     :: corrD
       real*8,intent(out)     :: corrT
       
       integer     :: nm
       integer     :: a,b,c
       integer     :: INFO                   ! Error code for dsygst subroutine.
       real*8      :: Qtmp(ngaus,ngaus)
       real*8      :: tScho(ngaus,ngaus)
       real*8      :: Ptmp(ngaus)
       real*8      :: tmp(ngaus)
       real*8      :: Qm1(nvdf)
       real*8      :: Qm2(nvdf)
       real*8      :: Qm3(nvdf)
 
       real*8,external :: ddot
 !     DEBUG
       integer     :: i,j
 
       corr=0.0D0
       do nm=1,nvdf
          tScho(:,:) = Scho(:,:,nm)
          if (nm == ex) then
             Ptmp(:) = Po(:,2,nm)
          else
             Ptmp(:) = Po(:,1,nm)
          end if
 
          tmp=0.0d0
          Qtmp(:,:) = 0d0
          Qtmp(:,:) = Q1(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
          Qm1(nm) = ddot(ngaus,Ptmp,1,tmp,1)
 
          tmp=0.0d0
          Qtmp(:,:) = 0d0
          Qtmp(:,:) = Q2(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
          Qm2(nm) = ddot(ngaus,Ptmp,1,tmp,1)
 
          tmp=0.0d0
          Qtmp(:,:) = 0d0
          Qtmp(:,:) = Q3(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
          Qm3(nm) = ddot(ngaus,Ptmp,1,tmp,1)
       end do
 
 !     2-MODE COUPLING POTENTIAL
       corrD=0.0d0
       do a=1,nvdf-1
          do b=a+1,nvdf
             corrD=corrD+&
             & tiij(a,b)*Qm2(a)*Qm1(b)/2.0d0+&
             & tjji(a,b)*Qm1(a)*Qm2(b)/2.0d0+&
             & uiiij(a,b)*Qm3(a)*Qm1(b)/6.0d0+&
             & ujjji(a,b)*Qm1(a)*Qm3(b)/6.0d0+&
             & uiijj(a,b)*Qm2(a)*Qm2(b)/4.0d0  
          end do
       end do
 
 !     3-MODE COUPLING POTENTIAL
       corrT=0d0
       do a=1,nvdf-2
          do b=a+1,nvdf-1
             do c=b+1,nvdf
                corrT=corrT+&
                 & tijk(a,b,c)*Qm1(a)*Qm1(b)*Qm1(c)+&
                 & uiijk(a,b,c)*Qm2(a)*Qm1(b)*Qm1(c)/2d0+&
                 & uijjk(a,b,c)*Qm1(a)*Qm2(b)*Qm1(c)/2d0+&
                 & uijkk(a,b,c)*Qm1(a)*Qm1(b)*Qm2(c)/2d0
             end do
          end do
       end do
 
       corr = corrD + 2d0*corrT
       end subroutine
 
 
       subroutine vscfcorr(Gmtrx,Scho,Po,nvdf,ngaus,corr)
 !     ------------------------------------------------------------------------ 
 !     COMPUTES VSCF ENERGY CORRECTION 
 !     ------------------------------------------------------------------------ 
       implicit none
       integer,intent(in)     :: nvdf
       integer,intent(in)     :: ngaus
       real*8,intent(in)      :: Gmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Scho(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Po(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(out)     :: corr
       
       integer     :: nm
       integer     :: INFO                   ! Error code for dsygst subroutine.
       real*8      :: Gtmp(ngaus,ngaus)
       real*8      :: tScho(ngaus,ngaus)
       real*8      :: Ptmp(ngaus)
       real*8      :: tmp(ngaus)
 
       real*8,external :: ddot
       
 !     DEBUG
       integer     :: i,j
       
       corr=0.0D0
       do nm=1,nvdf
          tmp=0.0d0
          Ptmp(:) = Po(:,1,nm)
          Gtmp(:,:) = Gmtrx(:,:,nm)
          tScho(:,:) = Scho(:,:,nm)
          call dsygst(1,'U',ngaus,Gtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymv('U',ngaus,1.0D0,Gtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
          corr = corr + ddot(ngaus,Ptmp,1,tmp,1)
 !        -----------------------------------------------------------
 !        DEBUG
 !        -----------------------------------------------------------
          if (nm==1) write(77,'(A)') '<veffD>vscf'
          write(77,'(D12.3)') ddot(ngaus,Ptmp,1,tmp,1)
 !        -----------------------------------------------------------
       end do
 
       end subroutine
 
 
       subroutine vscf(nmodes,eig,dy,nmcoup,gwidth,nqmatoms,nclatoms,&
                    &ngaus,qmcoords,clcoords,at_numbers,ndf,nvdf,&
                    &P,Po,Scho,Evscf,Evirt1,Emod,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,&
                    &hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk)
 
       implicit none
 
       integer,intent(in)  :: ndf                  ! Number of classical atoms
       integer,intent(in)  :: nvdf                 ! Number of classical atoms
       integer,intent(in)  :: ngaus                ! Number of classical atoms
       integer,intent(in)  :: nmcoup               ! Number of normal modes to couple in QFF.
       integer,intent(in)  :: nqmatoms             ! Number of QM atoms
       integer,intent(in)  :: nclatoms             ! Number of classical atoms
       integer,intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
       real*8, intent(in)  :: dy                   ! Step size factor for num derivatives.
       real*8, intent(in)  :: gwidth               ! Width factor for gaussian primitives.
       real*8, intent(in)  :: eig(ndf)             ! Eigenvalues of the hessian matrix.
       real*8, intent(in)  :: nmodes(ndf,ndf)      ! Eigenvectors of the hessian matrix.
       real*8, intent(in)  :: qmcoords(3,nqmatoms) ! QM atom coordinates
       real*8, intent(in)  :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
       real*8, intent(out) :: P(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
       real*8, intent(out) :: Po(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
       real*8, intent(out) :: Scho(ngaus,ngaus,nvdf)  ! Cholesky factor.
       real*8, intent(out) :: Evscf                ! VSCF energy.
       real*8, intent(out) :: Evirt1(nvdf)         ! VSCF 1st excited virt. states energies.
       real*8, intent(out) :: Emod(ngaus,nvdf)     ! Modal energies.
       real*8, intent(out) :: Q1(ngaus,ngaus,nvdf) !  <gi|Qa^1|gj> matrix
       real*8, intent(out) :: Q2(ngaus,ngaus,nvdf) !  <gi|Qa^2|gj> matrix
       real*8, intent(out) :: Q3(ngaus,ngaus,nvdf) !  <gi|Qa^3|gj> matrix
       real*8, intent(out) :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian (non-orthogonal basis)
       real*8, intent(out) :: GDmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8, intent(out) :: GTmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8, intent(in) :: hii(ndf)             ! Diagonal Hessian eigenvalues.
       real*8, intent(in) :: tiii(nvdf)           ! Diagonal cubic coupling terms
       real*8, intent(in) :: uiiii(nvdf)          ! Diagonal quartic coupling terms
       real*8, intent(in) :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
       real*8, intent(in) :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
       real*8, intent(in) :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8, intent(in) :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8, intent(in) :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8, intent(in) :: tijk(nvdf,nvdf,nvdf)  ! 3-mode quartic coupling trms.
       real*8, intent(in) :: uiijk(nvdf,nvdf,nvdf) ! 3-mode quartic coupling trms.
       real*8, intent(in) :: uijjk(nvdf,nvdf,nvdf) ! 3-mode quartic coupling trms.
       real*8, intent(in) :: uijkk(nvdf,nvdf,nvdf) ! 3-mode quartic coupling trms.
 !     LOCAL VARIABLES      
       integer     :: INFO                   ! Error code for dsygst subroutine.
       integer     :: iter                   ! Multipurpose indexes.
       integer     :: i,j,k                  ! Multipurpose indexes.
       integer     :: nm1,nm2,nm,at,ex       ! Multipurpose indexes.
       real*8      :: nmdes(ndf,ndf)         ! Eigenvectors of the hessian matrix.
       real*8      :: X0(ndf)                ! Initial geometry in cartesian coords and AU.
       real*8      :: X0m(ndf)               ! Initial geometry in mass-weighted cart coords (AU).
       real*8      :: Q0(nvdf)               ! Initial geomtetry in normal coords and AU.
       real*8      :: Mass(ndf)                 ! Sqrt(mass_i) diagonal matrix.
       real*8      :: L(ndf,nvdf)            ! Normal mode matrix.
       real*8      :: ll(ndf)                ! Dummy variable. Normal mode vector
       real*8      :: S(ngaus,ngaus,nvdf)    ! Overlap matrix.
       real*8      :: IScho(ngaus,ngaus,nvdf)! Inverse of Cholesky factor of S.
       real*8      :: Evscf_old              ! VSCF energy.
       real*8      :: fnd_trans(nvdf)        ! VSCF virtual energy.
       real*8      :: Gmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8      :: F(ngaus,ngaus,nvdf)    ! Effective potential matrix.
       real*8      :: tmp(ngaus)             ! Auxiliary variable for computing vef
       real*8      :: corr                   ! Auxiliary variable for computing vef
       real*8      :: corrD                  ! Auxiliary variable for computing vef
       real*8      :: corrT                  ! Auxiliary variable for computing vef
       real*8      :: Gtmp(ngaus,ngaus)
       real*8      :: tScho(ngaus,ngaus)
       real*8      :: Ptmp(ngaus)
       real*8      :: Ediag(nvdf+1)
       real*8      :: Evhf(nvdf+1)
       real*8      :: harm
 !     -----------------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------------
       real*8      :: Vdavg    !  <psi_0|Vc2MC|psi_0>
       real*8      :: Vtavg    !  <psi_0|Vc3MC|psi_0>
       real*8      :: veff2avg ! <psi_0|veff2MC|psi_0>
       real*8      :: veff3avg ! <psi_0|veff3MC|psi_0>
       real*8      :: vtmp
 !     -----------------------------------------------------------------
 !     FUNCTIONS
       real*8,external :: ddot 
 !     PARAMETERS
       real*8, parameter :: hartree2cm = 219474.63d0 ! cm-1/Ha
       include "qvmbia_param.f"
 
 !     ngaus=16
 !     ndf  = 3*nqmatoms
 !     nvdf = ndf-6
 !     dy=0.5d0
 !----------------------------------------------------------------------
 !     INITIAL GUESS
 !----------------------------------------------------------------------
 !     Eliminating rotational and translational degrees of freedom.
 !     DANGER
 !      do nm=7,ndf
 !         hii(nm-6)=eig(nm)
 !      end do
 
 !     Build Core Hamiltonian and Qi^n matrices.
       call buildoperators(hii,tiii,uiiii,gwidth,ngaus,nvdf,S,Q1,Q2,Q3,Hcore)
       
 !     Compute Cholesky decomposition and the inverse of the cholesky
 !     factor.
       call cholesky(S,nvdf,ngaus,Scho,IScho)
 
 
 !     Transform core hamiltonian using IScho and diagonalize. Return
 !     transformed back eigenvectors coeficients and eigenvalues.
       call solveGEV(Hcore,Scho,IScho,nvdf,ngaus,P,Po,Emod)
      
       Evscf=0.0d0
       Ediag=0.0d0
       do nm=1,nvdf
          Evscf = Evscf + Emod(1,nm)
       end do
       Ediag(1)=Evscf
 
       do ex=1,nvdf
          do nm=1,nvdf
             if (nm == ex) then
                 Evirt1(ex) = Evirt1(ex) + Emod(2,nm)
             else
                 Evirt1(ex) = Evirt1(ex) + Emod(1,nm)
             end if
             Ediag(ex+1)=Evirt1(ex)
          end do
       end do
       do nm=1,nvdf
          fnd_trans(nm) = Evirt1(nm) - Evscf
       end do
       fnd_trans = fnd_trans*hartree2cm
       write(77,'(A)') 'INITIAL GUESS ZERO POINT ENERGY'
       write(77,'(F14.6)') Evscf*hartree2cm
       write(77,'(A)') 'INITIAL GUESS FUNDAMENTAL TRANSITIONS'
       write(77,'(F14.6)') fnd_trans
 !     ------------------------------------------------------------------
             
 
 !----------------------------------------------------------------------
 !     VSCF ITERATIONS
 !----------------------------------------------------------------------
 !     BEGINING VSCF ITERATIONS
 
       iter=1
       write(77,'(A)') 'ITER     Evscf   MP1 CORRECTION'
       write(77,'(A)') '----------------------------------'
       DO ! Begin infinite loop.
 
 !        Build effective potential matrix G
          Gmtrx=0.0d0
          Emod=0.0d0
          call calcVeff(Po,Q1,Q2,Q3,Scho,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk,nvdf,ngaus,Gmtrx,GDmtrx,GTmtrx)
 
 !        Build fock matrix.
          F = Hcore + Gmtrx
 
 !        Transform core hamiltonian using IScho and diagonalize. Return
 !        transformed back eigenvectors coeficients and eigenvalues.
          Evscf_old=Evscf
          call solveGEV(F,Scho,IScho,nvdf,ngaus,P,Po,Emod)
 
 !        Compute VSCF energy for this iteration.
 
          Evscf=0.0D0
          do nm=1,nvdf
             Evscf=Evscf+Emod(1,nm)
          end do
 
 !        Compute VSCF correction for double counting.
          call vscfcorr2(0,Scho,Po,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
                   & tijk,uiijk,uijjk,uijkk,nvdf,ngaus,corr,corrD,corrT)
 !     -----------------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------------
          Vdavg=corrD
          Vtavg=corrT
 !     -----------------------------------------------------------------
 !        DANGER
 !         call vscfcorr(GDmtrx,Scho,Po,nvdf,ngaus,corr)
 !         corrD=corr
 !         call vscfcorr(GTmtrx,Scho,Po,nvdf,ngaus,corr)
 !         corrT=corr
 !         corr=0.5d0*corrD+2d0*corrT/3d0
 
 !        Write to file iteration information.
          write(77,'(I3,F13.8,2F15.8)') iter,Evscf,corrD,corrT
 
          Evscf = Evscf - corr
          if (iter >= 60) STOP('VSCF DID NOT CONVERGE: Too many iterations.')
          iter=iter+1
 
 !        TEST FOR CONVERGENCE.
          if (Abs(Evscf_old-Evscf) < 1.0D-6) EXIT
       END DO
 
       write(77,'(A)')      '------------------------------------'
       write(77,'(A,I3,A)') 'VSCF CONVERGED IN ',iter-1,' ITERATIONS'
       write(77,'(A)')      '------------------------------------'
 
 !     Compute singly excited virtual energies.
 !     DANGER
 !      call calcEvirt(Po,Scho,GDmtrx,GTmtrx,Emod,ngaus,nvdf,Evirt1)
       call calcEvirt2(Emod,Po,Scho,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
                  & tijk,uiijk,uijjk,uijkk,ngaus,nvdf,Evirt1)
 
 !     Fundamental Transitions
 !     -----------------------------------------------------------------
       Evhf(1)=Evscf
       do nm=1,nvdf
          fnd_trans(nm)=Evirt1(nm)-Evscf
          Evhf(nm+1)=Evirt1(nm)
       end do
 
 !     CONVERTING TO CM**-1
 !     -----------------------------------------------------------------
       fnd_trans=fnd_trans*hartree2cm
       Evscf=Evscf*hartree2cm
 !     -----------------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------------
       write(79,*) 
       write(79,'(A)') 'VSCF P MATRIX in vscf main'
       do i=1,ngaus
          write(79,'(99D11.3)') (P(i,j,3),j=1,ngaus)
       end do
       write(79,*) 
       write(79,'(A)') 'VSCF Po MATRIX in vscf main'
       do i=1,ngaus
          write(79,'(99D11.3)') (Po(i,j,3),j=1,ngaus)
       end do
 !     -----------------------------------------------------------------
 !     -----------------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------------
          vtmp=0d0
          veff2avg=0d0
          veff3avg=0d0
          call vscfcorr(GDmtrx,Scho,Po,nvdf,ngaus,vtmp)
          veff2avg=vtmp
          call vscfcorr(GTmtrx,Scho,Po,nvdf,ngaus,vtmp)
          veff3avg=vtmp
       write(77,'(A,F16.9)') '<psi0|veffd|psi0> = ',veff2avg 
       write(77,'(A,F16.9)') '<psi0|vcd|psi0> = ',Vdavg 
       write(77,'(A,F16.9)') '<psi0|vefft|psi0> = ',veff3avg 
       write(77,'(A,F16.9)') '<psi0|vct|psi0> = ',Vtavg 
 !     -----------------------------------------------------------------
 
 !     PRINT VIRTUAL TRANSITIONS
 !     -----------------------------------------------------------------
       write(77,'(A)') ' MODE      HARMONIC      DIAGONAL         VSCF '
       write(77,'(A)') '-----------------------------------------------'
       do nm=1,nvdf
          harm = sqrt(abs(hii(nm)))*hartree2cm
          write(77,'(I4,99F14.2)') nm,harm,       &
                                 & (Ediag(nm+1)-Ediag(1))*hartree2cm, &
                                 & (Evhf(nm+1)-Evhf(1))*hartree2cm
       end do
       write(77,'(A,99F14.2)') 'ZERO P. ENERGIES  ',  &
                               & Ediag(1)*hartree2cm,   &
                               & Evhf(1)*hartree2cm
 
       end subroutine
 
 
 !#######################################################################
 !     STATE-SPECIFIC VSCF (SSVSCF) SECTION
 !#######################################################################
 
       subroutine ssvscf2(ref,nmodes,eig,dy,nmcoup,gwidth,nqmatoms,nclatoms,&
                    &ngaus,qmcoords,clcoords,at_numbers,ndf,nvdf,&
                    &P,Po,Scho,Evscf,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Emod,&
                    &hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk,vscfstat)
 
       implicit none
 
       integer,intent(in)  :: ref(4)               ! Reference state          
       integer,intent(in)  :: ndf                  ! Number of classical atoms
       integer,intent(in)  :: nvdf                 ! Number of classical atoms
       integer,intent(in)  :: ngaus                ! Number of classical atoms
       integer,intent(in)  :: nmcoup               ! Number of normal modes to couple in QFF.
       integer,intent(in)  :: nqmatoms             ! Number of QM atoms
       integer,intent(in)  :: nclatoms             ! Number of classical atoms
       integer,intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
       integer, intent(out):: vscfstat
       real*8, intent(in)  :: dy                   ! Step size factor for num derivatives.
       real*8, intent(in)  :: gwidth               ! Width factor for gaussian primitives.
       real*8, intent(in)  :: eig(ndf)             ! Eigenvalues of the hessian matrix.
       real*8, intent(in)  :: nmodes(ndf,ndf)      ! Eigenvectors of the hessian matrix.
       real*8, intent(in)  :: qmcoords(3,nqmatoms) ! QM atom coordinates
       real*8, intent(in)  :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
       real*8, intent(in)  :: hii(ndf)             ! Diagonal Hessian eigenvalues.
       real*8, intent(in)  :: tiii(nvdf)           ! Diagonal cubic coupling terms
       real*8, intent(in)  :: uiiii(nvdf)          ! Diagonal quartic coupling terms
       real*8, intent(in)  :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
       real*8, intent(in)  :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
       real*8, intent(in)  :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8, intent(in)  :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8, intent(in)  :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling terms
       real*8, intent(in)  :: tijk(nvdf,nvdf,nvdf)  ! 3-mode quartic coupling trms.
       real*8, intent(in)  :: uiijk(nvdf,nvdf,nvdf) ! 3-mode quartic coupling trms.
       real*8, intent(in)  :: uijjk(nvdf,nvdf,nvdf) ! 3-mode quartic coupling trms.
       real*8, intent(in)  :: uijkk(nvdf,nvdf,nvdf) ! 3-mode quartic coupling trms.
       real*8, intent(out) :: GDmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8, intent(out) :: GTmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8, intent(out) :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8, intent(out) :: P(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
       real*8, intent(out) :: Po(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
       real*8, intent(out) :: Scho(ngaus,ngaus,nvdf)  ! Cholesky factor.
       real*8, intent(out) :: Evscf                ! VSCF energy.
       real*8, intent(out) :: Q1(ngaus,ngaus,nvdf) !  <gi|Qa^1|gj> matrix
       real*8, intent(out) :: Q2(ngaus,ngaus,nvdf) !  <gi|Qa^2|gj> matrix
       real*8, intent(out) :: Q3(ngaus,ngaus,nvdf) !  <gi|Qa^3|gj> matrix
       real*8, intent(out) :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian (non-orthogonal basis)
 !     LOCAL VARIABLES      
       integer     :: psi(nvdf)            ! Reference state          
       integer     :: INFO                   ! Error code for dsygst subroutine.
       integer     :: iter                   ! Multipurpose indexes.
       integer     :: i,j,k                  ! Multipurpose indexes.
       integer     :: nm1,nm2,nm,at,ex       ! Multipurpose indexes.
       integer     :: qn1,qn2
       real*8      :: nmdes(ndf,ndf)         ! Eigenvectors of the hessian matrix.
       real*8      :: X0(ndf)                ! Initial geometry in cartesian coords and AU.
       real*8      :: X0m(ndf)               ! Initial geometry in mass-weighted cart coords (AU).
       real*8      :: Q0(nvdf)               ! Initial geomtetry in normal coords and AU.
       real*8      :: Mass(ndf)                 ! Sqrt(mass_i) diagonal matrix.
       real*8      :: L(ndf,nvdf)            ! Normal mode matrix.
       real*8      :: ll(ndf)                ! Dummy variable. Normal mode vector
       real*8      :: S(ngaus,ngaus,nvdf)    ! Overlap matrix.
       real*8      :: IScho(ngaus,ngaus,nvdf)! Inverse of Cholesky factor of S.
       real*8      :: Evscf_old              ! VSCF energy.
       real*8      :: fnd_trans(nvdf)        ! VSCF virtual energy.
       real*8      :: Gmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8      :: F(ngaus,ngaus,nvdf)    ! Effective potential matrix.
       real*8      :: tmp(ngaus)             ! Auxiliary variable for computing vef
       real*8      :: corr                   ! Auxiliary variable for computing vef
       real*8      :: corrD                  ! Auxiliary variable for computing vef
       real*8      :: corrT                  ! Auxiliary variable for computing vef
       real*8      :: Gtmp(ngaus,ngaus)
       real*8      :: tScho(ngaus,ngaus)
       real*8      :: Ptmp(ngaus)
       real*8      :: Ediag(nvdf+1)
       real*8      :: Evhf(nvdf+1)
       real*8      :: harm
 !     -----------------------------------------------------------------
 !     DEBUG
 !     -----------------------------------------------------------------
       real*8      :: Vdavg    !  <psi_0|Vc2MC|psi_0>
       real*8      :: Vtavg    !  <psi_0|Vc3MC|psi_0>
       real*8      :: veff2avg ! <psi_0|veff2MC|psi_0>
       real*8      :: veff3avg ! <psi_0|veff3MC|psi_0>
       real*8      :: vtmp
 !     -----------------------------------------------------------------
 !     FUNCTIONS
       real*8,external :: ddot 
 !     PARAMETERS
       real*8, parameter :: hartree2cm = 219474.63d0 ! cm-1/Ha
       include "qvmbia_param.f"
 !----------------------------------------------------------------------
 !----------------------------------------------------------------------
 
 !     ngaus=16
 !     ndf  = 3*nqmatoms
 !     nvdf = ndf-6
 !     dy=0.5d0
 !     -----------------------------------------------------------------
 !     INITIAL GUESS
 !     -----------------------------------------------------------------
 !     Eliminating rotational and translational degrees of freedom.
 !     DANGER
 !      do nm=7,ndf
 !         hii(nm-6)=eig(nm)
 !      end do
       vscfstat=0
 
       nm1=ref(1)
       qn1=ref(2)
       nm2=ref(3)
       qn2=ref(4)
 
       psi=0.0d0
       if (nm1 /= 0) psi(nm1) = qn1
       if (nm2 /= 0) psi(nm2) = qn2
 
 !     Build Core Hamiltonian and Qi^n matrices.
       call buildoperators(hii,tiii,uiiii,gwidth,ngaus,nvdf,S,Q1,Q2,Q3,Hcore)
       
 !     Compute Cholesky decomposition and the inverse of the cholesky
 !     factor.
       call cholesky(S,nvdf,ngaus,Scho,IScho)
 
 
 !     Transform core hamiltonian using IScho and diagonalize. Return
 !     transformed back eigenvectors coeficients and eigenvalues.
       call solveGEV(Hcore,Scho,IScho,nvdf,ngaus,P,Po,Emod)
      
 
       do nm=1,nvdf
          if (nm == nm1) then
             Evscf = Evscf + Emod(qn1+1,nm)
          else if (nm == nm2) then
             Evscf = Evscf + Emod(qn2+1,nm)
          else
             Evscf = Evscf + Emod(1,nm)
          end if
       end do
 
       write(77,'(A,8I3)') 'INITIAL GUESS ENERGY FOR REFERENCE STATE',ref
       write(77,'(F14.6)') Evscf*hartree2cm
 
 
 !     -----------------------------------------------------------------
 !     VSCF ITERATIONS
 !     -----------------------------------------------------------------
 !     BEGINING VSCF ITERATIONS
 
       iter=1
       write(77,'(A)') 'ITER     Evscf   MP1 CORRECTION'
       write(77,'(A)') '----------------------------------'
       DO ! Begin infinite loop.
 
 !        Build effective potential matrix G
          Gmtrx=0.0d0
          Emod=0.0d0
          call ssVeff2(ref,Po,Q1,Q2,Q3,Scho,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk,nvdf,ngaus,Gmtrx,GDmtrx,GTmtrx)
 
 !        Build fock matrix.
          F = Hcore + Gmtrx
 
 !        Transform core hamiltonian using IScho and diagonalize. Return
 !        transformed back eigenvectors coeficients and eigenvalues.
          Evscf_old=Evscf
          call solveGEV(F,Scho,IScho,nvdf,ngaus,P,Po,Emod)
 
 !        Compute ssVSCF energy for ref state in this iteration.
 
          Evscf=0.0D0
          do nm=1,nvdf
             if (nm == nm1) then
                Evscf = Evscf + Emod(qn1+1,nm)
             else if (nm == nm2) then
                Evscf = Evscf + Emod(qn2+1,nm)
             else
                Evscf = Evscf + Emod(1,nm)
             end if
          end do
 
 
 !        Compute VSCF correction for double counting.
          call ssvscfcorr(ref,Scho,Po,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
                   & tijk,uiijk,uijjk,uijkk,nvdf,ngaus,corr,corrD,corrT)
 
 !        Write to file iteration information.
          write(77,'(I3,F13.8,2F15.8)') iter,Evscf,corrD,corrT
 
          Evscf = Evscf - corr
 !         if (iter >= 60) STOP('VSCF DID NOT CONVERGE: Too many iterations.')
          if (iter >= 60) then
             write(77,'(A)') 'VSCF DID NOT CONVERGE: Too many iterations.'
             vscfstat=1
             EXIT
          end if
          iter=iter+1
 
 !        TEST FOR CONVERGENCE.
          if (Abs(Evscf_old-Evscf) < 1.0D-9) EXIT
       END DO
 
       write(77,'(A)')      '------------------------------------'
       write(77,'(A,I3,A)') 'VSCF CONVERGED IN ',iter-1,' ITERATIONS'
       write(77,'(A)')      '------------------------------------'
 
       Evscf=Evscf*hartree2cm
       write(77,'(A)')      '------------------------------------'
       write(77,'(A,4I3)') 'ssVSCF ENERGY AT REF STATE',ref
       write(77,'(F16.2)') Evscf
       write(77,'(A)')      '------------------------------------'
 
       end subroutine
 
       subroutine ssVeff2(ref,Po,Q1,Q2,Q3,Scho,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk,nvdf,ngaus,Gmtrx,GDmtrx,GTmtrx)
 !     --------------------------------------------------------------
 !     COMPUTE VSCF EFFECTIVE POTENTIAL MATRIX "G"
 !     --------------------------------------------------------------
       implicit none
 
       integer,intent(in) :: ref(4)               ! Reference state          
       integer,intent(in) :: nvdf
       integer,intent(in) :: ngaus
       real*8,intent(in)  :: Po(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q1(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q2(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Q3(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: Scho(ngaus,ngaus,nvdf)
       real*8,intent(in)  :: tiij(nvdf,nvdf)
       real*8,intent(in)  :: tjji(nvdf,nvdf)
       real*8,intent(in)  :: uiiij(nvdf,nvdf)
       real*8,intent(in)  :: ujjji(nvdf,nvdf)
       real*8,intent(in)  :: uiijj(nvdf,nvdf)
       real*8,intent(in)  :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)  :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)  :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)  :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(out) :: Gmtrx(ngaus,ngaus,nvdf)
       real*8,intent(out) :: GDmtrx(ngaus,ngaus,nvdf)
       real*8,intent(out) :: GTmtrx(ngaus,ngaus,nvdf)
 
       integer :: nm,nm1,nm2,nm3,n1,n2,n3
       integer :: nnm1,nnm2,qn1,qn2
       integer :: INFO
       real*8  :: Coef1
       real*8  :: Coef2
       real*8  :: Coef3
       real*8  :: tmp(ngaus)
       real*8  :: Q1avg(nvdf)
       real*8  :: Q2avg(nvdf)
       real*8  :: Q3avg(nvdf)
       real*8  :: tScho(ngaus,ngaus)
       real*8  :: Qtmp(ngaus,ngaus)
       real*8  :: Pt(ngaus)
       real*8,external :: ddot
 !     --------------------------------------------------------------
 
       nnm1=ref(1)
       qn1=ref(2)
       nnm2=ref(3)
       qn2=ref(4)
 
 !     ------------------------------------------------------------------
 !     COMPUTING AVERAGES OF QM OPERATORS OVER MODALS
 !     ------------------------------------------------------------------
 !     Computing average normal coordinate operators over ground state modals
 !                           <phi_a|(Q_a)**n|phi_a>
       do nm=1,nvdf
 
          Pt=0.0d0
          if (nm == nnm1) then
             Pt=Po(:,qn1+1,nm)
          else if (nm == nnm2) then
             Pt=Po(:,qn2+1,nm)
          else
             Pt=Po(:,1,nm)
          end if
          tScho=Scho(:,:,nm)
 
          tmp=0.0d0
          Qtmp=Q1(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Pt,1,0.0D0,tmp,1)
          Q1avg(nm) = ddot(ngaus,Pt,1,tmp,1)
 
          tmp=0.0d0
          Qtmp=Q2(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Pt,1,0.0D0,tmp,1)
          Q2avg(nm) = ddot(ngaus,Pt,1,tmp,1)
 
          tmp=0.0d0
          Qtmp=Q3(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Pt,1,0.0D0,tmp,1)
          Q3avg(nm) = ddot(ngaus,Pt,1,tmp,1)
       end do
 
 !     Computing vscf effective potential matrix G
 !     ------------------------------------------------------------------
 !     BEGINING COMPUTATION OF 2 MODE COUPLING TERMS
 !     ------------------------------------------------------------------
       GDmtrx=0.0d0
       GTmtrx=0.0d0
       do nm1=1,nvdf
          do nm2=1,nvdf
             if (nm1 == nm2) CYCLE
             if (nm1 > nm2) then
                n1=nm2
                n2=nm1
 
                Coef1=0.0d0
                Coef2=0.0d0
                Coef3=0.0d0
 
                Coef1=tiij(n1,n2)*Q2avg(n1)/2d0+uiiij(n1,n2)*Q3avg(n1)/6d0
                Coef2=tjji(n1,n2)*Q1avg(n1)/2d0+uiijj(n1,n2)*Q2avg(n1)/4d0
                Coef3=ujjji(n1,n2)*Q1avg(n1)/6d0
 
                GDmtrx(:,:,nm1)=GDmtrx(:,:,n2)+&
                             &Coef1*Q1(:,:,n2)+&
                             &Coef2*Q2(:,:,n2)+&
                             &Coef3*Q3(:,:,n2)
             else
                n1=nm1
                n2=nm2
 
                Coef1=0.0d0
                Coef2=0.0d0
                Coef3=0.0d0
 
                Coef1=tjji(n1,n2)*Q2avg(n2)/2d0+ujjji(n1,n2)*Q3avg(n2)/6d0
                Coef2=tiij(n1,n2)*Q1avg(n2)/2d0+uiijj(n1,n2)*Q2avg(n2)/4d0
                Coef3=uiiij(n1,n2)*Q1avg(n2)/6d0
 
                GDmtrx(:,:,nm1)=GDmtrx(:,:,n1)+&
                             &Coef1*Q1(:,:,n1)+&
                             &Coef2*Q2(:,:,n1)+&
                             &Coef3*Q3(:,:,n1)
             end if
          end do
       end do
 
 !     ------------------------------------------------------------------
 !     BEGINING COMPUTATION OF 3 MODE COUPLING TERMS
 !     ------------------------------------------------------------------
       do nm1=1,nvdf
          do nm2=1,nvdf-1
             do nm3=nm2+1,nvdf
                if (nm1 == nm2) CYCLE
                if (nm1 == nm3) CYCLE
                if (nm1<nm2) then
                   n1=nm1
                   n2=nm2
                   n3=nm3
 
                   coef1=0d0 
                   coef2=0d0
    
                   coef1=tijk(n1,n2,n3)*Q1avg(n2)*Q1avg(n3)+&
                        &uijjk(n1,n2,n3)*Q2avg(n2)*Q1avg(n3)/2d0+&
                        &uijkk(n1,n2,n3)*Q1avg(n2)*Q2avg(n3)/2d0
    
                   coef2=uiijk(n1,n2,n3)*Q1avg(n2)*Q1avg(n3)/2d0
    
                   GTmtrx(:,:,n1)=GTmtrx(:,:,n1)+&
                               &coef1*Q1(:,:,n1)+&
                               &coef2*Q2(:,:,n1)
 
                else if (nm1>nm2 .AND. nm1<nm3) then
                   n1=nm2
                   n2=nm1
                   n3=nm3
 
                   coef1=0d0 
                   coef2=0d0
    
                   coef1=tijk(n1,n2,n3)*Q1avg(n1)*Q1avg(n3)+&
                        &uiijk(n1,n2,n3)*Q2avg(n1)*Q1avg(n3)/2d0+&
                        &uijkk(n1,n2,n3)*Q1avg(n1)*Q2avg(n3)/2d0
    
                   coef2=uijjk(n1,n2,n3)*Q1avg(n1)*Q1avg(n3)/2d0
    
                   GTmtrx(:,:,n2)=GTmtrx(:,:,n2)+&
                               &coef1*Q1(:,:,n2)+&
                               &coef2*Q2(:,:,n2)
                else if (nm1>nm3) then
                   n1=nm2
                   n2=nm3
                   n3=nm1
 
                   coef1=0d0 
                   coef2=0d0
    
                   coef1=tijk(n1,n2,n3)*Q1avg(n1)*Q1avg(n2)+&
                        &uiijk(n1,n2,n3)*Q2avg(n1)*Q1avg(n2)/2d0+&
                        &uijjk(n1,n2,n3)*Q1avg(n1)*Q2avg(n2)/2d0
    
                   coef2=uijkk(n1,n2,n3)*Q1avg(n1)*Q1avg(n2)/2d0
    
                   GTmtrx(:,:,n3)=GTmtrx(:,:,n3)+&
                               &coef1*Q1(:,:,n3)+&
                               &coef2*Q2(:,:,n3)
 
                end if
                
             end do
          end do
       end do 
 
       Gmtrx=GDmtrx+GTmtrx
                   
       end subroutine
 
       subroutine ssvscfcorr(ref,Scho,Po,Q1,Q2,Q3,tiij,tjji,uiiij,ujjji,uiijj,&
                   & tijk,uiijk,uijjk,uijkk,nvdf,ngaus,corr,corrD,corrT)
 !     ------------------------------------------------------------------------ 
 !     COMPUTES VSCF ENERGY CORRECTION
 !     CORR = <PSI|VcD|PSI> + 2 <PSI|VcT|PSI>
 !     WHERE VcD AND VcT ARE THE DOUBLES AND TRIPLES COUPLING POTENTIAL,
 !     RESPECTIVELY.
 !     ------------------------------------------------------------------------ 
       implicit none
 
       integer,intent(in)     :: ref(4)
       integer,intent(in)     :: nvdf
       integer,intent(in)     :: ngaus
       real*8,intent(in)      :: Scho(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Po(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Q1(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Q2(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: Q3(ngaus,ngaus,nvdf)! Effective potential matrix.
       real*8,intent(in)      :: tiij(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: tjji(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uiiij(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: ujjji(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uiijj(nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: tijk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uiijk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uijjk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(in)      :: uijkk(nvdf,nvdf,nvdf)     ! Effective potential matrix.
       real*8,intent(out)     :: corr
       real*8,intent(out)     :: corrD
       real*8,intent(out)     :: corrT
       
       integer     :: nm
       integer     :: nm1,nm2,qn1,qn2
       integer     :: a,b,c
       integer     :: INFO                   ! Error code for dsygst subroutine.
       real*8      :: Qtmp(ngaus,ngaus)
       real*8      :: tScho(ngaus,ngaus)
       real*8      :: Ptmp(ngaus)
       real*8      :: tmp(ngaus)
       real*8      :: Qm1(nvdf)
       real*8      :: Qm2(nvdf)
       real*8      :: Qm3(nvdf)
 
       real*8,external :: ddot
 !     DEBUG
       integer     :: i,j
 
       nm1=ref(1)
       qn1=ref(2)
       nm2=ref(3)
       qn2=ref(4)
 
       corr=0.0D0
       do nm=1,nvdf
          tScho(:,:) = Scho(:,:,nm)
          if (nm == nm1) then
             Ptmp=Po(:,qn1+1,nm)
          else if (nm == nm2) then
             Ptmp=Po(:,qn2+1,nm)
          else
             Ptmp=Po(:,1,nm)
          end if
 
          tmp=0.0d0
          Qtmp(:,:) = 0d0
          Qtmp(:,:) = Q1(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
          Qm1(nm) = ddot(ngaus,Ptmp,1,tmp,1)
 
          tmp=0.0d0
          Qtmp(:,:) = 0d0
          Qtmp(:,:) = Q2(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
          Qm2(nm) = ddot(ngaus,Ptmp,1,tmp,1)
 
          tmp=0.0d0
          Qtmp(:,:) = 0d0
          Qtmp(:,:) = Q3(:,:,nm)
          call dsygst(1,'U',ngaus,Qtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymv('U',ngaus,1.0D0,Qtmp,ngaus,Ptmp,1,0.0D0,tmp,1)
          Qm3(nm) = ddot(ngaus,Ptmp,1,tmp,1)
       end do
 
 !     2-MODE COUPLING POTENTIAL
       corrD=0.0d0
       do a=1,nvdf-1
          do b=a+1,nvdf
             corrD=corrD+&
             & tiij(a,b)*Qm2(a)*Qm1(b)/2.0d0+&
             & tjji(a,b)*Qm1(a)*Qm2(b)/2.0d0+&
             & uiiij(a,b)*Qm3(a)*Qm1(b)/6.0d0+&
             & ujjji(a,b)*Qm1(a)*Qm3(b)/6.0d0+&
             & uiijj(a,b)*Qm2(a)*Qm2(b)/4.0d0  
          end do
       end do
 
 !     3-MODE COUPLING POTENTIAL
       corrT=0d0
       do a=1,nvdf-2
          do b=a+1,nvdf-1
             do c=b+1,nvdf
                corrT=corrT+&
                 & tijk(a,b,c)*Qm1(a)*Qm1(b)*Qm1(c)+&
                 & uiijk(a,b,c)*Qm2(a)*Qm1(b)*Qm1(c)/2d0+&
                 & uijjk(a,b,c)*Qm1(a)*Qm2(b)*Qm1(c)/2d0+&
                 & uijkk(a,b,c)*Qm1(a)*Qm1(b)*Qm2(c)/2d0
             end do
          end do
       end do
 
       corr = corrD + 2d0*corrT
       end subroutine
 
 
 
!       subroutine ssvscfci(qumvia_nmc,qumvia_qff,nrst,ndf,nvdf,ngaus,nmcoup,nqmatoms,nclatoms,&
!                          &qmaxx1,qmaxx2,nconf,at_numbers,dy,gwidth,eig,&
!                          &nmodes,qmcoords,clcoords)
! !     -----------------------------------------------------------------
! !     THIS SUBROUTINE PERFORMS VIBRATIONAL SELF-CONSISTENT FIELD 
! !     FOLLOWED BY CONFIGURATION INTERACTION CALCULATION IN A 
! !     STATE-SPECIFIC FASHION.
! !     A separate VSCF and VCI computation is performed for each state
! !     of interest, by default the ground state and the first excited 
! !     states of each normal mode, using that state as the reference.
! !     The reference state is the vscf configuration over which the
! !     effective potential is computed.
! !     -----------------------------------------------------------------
!       implicit none
! 
!       integer,intent(in)  :: qumvia_nmc           ! # of coupled mode in Hci
!       integer,intent(in)  :: qumvia_qff           ! # keyword for qff read/calc
!       integer,intent(in)  :: nrst
!       integer,intent(in)  :: ndf                  ! Number of classical atoms
!       integer,intent(in)  :: nvdf                 ! Number of classical atoms
!       integer,intent(in)  :: ngaus                ! Number of classical atoms
!       integer,intent(in)  :: nmcoup               ! Number of normal modes to couple in QFF.
!       integer,intent(in)  :: nqmatoms             ! Number of QM atoms
!       integer,intent(in)  :: nclatoms             ! Number of classical atoms
!       integer,intent(in)  :: qmaxx1               ! Max excitation for singles.
!       integer,intent(in)  :: qmaxx2               ! Max excitation for doules.
!       integer,intent(in)  :: nconf                ! Dimension of CI basis set.
!       integer,intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
!       real*8,intent(in)   :: dy                   ! Step size factor for num derivatives.
!       real*8,intent(in)   :: gwidth               ! Width factor for gaussian primitives.
!       real*8,intent(in)   :: eig(ndf)             ! Eigenvalues of the hessian matrix.
!       real*8,intent(in)   :: nmodes(ndf,ndf)      ! Eigenvectors of the hessian matrix.
!       real*8,intent(in)   :: qmcoords(3,nqmatoms) ! QM atom coordinates
!       real*8,intent(in)   :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
! 
!       integer            :: ref(4)
!       integer            :: cnf
!       integer            :: bdim
!       integer            :: refstat(4,nrst)
!       integer            :: vscfstat
!       real*8             :: P(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
!       real*8             :: Po(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
!       real*8             :: Scho(ngaus,ngaus,nvdf)  ! Cholesky factor.
!       real*8             :: Evscf                ! VSCF energy.
!       real*8             :: Q1(ngaus,ngaus,nvdf) !  <gi|Qa^1|gj> matrix
!       real*8             :: Q2(ngaus,ngaus,nvdf) !  <gi|Qa^2|gj> matrix
!       real*8             :: Q3(ngaus,ngaus,nvdf) !  <gi|Qa^3|gj> matrix
!       real*8             :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian (non-orthogonal basis)
!       real*8             :: Essvhf(nrst)
!       real*8             :: Pss(ngaus,nrst)
!       real*8             :: Poss(ngaus,nrst)
!       real*8             :: ECI(nconf)         ! CI energy.
!       real*8             :: Cci(nconf,nconf)   ! CI eigenvectors.
!       real*8             :: GDmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
!       real*8             :: GTmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
!       real*8             :: Emod(ngaus,nvdf)       ! Modal energies.
!       real*8             :: hii(ndf)             ! Diagonal Hessian eigenvalues.
!       real*8             :: tiii(nvdf)           ! Diagonal cubic coupling terms
!       real*8             :: uiiii(nvdf)          ! Diagonal quartic coupling terms
!       real*8             :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
!       real*8             :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
!       real*8             :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling terms
!       real*8             :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling terms
!       real*8             :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling terms
!       real*8             :: tijk(nvdf,nvdf,nvdf)      ! 3-mode cubic coupling terms
!       real*8             :: uiijk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
!       real*8             :: uijjk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
!       real*8             :: uijkk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
! 
! !     -----------------------------------------------------------------
! !     GENERATING REFERENCE FONFIGURATIONS
!       refstat=0
!       call genConf2b(0,0,nvdf,nvdf+1,refstat)
! 
! !     QUARTIC FORCE FIELD
! !     ---------------------------------------------------------------
!       call selectqff(qumvia_qff,nmodes,eig,dy,nqmatoms,nclatoms,&
!                      &qmcoords,clcoords,at_numbers,ndf,nvdf,&
!                      &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj,&
!                      &tijk,uiijk,uijjk,uijkk)
! 
! !     COMPUTING VSCF/VCI OVER THESE CONFIGURATIONS      
!       do cnf=1,nrst
! !        Computing VSCF
!          ref(:) = refstat(:,cnf)
!          call ssvscf2(ref,nmodes,eig,dy,nmcoup,gwidth,nqmatoms,nclatoms,&
!                    &ngaus,qmcoords,clcoords,at_numbers,ndf,nvdf,&
!                    &P,Po,Scho,Evscf,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Emod,&
!                    &hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
!                    &tijk,uiijk,uijjk,uijkk,vscfstat)
!          if (vscfstat /=0) CYCLE
!          Essvhf(cnf)=Evscf
!          Pss(:,cnf)=P(:,cnf,cnf)
!          Poss(:,cnf)=Po(:,cnf,cnf)
! 
! !        Computing VCI
!          write(77,'(A,4I3)') 'COMPUTING VCI FOR REFERENCE STATE ',ref
!          bdim=qmaxx1+1
!          call virtCI2(qumvia_nmc,Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,Emod,&
!                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
!                    &nconf,ngaus,nvdf,qmaxx1,qmaxx2,bdim,ECI,Cci)
!       end do
! 
!       write(77,'(A)') 'ssVSCF ENERGIES AND TRANSITIONS (CM-1)'
!       do cnf=1,nrst
!          write(77,'(I3,2F16.2)') cnf-1, Essvhf(cnf), Essvhf(cnf)-Essvhf(1)
!       end do
! 
! 
!       end subroutine
 
 
 
 
 
 
 !#######################################################################
 !     VIRTUAL CONFIGURATION INTERACTION SECTION
 !#######################################################################
 
 
 
 
 
       subroutine build_CIop2(Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,ngaus,nvdf,nmods,&
                     &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 !     -----------------------------------------------------------------
 !     CHANGE OF BASIS FOR Q AND CORE HAMILTONIAN OPERATORS TO VSCF 
 !     MODALS BASIS SET.
 !     THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT P^T*Qa^n*P
 !     THE RESULTING MATRIX ELEMENTS ARE <modal_i|Qa^n|modal_j>
 !     WHERE I AND J INDICATE DIFFERENT (VIRTUAL) MODALS OF THE SAME
 !     NORMAL MODE "a". THE SAME FOR THE CORE HAMILTONIAN MATRIX.
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: ngaus
       integer,intent(in)  :: nvdf
       integer,intent(in)  :: nmods ! Dimension of CI operator matrices.
       real*8,intent(in)   :: Po(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Hcore(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Q1(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Q2(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Q3(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: GDmtrx(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: GTmtrx(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Scho(ngaus,ngaus,nvdf)
 
       real*8,intent(out)  :: Qm1(nmods,nmods,nvdf)
       real*8,intent(out)  :: Qm2(nmods,nmods,nvdf)
       real*8,intent(out)  :: Qm3(nmods,nmods,nvdf)
       real*8,intent(out)  :: Hmc(nmods,nmods,nvdf)
       real*8,intent(out)  :: GDm(nmods,nmods,nvdf)
       real*8,intent(out)  :: GTm(nmods,nmods,nvdf)
  
       integer      :: nm
       integer      :: INFO
       real*8       :: tScho(ngaus,ngaus)
       real*8       :: Itmp(ngaus,ngaus)
       real*8       :: Ptmp(ngaus,nmods)
       real*8       :: Mtmp(ngaus,nmods)
       real*8       :: Rtmp(nmods,nmods)
 !     -----------------------------------------------------------------
       do nm=1,nvdf
          Ptmp=Po(:,1:nmods,nm)
          tScho=Scho(:,:,nm)
 
          Mtmp=0.0d0
          Rtmp=0.0d0
          Itmp=GDmtrx(:,:,nm)
          call dsygst(1,'U',ngaus,Itmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymm('L','U',ngaus,nmods,1.0D0,Itmp,ngaus,&
               &Ptmp,ngaus,0.0D0,Mtmp,ngaus)
          call dgemm('T','N',nmods,nmods,ngaus,1.0D0,Ptmp,&
               &ngaus,Mtmp,ngaus,0.0D0,Rtmp,nmods)
          GDm(:,:,nm)=Rtmp
 
          Mtmp=0.0d0
          Rtmp=0.0d0
          Itmp=GTmtrx(:,:,nm)
          call dsygst(1,'U',ngaus,Itmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymm('L','U',ngaus,nmods,1.0D0,Itmp,ngaus,&
               &Ptmp,ngaus,0.0D0,Mtmp,ngaus)
          call dgemm('T','N',nmods,nmods,ngaus,1.0D0,Ptmp,&
               &ngaus,Mtmp,ngaus,0.0D0,Rtmp,nmods)
          GTm(:,:,nm)=Rtmp
 
          Mtmp=0.0d0
          Rtmp=0.0d0
          Itmp=Hcore(:,:,nm)
          call dsygst(1,'U',ngaus,Itmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymm('L','U',ngaus,nmods,1.0D0,Itmp,ngaus,&
               &Ptmp,ngaus,0.0D0,Mtmp,ngaus)
          call dgemm('T','N',nmods,nmods,ngaus,1.0D0,Ptmp,&
               &ngaus,Mtmp,ngaus,0.0D0,Rtmp,nmods)
          Hmc(:,:,nm)=Rtmp
 
          Mtmp=0.0d0
          Rtmp=0.0d0
          Itmp=Q1(:,:,nm)
          call dsygst(1,'U',ngaus,Itmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymm('L','U',ngaus,nmods,1.0D0,Itmp,ngaus,&
               &Ptmp,ngaus,0.0D0,Mtmp,ngaus)
          call dgemm('T','N',nmods,nmods,ngaus,1.0D0,Ptmp,&
               &ngaus,Mtmp,ngaus,0.0D0,Rtmp,nmods)
          Qm1(:,:,nm)=Rtmp
 
          Mtmp=0.0d0
          Rtmp=0.0d0
          Itmp=Q2(:,:,nm)
          call dsygst(1,'U',ngaus,Itmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymm('L','U',ngaus,nmods,1.0D0,Itmp,ngaus,&
               &Ptmp,ngaus,0.0D0,Mtmp,ngaus)
          call dgemm('T','N',nmods,nmods,ngaus,1.0D0,Ptmp,&
               &ngaus,Mtmp,ngaus,0.0D0,Rtmp,nmods)
          Qm2(:,:,nm)=Rtmp
 
          Mtmp=0.0d0
          Rtmp=0.0d0
          Itmp=Q3(:,:,nm)
          call dsygst(1,'U',ngaus,Itmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix trasnformation')
          call dsymm('L','U',ngaus,nmods,1.0D0,Itmp,ngaus,&
               &Ptmp,ngaus,0.0D0,Mtmp,ngaus)
          call dgemm('T','N',nmods,nmods,ngaus,1.0D0,Ptmp,&
               &ngaus,Mtmp,ngaus,0.0D0,Rtmp,nmods)
          Qm3(:,:,nm)=Rtmp
 
       end do
       end subroutine
 
       subroutine build_CIop3(Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,ngaus,nvdf,nmods,&
                     &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 !     -----------------------------------------------------------------
 !     CHANGE OF BASIS FOR Q AND CORE HAMILTONIAN OPERATORS TO VSCF 
 !     MODALS BASIS SET.
 !     THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT P^T*Qa^n*P
 !     THE RESULTING MATRIX ELEMENTS ARE <modal_i|Qa^n|modal_j>
 !     WHERE I AND J INDICATE DIFFERENT (VIRTUAL) MODALS OF THE SAME
 !     NORMAL MODE "a". THE SAME FOR THE CORE HAMILTONIAN MATRIX.
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: ngaus
       integer,intent(in)  :: nvdf
       integer,intent(in)  :: nmods ! Dimension of CI operator matrices.
       real*8,intent(in)   :: Po(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Hcore(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Q1(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Q2(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Q3(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: GDmtrx(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: GTmtrx(ngaus,ngaus,nvdf)
       real*8,intent(in)   :: Scho(ngaus,ngaus,nvdf)
 
       real*8,intent(out)  :: Qm1(nmods,nmods,nvdf)
       real*8,intent(out)  :: Qm2(nmods,nmods,nvdf)
       real*8,intent(out)  :: Qm3(nmods,nmods,nvdf)
       real*8,intent(out)  :: Hmc(nmods,nmods,nvdf)
       real*8,intent(out)  :: GDm(nmods,nmods,nvdf)
       real*8,intent(out)  :: GTm(nmods,nmods,nvdf)
  
       integer      :: nm
       integer      :: n1,n2
       integer      :: INFO
       real*8       :: tmp(ngaus)
       real*8       :: Ptmp1(ngaus)
       real*8       :: Ptmp2(ngaus)
       real*8       :: tScho(ngaus,ngaus)
       real*8       :: Q1tmp(ngaus,ngaus)
       real*8       :: Q2tmp(ngaus,ngaus)
       real*8       :: Q3tmp(ngaus,ngaus)
       real*8       :: GDtmp(ngaus,ngaus)
       real*8       :: GTtmp(ngaus,ngaus)
 
       real*8,external :: ddot
 !     DEBUG
       integer      :: i,j
 !     -----------------------------------------------------------------
       Qm1=0d0
       Qm2=0d0
       Qm3=0d0
       GDm=0d0
       GTm=0d0
       Hmc=0d0
       do nm=1,nvdf
 
          Q1tmp(:,:)=0d0
          Q2tmp(:,:)=0d0
          Q3tmp(:,:)=0d0
          GDtmp(:,:)=0d0
          GTtmp(:,:)=0d0
          tScho(:,:)=0d0
 
          tScho=Scho(:,:,nm)
          Q1tmp(:,:)=Q1(:,:,nm)
          Q2tmp(:,:)=Q2(:,:,nm)
          Q3tmp(:,:)=Q3(:,:,nm)
          GDtmp(:,:)=GDmtrx(:,:,nm)
          GTtmp(:,:)=GTmtrx(:,:,nm)
 
          call dsygst(1,'U',ngaus,Q1tmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsygst(1,'U',ngaus,Q2tmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsygst(1,'U',ngaus,Q3tmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsygst(1,'U',ngaus,GDtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
          call dsygst(1,'U',ngaus,GTtmp,ngaus,tScho,ngaus,INFO)
          if (INFO /= 0) STOP('Error during Matrix transformation')
 
          do n1=1,nmods
             do n2=n1,nmods
                Ptmp1=Po(:,n1,nm)
                Ptmp2=Po(:,n2,nm)
       
                tmp=0.0d0
                call dsymv('U',ngaus,1.0D0,Q1tmp,ngaus,Ptmp1,1,0.0D0,tmp,1)
                Qm1(n1,n2,nm) = ddot(ngaus,Ptmp2,1,tmp,1)
                if (n1 /= n2) Qm1(n2,n1,nm) = Qm1(n1,n2,nm)
       
                tmp=0.0d0
                call dsymv('U',ngaus,1.0D0,Q2tmp,ngaus,Ptmp1,1,0.0D0,tmp,1)
                Qm2(n1,n2,nm) = ddot(ngaus,Ptmp2,1,tmp,1)
                if (n1 /= n2) Qm2(n2,n1,nm) = Qm2(n1,n2,nm)
       
                tmp=0.0d0
                call dsymv('U',ngaus,1.0D0,Q3tmp,ngaus,Ptmp1,1,0.0D0,tmp,1)
                Qm3(n1,n2,nm) = ddot(ngaus,Ptmp2,1,tmp,1)
                if (n1 /= n2) Qm3(n2,n1,nm) = Qm3(n1,n2,nm)
       
                tmp=0.0d0
                call dsymv('U',ngaus,1.0D0,GDtmp,ngaus,Ptmp1,1,0.0D0,tmp,1)
                GDm(n1,n2,nm) = ddot(ngaus,Ptmp2,1,tmp,1)
                if (n1 /= n2) GDm(n2,n1,nm) = GDm(n1,n2,nm)
       
                tmp=0.0d0
                call dsymv('U',ngaus,1.0D0,GTtmp,ngaus,Ptmp1,1,0.0D0,tmp,1)
                GTm(n1,n2,nm) = ddot(ngaus,Ptmp2,1,tmp,1)
                if (n1 /= n2) GTm(n2,n1,nm) = GTm(n1,n2,nm)
             end do
          end do
       end do
       end subroutine
 
       subroutine genConf2b(qmaxx1,qmaxx2,nvdf,nconf,configs)
 !     -----------------------------------------------------------------
 !     THIS SUBROUTINE GENERATES ALL CI SINGLE AND DOUBLE CONFIGURATIONS
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: qmaxx1               ! Max excitation singles
       integer,intent(in)  :: qmaxx2               ! Max excitation doubles
       integer,intent(in)  :: nvdf                ! No of vibrational dof.
       integer,intent(in)  :: nconf               ! CI basis set size.
       integer,intent(out) :: configs(4,nconf)    ! Configurations.
       integer    :: nc
       integer    :: nm
       integer    :: qn
       integer    :: nm1
       integer    :: nm2
       integer    :: qn1
       integer    :: qn2
 !     -----------------------------------------------------------------
 !     Total number of configurations
 !     nconf = 1 + nvdf*qmaxx1 + nvdf*(nvdf-1)*qmaxx2*qmaxx2/2
       nc=1
 !     GROUND STATE
 !     -----------------------------------------------------------------
       configs(:,nc)=0
       nc=nc+1
 
 !     FUNDAMENTALS
 !     -----------------------------------------------------------------
       do nm=1,nvdf
          configs(:,nc) = (/ nm, 1, 0, 0 /)
          nc=nc+1
       end do
 
       if (qmaxx1 == 0) RETURN
 !     OVERTONES
 !     -----------------------------------------------------------------
       do nm=1,nvdf
          do qn=2,qmaxx1
             configs(:,nc) = (/ nm, qn, 0, 0 /)
             nc=nc+1
          end do
       end do
 
       if (qmaxx2 == 0) RETURN
 
 !     DOUBLES
 !     -----------------------------------------------------------------
       do nm1=1,nvdf-1
          do nm2=nm1+1,nvdf
             do qn1=1,qmaxx2
                do qn2=1,qmaxx2
                   configs(:,nc)=(/ nm1, qn1, nm2, qn2 /)
                   nc=nc+1
                end do
             end do
          end do
       end do
 
       end subroutine
 
       subroutine Hci_diag2(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,tiij,tjji,uiiij,ujjji,uiijj,&
                      &tijk,uiijk,uijjk,uijkk,nvdf,ngaus,nmods,nconf,configs,Hci)
 !     -----------------------------------------------------------------
 !     This subroutine computes the hamiltonian DIAGONAL matrix elements
 !     over VSCF states Psi_a. 
 !                      Haa = < Psi_a | H | Psi_a >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc ! # of coupled modes.
       integer,intent(in)   :: nvdf   ! # of vibrational deg of freedom.
       integer,intent(in)   :: nmods  ! # of selected VSCF virtual states.
       integer,intent(in)   :: nconf  ! Total number of configurations
       integer,intent(in)   :: ngaus  ! Dimension of gaussian basis set
       integer,intent(in)   :: configs(4,nconf)  ! List of configurations.
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf) !Q operator in modals basis
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf) !Q^2 operator in modals basis
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf) !Q^3 operator in modals basis
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8,intent(in)    :: tiij(nvdf,nvdf)   ! Coupling potential parameters
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(inout) :: Hci(nconf+nconf*(nconf-1)/2) ! VCI hamiltonian matrix in packed storage format.
 
       integer :: i,a,b,c
       integer :: nm1,nm2
       integer :: qn1,qn2
       integer :: cnf
       integer :: n1,n2,n3 ! Various indices.
       integer :: psi(nvdf) ! VSCF wavefunction. Array of vib. quantum numbers.
       real*8 :: Vc2       ! Auxiliary variable for Hci calculation.
       real*8 :: Vc3       ! Auxiliary variable for Hci calculation.
       real*8 :: ehf       ! Auxiliary variable for Hci calculation.
       real*8 :: veff      ! Auxiliary variable for Hci calculation.
       real*8 :: veff2     ! Auxiliary variable for Hci calculation.
       real*8 :: veff3     ! Auxiliary variable for Hci calculation.
       real*8 :: Ehf0      ! DEBUG
       real*8 :: corr      ! DEBUG
 
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
 !     -----------------------------------------------------------------
       
       do i=1,nconf
          cnf=i*(i+1)/2         ! Label for packed matrix diagonal terms.
          nm1=configs(1,i)
          qn1=configs(2,i)
          nm2=configs(3,i)
          qn2=configs(4,i)
 
          psi=0
          if (nm1 > 0) psi(nm1)=qn1
          if (nm2 > 0) psi(nm2)=qn2
 
 !        VSCF energy.
          ehf=0d0
          do a=1,nvdf
             ehf = ehf + Emod(psi(a)+1,a)
          end do
 
 !        Average on effective potential.
          veff2=0d0
          do a=1,nvdf
             veff2 = veff2 + GDm(psi(a)+1,psi(a)+1,a)
          end do
          veff3=0d0
          do a=1,nvdf
             veff3 = veff3 + GTm(psi(a)+1,psi(a)+1,a)
          end do
          veff = veff2 + veff3
 
 !        2-MODE COUPLING POTENTIAL 
          Vc2=0.0d0
          do a=1,nvdf-1
             do b=a+1,nvdf
                n1=psi(a)+1
                n2=psi(b)+1
                Vc2=Vc2+&
                & tiij(a,b)*Qm2(n1,n1,a)*Qm1(n2,n2,b)/2.0d0+&
                & tjji(a,b)*Qm1(n1,n1,a)*Qm2(n2,n2,b)/2.0d0+&
                & uiiij(a,b)*Qm3(n1,n1,a)*Qm1(n2,n2,b)/6.0d0+&
                & ujjji(a,b)*Qm1(n1,n1,a)*Qm3(n2,n2,b)/6.0d0+&
                & uiijj(a,b)*Qm2(n1,n1,a)*Qm2(n2,n2,b)/4.0d0  
             end do
          end do
 
 !        3-MODE COUPLING POTENTIAL
          Vc3=0d0
          IF (qumvia_nmc == 3) THEN
          do a=1,nvdf-2
             do b=a+1,nvdf-1
                do c=b+1,nvdf
                   n1=psi(a)+1
                   n2=psi(b)+1
                   n3=psi(c)+1
 
                   Vc3=Vc3+&
                    & tijk(a,b,c)*Qm1(n1,n1,a)*Qm1(n2,n2,b)*Qm1(n3,n3,c)+&
                    & uiijk(a,b,c)*Qm2(n1,n1,a)*Qm1(n2,n2,b)*Qm1(n3,n3,c)/2d0+&
                    & uijjk(a,b,c)*Qm1(n1,n1,a)*Qm2(n2,n2,b)*Qm1(n3,n3,c)/2d0+&
                    & uijkk(a,b,c)*Qm1(n1,n1,a)*Qm1(n2,n2,b)*Qm2(n3,n3,c)/2d0
 
                end do
             end do
          end do
          END IF
 
 !        COMPUTING FINAL MATRIX ELEMENT
 !         Hci(cnf)=ehf-veff+Vc2+Vc3
 !        GAMESS FORM OF Hci DIAGONAL ELEMENTS
          Hci(cnf)=ehf-Vc2-2.0d0*Vc3
 !        ---------------------------------------------------------------
 !        DEBUG DANGER
 !        ---------------------------------------------------------------
         if (cnf < 30) then
          if (cnf==1) Ehf0=(ehf-Vc2-2d0*Vc3)*h2cm
          if (cnf==1) corr=Vc2+2d0*Vc3
          if (cnf ==1) write(77,'(A)') 'CNF        PSI       Evscf          <PSI|Vc|PSI>   <PSI|vef|PSI>   Evscf-corr' 
          write(77,'(I6,A,3I2,4F16.9,2F13.2)') cnf,'  ',psi,ehf,Vc2,Vc3,veff2,(ehf-corr)*h2cm,(ehf-corr)*h2cm-Ehf0
 !         if (cnf ==1) write(77,'(A)') '<veffD> CI'
 !         if (cnf==1) write(77,'(D12.3)') (GDm(1,1,a),a=1,nvdf)
         end if
 !        ---------------------------------------------------------------
          
       end do
       
 
       end subroutine
 
       subroutine calcHterm1b(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,tiij,tjji,uiiij,ujjji,uiijj,&
                    &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,Hci,cnf,psi1,psi2,dif)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY ONE MODAL.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf,cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(inout) :: Hci(nconf+nconf*(nconf-1)/2)
 
       integer :: i,j,k
       integer :: a,b
       integer :: n1,n2
       integer :: ni1,ni2
       integer :: nj1,nj2
       integer :: nk1,nk2
       real*8 :: Vc2
       real*8 :: Vc3
       real*8 :: veff      ! Auxiliary variable for Hci calculation.
 !     -----------------------------------------------------------------
 
 
 !     Computing VSCF effective potential.
       veff=0d0
       n1=psi1(dif(1))+1
       n2=psi2(dif(1))+1
       veff=GDm(n1,n2,dif(1))+GTm(n1,n2,dif(1))
 
 !     Computing Coupling Potential
       Vc2 = 0.0d0
       do a=1,nvdf
          if (a == dif(1)) CYCLE
 
          if (a<dif(1)) then
             i=a
             j=dif(1)
             ni1=psi1(a)+1
             ni2=ni1
             nj1=psi1(dif(1))+1
             nj2=psi2(dif(1))+1
          else
             i=dif(1)
             j=a
             ni1=psi1(dif(1))+1
             ni2=psi2(dif(1))+1
             nj1=psi1(a)+1
             nj2=nj1
          end if 
              
          Vc2 = Vc2 + &
               & tiij(i,j)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)/2.0d0 + &
               & tjji(i,j)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)/2.0d0 + &
               & uiiij(i,j)*Qm3(ni1,ni2,i)*Qm1(nj1,nj2,j)/6.0d0+ &
               & ujjji(i,j)*Qm1(ni1,ni2,i)*Qm3(nj1,nj2,j)/6.0d0+ &
               & uiijj(i,j)*Qm2(ni1,ni2,i)*Qm2(nj1,nj2,j)/4.0d0
       end do
 
       Vc3 = 0d0
       IF (qumvia_nmc == 3) THEN
       do a=1,nvdf-1
          do b=a+1,nvdf
 
             if (a == dif(1)) CYCLE
             if (b == dif(1)) CYCLE
 
             if (dif(1) < a) then
                i=dif(1)
                j=a
                k=b
 
                ni1=psi1(dif(1))+1
                ni2=psi2(dif(1))+1
                nj1=psi1(a)+1
                nj2=psi2(a)+1
                nk1=psi1(b)+1
                nk2=psi2(b)+1
             else if (a<dif(1) .AND. dif(1)<b) then
                i=a
                j=dif(1)
                k=b
 
                ni1=psi1(a)+1
                ni2=psi2(a)+1
                nj1=psi1(dif(1))+1
                nj2=psi2(dif(1))+1
                nk1=psi1(b)+1
                nk2=psi2(b)+1
             else if (b<dif(1)) then
                i=a
                j=b
                k=dif(1)
 
                ni1=psi1(a)+1
                ni2=psi2(a)+1
                nj1=psi1(b)+1
                nj2=psi2(b)+1
                nk1=psi1(dif(1))+1
                nk2=psi2(dif(1))+1
             end if
 
             Vc3=Vc3+&
              & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)+ &
              & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0+ &
              & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0+ &
              & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
          end do
       end do 
       END IF
 
 
 !     COMPUTING HAMILTONIAN MATRIX ELEMENT
 !     DANGER DANGER!   
       Hci(cnf)=Vc2 + Vc3 - veff
 !      Hci(cnf)=Hmc(n1,n2,dif(1)) + Vc2 + Vc3 - veff
 
       end subroutine
 
       subroutine calcHterm2b(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY TWO MODALS.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf
       integer,intent(in)   :: nmods
       integer,intent(in)   :: nconf
       integer,intent(in)   :: cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(inout) :: Hci(nconf+nconf*(nconf-1)/2)
 
       integer :: i,j,k
       integer :: a,b,c
       integer :: ni1,ni2,nj1,nj2,nk1,nk2
       real*8  :: Vc2,Vc3
 !     -----------------------------------------------------------------
 
 
 
       if (dif(1) < dif(2)) then
          i=dif(1)
          j=dif(2)
       else
          i=dif(2)
          j=dif(1)
       end if 
 
       ni1=psi1(i)+1
       ni2=psi2(i)+1
       nj1=psi1(j)+1
       nj2=psi2(j)+1
           
       Vc2 = 0d0
       Vc2 =  &
          & tiij(i,j)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)/2.0d0+&
          & tjji(i,j)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)/2.0d0+&
          & uiiij(i,j)*Qm3(ni1,ni2,i)*Qm1(nj1,nj2,j)/6.0d0+&
          & ujjji(i,j)*Qm1(ni1,ni2,i)*Qm3(nj1,nj2,j)/6.0d0+&
          & uiijj(i,j)*Qm2(ni1,ni2,i)*Qm2(nj1,nj2,j)/4.0d0
 
       Vc3 = 0d0
       IF (qumvia_nmc == 3) THEN
       do a=1,nvdf
 
          if (a == dif(1)) CYCLE
          if (a == dif(2)) CYCLE
 
          if (a < dif(1) .AND. dif(1) < dif(2)) then
             i=a
             j=dif(1)
             k=dif(2)
          else if (dif(1) < a .AND. a < dif(2)) then
             i=dif(1)
             j=a
             k=dif(2)
          else if ( dif(2) < a ) then
             i=dif(1)
             j=dif(2)
             k=a
          end if
 
          ni1=psi1(i)+1
          ni2=psi2(i)+1
          nj1=psi1(j)+1
          nj2=psi2(j)+1
          nk1=psi1(k)+1
          nk2=psi2(k)+1
 
          Vc3=Vc3+ &
           & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k) + &
           & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
       end do
       END IF 
    
 !     VCI MATRIX ELEMENT
       Hci(cnf) = Vc2 + Vc3
 
       end subroutine
 
       subroutine calcHterm3(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY TWO MODALS.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf
       integer,intent(in)   :: nmods
       integer,intent(in)   :: nconf
       integer,intent(in)   :: cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(inout) :: Hci(nconf+nconf*(nconf-1)/2)
 
       integer :: i,j,k
       integer :: ni1,ni2,nj1,nj2,nk1,nk2
 !     -----------------------------------------------------------------
 
 
       IF (qumvia_nmc == 3) THEN
 !        COMPUTING 3-MODE COUPLING TERMS
          if (dif(1) > dif(2)) STOP ('dif(1) > dif(2) in calcHterm3')
          if (dif(2) > dif(3)) STOP ('dif(2) > dif(3) in calcHterm3')
          if (dif(1) > dif(3)) STOP ('dif(1) > dif(3) in calcHterm3')
    
          i=dif(1)
          j=dif(2)
          k=dif(3)
    
          ni1=psi1(i)+1
          ni2=psi2(i)+1
          nj1=psi1(j)+1
          nj2=psi2(j)+1
          nk1=psi1(k)+1
          nk2=psi2(k)+1
    
          Hci(cnf)= &
           & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k) + &
           & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
       ELSE
          Hci(cnf) = 0d0
       END IF
 
       end subroutine
 
       subroutine Hci_Offdiag2(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,tiij,tjji,uiiij,ujjji,&
                    &uiijj,tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,configs,Hci)
 !     -----------------------------------------------------------------
 !     Computes Off-diagonal Hamiltonian matrix elements 
 !                      < K | H | L >
 !     where | K > and | L > may differ in 1 (case 1), 2 (case 2) and 
 !     more than 2 (case 3) modals. This subroutine runs over all off-
 !     diagonal matrix elements, determines to which case each belong
 !     and then computes the integral using subroutines calcHterm1 
 !     for case 1 matrix elements, calcHterm2 for case 2 ME or
 !     sets Hci(cnf)=0.0 for case 3.
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf
       integer,intent(in)   :: configs(4,nconf)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(inout) :: Hci(nconf+nconf*(nconf-1)/2)
 
       integer    :: a
       integer    :: b
       integer    :: i
       integer    :: j
       integer    :: cnf
       integer    :: nm1
       integer    :: nm2
       integer    :: nm3
       integer    :: nm4
       integer    :: qn1
       integer    :: qn2
       integer    :: qn3
       integer    :: qn4
       integer    :: ndif
       integer    :: dif(3)
       integer    :: psi1(nvdf)
       integer    :: psi2(nvdf)
 !     -----------------------------------------------------------------
 
       do i=1,nconf-1
          do j=i+1,nconf
             cnf=i+j*(j-1)/2
             nm1=configs(1,i)
             qn1=configs(2,i)
             nm2=configs(3,i)
             qn2=configs(4,i)
 
             nm3=configs(1,j)
             qn3=configs(2,j)
             nm4=configs(3,j)
             qn4=configs(4,j)
 
             psi1=0
             psi2=0
 
             if (nm1 > 0) psi1(nm1)=qn1
             if (nm2 > 0) psi1(nm2)=qn2
             if (nm3 > 0) psi2(nm3)=qn3
             if (nm4 > 0) psi2(nm4)=qn4
 
             dif=0
             ndif=0
             do a=1,nvdf
                if ( psi1(a) /= psi2(a) ) then
                   ndif=ndif+1
                   if (ndif >= 4) EXIT 
                   dif(ndif)=a
                end if
             end do
 
             select case (ndif)
                case (1)
                   call calcHterm1b(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,Hci,cnf,psi1,psi2,dif)
                case (2)
                   call calcHterm2b(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
                case (3)
                   call calcHterm3(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
                case default
                   Hci(cnf)=0.0d0
             end select
 
 !            if (Hci(cnf) < 1d-8) Hci(cnf)=0d0
 
          end do
       end do
       end subroutine
 
 
       subroutine virtCI2(qumvia_nmc,Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,Emod,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nconf,ngaus,nvdf,qmaxx1,qmaxx2,bdim,ECI,Cci)
 
 !     -----------------------------------------------------------------
 !     SINGLE REFERENCE VIRTUAL CONFIGURATION INTERACTION
 !     This subroutine is the main program for virtual configuration
 !     interaction (VCI). 
 !     First a matrix representation of the hamiltonian in the VSCF 
 !     virtual wavefunctions (configurations) is built. Then it attempts
 !     to diagonalize it.
 !     To save memory, the hamiltonian is stored in packed format for
 !     symmetric matrices
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: qumvia_nmc ! Number of coupled modes.
       integer,intent(in)  :: ngaus ! Dimension of DGBS
       integer,intent(in)  :: nconf ! Dimension of CI basis set (VSCF configurations)
       integer,intent(in)  :: nvdf  ! Number of vib. deg. of freedom.
       integer,intent(in)  :: qmaxx1 ! Max quantum number allowed for sigles.
       integer,intent(in)  :: qmaxx2 ! Max quantum number allowed for doubles.
       integer,intent(in)  :: bdim  ! Dimension of operator matrices (qmaxx1+1 by default)
       real*8,intent(in)   :: Emod(ngaus,nvdf) ! VSCF modal energies.
       real*8,intent(in)   :: Po(ngaus,ngaus,nvdf) ! VSCF coeficients.
       real*8,intent(in)   :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian in DGB.
       real*8,intent(in)   :: Scho(ngaus,ngaus,nvdf) ! Core Hamiltonian in DGB.
       real*8,intent(in)   :: Q1(ngaus,ngaus,nvdf) ! Q^1 operator in DGB
       real*8,intent(in)   :: Q2(ngaus,ngaus,nvdf) ! Q^2 operator in DGB
       real*8,intent(in)   :: Q3(ngaus,ngaus,nvdf) ! Q^3 operator in DGB
       real*8,intent(in)   :: GDmtrx(ngaus,ngaus,nvdf) ! Effective 2mc potential
       real*8,intent(in)   :: GTmtrx(ngaus,ngaus,nvdf) ! Effective 3mc potential
       real*8,intent(in)   :: tiij(nvdf,nvdf)      !! 2-mode
       real*8,intent(in)   :: tjji(nvdf,nvdf)      !! Coupling
       real*8,intent(in)   :: uiiij(nvdf,nvdf)     !! potential 
       real*8,intent(in)   :: ujjji(nvdf,nvdf)     !! paramters
       real*8,intent(in)   :: uiijj(nvdf,nvdf)     !!
       real*8,intent(in)   :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(out)  :: ECI(nconf)         ! CI energy.
       real*8,intent(out)  :: Cci(nconf,nconf)   ! CI eigenvectors.
 
       integer :: i,j
       integer :: nmods
       integer :: ev
       integer :: cnf
       integer :: configs(4,nconf)
       real*8  :: Qm1(bdim,bdim,nvdf)
       real*8  :: Qm2(bdim,bdim,nvdf)
       real*8  :: Qm3(bdim,bdim,nvdf)
       real*8  :: Hmc(bdim,bdim,nvdf)
       real*8  :: GDm(bdim,bdim,nvdf)
       real*8  :: GTm(bdim,bdim,nvdf)
       real*8  :: Hci(nconf*(nconf+1)/2)
       real*8  :: Hprint(nconf,nconf)
       real*8  :: tol
       real*8  :: Evhf(nconf)
 !     DEBUG
       real*8  :: CIH(nconf,nconf)
 
 !     Workspace for diagonalization. 
       real*8,  dimension(:), allocatable :: WORK
       integer, dimension(:), allocatable :: IWORK,IFAIL
       integer :: INFO
 
       real*8,external  :: dnrm2
       real*8,external  :: dlamch
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
 !     -----------------------------------------------------------------
       nmods=qmaxx1+1
 
       write(77,'(A)') '------------'
       write(77,'(A)') 'STARTING VCI'
       write(77,'(A)') '------------'
 
 !     CHANGE OF BASIS FROM GAUSSIAN TO VSCF MODALS FOR Q AND Hcore
 !     OPERATORS.
       write(*,'(A)') 'BUILDING VCI OPERATORS   '
       call build_CIop3(Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,ngaus,nvdf,nmods,&
                     &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 !      call build_CIop2(Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,ngaus,nvdf,nmods,&
 !                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 
 !     Build the list of configurations to include in CI.
 !     nconf= 1 + nvdf*qmaxx1 + nvdf*(nvdf-1)*qmaxx2*qmaxx2/2
       call genConf2b(qmaxx1,qmaxx2,nvdf,nconf,configs)
 
 !     Building Diagonal Elements of CI Hamiltonian Matrix.
       Hci=0.0d0
       call Hci_diag2(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,tiij,tjji,uiiij,ujjji,uiijj,&
                      &tijk,uiijk,uijjk,uijkk,nvdf,ngaus,nmods,nconf,configs,Hci)
 
 !     Building off-diagonal elements of CI Hamiltonian matrix.
       call Hci_Offdiag2(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,tiij,tjji,uiiij,ujjji,&
                    &uiijj,tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,configs,Hci)
 !     ------------------------------------------------------------------
 !     DEBUG
 !     ------------------------------------------------------------------
       CIH=0d0
       do i=1,nconf
          do j=i,nconf
             cnf = i+j*(j-1)/2
             CIH(i,j) = Hci(cnf)
          end do
       end do
       do i=1,15
          write(77,'(999F12.7)') (CIH(j,i),j=1,15)
       end do
 !     ------------------------------------------------------------------
 
       do i=1,nconf
          Evhf(i)=Hci(i+i*(i-1)/2)
       end do
 
 !     Diagonalizing Hamiltonian
       write(*,'(A)') 'DIAGONALIZING VCI HAMILTONIAN'
       tol=2*dlamch('S')
       allocate ( WORK(8*nconf), IWORK(5*nconf) ,IFAIL(nconf))
       call dspevx('V','A','U',nconf,Hci,0.0d0,0.05d0,1,10*nvdf,&
              &tol,10*nvdf,Eci,Cci,nconf,WORK,IWORK,IFAIL,INFO)
 !     call dspevx('V','I','U',nconf,Hci,0.0d0,0.0d0,1,10*nvdf,&
 !            &0.0d0,10*nvdf,Eci,Cci,nconf,WORK,IWORK,IFAIL,INFO)
      
       if (INFO == 0) then
          do ev=1,nconf
 !         do ev=1,10*nvdf
             write(77,'(A)') '---------------------------------------------'
             write(77,'(A,I6,A,F18.8,F11.2)')'EVEC  ',ev,'Eci=',Eci(ev),(Eci(ev)-Eci(1))*h2cm
             write(77,'(A)') 'COEFICIENTS'
             do cnf=1,nconf
                if (abs(Cci(cnf,ev))>=0.20) then
                   write(77,'(I6,F12.5)') cnf, Cci(cnf,ev)
                end if 
             end do
          end do
          write(77,'(A)') 'No   CONFIGURATIONS         VSCF ENERGY VSCF TRANS.   VCI ENERGY  VCI TRANS  '
          do i=1,nconf
             write(77,'(I4,A,4I3,F22.8,F11.2,F14.8,F11.2)') i,' ',configs(:,i),&
                  &Evhf(i),(Evhf(i)-Evhf(1))*h2cm,  &
                  &Eci(i),(Eci(i)-Eci(1))*h2cm
          end do
       else 
          print*,'ERROR DURING CI HAMILTONIAN DIAGONALIZATION'
          print*,'INFO=',INFO
          print*,'IFAIL=',IFAIL
       end if
       deallocate ( WORK,IWORK,IFAIL )
 
       end subroutine
 
 !#######################################################################
 !     CONFIGURATION SELECTION VCI SECTION
 !#######################################################################

!       subroutine ssvscf_csvci(qumvia_nmc,qumvia_qff,ethresh,resthresh,selcut1,&
!                          &selcut2,ndf,nvdf,ngaus,nmcoup,nqmatoms,nclatoms,&
!                          &qmaxx1,qmaxx2,qmaxx3,qmaxx4,nconf,at_numbers,dy,gwidth,&
!                          &eig,nmodes,qmcoords,clcoords,naddref,nrst)
! !     -----------------------------------------------------------------
! !     THIS SUBROUTINE PERFORMS VIBRATIONAL SELF-CONSISTENT FIELD 
! !     FOLLOWED BY CONFIGURATION INTERACTION CALCULATION IN A 
! !     STATE-SPECIFIC FASHION.
! !     A separate VSCF and VCI computation is performed for each state
! !     of interest, by default the ground state and the first excited 
! !     states of each normal mode, using that state as the reference.
! !     The reference state is the vscf configuration over which the
! !     effective potential is computed.
! !     -----------------------------------------------------------------
!       implicit none
! 
!       integer,intent(in)  :: nrst
!       integer,intent(in)  :: qumvia_nmc           ! # of coupled mode in Hci
!       integer,intent(in)  :: qumvia_qff           ! # keyword for qff read/calc
!       integer,intent(in)  :: naddref              ! # of additional ref confs.
!       integer,intent(in)  :: ndf                  ! Number of classical atoms
!       integer,intent(in)  :: nvdf                 ! Number of classical atoms
!       integer,intent(in)  :: ngaus                ! Number of classical atoms
!       integer,intent(in)  :: nmcoup               ! Number of normal modes to couple in QFF.
!       integer,intent(in)  :: nqmatoms             ! Number of QM atoms
!       integer,intent(in)  :: nclatoms             ! Number of classical atoms
!       integer,intent(in)  :: qmaxx1               ! Max excitation for singles.
!       integer,intent(in)  :: qmaxx2               ! Max excitation for doules.
!       integer,intent(in)  :: qmaxx3               ! Max excitation for triples
!       integer,intent(in)  :: qmaxx4               ! Max excitation for quadruples
!       integer,intent(in)  :: nconf                ! Dimension of CI basis set.
!       integer,intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
!       real*8,intent(in)   :: dy                   ! Step size factor for num derivatives.
!       real*8,intent(in)   :: resthresh            ! Step size factor for num derivatives.
!       real*8,intent(in)   :: ethresh              ! Step size factor for num derivatives.
!       real*8,intent(in)   :: selcut1              ! Step size factor for num derivatives.
!       real*8,intent(in)   :: selcut2              ! Step size factor for num derivatives.
!       real*8,intent(in)   :: gwidth               ! Width factor for gaussian primitives.
!       real*8,intent(in)   :: eig(ndf)             ! Eigenvalues of the hessian matrix.
!       real*8,intent(in)   :: nmodes(ndf,ndf)      ! Eigenvectors of the hessian matrix.
!       real*8,intent(in)   :: qmcoords(3,nqmatoms) ! QM atom coordinates
!       real*8,intent(in)   :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
! 
!       integer            :: i
!       integer            :: ref(8)
!       integer            :: cnf
!       integer            :: bdim
!       integer            :: trash
!       integer            :: vscfstat
!       real*8             :: P(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
!       real*8             :: Po(ngaus,ngaus,nvdf)  ! Modals coeficient matrices.
!       real*8             :: Scho(ngaus,ngaus,nvdf)  ! Cholesky factor.
!       real*8             :: Evscf                ! VSCF energy.
!       real*8             :: Q1(ngaus,ngaus,nvdf) !  <gi|Qa^1|gj> matrix
!       real*8             :: Q2(ngaus,ngaus,nvdf) !  <gi|Qa^2|gj> matrix
!       real*8             :: Q3(ngaus,ngaus,nvdf) !  <gi|Qa^3|gj> matrix
!       real*8             :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian (non-orthogonal basis)
!       real*8             :: Essvhf(nrst)
!       real*8             :: Pss(ngaus,nrst)
!       real*8             :: Poss(ngaus,nrst)
!       real*8             :: GDmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
!       real*8             :: GTmtrx(ngaus,ngaus,nvdf)! Effective potential matrix.
!       real*8             :: Emod(ngaus,nvdf)       ! Modal energies.
!       real*8             :: hii(nvdf)             ! Diagonal Hessian eigenvalues.
!       real*8             :: tiii(nvdf)           ! Diagonal cubic coupling terms
!       real*8             :: uiiii(nvdf)          ! Diagonal quartic coupling terms
!       real*8             :: tiij(nvdf,nvdf)      ! 2-mode cubic coupling terms
!       real*8             :: tjji(nvdf,nvdf)      ! 2-mode cubic coupling terms
!       real*8             :: uiiij(nvdf,nvdf)     ! 2-mode quartic coupling terms
!       real*8             :: ujjji(nvdf,nvdf)     ! 2-mode quartic coupling terms
!       real*8             :: uiijj(nvdf,nvdf)     ! 2-mode quartic coupling terms
!       real*8             :: tijk(nvdf,nvdf,nvdf)      ! 3-mode cubic coupling terms
!       real*8             :: uiijk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
!       real*8             :: uijjk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
!       real*8             :: uijkk(nvdf,nvdf,nvdf)     ! 3-mode quartic coupling terms
!       real*8             :: Evci(nrst)
!       real*8             :: Eref
!       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
! 
! !     Allocatable variables.
!       integer,allocatable,dimension(:,:) :: addrefs
!       integer,allocatable,dimension(:,:) :: fundrefs
!       integer,allocatable,dimension(:,:) :: refstat
! 
! !     -----------------------------------------------------------------
! !     GENERATING REFERENCE FONFIGURATIONS
!      
!       write(77,'(A,I3)') 'NUMBER OF ADDITIONAL REFERENCE CONFIGS:',naddref
!       if (naddref > 0) then 
!          allocate ( addrefs(8,naddref),fundrefs(8,nvdf+1) )
!          addrefs=0
!          fundrefs=0
!          call readaddref(naddref,addrefs)
!          call genConf3(0,0,0,0,hii,ethresh,nvdf,nvdf+1,trash,fundrefs)
!          allocate ( refstat(8,nrst) )
!          refstat=0
!          refstat(:,1:nvdf+1) = fundrefs
!          refstat(:,nvdf+2:nrst) = addrefs
!          deallocate (addrefs,fundrefs)
!       else
!          allocate ( refstat(8,nrst) )
!          refstat=0
!          call genConf3(0,0,0,0,hii,ethresh,nvdf,nvdf+1,trash,refstat)
!       end if
!       write(77,'(A)') 'REFERENCE CONFIGURATIONS'
!       do i=1,nrst+naddref
!          write(77,'(8I3)') refstat(:,i)
!       end do
! 
! 
! !     QUARTIC FORCE FIELD
! !     ---------------------------------------------------------------
!       call selectqff(qumvia_qff,nmodes,eig,dy,nqmatoms,nclatoms,&
!                      &qmcoords,clcoords,at_numbers,ndf,nvdf,&
!                      &hii,tiii,tiij,tjji,uiiii,uiiij,ujjji,uiijj,&
!                      &tijk,uiijk,uijjk,uijkk)
! 
! !     COMPUTING VSCF/VCI OVER THESE CONFIGURATIONS      
!       do cnf=1,nrst
! !        Computing VSCF
!          ref(:) = refstat(:,cnf)
!          call ssvscf2(ref,nmodes,eig,dy,nmcoup,gwidth,nqmatoms,nclatoms,&
!                    &ngaus,qmcoords,clcoords,at_numbers,ndf,nvdf,&
!                    &P,Po,Scho,Evscf,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Emod,&
!                    &hii,tiii,uiiii,tiij,tjji,uiiij,ujjji,uiijj,&
!                    &tijk,uiijk,uijjk,uijkk,vscfstat)
!          if (vscfstat /=0) CYCLE
!          Essvhf(cnf)=Evscf
!          Pss(:,cnf)=P(:,cnf,cnf)
!          Poss(:,cnf)=Po(:,cnf,cnf)
! 
! !        Computing VCI
!          write(77,'(A,8I3)') 'COMPUTING VCI FOR REFERENCE STATE ',ref
!          bdim=qmaxx1+1
!          call csVCI(ref,qumvia_nmc,ethresh,resthresh,selcut1,selcut2,&
!                    &Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,Emod,&
!                    &hii,tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
!                    &nconf,ngaus,nvdf,qmaxx1,qmaxx2,qmaxx3,qmaxx4,bdim,Eref)
!          Evci(cnf)=Eref*h2cm
!       end do
! 
!       write(77,'(A)') 'ssVSCF AND csVCI ENERGIES AND TRANSITIONS (CM-1)'
!       write(77,'(A)') 'No.       E(ssVSCF)      dE(ssVSCF)        E(csVCI)        dE(csVCI)'
!       do cnf=1,nrst
!          write(77,'(I3,4F16.2)') cnf-1, Essvhf(cnf), Essvhf(cnf)-Essvhf(1),&
!                                 & Evci(cnf),Evci(cnf)-Evci(1)
!       end do
! 
! 
!       end subroutine
 
       subroutine readaddref(naddref,addrefs)
 !     ------------------------------------------------------------------
 !     THIS SUBROUTINE READS QUARTIC FORCE FIELD FROM FILE NAMED qff.qba
 !     ------------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: naddref
       integer,intent(out) :: addrefs(8,naddref)
 
       integer   :: i,j,k
       integer   :: nm1,nm2,nm3,nm4
       integer   :: qn1,qn2,qn3,qn4
       integer   :: openstatus
       integer   :: closestatus
       character*6 :: chaff
       real*8,parameter  :: amu2au = 1d0/0.0005485799111d0
 !     ------------------------------------------------------------------
 
       open(UNIT=15, FILE='addref.qva', ACTION='READ', IOSTAT=openstatus)
       if (openstatus /= 0) then
          write(77,'(A)') 'COULD NOT OPEN FILE addref.qva'
          STOP
       end if
 
       do i=1,naddref
          read(15,*) addrefs(:,i)
       end do
 
       close(15,iostat=closestatus)
       if (closestatus/=0) then
          write(77,'(A)') 'ERROR: COULD NOT CLOSE FILE qff.gqff'
          STOP
       end if 
 
       write(77,'(A)') '----------------------------------------------'
       write(77,'(A)') ' ADDITIONAL REFERENCE CONFIGURATIONS READ     '
       write(77,'(A)') '----------------------------------------------'
       do i=1,naddref
          write(77,'(8I3)') (addrefs(j,i),j=1,8)
       end do
 
       end subroutine
 
       subroutine csVCI(ref,qumvia_nmc,ethresh,resthresh,selcut1,selcut2,&
                    &Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,Emod,&
                    &hii,tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nconf,ngaus,nvdf,qmaxx1,qmaxx2,qmaxx3,qmaxx4,bdim,Eref)
 
 !     -----------------------------------------------------------------
 !     SINGLE REFERENCE VIRTUAL CONFIGURATION INTERACTION
 !     This subroutine is the main program for virtual configuration
 !     interaction (VCI). 
 !     First a matrix representation of the hamiltonian in the VSCF 
 !     virtual wavefunctions (configurations) is built. Then it attempts
 !     to diagonalize it.
 !     To save memory, the hamiltonian is stored in packed format for
 !     symmetric matrices
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: qumvia_nmc ! Number of coupled modes.
       integer,intent(in)  :: ref(8) ! Number of coupled modes.
       integer,intent(in)  :: ngaus ! Dimension of DGBS
       integer,intent(in)  :: nconf ! Dimension of CI basis set (VSCF configurations)
       integer,intent(in)  :: nvdf  ! Number of vib. deg. of freedom.
       integer,intent(in)  :: qmaxx1 ! Max quantum number allowed for sigles.
       integer,intent(in)  :: qmaxx2 ! Max quantum number allowed for doubles.
       integer,intent(in)  :: qmaxx3 ! Max quantum number allowed for triples.
       integer,intent(in)  :: qmaxx4 ! Max quantum number allowed for quadruples
       integer,intent(in)  :: bdim  ! Dimension of operator matrices (qmaxx1+1 by default)
       real*8,intent(in)   :: ethresh    ! Energy threshold for CS (cm-1)
       real*8,intent(in)   :: resthresh  ! Resonance threshold
       real*8,intent(in)   :: Emod(ngaus,nvdf) ! VSCF modal energies.
       real*8,intent(in)   :: Po(ngaus,ngaus,nvdf) ! VSCF coeficients.
       real*8,intent(in)   :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian in DGB.
       real*8,intent(in)   :: Scho(ngaus,ngaus,nvdf) ! Core Hamiltonian in DGB.
       real*8,intent(in)   :: Q1(ngaus,ngaus,nvdf) ! Q^1 operator in DGB
       real*8,intent(in)   :: Q2(ngaus,ngaus,nvdf) ! Q^2 operator in DGB
       real*8,intent(in)   :: Q3(ngaus,ngaus,nvdf) ! Q^3 operator in DGB
       real*8,intent(in)   :: GDmtrx(ngaus,ngaus,nvdf) ! Effective 2mc potential
       real*8,intent(in)   :: GTmtrx(ngaus,ngaus,nvdf) ! Effective 3mc potential
       real*8,intent(in)   :: hii(nvdf)
       real*8,intent(in)   :: tiij(nvdf,nvdf)      !! 2-mode
       real*8,intent(in)   :: tjji(nvdf,nvdf)      !! Coupling
       real*8,intent(in)   :: uiiij(nvdf,nvdf)     !! potential 
       real*8,intent(in)   :: ujjji(nvdf,nvdf)     !! paramters
       real*8,intent(in)   :: uiijj(nvdf,nvdf)     !!
       real*8,intent(in)   :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: selcut1
       real*8,intent(in)   :: selcut2
       real*8,intent(out)  :: Eref
 
       integer :: i,j
       integer :: ns
       integer :: nmods
       integer :: nsconf
       integer :: nsconf2
       integer :: ntrsh
       integer :: nref
       integer :: ev
       integer :: cnf
       integer :: err
       integer :: ml(1)
       real*8  :: Qm1(bdim,bdim,nvdf)
       real*8  :: Qm2(bdim,bdim,nvdf)
       real*8  :: Qm3(bdim,bdim,nvdf)
       real*8  :: Hmc(bdim,bdim,nvdf)
       real*8  :: GDm(bdim,bdim,nvdf)
       real*8  :: GTm(bdim,bdim,nvdf)
 !      real*8  :: Hprint(nconf,nconf)
       real*8  :: tol
 !     DEBUG
 !      real*8  :: CIH(nconf,nconf)
 
 !     Alocatable variables.
       real*8,dimension(:),allocatable   :: Hci
       real*8,dimension(:),allocatable   :: diag
       real*8,dimension(:),allocatable   :: ECI
       real*8,dimension(:),allocatable   :: Evhf
       real*8,dimension(:),allocatable   :: vec
       real*8,dimension(:,:),allocatable :: Cci
 
       integer,dimension(:),allocatable    :: sort
       integer,dimension(:),allocatable    :: gsel
       integer,dimension(:,:), allocatable :: configs
       integer,dimension(:,:), allocatable :: sconfigs
 
 !     Workspace for diagonalization. 
       real*8,  dimension(:), allocatable :: WORK
       integer, dimension(:), allocatable :: IWORK,IFAIL
       integer :: INFO
 
       real*8,external  :: dnrm2
       real*8,external  :: dlamch
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
 !     -----------------------------------------------------------------
       nmods=qmaxx1+1 ! Dimension of CI operator matrices.
 
       write(77,'(A)') '--------------'
       write(77,'(A)') 'STARTING csVCI'
       write(77,'(A)') '--------------'
 
 !     CHANGE OF BASIS FROM GAUSSIAN TO VSCF MODALS FOR Q AND Hcore
 !     OPERATORS.
       call build_CIop3(Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,ngaus,nvdf,nmods,&
                     &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
       write(77,'(A)') 'FINNISHED VCI OPERATORS   '
 
 !     Build the list of configurations to include in CI.
 !     nconf=1+
 !          +nvdf*qmaxx1+
 !          +(1/2)*nvdf*(nvdf-1)*qmaxx2**2+
 !          +(1/6)*nvdf*(nvdf-1)*(nvdf-2)*qmaxx3**3
       allocate(configs(8,nconf))
 
 
       call genConf3(qmaxx1,qmaxx2,qmaxx3,qmaxx4,hii,ethresh,nvdf,nconf,nsconf,configs)
 
 
       allocate(sconfigs(8,nsconf))
       sconfigs = configs(:,1:nsconf)
 
 
       if (allocated(configs)) deallocate(configs,stat=err)
       if (err /= 0) STOP ('ALLOCATION ERROR')
 
       allocate(configs(8,nsconf))
       configs = sconfigs
 
       if (allocated(sconfigs)) deallocate(sconfigs,stat=err)
       if (err/=0) STOP('ALLOCATION ERROR: deallocating sconfigs')
       
       write(77,'(A,I15)') 'Initial configurations=',nconf
       write(77,'(A,I15)') 'Configs after HO energy threshold =',nsconf
 
 !     Building Diagonal Elements of CI Hamiltonian Matrix and imposing the
 !     energy threshold.
       allocate (diag(nsconf))
       call energythresh(ref,qumvia_nmc,ethresh,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,&
                      &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                      &nvdf,ngaus,nmods,nsconf,configs,diag,ntrsh,nref)
 
       write(77,'(A,I15)') 'Configs after diagonal energy threshold =',ntrsh
 
 !     Recursive configuration selection algorithm.
       allocate (gsel(ntrsh))
       call configsel(ref,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,&
                      &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                      &nvdf,ngaus,nmods,nsconf,configs,diag,ntrsh,nref,&
                      &gsel,selcut1,selcut2,nsconf2)
 
       write(77,'(A,I15)') 'Configs after recursive CS =',nsconf2
 
 !     Building diagonal Hamiltonian and new list of (selected) configurations.
       allocate (sconfigs(8,nsconf2))
       allocate (Hci(nsconf2+nsconf2*(nsconf2-1)/2))
       ns=1
       Hci=0d0
       do i=1,ntrsh
          if (gsel(i)==1) then
             Hci(ns*(ns+1)/2)=diag(i)
             sconfigs(:,ns)=configs(:,i)
             if (i==nref) nref=ns
             ns=ns+1
          end if
       end do
       deallocate(configs, gsel, diag)
 
       allocate (Evhf(nsconf2))
       do i=1,nsconf2
          Evhf(i)=Hci(i+i*(i-1)/2)
       end do
 
 !     Building off-diagonal elements of CI Hamiltonian matrix.
       call csvci_mtrx(qumvia_nmc,nvdf,nmods,nsconf2,sconfigs,Hci,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 
 
 !     Diagonalizing Hamiltonian
       write(*,'(A)') 'DIAGONALIZING VCI HAMILTONIAN'
       tol=2*dlamch('S')
       allocate(Eci(nsconf2),Cci(nsconf2,nsconf2))
       allocate ( WORK(8*nsconf), IWORK(5*nsconf) ,IFAIL(nsconf))
       call dspevx('V','A','U',nsconf2,Hci,0.0d0,0.05d0,1,10*nvdf,&
              &tol,10*nvdf,Eci,Cci,nsconf2,WORK,IWORK,IFAIL,INFO)
 !     call dspevx('V','I','U',nconf,Hci,0.0d0,0.0d0,1,10*nvdf,&
 !            &0.0d0,10*nvdf,Eci,Cci,nconf,WORK,IWORK,IFAIL,INFO)
 
       allocate (sort(nsconf2),vec(nsconf2))
       if (INFO == 0) then
          do ev=1,nsconf2
 !         do ev=1,10*nvdf
             ml=maxloc(abs(Cci(:,ev)))
             if (ml(1)==nref) Eref=Eci(ev)
             write(77,'(A)') '---------------------------------------------'
             write(77,'(A,I6,A,F18.8,F11.2)')'EVEC  ',ev,'Eci=',Eci(ev),(Eci(ev)-Eci(1))*h2cm
             vec = Abs(Cci(:,ev))
             call quick_sort(vec,sort,nsconf2)
             write(77,'(A,8I3)') 'MAIN CONFIGURATION:',sconfigs(:,sort(nsconf2))
             write(77,'(A)') 'COEFICIENTS'
             do cnf=1,nsconf2
                if (abs(Cci(sort(nsconf2-cnf+1),ev))<=0.05) EXIT
                write(77,'(8I3,F12.5)') sconfigs(:,sort(nsconf2-cnf+1)), Cci(sort(nsconf2-cnf+1),ev)
             end do
          end do
          write(77,'(A)') 'No   CONFIGURATIONS         VSCF ENERGY VSCF TRANS.   VCI ENERGY  VCI TRANS  '
 !         do i=1,nsconf2
 !            write(77,'(I4,A,6I3,F22.8,F11.2,F14.8,F11.2)') i,' ',sconfigs(:,i),&
 !                 &Evhf(i),(Evhf(i)-Evhf(1))*h2cm,  &
 !                 &Eci(i),(Eci(i)-Eci(1))*h2cm
 !         end do
       else 
          print*,'ERROR DURING CI HAMILTONIAN DIAGONALIZATION'
          print*,'INFO=',INFO
          print*,'IFAIL=',IFAIL
       end if
       deallocate ( WORK,IWORK,IFAIL,Eci,Cci,Evhf,Hci )
 
       end subroutine
 
       subroutine genConf3(qmaxx1,qmaxx2,qmaxx3,qmaxx4,hii,ethresh,nvdf,nconf,nsconf,configs)
 !     -----------------------------------------------------------------
 !     THIS SUBROUTINE GENERATES ALL CI SINGLE AND DOUBLE CONFIGURATIONS
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: qmaxx1               ! Max excitation singles
       integer,intent(in)  :: qmaxx2               ! Max excitation doubles
       integer,intent(in)  :: qmaxx3               ! Max excitation triples
       integer,intent(in)  :: qmaxx4               ! Max excitation quadruples
       integer,intent(in)  :: nvdf                ! No of vibrational dof.
       integer,intent(in)  :: nconf               ! CI basis set size.
       integer,intent(out) :: nsconf              ! CI selected basis set size.
       integer,intent(out) :: configs(8,nconf)    ! Configurations.
       real*8,intent(in)   :: ethresh
       real*8,intent(in)   :: hii(nvdf)
 
       integer    :: i
       integer    :: nc
       integer    :: nm
       integer    :: qn
       integer    :: nm1
       integer    :: nm2
       integer    :: nm3
       integer    :: nm4
       integer    :: qn1
       integer    :: qn2
       integer    :: qn3
       integer    :: qn4
       real*8     :: zpe
       real*8     :: ecut
       real*8     :: hoe
       real*8     :: hof(nvdf)
 
       real*8, parameter :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 !     Compute Harmoni Oscillator frequencies in cm-1
       hof=sqrt(abs(hii))*h2cm  
 
 !     Compute Harmonic Oscillator zero point energy
       zpe=0d0
       do i=1,nvdf
          zpe=zpe+hof(i)/2d0
       end do
       ecut=zpe+ethresh+5000d0
       write(77,'(A)') 'HARMONIC OSCILLATOR TRANSITIONS'
       write(77,'(99F16.8)') hof
       write(77,'(A)') 'zpe         ethresh         ecut'
       write(77,'(99F16.2)') zpe, ethresh, ecut
 
 !     Total number of configurations
 !     nconf = 1 + nvdf*qmaxx1 + nvdf*(nvdf-1)*qmaxx2*qmaxx2/2 +
 !     nvdf*(nvdf-1)*(nvdf-2)*qmaxx3**3/6
 
 
 !     GROUND STATE
 !     -----------------------------------------------------------------
       nc=1
       configs(:,nc)=0
 
 !     FUNDAMENTALS
 !     -----------------------------------------------------------------
       do nm=1,nvdf
          nc=nc+1
          configs(:,nc) = (/nm,1,0,0,0,0,0,0/)
       end do
 
       if (qmaxx1 == 0) then
          nsconf=nc
          RETURN
       end if
 
 !     OVERTONES
 !     -----------------------------------------------------------------
       do nm=1,nvdf
          do qn=2,qmaxx1
             hoe=zpe+qn*hof(nm)
             if (hoe < ecut) then
 !               write(77,'(F16.8)') hoe
                nc=nc+1
                configs(:,nc) = (/nm,qn,0,0,0,0,0,0/)
             end if
          end do
       end do
 
       if (qmaxx2 == 0) then
          nsconf=nc
          RETURN
       end if
 
 !     DOUBLES
 !     -----------------------------------------------------------------
       do nm1=1,nvdf-1
          do nm2=nm1+1,nvdf
             do qn1=1,qmaxx2
                do qn2=1,qmaxx2
                   hoe=zpe+qn1*hof(nm1)+qn2*hof(nm2)
                   if (hoe < ecut) then
 !                     write(77,'(F16.8)') hoe
                      nc=nc+1
                      configs(:,nc)=(/nm1,qn1,nm2,qn2,0,0,0,0/)
                   end if
                end do
             end do
          end do
       end do
 
 !     TRIPLES
 !     -----------------------------------------------------------------
       do nm1=1,nvdf-2
          do nm2=nm1+1,nvdf-1
             do nm3=nm2+1,nvdf
                do qn1=1,qmaxx3
                   do qn2=1,qmaxx3
                      do qn3=1,qmaxx3
                         hoe=zpe+qn1*hof(nm1)+qn2*hof(nm2)+qn3*hof(nm3)
                         if (hoe < ecut) then
 !                           write(77,'(F16.8)') hoe
                            nc=nc+1
                            configs(:,nc)=(/nm1,qn1,nm2,qn2,nm3,qn3,0,0/)
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do
 
       if (nvdf<4) then
          nsconf=nc
          RETURN
       end if
          
 
 !     QUADRUPLES
 !     -----------------------------------------------------------------
       do nm1=1,nvdf-3
          do nm2=nm1+1,nvdf-2
             do nm3=nm2+1,nvdf-1
                do nm4=nm3+1,nvdf
                   do qn1=1,qmaxx4
                      do qn2=1,qmaxx4
                         do qn3=1,qmaxx4
                            do qn4=1,qmaxx4
                               hoe=zpe+qn1*hof(nm1)+qn2*hof(nm2)+qn3*hof(nm3)+qn4*hof(nm4)
                               if (hoe < ecut) then
                                  nc=nc+1
                                  configs(:,nc)=(/nm1,qn1,nm2,qn2,nm3,qn3,nm4,qn4/)
                               end if
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
 
 
 !     TOTAL NUMBER OF SELECTED CONFIGURATIONS
 !     -----------------------------------------------------------------
       nsconf=nc
 
 
       end subroutine
 
 
 
 
 
 
       subroutine energythresh(ref,qumvia_nmc,ethresh,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,&
                      &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                      &nvdf,ngaus,nmods,nconf,configs,diag,ntrsh,nref)
 !     -----------------------------------------------------------------
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc ! # of coupled modes.
       integer,intent(in)   :: nvdf   ! # of vibrational deg of freedom.
       integer,intent(in)   :: nmods  ! # of selected VSCF virtual states.
       integer,intent(in)   :: nconf  ! Total number of configurations
       integer,intent(in)   :: ngaus  ! Dimension of gaussian basis set
       integer,intent(in)   :: ref(8)  ! Reference configuration
       integer,intent(out)  :: ntrsh
       integer,intent(out)  :: nref
       integer,intent(inout):: configs(8,nconf)  ! List of configurations.
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf) !Q operator in modals basis
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf) !Q^2 operator in modals basis
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf) !Q^3 operator in modals basis
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8,intent(in)    :: tiij(nvdf,nvdf)   ! Coupling potential parameters
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: ethresh ! Energy threshold (rel 2 grnd,cm-1)
       real*8,intent(out)   :: diag(nconf)
 
       integer :: i,a,b,c
       integer :: nm1,nm2,nm3,nm4
       integer :: qn1,qn2,qn3,qn4
       integer :: cnf
       integer :: sref
       integer :: n1,n2,n3 ! Various indices.
       integer :: psi(nvdf) ! VSCF wavefunction. Array of vib. quantum numbers.
       integer :: sort(nconf)
       real*8 :: Vc2       ! Auxiliary variable for Hci calculation.
       real*8 :: Vc3       ! Auxiliary variable for Hci calculation.
       real*8 :: ehf       ! Auxiliary variable for Hci calculation.
       real*8 :: veff      ! Auxiliary variable for Hci calculation.
       real*8 :: veff2     ! Auxiliary variable for Hci calculation.
       real*8 :: veff3     ! Auxiliary variable for Hci calculation.
       real*8 :: absethresh
       real*8 :: Ehf0      ! DEBUG
       real*8 :: corr      ! DEBUG
 
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
 !     -----------------------------------------------------------------
       
       write(77,'(A)') 'BEGINNING ENERGYTHRESH'
       do cnf=1,nconf
          nm1=configs(1,cnf)
          qn1=configs(2,cnf)
          nm2=configs(3,cnf)
          qn2=configs(4,cnf)
          nm3=configs(5,cnf)
          qn3=configs(6,cnf)
          nm4=configs(7,cnf)
          qn4=configs(8,cnf)
 
          psi=0
          if (nm1 > 0) psi(nm1)=qn1
          if (nm2 > 0) psi(nm2)=qn2
          if (nm3 > 0) psi(nm3)=qn3
          if (nm4 > 0) psi(nm4)=qn4
 
 !        VSCF energy.
          ehf=0d0
          do a=1,nvdf
             ehf = ehf + Emod(psi(a)+1,a)
          end do
 
 !        Average on effective potential.
          veff2=0d0
          do a=1,nvdf
             veff2 = veff2 + GDm(psi(a)+1,psi(a)+1,a)
          end do
          veff3=0d0
          do a=1,nvdf
             veff3 = veff3 + GTm(psi(a)+1,psi(a)+1,a)
          end do
          veff = veff2 + veff3
 
 !        2-MODE COUPLING POTENTIAL 
          Vc2=0.0d0
          do a=1,nvdf-1
             do b=a+1,nvdf
                n1=psi(a)+1
                n2=psi(b)+1
                Vc2=Vc2+&
                & tiij(a,b)*Qm2(n1,n1,a)*Qm1(n2,n2,b)/2.0d0+&
                & tjji(a,b)*Qm1(n1,n1,a)*Qm2(n2,n2,b)/2.0d0+&
                & uiiij(a,b)*Qm3(n1,n1,a)*Qm1(n2,n2,b)/6.0d0+&
                & ujjji(a,b)*Qm1(n1,n1,a)*Qm3(n2,n2,b)/6.0d0+&
                & uiijj(a,b)*Qm2(n1,n1,a)*Qm2(n2,n2,b)/4.0d0  
             end do
          end do
 
 !        3-MODE COUPLING POTENTIAL
          Vc3=0d0
          IF (qumvia_nmc == 3) THEN
          do a=1,nvdf-2
             do b=a+1,nvdf-1
                do c=b+1,nvdf
                   n1=psi(a)+1
                   n2=psi(b)+1
                   n3=psi(c)+1
 
                   Vc3=Vc3+&
                    & tijk(a,b,c)*Qm1(n1,n1,a)*Qm1(n2,n2,b)*Qm1(n3,n3,c)+&
                    & uiijk(a,b,c)*Qm2(n1,n1,a)*Qm1(n2,n2,b)*Qm1(n3,n3,c)/2d0+&
                    & uijjk(a,b,c)*Qm1(n1,n1,a)*Qm2(n2,n2,b)*Qm1(n3,n3,c)/2d0+&
                    & uijkk(a,b,c)*Qm1(n1,n1,a)*Qm1(n2,n2,b)*Qm2(n3,n3,c)/2d0
 
                end do
             end do
          end do
          END IF
 
 !        DANGER
          diag(cnf)=ehf-Vc2-2.0d0*Vc3
 !         diag(cnf)=ehf-veff+Vc2+Vc3
          
       end do
       
 !     DEBUG------------------------------------------------
 !      write(77,'(A)') 'DIAGONAL HAMILTONIAN'
 !      do i=1,nconf
 !         write(77,'(F16.8)') diag(i)
 !      end do
 !     -----------------------------------------------------
 
       sort=0
       call quick_sort(diag,sort,nconf)
       configs = configs(:,sort)
       
 !     DEBUG------------------------------------------------
 !      write(77,'(A)') 'SORTED DIAGONAL HAMILTONIAN'
 !      do i=1,nconf
 !         write(77,'(F16.8)') diag(i)
 !      end do
 !
 !      write(77,'(A)') 'SORTED CONFIGS'
 !      do i=1,nconf
 !         write(77,'(6I2)') configs(:,i)
 !      end do
 !     -----------------------------------------------------
 
 !     Computing energy threshold.
       do i=1,nconf
          if ( equal(ref, configs(:,i)) ) then
             nref=i
             absethresh=diag(nref)+ethresh/h2cm
             exit
          end if
       end do
 
       write(77,'(A,I14)') 'REFERENCE CONFIG:',nref
       write(77,'(A,F16.8)') 'REFERENCE ENERGY:',diag(nref)
       write(77,'(A,F16.8)') 'ENERGY THRESHOLD IN CM-1:',ethresh
       write(77,'(A,F16.8)') 'ENERGY THRESHOLD IN AU:',ethresh/h2cm
       write(77,'(A,F16.8)') 'Abs ENERGY THRESHOLD IN AU:',absethresh
 
 !     Determining the number of configurations selected by energy threshold.
       if ( diag(nconf) < absethresh) then
          ntrsh=nconf
       else
          do i=1,nconf
             if (diag(i) > absethresh) then
                ntrsh=i-1
                exit
             end if
          end do
       end if
 
       write(77,'(A,I14)') 'UPPER LIMIT CONFIG:',ntrsh
 
       end subroutine
 
 
       logical function equal( array1, array2 )
       integer, dimension(:), intent(in) :: array1, array2
       integer :: i
       
       equal =size(array1) == size(array2)
       if ( equal ) then
       do i = 1,size(array1)
       equal = array1(i) == array2(i)
       if ( .not. equal )exit
       enddo
       endif
       return
       end function equal
 
 
       subroutine configsel(ref,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,&
                      &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                      &nvdf,ngaus,nmods,nconf,configs,diag,ntrsh,nref,&
                      &gsel,selcut1,selcut2,nsconf)
 !     -----------------------------------------------------------------
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc ! # of coupled modes.
       integer,intent(in)   :: nvdf   ! # of vibrational deg of freedom.
       integer,intent(in)   :: nmods  ! # of selected VSCF virtual states.
       integer,intent(in)   :: nconf  ! Total number of configurations
       integer,intent(in)   :: ngaus  ! Dimension of gaussian basis set
       integer,intent(in)   :: configs(8,nconf)  ! List of configurations.
       integer,intent(in)   :: ref(8)  ! Reference configuration
       integer,intent(in)   :: ntrsh
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf) !Q operator in modals basis
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf) !Q^2 operator in modals basis
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf) !Q^3 operator in modals basis
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8,intent(in)    :: tiij(nvdf,nvdf)   ! Coupling potential parameters
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: diag(nconf)
       real*8,intent(in)    :: selcut1
       real*8,intent(in)    :: selcut2
       integer,intent(out)  :: nsconf
       integer,intent(out)  :: gsel(ntrsh) ! General selections
       integer,intent(inout):: nref
 
       integer :: i
       integer :: sref
       integer :: sel1(ntrsh)
       integer :: sel2(ntrsh)
 !     -----------------------------------------------------------------
 
       write(77,'(A)') 'STARTING CONFIGSEL'
       gsel = 0
       sel1 = 0
       sel2 = 0
       gsel(nref) = 1
 
       nsconf=1
       call pcconfsel(nref,ntrsh,diag,gsel,sel1,nsconf,selcut1,qumvia_nmc,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,configs)
 
       do i=1,ntrsh
          if (sel1(i) == 0) CYCLE
          sref=i
          call pcconfsel(sref,ntrsh,diag,gsel,sel2,nsconf,selcut2,qumvia_nmc,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,configs)
       end do
 
       nsconf=0
       do i=1,ntrsh
          if (gsel(i)==1) nsconf=nsconf+1
          if (i==nref) nref=i
       end do
 
       end subroutine
 
 
 
 
       subroutine pcconfsel(nref,ntrsh,diag,gsel,rsel,nsconf,selcut,qumvia_nmc,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,configs)
 !     -----------------------------------------------------------------
 !     Computes Off-diagonal Hamiltonian matrix elements 
 !                      < K | H | L >
 !     where | K > and | L > may differ in 1 (case 1), 2 (case 2) and 
 !     more than 2 (case 3) modals. This subroutine runs over all off-
 !     diagonal matrix elements, determines to which case each belong
 !     and then computes the integral using subroutines calcHterm1 
 !     for case 1 matrix elements, calcHterm2 for case 2 ME or
 !     sets Hci(cnf)=0.0 for case 3.
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: nref,ntrsh
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf
       integer,intent(in)   :: configs(8,nconf)
       integer,intent(out)  :: rsel(ntrsh) ! Selections in this round
       integer,intent(inout):: gsel(ntrsh) ! General selections
       integer,intent(inout):: nsconf ! # of selected configurations
       real*8,intent(in)    :: selcut      ! Energy cutoff for selection
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: diag(nconf)
 
       integer    :: a
       integer    :: b
       integer    :: i
       integer    :: j
       integer    :: r
       integer    :: cnf
       integer    :: nm1
       integer    :: nm2
       integer    :: nm3
       integer    :: nm4
       integer    :: nm5
       integer    :: nm6
       integer    :: nm7
       integer    :: nm8
       integer    :: qn1
       integer    :: qn2
       integer    :: qn3
       integer    :: qn4
       integer    :: qn5
       integer    :: qn6
       integer    :: qn7
       integer    :: qn8
       integer    :: ndif
       integer    :: dif(3)
       integer    :: psi1(nvdf)
       integer    :: psi2(nvdf)
       real*8     :: CSPCE
 !     -----------------------------------------------------------------
 
       do r=1,ntrsh
          if (gsel(r) == 1) CYCLE
          if (r<nref) then
             i=r
             j=nref
          else 
             i=nref
             j=r
          end if
 
          nm1=configs(1,i)
          qn1=configs(2,i)
          nm2=configs(3,i)
          qn2=configs(4,i)
          nm3=configs(5,i)
          qn3=configs(6,i)
          nm4=configs(7,i)
          qn4=configs(8,i)
 
          nm5=configs(1,j)
          qn5=configs(2,j)
          nm6=configs(3,j)
          qn6=configs(4,j)
          nm7=configs(5,j)
          qn7=configs(6,j)
          nm8=configs(7,j)
          qn8=configs(8,j)
 
          psi1=0
          psi2=0
 
          if (nm1 > 0) psi1(nm1)=qn1
          if (nm2 > 0) psi1(nm2)=qn2
          if (nm3 > 0) psi1(nm3)=qn3
          if (nm4 > 0) psi1(nm4)=qn4
 
          if (nm5 > 0) psi2(nm5)=qn5
          if (nm6 > 0) psi2(nm6)=qn6
          if (nm7 > 0) psi2(nm7)=qn7
          if (nm8 > 0) psi2(nm8)=qn8
 
          dif=0
          ndif=0
          do a=1,nvdf
             if ( psi1(a) /= psi2(a) ) then
                ndif=ndif+1
                if (ndif >= 4) EXIT 
                dif(ndif)=a
             end if
          end do
 
          select case (ndif)
             case (1)
               call cspaircorr1(nref,r,diag,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
               &nvdf,nmods,nconf,CSPCE,cnf,psi1,psi2,dif)
             case (2)
               call cspaircorr2(nref,r,diag,qumvia_nmc,Qm1,Qm2,Qm3,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
               &nvdf,nmods,nconf,CSPCE,cnf,psi1,psi2,dif)
             case (3)
               call cspaircorr3(nref,r,diag,qumvia_nmc,Qm1,Qm2,Qm3,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
               &nvdf,nmods,nconf,CSPCE,cnf,psi1,psi2,dif)
             case default
               CSPCE=0D0
          end select
 
 !        DEBUG--------------------------------
 !         write(77,'(A,F16.8)') 'CSPE=',CSPCE
 !        -------------------------------------
 
 !        If pair correlation energy is greater than cutoff add this
 !        configuration to the selection.
          if (abs(CSPCE) > selcut) then
             gsel(r) = 1
             rsel(r) = 1
             nsconf=nsconf+1
          end if
 
       end do
       end subroutine
 
 
       subroutine cspaircorr1(nref,r,diag,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,CSPCE,cnf,psi1,psi2,dif)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY ONE MODAL.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: nref
       integer,intent(in)   :: r
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf,cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: diag(nconf)
       real*8,intent(out)   :: CSPCE
 
       integer :: i,j,k
       integer :: a,b
       integer :: n1,n2
       integer :: ni1,ni2
       integer :: nj1,nj2
       integer :: nk1,nk2
       real*8 :: Vc2
       real*8 :: Vc3
       real*8 :: veff      ! Auxiliary variable for Hci calculation.
       real*8 :: dE
       real*8, PARAMETER :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 
 
 !     Computing VSCF effective potential.
       veff=0d0
       n1=psi1(dif(1))+1
       n2=psi2(dif(1))+1
       veff=GDm(n1,n2,dif(1))+GTm(n1,n2,dif(1))
 
 !     Computing Coupling Potential
       Vc2 = 0.0d0
       do a=1,nvdf
          if (a == dif(1)) CYCLE
 
          if (a<dif(1)) then
             i=a
             j=dif(1)
             ni1=psi1(a)+1
             ni2=ni1
             nj1=psi1(dif(1))+1
             nj2=psi2(dif(1))+1
          else
             i=dif(1)
             j=a
             ni1=psi1(dif(1))+1
             ni2=psi2(dif(1))+1
             nj1=psi1(a)+1
             nj2=nj1
          end if 
              
          Vc2 = Vc2 + &
               & tiij(i,j)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)/2.0d0 + &
               & tjji(i,j)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)/2.0d0 + &
               & uiiij(i,j)*Qm3(ni1,ni2,i)*Qm1(nj1,nj2,j)/6.0d0+ &
               & ujjji(i,j)*Qm1(ni1,ni2,i)*Qm3(nj1,nj2,j)/6.0d0+ &
               & uiijj(i,j)*Qm2(ni1,ni2,i)*Qm2(nj1,nj2,j)/4.0d0
       end do
 
       Vc3 = 0d0
       IF (qumvia_nmc == 3) THEN
       do a=1,nvdf-1
          do b=a+1,nvdf
 
             if (a == dif(1)) CYCLE
             if (b == dif(1)) CYCLE
 
             if (dif(1) < a) then
                i=dif(1)
                j=a
                k=b
 
                ni1=psi1(dif(1))+1
                ni2=psi2(dif(1))+1
                nj1=psi1(a)+1
                nj2=psi2(a)+1
                nk1=psi1(b)+1
                nk2=psi2(b)+1
             else if (a<dif(1) .AND. dif(1)<b) then
                i=a
                j=dif(1)
                k=b
 
                ni1=psi1(a)+1
                ni2=psi2(a)+1
                nj1=psi1(dif(1))+1
                nj2=psi2(dif(1))+1
                nk1=psi1(b)+1
                nk2=psi2(b)+1
             else if (b<dif(1)) then
                i=a
                j=b
                k=dif(1)
 
                ni1=psi1(a)+1
                ni2=psi2(a)+1
                nj1=psi1(b)+1
                nj2=psi2(b)+1
                nk1=psi1(dif(1))+1
                nk2=psi2(dif(1))+1
             end if
 
             Vc3=Vc3+&
              & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)+ &
              & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0+ &
              & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0+ &
              & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
          end do
       end do 
       END IF
 
 
 !     COMPUTING HAMILTONIAN MATRIX ELEMENT
       dE=(diag(r)-diag(nref))*h2cm
       CSPCE=(Vc2 + Vc3 - veff)*h2cm
       CSPCE=CSPCE**2/dE
 
       end subroutine
 
       subroutine cspaircorr2(nref,r,diag,qumvia_nmc,Qm1,Qm2,Qm3,&
                   &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                   &nvdf,nmods,nconf,CSPCE,cnf,psi1,psi2,dif)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY TWO MODALS.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: nref
       integer,intent(in)   :: r
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf
       integer,intent(in)   :: nmods
       integer,intent(in)   :: nconf
       integer,intent(in)   :: cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: diag(nconf)
       real*8,intent(out)   :: CSPCE
 
       integer :: i,j,k
       integer :: a,b,c
       integer :: ni1,ni2,nj1,nj2,nk1,nk2
       real*8  :: Vc2,Vc3
       real*8  :: dE
       real*8, PARAMETER :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 
 
 
       if (dif(1) < dif(2)) then
          i=dif(1)
          j=dif(2)
          ni1=psi1(dif(1))+1
          ni2=psi2(dif(1))+1
          nj1=psi1(dif(2))+1
          nj2=psi2(dif(2))+1
       else
          i=dif(2)
          j=dif(1)
          ni1=psi1(dif(2))+1
          ni2=psi2(dif(2))+1
          nj1=psi1(dif(1))+1
          nj2=psi2(dif(1))+1
       end if 
           
       Vc2 = 0d0
       Vc2 =  &
          & tiij(i,j)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)/2.0d0+&
          & tjji(i,j)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)/2.0d0+&
          & uiiij(i,j)*Qm3(ni1,ni2,i)*Qm1(nj1,nj2,j)/6.0d0+&
          & ujjji(i,j)*Qm1(ni1,ni2,i)*Qm3(nj1,nj2,j)/6.0d0+&
          & uiijj(i,j)*Qm2(ni1,ni2,i)*Qm2(nj1,nj2,j)/4.0d0
 
       Vc3 = 0d0
       IF (qumvia_nmc == 3) THEN
       do a=1,nvdf
 
          if (a == dif(1)) CYCLE
          if (a == dif(2)) CYCLE
 
          if (a < dif(1) .AND. dif(1) < dif(2)) then
             i=a
             j=dif(1)
             k=dif(2)
 
             ni1=psi1(i)+1
             ni2=psi2(i)+1
             nj1=psi1(j)+1
             nj2=psi2(j)+1
             nk1=psi1(k)+1
             nk2=psi2(k)+1
          else if (dif(1) < a .AND. a < dif(2)) then
             i=dif(1)
             j=a
             k=dif(2)
 
             ni1=psi1(i)+1
             ni2=psi2(i)+1
             nj1=psi1(j)+1
             nj2=psi2(j)+1
             nk1=psi1(k)+1
             nk2=psi2(k)+1
          else if ( dif(2) < a ) then
             i=dif(1)
             j=dif(2)
             k=a
 
             ni1=psi1(i)+1
             ni2=psi2(i)+1
             nj1=psi1(j)+1
             nj2=psi2(j)+1
             nk1=psi1(k)+1
             nk2=psi2(k)+1
          end if
 
          Vc3=Vc3+ &
           & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k) + &
           & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
       end do
       END IF 
    
 !     VCI MATRIX ELEMENT
       dE=(diag(r)-diag(nref))*h2cm
       CSPCE = (Vc2 + Vc3)*h2cm
       CSPCE = CSPCE**2/dE
 
       end subroutine
 
 
 
 
       subroutine cspaircorr3(nref,r,diag,qumvia_nmc,Qm1,Qm2,Qm3,&
                   &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                   &nvdf,nmods,nconf,CSPCE,cnf,psi1,psi2,dif)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY TWO MODALS.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: nref
       integer,intent(in)   :: r
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf
       integer,intent(in)   :: nmods
       integer,intent(in)   :: nconf
       integer,intent(in)   :: cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: diag(nconf)
       real*8,intent(out)   :: CSPCE
 
       integer :: i,j,k
       integer :: ni1,ni2,nj1,nj2,nk1,nk2
       real*8  :: dE
       real*8, PARAMETER :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 
 
       IF (qumvia_nmc == 3) THEN
 !        COMPUTING 3-MODE COUPLING TERMS
          if (dif(1) > dif(2)) STOP ('dif(1) > dif(2) in calcHterm3')
          if (dif(2) > dif(3)) STOP ('dif(2) > dif(3) in calcHterm3')
          if (dif(1) > dif(3)) STOP ('dif(1) > dif(3) in calcHterm3')
    
          i=dif(1)
          j=dif(2)
          k=dif(3)
    
          ni1=psi1(i)+1
          ni2=psi2(i)+1
          nj1=psi1(j)+1
          nj2=psi2(j)+1
          nk1=psi1(k)+1
          nk2=psi2(k)+1
    
          CSPCE= &
           & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k) + &
           & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
          dE=(diag(r)-diag(nref))*h2cm
          CSPCE=(CSPCE*h2cm)**2/dE
       ELSE
          CSPCE = 0d0
       END IF
 
       end subroutine
 
       subroutine csvci_mtrx(qumvia_nmc,nvdf,nmods,nconf,configs,Hci,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 !     -----------------------------------------------------------------
 !     Computes Off-diagonal Hamiltonian matrix elements 
 !                      < K | H | L >
 !     where | K > and | L > may differ in 1 (case 1), 2 (case 2) and 
 !     more than 2 (case 3) modals. This subroutine runs over all off-
 !     diagonal matrix elements, determines to which case each belong
 !     and then computes the integral using subroutines calcHterm1 
 !     for case 1 matrix elements, calcHterm2 for case 2 ME or
 !     sets Hci(cnf)=0.0 for case 3.
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf
       integer,intent(in)   :: configs(8,nconf)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(inout) :: Hci(nconf+nconf*(nconf-1)/2)
 
       integer    :: a
       integer    :: b
       integer    :: i
       integer    :: j
       integer    :: cnf
       integer    :: nm1
       integer    :: nm2
       integer    :: nm3
       integer    :: nm4
       integer    :: nm5
       integer    :: nm6
       integer    :: nm7
       integer    :: nm8
       integer    :: qn1
       integer    :: qn2
       integer    :: qn3
       integer    :: qn4
       integer    :: qn5
       integer    :: qn6
       integer    :: qn7
       integer    :: qn8
       integer    :: ndif
       integer    :: dif(3)
       integer    :: psi1(nvdf)
       integer    :: psi2(nvdf)
 !     -----------------------------------------------------------------
 
       do i=1,nconf-1
          do j=i+1,nconf
             cnf=i+j*(j-1)/2
             nm1=configs(1,i)
             qn1=configs(2,i)
             nm2=configs(3,i)
             qn2=configs(4,i)
             nm3=configs(5,i)
             qn3=configs(6,i)
             nm4=configs(7,i)
             qn4=configs(8,i)
 
             nm5=configs(1,j)
             qn5=configs(2,j)
             nm6=configs(3,j)
             qn6=configs(4,j)
             nm7=configs(5,j)
             qn7=configs(6,j)
             nm8=configs(7,j)
             qn8=configs(8,j)
 
             psi1=0
             psi2=0
 
             if (nm1 > 0) psi1(nm1)=qn1
             if (nm2 > 0) psi1(nm2)=qn2
             if (nm3 > 0) psi1(nm3)=qn3
             if (nm4 > 0) psi1(nm4)=qn4
 
             if (nm5 > 0) psi2(nm5)=qn5
             if (nm6 > 0) psi2(nm6)=qn6
             if (nm7 > 0) psi2(nm7)=qn7
             if (nm8 > 0) psi2(nm8)=qn8
 
             dif=0
             ndif=0
             do a=1,nvdf
                if ( psi1(a) /= psi2(a) ) then
                   ndif=ndif+1
                   if (ndif >= 4) EXIT 
                   dif(ndif)=a
                end if
             end do
 
             select case (ndif)
                case (1)
                   call calcHterm1b(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,Hci,cnf,psi1,psi2,dif)
                case (2)
                   call calcHterm2b(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
                case (3)
                   call calcHterm3(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
                case default
                   Hci(cnf)=0.0d0
             end select
 
             if (Hci(cnf) < 1d-6) Hci(cnf)=0d0
 
          end do
       end do
       end subroutine
 
 
 
 
 
 !#######################################################################
 !     CONFIGURATION SELECTION VCI version 2 SECTION
 !#######################################################################
 
 
       subroutine csVCI2(ref,qumvia_nmc,ethresh,resthresh,selcut1,selcut2,&
                    &Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,Emod,&
                    &hii,tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nconf,ngaus,nvdf,qmaxx1,qmaxx2,qmaxx3,qmaxx4,bdim,Eref)
 
 !     -----------------------------------------------------------------
 !     SINGLE REFERENCE VIRTUAL CONFIGURATION INTERACTION
 !     This subroutine is the main program for virtual configuration
 !     interaction (VCI). 
 !     First a matrix representation of the hamiltonian in the VSCF 
 !     virtual wavefunctions (configurations) is built. Then it attempts
 !     to diagonalize it.
 !     To save memory, the hamiltonian is stored in packed format for
 !     symmetric matrices
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)  :: qumvia_nmc ! Number of coupled modes.
       integer,intent(in)  :: ref(8) ! Reference configuration.
       integer,intent(in)  :: ngaus ! Dimension of DGBS
       integer,intent(in)  :: nconf ! Dimension of CI basis set (VSCF configurations)
       integer,intent(in)  :: nvdf  ! Number of vib. deg. of freedom.
       integer,intent(in)  :: qmaxx1 ! Max quantum number allowed for sigles.
       integer,intent(in)  :: qmaxx2 ! Max quantum number allowed for doubles.
       integer,intent(in)  :: qmaxx3 ! Max quantum number allowed for triples.
       integer,intent(in)  :: qmaxx4 ! Max quantum number allowed for quadruples
       integer,intent(in)  :: bdim  ! Dimension of operator matrices (qmaxx1+1 by default)
       real*8,intent(in)   :: ethresh    ! Energy threshold for CS (cm-1)
       real*8,intent(in)   :: resthresh  ! Resonance threshold
       real*8,intent(in)   :: Emod(ngaus,nvdf) ! VSCF modal energies.
       real*8,intent(in)   :: Po(ngaus,ngaus,nvdf) ! VSCF coeficients.
       real*8,intent(in)   :: Hcore(ngaus,ngaus,nvdf) ! Core Hamiltonian in DGB.
       real*8,intent(in)   :: Scho(ngaus,ngaus,nvdf) ! Core Hamiltonian in DGB.
       real*8,intent(in)   :: Q1(ngaus,ngaus,nvdf) ! Q^1 operator in DGB
       real*8,intent(in)   :: Q2(ngaus,ngaus,nvdf) ! Q^2 operator in DGB
       real*8,intent(in)   :: Q3(ngaus,ngaus,nvdf) ! Q^3 operator in DGB
       real*8,intent(in)   :: GDmtrx(ngaus,ngaus,nvdf) ! Effective 2mc potential
       real*8,intent(in)   :: GTmtrx(ngaus,ngaus,nvdf) ! Effective 3mc potential
       real*8,intent(in)   :: hii(nvdf)
       real*8,intent(in)   :: tiij(nvdf,nvdf)      !! 2-mode
       real*8,intent(in)   :: tjji(nvdf,nvdf)      !! Coupling
       real*8,intent(in)   :: uiiij(nvdf,nvdf)     !! potential 
       real*8,intent(in)   :: ujjji(nvdf,nvdf)     !! paramters
       real*8,intent(in)   :: uiijj(nvdf,nvdf)     !!
       real*8,intent(in)   :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)   :: selcut1
       real*8,intent(in)   :: selcut2
       real*8,intent(out)  :: Eref
 
       integer :: csdepth ! Number of CS itterations.
       integer :: stopcnf ! Number of selected confs in prev round.
       integer :: i,j
       integer :: ns
       integer :: nmods
       integer :: nsc     ! number of selected confs
       integer :: ntrsh
       integer :: nref
       integer :: ev
       integer :: cnf
       integer :: err
       integer :: ml(1)
       real*8  :: Qm1(bdim,bdim,nvdf)
       real*8  :: Qm2(bdim,bdim,nvdf)
       real*8  :: Qm3(bdim,bdim,nvdf)
       real*8  :: Hmc(bdim,bdim,nvdf)
       real*8  :: GDm(bdim,bdim,nvdf)
       real*8  :: GTm(bdim,bdim,nvdf)
       real*8  :: tol
       real*8  :: cut
 
 !     Alocatable variables.
       real*8,dimension(:),allocatable   :: Hci
       real*8,dimension(:),allocatable   :: diag
       real*8,dimension(:),allocatable   :: sdiag
       real*8,dimension(:),allocatable   :: ECI
       real*8,dimension(:),allocatable   :: Evhf
       real*8,dimension(:),allocatable   :: vec
       real*8,dimension(:,:),allocatable :: Cci
 
       integer,dimension(:),allocatable    :: sort
       integer,dimension(:),allocatable    :: gsel
       integer,dimension(:,:), allocatable :: configs
       integer,dimension(:,:), allocatable :: sconfigs
 
 !     Workspace for diagonalization. 
       real*8,  dimension(:), allocatable :: WORK
       integer, dimension(:), allocatable :: IWORK,IFAIL
       integer :: INFO
 
       real*8,external  :: dnrm2
       real*8,external  :: dlamch
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
 !     -----------------------------------------------------------------
       nmods=qmaxx1+1 ! Dimension of CI operator matrices.
 
       write(77,'(A)') '--------------'
       write(77,'(A)') 'STARTING csVCI'
       write(77,'(A)') '--------------'
 
 !     CHANGE OF BASIS FROM GAUSSIAN TO VSCF MODALS FOR Q AND Hcore OPERATORS.
 !     -----------------------------------------------------------------
       call build_CIop3(Po,Q1,Q2,Q3,Hcore,GDmtrx,GTmtrx,Scho,ngaus,nvdf,nmods,&
                     &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
       write(77,'(A)') 'FINNISHED VCI OPERATORS   '
 
 
 !     CONFIGURATION SELECTION
 !     -----------------------------------------------------------------
       allocate(configs(8,nconf),diag(nconf),stat=err)
       if (err /= 0) STOP ('ALLOCATION ERROR: Configuration selection')
 
       nsc=1
       csdepth=2
       configs=0
       do i=1,csdepth
          if (i==1) cut=selcut1
          if (i==2) cut=5*selcut1
          call confsel(qmaxx1,qmaxx2,qmaxx3,qmaxx4,nvdf,nconf,configs,ngaus,&
               & nmods,qumvia_nmc,ref,nsc,cut,diag,ethresh,Emod,&
               & Qm1,Qm2,Qm3,Hmc,GDm,GTm,hii,tiij,tjji,uiiij,ujjji,uiijj,&
               & tijk,uiijk,uijjk,uijkk)
       end do
       write(77,'(A,I15)') 'Configs after RECURSIVE CS =',nsc
 
 
 
 !     REDUCING SIZE OF ARRAYS
 !     -----------------------------------------------------------------
       allocate(sconfigs(8,nsc),sdiag(nsc),stat=err)
       if (err /= 0) STOP ('ALLOCATION ERROR: sconfigs')
 
       sconfigs = configs(:,1:nsc)
       sdiag = diag(1:nsc)
 
       if (allocated(configs)) deallocate(configs,diag,stat=err)
       if (err /= 0) STOP ('DEALLOCATION ERROR: configs,diag')
       allocate(configs(8,nsc),diag(nsc),stat=err)
       if (err /= 0) STOP ('ALLOCATION ERROR: configs,diag')
 
       configs = sconfigs
       diag = sdiag
 
       if (allocated(sconfigs)) deallocate(sconfigs,sdiag,stat=err)
       if (err/=0) STOP('ALLOCATION ERROR: deallocating sconfigs/sdiag')
       
 !      write(77,'(A)') 'UNSORTED CONFIGS'
 !      do i=1,nsc
 !         write(77,'(8I3,D15.6)') configs(:,i),diag(i)
 !      end do
 
 !     SORTING BY DIAGONAL ENERGY
 !     -----------------------------------------------------------------
       write(77,'(A)') 'SORTING CONFIGURATIONS BY DIAGONAL ENERGY'
       allocate(sort(nsc))
       sort=0
       call quick_sort(diag,sort,nsc)
       configs = configs(:,sort)
       do i=1,nsc
          if (sort(i)==1) then
             nref=i
             exit
          end if
       end do
 !      write(77,'(A)') 'SORT'
 !      write(77,'(99999I4)') sort
 
 !      write(77,'(A)') 'SORTED CONFIGS'
 !      do i=1,nsc
 !         write(77,'(8I3,D15.6)') configs(:,i),diag(i)
 !      end do
 
 !     BUILDING DIAGONAL HAMILTONIAN
 !     -----------------------------------------------------------------
       write(77,'(A)') 'BUILDING DIAGONAL HAMILTONIAN'
       allocate (Hci(nsc+nsc*(nsc-1)/2))
       Hci=0d0
       do ns=1,nsc
          Hci(ns*(ns+1)/2)=diag(ns)
       end do
       allocate (Evhf(nsc))
       Evhf=diag
       deallocate(diag,stat=err)
       if (err/=0) STOP('DEALLOCATION ERROR: Building diag hamiltonian')
 
 
 !     BUILDING OFF-DIAGONAL HAMILTONIAN 
 !     -----------------------------------------------------------------
       write(77,'(A)') 'BUILDING OFF-DIAGONAL HAMILTONIAN'
       call csvci_mtrx2(qumvia_nmc,nvdf,nmods,nsc,configs,Hci,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 
 
 !     DIAGONALIZING HAMILTONIAN
 !     -----------------------------------------------------------------
       write(77,'(A)') 'DIAGONALIZING VCI HAMILTONIAN'
       tol=2*dlamch('S')
       allocate(Eci(nsc),Cci(nsc,nsc))
       allocate ( WORK(8*nsc), IWORK(5*nsc) ,IFAIL(nsc))
       call dspevx('V','A','U',nsc,Hci,0.0d0,0.05d0,1,10*nvdf,&
              &tol,10*nvdf,Eci,Cci,nsc,WORK,IWORK,IFAIL,INFO)
 !     call dspevx('V','I','U',nconf,Hci,0.0d0,0.0d0,1,10*nvdf,&
 !            &0.0d0,10*nvdf,Eci,Cci,nconf,WORK,IWORK,IFAIL,INFO)
 
       allocate (vec(nsc))
       if (INFO == 0) then
          do ev=1,nsc
 !         do ev=1,10*nvdf
             ml=maxloc(abs(Cci(:,ev)))
             if (ml(1)==nref) Eref=Eci(ev)
             write(77,'(A)') '---------------------------------------------'
             write(77,'(A,I6,A,F18.8,F11.2)')'EVEC  ',ev,'Eci=',Eci(ev),(Eci(ev)-Eci(1))*h2cm
             vec = Abs(Cci(:,ev))
             call quick_sort(vec,sort,nsc)
             write(77,'(A,8I3)') 'MAIN CONFIGURATION:',configs(:,sort(nsc))
             write(77,'(A)') 'COEFICIENTS'
             do cnf=1,nsc
                if (abs(Cci(sort(nsc-cnf+1),ev))<=0.05) EXIT
                write(77,'(8I3,F12.5)') configs(:,sort(nsc-cnf+1)), Cci(sort(nsc-cnf+1),ev)
             end do
          end do
 !         write(77,'(A)') 'No   CONFIGURATIONS         VSCF ENERGY VSCF TRANS.   VCI ENERGY  VCI TRANS  '
 !         do i=1,nsc
 !            write(77,'(I4,A,6I3,F22.8,F11.2,F14.8,F11.2)') i,' ',sconfigs(:,i),&
 !                 &Evhf(i),(Evhf(i)-Evhf(1))*h2cm,  &
 !                 &Eci(i),(Eci(i)-Eci(1))*h2cm
 !         end do
       else 
          print*,'ERROR DURING CI HAMILTONIAN DIAGONALIZATION'
          print*,'INFO=',INFO
          print*,'IFAIL=',IFAIL
       end if
       deallocate ( WORK,IWORK,IFAIL,Eci,Cci,Evhf,Hci )
 
       end subroutine
 
 
 
 
       subroutine confsel(qmaxx1,qmaxx2,qmaxx3,qmaxx4,nvdf,nconf,configs,ngaus,&
                         &nmods,qumvia_nmc,ref,stopcnf,pcecut,diag,ethresh,Emod,     &
                         &Qm1,Qm2,Qm3,Hmc,GDm,GTm,hii,tiij,tjji,uiiij,ujjji,uiijj,   &
                         &tijk,uiijk,uijjk,uijkk)
 !     -----------------------------------------------------------------
 !     THIS SUBROUTINE GENERATES ALL CI SINGLE AND DOUBLE CONFIGURATIONS
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qmaxx1               ! Max excitation singles
       integer,intent(in)   :: qmaxx2               ! Max excitation doubles
       integer,intent(in)   :: qmaxx3               ! Max excitation triples
       integer,intent(in)   :: qmaxx4               ! Max excitation quadruples
       integer,intent(in)   :: nvdf                ! No of vibrational dof.
       integer,intent(in)   :: nconf               ! CI basis set size.
       integer,intent(in)   :: nmods
       integer,intent(in)   :: ngaus
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: ref(8)
       integer,intent(inout):: stopcnf
       integer,intent(inout):: configs(8,nconf)
 
       real*8,intent(in)    :: pcecut
       real*8,intent(inout) :: diag(nconf)
 
       real*8,intent(in)    :: ethresh
       real*8,intent(in)    :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
 
       real*8,intent(in)    :: hii(nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
 
       integer    :: i,j
       integer    :: nc
       integer    :: nm
       integer    :: qn
       integer    :: nm1
       integer    :: nm2
       integer    :: nm3
       integer    :: nm4
       integer    :: qn1
       integer    :: qn2
       integer    :: qn3
       integer    :: qn4
       integer    :: conf(8)
       real*8     :: zpe
       real*8     :: ecut
       real*8     :: hoe
       real*8     :: hof(nvdf)
       real*8     :: Ediag
       real*8     :: Erefho
       real*8     :: Edho
 
       real*8, parameter :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 
 !     Compute Harmoni Oscillator frequencies in cm-1
       hof=sqrt(abs(hii))*h2cm  
 
 !     Compute Harmonic Oscillator zero point energy
       zpe=0d0
       do i=1,nvdf
          zpe=zpe+hof(i)/2d0
       end do
       write(77,'(A,D15.5)') 'ZERO POINT ENERGY = ',zpe
 
 !     -----------------------------------------------------------------
 !     VSCF REFERENCE STATE
 !     -----------------------------------------------------------------
       nc=0
       if (stopcnf == 1) then
          conf=ref
          Ediag=0d0
          Erefho=0d0
          call calcHOE(conf,nvdf,zpe,hof,Erefho)
          call diagCIme(conf,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,&
                   &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                   &nvdf,ngaus,nmods,nconf,Ediag)
          configs(:,1)=ref
          diag(1)=Ediag
          ecut=Erefho+ethresh
       else
          conf=ref
          Erefho=0d0
          call calcHOE(conf,nvdf,zpe,hof,Erefho)
          ecut=Erefho+ethresh
       end if
       write(77,'(A,D15.6)') 'ENERGY CUTOFF FOR CS IS ',ecut
 
 !     GROUND STATE
 !     -----------------------------------------------------------------
       conf=0
       if ( .NOT. equal(ref,conf) ) then
          call calcHOE(conf,nvdf,zpe,hof,Edho)
          call tester(nc,nconf,conf,stopcnf,pcecut,hof,zpe,configs,diag,ethresh,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,Edho,&
               &qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,nvdf,nmods,ngaus)
       end if
 
 !     FUNDAMENTALS
 !     -----------------------------------------------------------------
       do nm=1,nvdf
          conf = (/nm,1,0,0,0,0,0,0/)
          if ( equal(ref,conf) ) CYCLE
          call calcHOE(conf,nvdf,zpe,hof,Edho)
          call tester(nc,nconf,conf,stopcnf,pcecut,hof,zpe,configs,diag,ethresh,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,Edho,&
               &qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,nvdf,nmods,ngaus)
       end do
 
 !     OVERTONES
 !     -----------------------------------------------------------------
       do nm=1,nvdf
          do qn=2,qmaxx1
             conf = (/nm,qn,0,0,0,0,0,0/)
             if ( equal(ref,conf) ) CYCLE
             call calcHOE(conf,nvdf,zpe,hof,Edho)
 !            write(77,'(A)') 'OVERTONES HO ENERGIES'
 !            write(77,'(8I4,D15.6)') conf,Edho
             if ( Edho > ecut ) CYCLE
             call tester(nc,nconf,conf,stopcnf,pcecut,hof,zpe,configs,diag,ethresh,&
                  &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,Edho,&
                  &qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,nvdf,nmods,ngaus)
          end do
       end do
 
 
 !     DOUBLES
 !     -----------------------------------------------------------------
       do nm1=1,nvdf-1
          do nm2=nm1+1,nvdf
             do qn1=1,qmaxx2
                do qn2=1,qmaxx2
                   conf=(/nm1,qn1,nm2,qn2,0,0,0,0/)
                   if ( equal(ref,conf) ) CYCLE
                   call calcHOE(conf,nvdf,zpe,hof,Edho)
                   if ( Edho > ecut ) CYCLE
                   call tester(nc,nconf,conf,stopcnf,pcecut,hof,zpe,configs,diag,ethresh,&
                        &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,Edho,&
                        &qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,nvdf,nmods,ngaus)
                end do
             end do
          end do
       end do
 
 !     TRIPLES
 !     -----------------------------------------------------------------
       do nm1=1,nvdf-2
          do nm2=nm1+1,nvdf-1
             do nm3=nm2+1,nvdf
                do qn1=1,qmaxx3
                   do qn2=1,qmaxx3
                      do qn3=1,qmaxx3
                         conf=(/nm1,qn1,nm2,qn2,nm3,qn3,0,0/)
                         if ( equal(ref,conf) ) CYCLE
                         call calcHOE(conf,nvdf,zpe,hof,Edho)
                         if ( Edho > ecut ) CYCLE
                         call tester(nc,nconf,conf,stopcnf,pcecut,hof,zpe,configs,diag,ethresh,&
                              &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,Edho,&
                              &qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,nvdf,nmods,ngaus)
                      end do
                   end do
                end do
             end do
          end do
       end do
 
       if (nvdf<4) then
          stopcnf=stopcnf+nc
          RETURN
       end if
          
 
 !     QUADRUPLES
 !     -----------------------------------------------------------------
       do nm1=1,nvdf-3
          do nm2=nm1+1,nvdf-2
             do nm3=nm2+1,nvdf-1
                do nm4=nm3+1,nvdf
                   do qn1=1,qmaxx4
                      do qn2=1,qmaxx4
                         do qn3=1,qmaxx4
                            do qn4=1,qmaxx4
                               conf=(/nm1,qn1,nm2,qn2,nm3,qn3,nm4,qn4/)
                               if ( equal(ref,conf) ) CYCLE
                               call calcHOE(conf,nvdf,zpe,hof,Edho)
                               if ( Edho > ecut ) CYCLE
                               call tester(nc,nconf,conf,stopcnf,pcecut,hof,zpe,configs,diag,ethresh,&
                                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,Edho,&
                                    &qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,nvdf,nmods,ngaus)
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
 
 
 !     TOTAL NUMBER OF SELECTED CONFIGURATIONS
 !     -----------------------------------------------------------------
       stopcnf=stopcnf+nc
 
 !     WRITE UPDATE
 !     -----------------------------------------------------------------
       write(77,'(I15,A)') nc,' SELECTED CONFIGURATIONS IN THIS ROUND OF RSC'
 
       end subroutine
 
 
 
 
       subroutine tester(nc,nconf,conf,stopcnf,pcecut,hof,zpe,configs,diag,ethresh,&
                        &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,Edho,&
                        &qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,nvdf,nmods,ngaus)
 !     -----------------------------------------------------------------
 !     DECIDES IF A CONFIGURATION conf IS ACCEPTED OR REJECTED
 !     If accepted, adds it to the selected configurations list "configs'
 !     Its energy is also added to the diagonal energies list "diag".
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: ngaus,nvdf,nmods,nconf
       integer,intent(in)   :: stopcnf
       integer,intent(in)   :: conf(8)
       integer,intent(inout):: nc
       integer,intent(inout):: configs(8,nconf)
 
       real*8,intent(in)    :: zpe
       real*8,intent(in)    :: pcecut
       real*8,intent(in)    :: ethresh
       real*8,intent(in)    :: Edho
       real*8,intent(in)    :: hof(nvdf)
       real*8,intent(inout) :: diag(nconf)
 
       real*8,intent(in)    :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
 
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
 
       integer              :: i
       integer              :: rcnf(8)
       real*8               :: ecut
       real*8               :: Erho
       real*8               :: Ediag
       real*8               :: PCE
       real*8               :: MEL
       real*8               :: dE
       real*8, parameter    :: h2cm = 219474.63d0 ! cm-1/Ha
 
 !     -----------------------------------------------------------------
       do i=1,stopcnf+nc
          rcnf=configs(:,i)
          if ( equal(rcnf,conf) ) RETURN
       end do
       do i=1,stopcnf
 !         if ( equal(configs(:,i),conf) ) CYCLE
          PCE=0d0
          MEL=0d0
          rcnf=configs(:,i)
          call calcHOE(rcnf,nvdf,zpe,hof,Erho)
          if (abs(Erho-Edho) > 5000d0) CYCLE
          call calcPCE(rcnf,conf,Erho,Edho,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                      &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                      &nvdf,nmods,nconf,configs,PCE,MEL)
          if (abs(PCE) > pcecut/1000d0) then
             call diagCIme(conf,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,&
                   &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                   &nvdf,ngaus,nmods,nconf,Ediag)
             dE=(diag(i)-Ediag)*h2cm
             PCE=MEL**2/dE
             if (abs(PCE) > pcecut) then
                nc=nc+1
                configs(:,stopcnf+nc)=conf
                diag(stopcnf+nc)=Ediag
                EXIT
             end if
          end if
       end do
 
       end subroutine
 
 
 
 
       subroutine calcHOE(conf,nvdf,zpe,hof,hoe)
 !     -----------------------------------------------------------------
 !     COMPUTES HARMONIC OSCILLATOR ENERGY OF A CONFIGURATION conf
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: nvdf
       integer,intent(in)   :: conf(8)
       real*8,intent(in)    :: zpe
       real*8,intent(in)    :: hof(nvdf)
       real*8,intent(out)   :: hoe
 
       integer     :: i,j
       integer     :: nm,qn
 !     -----------------------------------------------------------------
       hoe=zpe
       do i=1,4
          nm=conf((i-1)*2+1)
          qn=conf(2*i)
          hoe = hoe + qn*hof(nm)
       end do
 
       end subroutine
 
 
 
       subroutine diagCIme(conf,qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,Emod,&
                      &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                      &nvdf,ngaus,nmods,nconf,Ediag)
 !     -----------------------------------------------------------------
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc ! # of coupled modes.
       integer,intent(in)   :: nvdf   ! # of vibrational deg of freedom.
       integer,intent(in)   :: nmods  ! # of selected VSCF virtual states.
       integer,intent(in)   :: nconf  ! Total number of configurations
       integer,intent(in)   :: ngaus  ! Dimension of gaussian basis set
       integer,intent(in)   :: conf(8)  ! Reference configuration
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf) !Q operator in modals basis
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf) !Q^2 operator in modals basis
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf) !Q^3 operator in modals basis
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf) !Core hamiltonian in modals basis.
       real*8,intent(in)    :: Emod(ngaus,nvdf)       ! Modal energies.
       real*8,intent(in)    :: tiij(nvdf,nvdf)   ! Coupling potential parameters
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(out)   :: Ediag
 
       integer :: i,a,b,c
       integer :: nm1,nm2,nm3,nm4
       integer :: qn1,qn2,qn3,qn4
       integer :: cnf
       integer :: sref
       integer :: n1,n2,n3 ! Various indices.
       integer :: psi(nvdf) ! VSCF wavefunction. Array of vib. quantum numbers.
       integer :: sort(nconf)
       real*8 :: Vc2       ! Auxiliary variable for Hci calculation.
       real*8 :: Vc3       ! Auxiliary variable for Hci calculation.
       real*8 :: ehf       ! Auxiliary variable for Hci calculation.
       real*8 :: Hcore     ! Auxiliary variable for Hci calculation.
       real*8 :: veff      ! Auxiliary variable for Hci calculation.
       real*8 :: veff2     ! Auxiliary variable for Hci calculation.
       real*8 :: veff3     ! Auxiliary variable for Hci calculation.
       real*8 :: absethresh
       real*8 :: Ehf0      ! DEBUG
       real*8 :: corr      ! DEBUG
 
       real*8, parameter :: h2cm = 219474.63d0 ! cm-1/Ha
 !     -----------------------------------------------------------------
       
       nm1=conf(1)
       qn1=conf(2)
       nm2=conf(3)
       qn2=conf(4)
       nm3=conf(5)
       qn3=conf(6)
       nm4=conf(7)
       qn4=conf(8)
 
       psi=0
       if (nm1 > 0) psi(nm1)=qn1
       if (nm2 > 0) psi(nm2)=qn2
       if (nm3 > 0) psi(nm3)=qn3
       if (nm4 > 0) psi(nm4)=qn4
 
 !     VSCF energy.
       ehf=0d0
       do a=1,nvdf
          ehf = ehf + Emod(psi(a)+1,a)
       end do
 
 !     Core hamiltonian
       Hcore=0d0
       do a=1,nvdf
          Hcore= Hcore + Hmc(psi(a)+1,psi(a)+1,a)
       end do
 
 !     Average on effective potential.
       veff2=0d0
       do a=1,nvdf
          veff2 = veff2 + GDm(psi(a)+1,psi(a)+1,a)
       end do
       veff3=0d0
       do a=1,nvdf
          veff3 = veff3 + GTm(psi(a)+1,psi(a)+1,a)
       end do
       veff = veff2 + veff3
 
 !     2-MODE COUPLING POTENTIAL 
       Vc2=0.0d0
       do a=1,nvdf-1
          do b=a+1,nvdf
             n1=psi(a)+1
             n2=psi(b)+1
             Vc2=Vc2+&
             & tiij(a,b)*Qm2(n1,n1,a)*Qm1(n2,n2,b)/2.0d0+&
             & tjji(a,b)*Qm1(n1,n1,a)*Qm2(n2,n2,b)/2.0d0+&
             & uiiij(a,b)*Qm3(n1,n1,a)*Qm1(n2,n2,b)/6.0d0+&
             & ujjji(a,b)*Qm1(n1,n1,a)*Qm3(n2,n2,b)/6.0d0+&
             & uiijj(a,b)*Qm2(n1,n1,a)*Qm2(n2,n2,b)/4.0d0  
          end do
       end do
 
 !     3-MODE COUPLING POTENTIAL
       Vc3=0d0
       IF (qumvia_nmc == 3) THEN
       do a=1,nvdf-2
          do b=a+1,nvdf-1
             do c=b+1,nvdf
                n1=psi(a)+1
                n2=psi(b)+1
                n3=psi(c)+1
 
                Vc3=Vc3+&
                 & tijk(a,b,c)*Qm1(n1,n1,a)*Qm1(n2,n2,b)*Qm1(n3,n3,c)+&
                 & uiijk(a,b,c)*Qm2(n1,n1,a)*Qm1(n2,n2,b)*Qm1(n3,n3,c)/2d0+&
                 & uijjk(a,b,c)*Qm1(n1,n1,a)*Qm2(n2,n2,b)*Qm1(n3,n3,c)/2d0+&
                 & uijkk(a,b,c)*Qm1(n1,n1,a)*Qm1(n2,n2,b)*Qm2(n3,n3,c)/2d0
 
             end do
          end do
       end do
       END IF
 
 !     DANGER
 !      Ediag=Hcore+Vc2+Vc3
 !      Ediag=ehf-Vc2-2.0d0*Vc3
       Ediag=ehf-veff+Vc2+Vc3
       
    
       end subroutine
 
 
 
 
 
       subroutine calcPCE(conf1,conf2,E1,E2,qumvia_nmc,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,configs,CSPCE,MEL)
 !     -----------------------------------------------------------------
 !     Computes Off-diagonal Hamiltonian matrix elements 
 !                      < K | H | L >
 !     where | K > and | L > may differ in 1 (case 1), 2 (case 2) and 
 !     more than 2 (case 3) modals. This subroutine runs over all off-
 !     diagonal matrix elements, determines to which case each belong
 !     and then computes the integral using subroutines calcHterm1 
 !     for case 1 matrix elements, calcHterm2 for case 2 ME or
 !     sets Hci(cnf)=0.0 for case 3.
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: conf1(8)
       integer,intent(in)   :: conf2(8)
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf
       integer,intent(in)   :: configs(8,nconf)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: E1,E2
       real*8,intent(out)   :: CSPCE
       real*8,intent(out)   :: MEL
 
       integer    :: a
       integer    :: b
       integer    :: i
       integer    :: j
       integer    :: r
       integer    :: cnf
       integer    :: nm1
       integer    :: nm2
       integer    :: nm3
       integer    :: nm4
       integer    :: nm5
       integer    :: nm6
       integer    :: nm7
       integer    :: nm8
       integer    :: qn1
       integer    :: qn2
       integer    :: qn3
       integer    :: qn4
       integer    :: qn5
       integer    :: qn6
       integer    :: qn7
       integer    :: qn8
       integer    :: ndif
       integer    :: dif(3)
       integer    :: psi1(nvdf)
       integer    :: psi2(nvdf)
 !     -----------------------------------------------------------------
 
          nm1=conf1(1)
          qn1=conf1(2)
          nm2=conf1(3)
          qn2=conf1(4)
          nm3=conf1(5)
          qn3=conf1(6)
          nm4=conf1(7)
          qn4=conf1(8)
 
          nm5=conf2(1)
          qn5=conf2(2)
          nm6=conf2(3)
          qn6=conf2(4)
          nm7=conf2(5)
          qn7=conf2(6)
          nm8=conf2(7)
          qn8=conf2(8)
 
          psi1=0
          psi2=0
 
          if (nm1 > 0) psi1(nm1)=qn1
          if (nm2 > 0) psi1(nm2)=qn2
          if (nm3 > 0) psi1(nm3)=qn3
          if (nm4 > 0) psi1(nm4)=qn4
 
          if (nm5 > 0) psi2(nm5)=qn5
          if (nm6 > 0) psi2(nm6)=qn6
          if (nm7 > 0) psi2(nm7)=qn7
          if (nm8 > 0) psi2(nm8)=qn8
 
          dif=0
          ndif=0
          do a=1,nvdf
             if ( psi1(a) /= psi2(a) ) then
                ndif=ndif+1
                if (ndif >= 4) EXIT 
                dif(ndif)=a
             end if
          end do
 
          select case (ndif)
             case (1)
               call CSPCEcase1(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
               &nvdf,nmods,nconf,CSPCE,MEL,cnf,psi1,psi2,E1,E2,dif)
             case (2)
               call CSPCEcase2(qumvia_nmc,Qm1,Qm2,Qm3,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
               &nvdf,nmods,nconf,CSPCE,MEL,cnf,psi1,psi2,dif,E1,E2)
             case (3)
               call CSPCEcase3(qumvia_nmc,Qm1,Qm2,Qm3,&
               &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
               &nvdf,nmods,nconf,CSPCE,MEL,cnf,psi1,psi2,dif,E1,E2)
             case default
               CSPCE=0D0
               MEL=0D0
          end select
 
 !        DEBUG--------------------------------
 !         write(77,'(A,F16.8)') 'CSPE=',CSPCE
 !        -------------------------------------
 
       end subroutine
 
 
       subroutine CSPCEcase1(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,CSPCE,MEL,cnf,psi1,psi2,E1,E2,dif)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY ONE MODAL.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf,cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: E1,E2
       real*8,intent(out)   :: CSPCE
       real*8,intent(out)   :: MEL
 
       integer :: i,j,k
       integer :: a,b
       integer :: n1,n2
       integer :: ni1,ni2
       integer :: nj1,nj2
       integer :: nk1,nk2
       real*8 :: Vc2
       real*8 :: Vc3
       real*8 :: veff      ! Auxiliary variable for Hci calculation.
       real*8 :: dE
       real*8, PARAMETER :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 
 
 !     Computing VSCF effective potential.
       veff=0d0
       n1=psi1(dif(1))+1
       n2=psi2(dif(1))+1
       veff=GDm(n1,n2,dif(1))+GTm(n1,n2,dif(1))
 
 !     Computing Coupling Potential
       Vc2 = 0.0d0
       do a=1,nvdf
          if (a == dif(1)) CYCLE
 
          if (a<dif(1)) then
             i=a
             j=dif(1)
             ni1=psi1(a)+1
             ni2=ni1
             nj1=psi1(dif(1))+1
             nj2=psi2(dif(1))+1
          else
             i=dif(1)
             j=a
             ni1=psi1(dif(1))+1
             ni2=psi2(dif(1))+1
             nj1=psi1(a)+1
             nj2=nj1
          end if 
              
          Vc2 = Vc2 + &
               & tiij(i,j)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)/2.0d0 + &
               & tjji(i,j)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)/2.0d0 + &
               & uiiij(i,j)*Qm3(ni1,ni2,i)*Qm1(nj1,nj2,j)/6.0d0+ &
               & ujjji(i,j)*Qm1(ni1,ni2,i)*Qm3(nj1,nj2,j)/6.0d0+ &
               & uiijj(i,j)*Qm2(ni1,ni2,i)*Qm2(nj1,nj2,j)/4.0d0
       end do
 
       Vc3 = 0d0
       IF (qumvia_nmc == 3) THEN
       do a=1,nvdf-1
          do b=a+1,nvdf
 
             if (a == dif(1)) CYCLE
             if (b == dif(1)) CYCLE
 
             if (dif(1) < a) then
                i=dif(1)
                j=a
                k=b
 
                ni1=psi1(dif(1))+1
                ni2=psi2(dif(1))+1
                nj1=psi1(a)+1
                nj2=psi2(a)+1
                nk1=psi1(b)+1
                nk2=psi2(b)+1
             else if (a<dif(1) .AND. dif(1)<b) then
                i=a
                j=dif(1)
                k=b
 
                ni1=psi1(a)+1
                ni2=psi2(a)+1
                nj1=psi1(dif(1))+1
                nj2=psi2(dif(1))+1
                nk1=psi1(b)+1
                nk2=psi2(b)+1
             else if (b<dif(1)) then
                i=a
                j=b
                k=dif(1)
 
                ni1=psi1(a)+1
                ni2=psi2(a)+1
                nj1=psi1(b)+1
                nj2=psi2(b)+1
                nk1=psi1(dif(1))+1
                nk2=psi2(dif(1))+1
             end if
 
             Vc3=Vc3+&
              & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)+ &
              & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0+ &
              & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0+ &
              & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
          end do
       end do 
       END IF
 
 
 !     COMPUTING HAMILTONIAN MATRIX ELEMENT
 !      dE=(E1-E2)*h2cm
       dE=(E1-E2)
 !     DANGER
       MEL=(Vc2 + Vc3 - veff)*h2cm
 !      MEL=(Vc2 + Vc3)*h2cm
       CSPCE=MEL**2/dE
 
       end subroutine
 
       subroutine CSPCEcase2(qumvia_nmc,Qm1,Qm2,Qm3,&
                   &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                   &nvdf,nmods,nconf,CSPCE,MEL,cnf,psi1,psi2,dif,E1,E2)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY TWO MODALS.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf
       integer,intent(in)   :: nmods
       integer,intent(in)   :: nconf
       integer,intent(in)   :: cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: E1,E2
       real*8,intent(out)   :: CSPCE
       real*8,intent(out)   :: MEL
 
       integer :: i,j,k
       integer :: a,b,c
       integer :: ni1,ni2,nj1,nj2,nk1,nk2
       real*8  :: Vc2,Vc3
       real*8  :: dE
       real*8, PARAMETER :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 
 
 
       if (dif(1) < dif(2)) then
          i=dif(1)
          j=dif(2)
          ni1=psi1(dif(1))+1
          ni2=psi2(dif(1))+1
          nj1=psi1(dif(2))+1
          nj2=psi2(dif(2))+1
       else
          i=dif(2)
          j=dif(1)
          ni1=psi1(dif(2))+1
          ni2=psi2(dif(2))+1
          nj1=psi1(dif(1))+1
          nj2=psi2(dif(1))+1
       end if 
           
       Vc2 = 0d0
       Vc2 =  &
          & tiij(i,j)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)/2.0d0+&
          & tjji(i,j)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)/2.0d0+&
          & uiiij(i,j)*Qm3(ni1,ni2,i)*Qm1(nj1,nj2,j)/6.0d0+&
          & ujjji(i,j)*Qm1(ni1,ni2,i)*Qm3(nj1,nj2,j)/6.0d0+&
          & uiijj(i,j)*Qm2(ni1,ni2,i)*Qm2(nj1,nj2,j)/4.0d0
 
       Vc3 = 0d0
       IF (qumvia_nmc == 3) THEN
       do a=1,nvdf
 
          if (a == dif(1)) CYCLE
          if (a == dif(2)) CYCLE
 
          if (a < dif(1) .AND. dif(1) < dif(2)) then
             i=a
             j=dif(1)
             k=dif(2)
 
             ni1=psi1(i)+1
             ni2=psi2(i)+1
             nj1=psi1(j)+1
             nj2=psi2(j)+1
             nk1=psi1(k)+1
             nk2=psi2(k)+1
          else if (dif(1) < a .AND. a < dif(2)) then
             i=dif(1)
             j=a
             k=dif(2)
 
             ni1=psi1(i)+1
             ni2=psi2(i)+1
             nj1=psi1(j)+1
             nj2=psi2(j)+1
             nk1=psi1(k)+1
             nk2=psi2(k)+1
          else if ( dif(2) < a ) then
             i=dif(1)
             j=dif(2)
             k=a
 
             ni1=psi1(i)+1
             ni2=psi2(i)+1
             nj1=psi1(j)+1
             nj2=psi2(j)+1
             nk1=psi1(k)+1
             nk2=psi2(k)+1
          end if
 
          Vc3=Vc3+ &
           & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k) + &
           & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
       end do
       END IF 
    
 !     VCI MATRIX ELEMENT
 !      dE=(E1-E2)*h2cm
       dE=(E1-E2)
       MEL = (Vc2 + Vc3)*h2cm
       CSPCE = MEL**2/dE
 
       end subroutine
 
 
 
 
       subroutine CSPCEcase3(qumvia_nmc,Qm1,Qm2,Qm3,&
                   &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                   &nvdf,nmods,nconf,CSPCE,MEL,cnf,psi1,psi2,dif,E1,E2)
 !     -----------------------------------------------------------------
 !     Given 2 configurations |K> and |L>, this subroutine computes
 !     the corresponding Hamiltonian matrix element suposing that
 !     |K > and |L> differ in ONLY TWO MODALS.
 !                          < K | H | L >
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf
       integer,intent(in)   :: nmods
       integer,intent(in)   :: nconf
       integer,intent(in)   :: cnf
       integer,intent(in)   :: psi1(nvdf)
       integer,intent(in)   :: psi2(nvdf)
       integer,intent(in)   :: dif(3)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: E1,E2
       real*8,intent(out)   :: CSPCE
       real*8,intent(out)   :: MEL
 
       integer :: i,j,k
       integer :: ni1,ni2,nj1,nj2,nk1,nk2
       real*8  :: dE
       real*8, PARAMETER :: h2cm=219475.64d0  ! Convert Hartree to cm-1
 !     -----------------------------------------------------------------
 
 
       IF (qumvia_nmc == 3) THEN
 !        COMPUTING 3-MODE COUPLING TERMS
          if (dif(1) > dif(2)) STOP ('dif(1) > dif(2) in calcHterm3')
          if (dif(2) > dif(3)) STOP ('dif(2) > dif(3) in calcHterm3')
          if (dif(1) > dif(3)) STOP ('dif(1) > dif(3) in calcHterm3')
    
          i=dif(1)
          j=dif(2)
          k=dif(3)
    
          ni1=psi1(i)+1
          ni2=psi2(i)+1
          nj1=psi1(j)+1
          nj2=psi2(j)+1
          nk1=psi1(k)+1
          nk2=psi2(k)+1
    
          MEL= &
           & tijk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k) + &
           & uiijk(i,j,k)*Qm2(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijjk(i,j,k)*Qm1(ni1,ni2,i)*Qm2(nj1,nj2,j)*Qm1(nk1,nk2,k)/2d0 + &
           & uijkk(i,j,k)*Qm1(ni1,ni2,i)*Qm1(nj1,nj2,j)*Qm2(nk1,nk2,k)/2d0
 
          dE=(E1-E2)
          MEL=MEL*h2cm
          CSPCE=MEL**2/dE
       ELSE
          CSPCE = 0d0
       END IF
 
       end subroutine
 
       subroutine csvci_mtrx2(qumvia_nmc,nvdf,nmods,nconf,configs,Hci,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &Qm1,Qm2,Qm3,Hmc,GDm,GTm)
 !     -----------------------------------------------------------------
 !     Computes Off-diagonal Hamiltonian matrix elements 
 !                      < K | H | L >
 !     where | K > and | L > may differ in 1 (case 1), 2 (case 2) and 
 !     more than 2 (case 3) modals. This subroutine runs over all off-
 !     diagonal matrix elements, determines to which case each belong
 !     and then computes the integral using subroutines calcHterm1 
 !     for case 1 matrix elements, calcHterm2 for case 2 ME or
 !     sets Hci(cnf)=0.0 for case 3.
 !     -----------------------------------------------------------------
       implicit none
 
       integer,intent(in)   :: qumvia_nmc
       integer,intent(in)   :: nvdf,nmods,nconf
       integer,intent(in)   :: configs(8,nconf)
       real*8,intent(in)    :: Qm1(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm2(nmods,nmods,nvdf)
       real*8,intent(in)    :: Qm3(nmods,nmods,nvdf)
       real*8,intent(in)    :: Hmc(nmods,nmods,nvdf)
       real*8,intent(in)    :: GDm(nmods,nmods,nvdf)
       real*8,intent(in)    :: GTm(nmods,nmods,nvdf)
       real*8,intent(in)    :: tiij(nvdf,nvdf)
       real*8,intent(in)    :: tjji(nvdf,nvdf)
       real*8,intent(in)    :: uiiij(nvdf,nvdf)
       real*8,intent(in)    :: ujjji(nvdf,nvdf)
       real*8,intent(in)    :: uiijj(nvdf,nvdf)
       real*8,intent(in)    :: tijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uiijk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijjk(nvdf,nvdf,nvdf)
       real*8,intent(in)    :: uijkk(nvdf,nvdf,nvdf)
       real*8,intent(inout) :: Hci(nconf+nconf*(nconf-1)/2)
 
       integer    :: a
       integer    :: b
       integer    :: i
       integer    :: j
       integer    :: cnf
       integer    :: nm1
       integer    :: nm2
       integer    :: nm3
       integer    :: nm4
       integer    :: nm5
       integer    :: nm6
       integer    :: nm7
       integer    :: nm8
       integer    :: qn1
       integer    :: qn2
       integer    :: qn3
       integer    :: qn4
       integer    :: qn5
       integer    :: qn6
       integer    :: qn7
       integer    :: qn8
       integer    :: ndif
       integer    :: dif(3)
       integer    :: psi1(nvdf)
       integer    :: psi2(nvdf)
 !     -----------------------------------------------------------------
 
       do i=1,nconf-1
          do j=i+1,nconf
             cnf=i+j*(j-1)/2
             nm1=configs(1,i)
             qn1=configs(2,i)
             nm2=configs(3,i)
             qn2=configs(4,i)
             nm3=configs(5,i)
             qn3=configs(6,i)
             nm4=configs(7,i)
             qn4=configs(8,i)
 
             nm5=configs(1,j)
             qn5=configs(2,j)
             nm6=configs(3,j)
             qn6=configs(4,j)
             nm7=configs(5,j)
             qn7=configs(6,j)
             nm8=configs(7,j)
             qn8=configs(8,j)
 
             psi1=0
             psi2=0
 
             if (nm1 > 0) psi1(nm1)=qn1
             if (nm2 > 0) psi1(nm2)=qn2
             if (nm3 > 0) psi1(nm3)=qn3
             if (nm4 > 0) psi1(nm4)=qn4
 
             if (nm5 > 0) psi2(nm5)=qn5
             if (nm6 > 0) psi2(nm6)=qn6
             if (nm7 > 0) psi2(nm7)=qn7
             if (nm8 > 0) psi2(nm8)=qn8
 
             dif=0
             ndif=0
             do a=1,nvdf
                if ( psi1(a) /= psi2(a) ) then
                   ndif=ndif+1
                   if (ndif >= 4) EXIT 
                   dif(ndif)=a
                end if
             end do
 
             select case (ndif)
                case (1)
                   call calcHterm1b(qumvia_nmc,Qm1,Qm2,Qm3,Hmc,GDm,GTm,&
                    &tiij,tjji,uiiij,ujjji,uiijj,tijk,uiijk,uijjk,uijkk,&
                    &nvdf,nmods,nconf,Hci,cnf,psi1,psi2,dif)
                case (2)
                   call calcHterm2b(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
                case (3)
                   call calcHterm3(qumvia_nmc,Qm1,Qm2,Qm3,tiij,tjji,uiiij,ujjji,uiijj,&
                   &tijk,uiijk,uijjk,uijkk,nvdf,nmods,nconf,&
                   &Hci,cnf,psi1,psi2,dif)
                case default
                   Hci(cnf)=0.0d0
             end select
 
             if (Abs(Hci(cnf)) < 1d-9) Hci(cnf)=0d0
 
          end do
       end do
       end subroutine
 
 
 
 
 
 
 
 
 
 
 
 
 
 !#######################################################################
 !     QUICKSORT SECTION 
 !#######################################################################
 
 RECURSIVE SUBROUTINE quick_sort(list, order, as)
 
 ! Quick sort routine from:
 ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
 ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
 ! Modified by Alan Miller to include an associated integer array which gives
 ! the positions of the elements in the original order.
 ! Modified by Diego Alonso de Armino for the use of double precision arrays.
 
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: as
 REAL*8, DIMENSION (as), INTENT(IN OUT)  :: list
 INTEGER, DIMENSION (as), INTENT(OUT)  :: order
 
 ! Local variable
 INTEGER :: i
 
 DO i = 1, SIZE(list)
   order(i) = i
 END DO
 
 CALL quick_sort_1(1, SIZE(list))
 
 CONTAINS
 
 RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
 
 INTEGER, INTENT(IN) :: left_end, right_end
 
 !     Local variables
 INTEGER             :: i, j, itemp
 REAL*8              :: reference, temp
 INTEGER, PARAMETER  :: max_simple_sort_size = 6
 
 IF (right_end < left_end + max_simple_sort_size) THEN
   ! Use interchange sort for small lists
   CALL interchange_sort(left_end, right_end)
 
 ELSE
   ! Use partition ("quick") sort
   reference = list((left_end + right_end)/2)
   i = left_end - 1; j = right_end + 1
 
   DO
     ! Scan list from left end until element >= reference is found
     DO
       i = i + 1
       IF (list(i) >= reference) EXIT
     END DO
     ! Scan list from right end until element <= reference is found
     DO
       j = j - 1
       IF (list(j) <= reference) EXIT
     END DO
 
 
     IF (i < j) THEN
       ! Swap two out-of-order elements
       temp = list(i); list(i) = list(j); list(j) = temp
       itemp = order(i); order(i) = order(j); order(j) = itemp
     ELSE IF (i == j) THEN
       i = i + 1
       EXIT
     ELSE
       EXIT
     END IF
   END DO
 
   IF (left_end < j) CALL quick_sort_1(left_end, j)
   IF (i < right_end) CALL quick_sort_1(i, right_end)
 END IF
 
 END SUBROUTINE quick_sort_1
 
 
 SUBROUTINE interchange_sort(left_end, right_end)
 
 INTEGER, INTENT(IN) :: left_end, right_end
 
 !     Local variables
 INTEGER             :: i, j, itemp
 REAL*8              :: temp
 
 DO i = left_end, right_end - 1
   DO j = i+1, right_end
     IF (list(i) > list(j)) THEN
       temp = list(i); list(i) = list(j); list(j) = temp
       itemp = order(i); order(i) = order(j); order(j) = itemp
     END IF
   END DO
 END DO
 
 END SUBROUTINE interchange_sort
 
 END SUBROUTINE quick_sort
 
 end module qvamod_common
 
  
