subroutine steepest_descent(nat,Xopt,at_numbers)
   use opt_data_mod, only :  X0, dX, Xold, lambda, max_opt_steps
   implicit none

!  -----------------------------------------------------------------------------
   integer,            intent(in)     :: nat
   integer,            intent(in)     :: at_numbers(nat)
   double precision,   intent(out)    :: Xopt(3,nat)
!  -----------------------------------------------------------------------------
   logical                            :: converged
   integer                            :: nclatoms=0
   integer                            :: i,j,k
   integer                            :: igmax
   double precision                   :: escf,eold
   double precision                   :: rms, gnorm
   double precision                   :: svec(3,nat)
   double precision                   :: dipxyz(3)
   double precision                   :: dxyzcl(3,nat)
   double precision                   :: clcoords(4,nat)
   double precision                   :: gradvec(3*nat)
!  -----------------------------------------------------------------------------
   double precision,parameter         :: maxg_thresh=4d-2
   double precision,parameter         :: rms_thresh=2d-2
!  -----------------------------------------------------------------------------
   double precision,external          :: DNRM2
   double precision,external          :: IDAMAX
!  -----------------------------------------------------------------------------

   converged=.FALSE.
   lambda=2.5D-2

!  Check if energy rises and update lambda
   if (stepno > 1) then
      if (i > 1 .AND. escf < eold) lambda = lambda*1.2d0
      if (i > 1 .AND. escf > eold) then
         lambda = lambda*0.5d0
         X0 = Xold
         CYCLE
      end if
   end if

!  Change gradient into vector format.
   do j=1,nat
      do k=1,3
         gradvec(3*(j-1)+k)=dX(k,j)
      end do
   end do

!  Check convergence
   rms   = DNRM2(3*nat,gradvec,1)
   rms   = rms/SQRT(3*DBLE(nat))
   igmax = idamax(3*nat,gradvec,1)

   if(gradvec(igmax) < maxg_thresh .AND. rms < rms_thresh) converged = .TRUE.

!  Move to new geometry
   if ( converged == .FALSE. ) then
      gnorm =   dnrm2(3*nat,gradvec,1)
      svec  = - dX/gnorm
      Xold  =   X0
      X0    =   X0+lambda*svec
      eold  =   escf
   else
      Xopt  =  X0
      write(77,'(A)') "Geometry Converged!"
      do j=1,nat
         write(77,'(I3,3D15.6)') at_numbers(j),(Xopt(k,j),k=1,3)
      end do
      exit
   end if

end subroutine steepest_descent