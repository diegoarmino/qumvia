subroutine conjugate_gradient(nat,Xopt,at_numbers)
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
   double precision,external          :: DDOT
   double precision,external          :: IDAMAX
!  -----------------------------------------------------------------------------

   converged=.FALSE.
   lambda=2.5D-2

!  Check convergence
   rms   = DNRM2(3*nat,dX_vec,1)
   rms   = rms/SQRT(3*DBLE(nat))
   igmax = idamax(3*nat,dX_vec,1)

   if(dX_vec(igmax) < maxg_thresh .AND. rms < rms_thresh) converged = .TRUE.

!  Computing direction vector
   if ( converged == .FALSE. ) then
      if (stepno > 1) then
!        Change gradient into vector format.
         do j=1,nat
         do k=1,3
            dX_vec(3*(j-1)+k)   =  dX(k,j)
         end do
         end do

!        Polak-Riviere algorithm
         dXdif = dX_vec - dXold_vec  ! Difference in gi+1 - gi for Polack-Riviere method.
         num = ddot(3*nat,dXdif_vec,1,gradvec,1)
         num = num**2
         den = dnrm(3*nat,dXold_vec,1)
         den = den**2
         gamma = num/den
         vnew = -dX + gamma*vold
      else
         vnew = -lambda*dX
      end if

!     Moving coordinates.
!     --------------------------------------------------------------------------
      X0 = X0 + vnew
!     --------------------------------------------------------------------------

!     Storing values of direction and gradient for next step.
      vold  = vnew
      dXold_vec = dX_vec
   else
      Xopt  =  X0
      write(77,'(A)') "Geometry Converged!"
      do j=1,nat
         write(77,'(I3,3D18.8)') at_numbers(j),(Xopt(k,j),k=1,3)
      end do
      exit
   end if

end subroutine conjugate_gradient