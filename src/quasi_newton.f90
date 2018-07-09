subroutine quasi_newton(nat,at_numbers)
   use opt_data_mod, only : dXold_vec, dXnew_vec, Xold_vec, Xnew_vec, Hold, &
                            lambda, max_opt_steps
   implicit none

!  -----------------------------------------------------------------------------
   integer,            intent(in)     :: nat
   integer,            intent(in)     :: at_numbers(nat)
   double precision,   intent(out)    :: Xopt(3,nat)
!  -----------------------------------------------------------------------------
   logical                            :: converged
   integer                            :: i,j,k
   integer                            :: igmax
   double precision                   :: rms, gnorm
   double precision                   :: vnew(3,nat)
   double precision                   :: dr(3*nat) ! Xnew-Xold vector format
   double precision                   :: dg(3*nat) ! dXnew-dXold vector format
   double precision                   :: dgdr      ! scalar product dg.dr
   double precision                   :: num1(3*nat,3*nat) ! dXnew-dXold vector format
   double precision                   :: Tmp(3*nat,3*nat) ! Temporary matrix
   double precision                   :: Hnew(3*nat,3*nat) ! dXnew-dXold vector format
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
!        Change gradient and position into vector format.
         do j=1,nat
         do k=1,3
            Xnew_vec(3*(j-1)+k) =  Xnew(k,j)
            dXnew_vec(3*(j-1)+k)   =  dXnew(k,j)
         end do
         end do

         dr = Xnew_vec  - Xold_vec
         dg = dXnew_vec - dXold_vec

!        Computing scalar products of denominators
         dgdr = ddot(3*nat,dg,1,dr,1)
         dgdr = 1d0/dgdr

!        Computing external products of numerators divided by denominators.
         num1=0d0
         dger(3*nat,3*nat,dgdr,dr,1,dg,1,num1,3*nat)

!        Computing (I-num1) where I=identity matrix.
         do i=1,3*nat
            num1(i,i) = 1d0 - num1(i,i)
         end do

!        Computing first term of the formula (I-num1) Hk (I-num1)^T
!        First term
         Tmp   = 0d0
         dgemm('N','N',3*nat,3*nat,3*nat,1d0,num1,3*nat,Hold,3*nat,1d0,Tmp,3*nat)
         dgemm('N','T',3*nat,3*nat,3*nat,1d0,Tmp,3*nat,num1,3*nat,1d0,Hnew,3*nat)

!        Computing second term of the formula dr.dr^T/(dg.dr) and adding to Hnew
         dger(3*nat,3*nat,dgdr,dr,1,dr,1,Hnew,3*nat)

!        Computing geometry updater v = Hnew.dg
         vnew_vec=0d0
         dgemv('N',3*nat,3*nat,1d0,Hnew,3*nat,dg,1,1d0,vnew_vec,1)
         do j=1,nat
         do k=1,3
            vnew(k,j)=vnew_vec(3*(j-1)+k)
         end do
         end do

      else
         vnew = -lambda*dX
         Hold = 0d0
         do i=1,3*nat
            Hold(i,i) = 1d0
         end do
      end if

!     Moving coordinates.
!     --------------------------------------------------------------------------
      Xnew = Xnew + vnew
!     --------------------------------------------------------------------------

!     Storing values of position, gradient and Hessian approximation for next step.
      dXold_vec = dXnew_vec
      Xold_vec  = Xnew_vec
      Hold      = Hnew

   else
      Xopt  =  Xnew
      write(77,'(A)') "Geometry Converged!"
      do j=1,nat
         write(77,'(I3,3D18.8)') at_numbers(j),(Xopt(k,j),k=1,3)
      end do
      exit
   end if

end subroutine conjugate_gradient
