subroutine conjugate_gradient()
   use opt_data_mod, only :  Xnew, dXnew, Xold, dXold, lambda, max_opt_steps, &
                             vold, nat, dXnew_vec, dXold_vec, stepno, gamma   
   implicit none


!  -----------------------------------------------------------------------------
   integer                            :: i,j,k
   double precision                   :: dXdif_vec(3*nat)
   double precision                   :: vnew(3,nat)
   double precision                   :: num, den
!  -----------------------------------------------------------------------------
   double precision,external          :: DDOT
!  -----------------------------------------------------------------------------

   lambda=2.5D-2
   gamma = 0d0
!  Change gradient into vector format.
   do j=1,nat
   do k=1,3
     dXnew_vec(3*(j-1)+k)   =  dXnew(k,j)
   end do
   end do

!  Computing direction vector
   if (stepno > 1) then
!     Polak-Riviere algorithm
      dXdif_vec = dXnew_vec - dXold_vec  ! Difference in gi+1 - gi for Polack-Riviere method.
      num = ddot(3*nat,dXdif_vec,1,dXnew_vec,1)
      den = ddot(3*nat,dXold_vec,1,dXold_vec,1)
      gamma = num/SIGN(MAX(DABS(den),1D-18),den) ! Preventing division by zero.
      if (ABS(den) < 1d-18) write(*,*) "Warning: norm of G smaller than 1d-18!"
      vnew = -dXnew + gamma*vold
   else
      vnew = -lambda*dXnew
   end if

!  Moving coordinates.
!  --------------------------------------------------------------------------
   Xnew = Xnew + vnew
!  --------------------------------------------------------------------------

!  Storing values of direction and gradient for next step.
   vold  = vnew
   dXold_vec = dXnew_vec

end subroutine conjugate_gradient
