subroutine steepest_descent()
   use opt_data_mod, only :  Xnew, dXnew, Xold, dXold, lambda, max_opt_steps,  &
                             stepno, escf, eold, dXnew_vec, nat
   implicit none

!  -----------------------------------------------------------------------------
   integer                            :: i,j,k
   double precision                   :: svec(3,nat)
   double precision                   :: gnorm
!  -----------------------------------------------------------------------------
   double precision,external          :: DNRM2
!  -----------------------------------------------------------------------------

!  Check if energy rises and update lambda
   if (stepno > 1) then
      if (i > 1 .AND. escf <= eold) lambda = lambda*1.2d0
      if (i > 1 .AND. escf > eold) then
         lambda = lambda*0.5d0
         Xnew   = Xold
         dXnew  = dXold
         escf   = eold
      end if
   end if

!  Change gradient into vector format.
   do j=1,nat
   do k=1,3
      dXnew_vec(3*(j-1)+k)=dXnew(k,j)
   end do
   end do

!  Move to new geometry
   gnorm =   dnrm2(3*nat,dXnew_vec,1)
   svec  = - dXnew/gnorm
   Xold  =   Xnew
   Xnew  =   Xnew+lambda*svec
   eold  =   escf

end subroutine steepest_descent
