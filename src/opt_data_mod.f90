module opt_data_mod
   private
   public ::
   implicit none
   contains

   integer                  :: nat
   integer                  :: max_opt_steps
   double precision         :: X0(3,nqmatoms) ! Initial geometry
   double precision         :: dX(3,nqmatoms) ! Gradient
   double precision         :: Xk(3,nqmatoms) ! Position in step k
   double precision         :: lambda



end module
