module opt_data_mod
   implicit none

   integer                  :: nat
   integer                  :: stepno
   integer                  :: max_opt_steps
   double precision         :: lambda
   double precision         :: gamma

   double precision, allocatable, dimension(:,:)     :: X0 ! Initial geometry
   double precision, allocatable, dimension(:,:)     :: dX ! Gradient
   double precision, allocatable, dimension(:,:)     :: Xold ! Position in step k

end module
