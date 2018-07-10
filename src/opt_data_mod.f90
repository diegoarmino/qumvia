module opt_data_mod
   implicit none

   integer                  :: nat
   integer                  :: stepno
   integer                  :: max_opt_steps
   double precision         :: lambda
   double precision         :: escf,eold
   double precision         :: gamma

   double precision, allocatable, dimension(:)     :: dXold_vec   ! Gradient in step k (column vector format)
   double precision, allocatable, dimension(:)     :: dXnew_vec   ! Gradient in step k+1 (column vector format)
   double precision, allocatable, dimension(:)     :: Xold_vec    ! Position in step k (column vector format)
   double precision, allocatable, dimension(:)     :: Xnew_vec    ! Position in step k+1 (column vector format)
   double precision, allocatable, dimension(:,:)   :: Hold        ! Hessian matrix estimation for step k
   double precision, allocatable, dimension(:,:)   :: Xnew        ! Position for step k+1 in lio format.
   double precision, allocatable, dimension(:,:)   :: Xold        ! Position for step k in lio format.
   double precision, allocatable, dimension(:,:)   :: dXnew       ! Gradient for step k+1 in Lio format.
   double precision, allocatable, dimension(:,:)   :: dXold       ! Gradient for step k in Lio format.
   double precision, allocatable, dimension(:,:)   :: vold        ! Position for step k+1 in lio format.

end module
