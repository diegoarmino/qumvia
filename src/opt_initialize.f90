subroutine opt_initialize(qva_nml)
   use opt_data_mod, only: Xnew,dXnew,vold,Xold_vec,Xnew_vec,dXnew_vec, &
                           dXold_vec, Hold
   use qvamod_common, only: qva_nml_type

   implicit none
   type(qva_nml_type), intent(in)  :: qva_nml  ! QUMVIA namelist parameters


   allocate( Xnew(3,nat) )
   allocate( dXnew(3,nat) )

   select case (qva_nml%opt_type)
     case ('SD')
        lambda=2.5D-2
     case ('CG')
        allocate( vold(3,nat)   )
     case ('QN')
        allocate( Xold_vec(3*nat)   )
        allocate( Xnew_vec(3*nat)   )
        allocate( dXnew_vec(3*nat)  )
        allocate( dXold_vec(3*nat)  )
        allocate( Hold(3*nat,3*nat) )
        Xold_vec  = 0d0
        Xnew_vec  = 0d0
        dXold_vec = 0d0
        dXnew_vec = 0d0
        Hold      = 0d0
     case default
        allocate( Xold_vec(3*nat)   )
        allocate( Xnew_vec(3*nat)   )
        allocate( dXnew_vec(3*nat)  )
        allocate( dXold_vec(3*nat)  )
        allocate( Hold(3*nat,3*nat) )
        Xold_vec  = 0d0
        Xnew_vec  = 0d0
        dXold_vec = 0d0
        dXnew_vec = 0d0
        Hold      = 0d0
   end select

end subroutine opt_initialize
