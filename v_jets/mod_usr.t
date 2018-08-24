module mod_usr
  use mod_mhd
  implicit none
  double precision :: rhoj, eta, vj
  double precision :: amp, B_strength, jet_time, alpha_val, tilt_pc
contains

  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid

    call mhd_activate()
    call params_read(par_files)

  end subroutine usr_init

!> Read parameters from a file
  subroutine params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /my_parameters/ amp,B_strength,jet_time,alpha_val,tilt_pc  
  do n = 1, size(files)
    open(unitpar, file=trim(files(n)), status="old")
    read(unitpar, my_parameters, end=113)
113 close(unitpar)
      end do
  end subroutine params_read


  subroutine initglobaldata_usr
    
    mhd_gamma=1.4d0
    rhoj=mhd_gamma
    eta=100.d0
    vj=800.d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    {^IFONED call mpistop("This is a multi-D MHD problem") }

    where(dabs(x(ix^S,1))<0.05d0.and.x(ix^S,2)<0.00d0)
       w(ix^S,rho_)=rhoj
       w(ix^S,mom(1))=0.0d0
       w(ix^S,mom(2))=rhoj*vj
       w(ix^S,e_)=one/(mhd_gamma-one)+0.5d0*rhoj*vj**2.0d0+0.5d0*B_strength**2.0d0
       w(ix^S,mag(1))=0.0d0
       w(ix^S,mag(2))=B_strength
    else where
       w(ix^S,rho_) = rhoj/eta
       w(ix^S,e_) = one/(mhd_gamma-one)+0.5d0*B_strength**2.0d0
       w(ix^S,mom(1)) = 0.0d0
       w(ix^S,mom(2)) = 0.0d0
       w(ix^S,mag(1))=0.0d0
       w(ix^S,mag(2))=B_strength
    end where

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer :: ixI^L, ix2

    ixImin^DD=ixBmin^DD;
    ixImax^DD=ixBmin^D-1+nghostcells^D%ixImax^DD=ixBmax^DD;
    ! Outflow:
    do ix2=ixImin2,ixImax2
       w(ixImin1:ixImax1,ix2,rho_) = w(ixImin1:ixImax1,ixImax2+1,rho_) 
       w(ixImin1:ixImax1,ix2,e_)   = w(ixImin1:ixImax1,ixImax2+1,e_) 
       w(ixImin1:ixImax1,ix2,mom(1))  = w(ixImin1:ixImax1,ixImax2+1,mom(1))
       w(ixImin1:ixImax1,ix2,mom(2))  = w(ixImin1:ixImax1,ixImax2+1,mom(2))
       w(ixImin1:ixImax1,ix2,mag(1))  = w(ixImin1:ixImax1,ixImax2+1,mag(1))
       w(ixImin1:ixImax1,ix2,mag(2))  = w(ixImin1:ixImax1,ixImax2+1,mag(2))
    end do
    where(dabs(x(ixI^S,1))<0.05d0)
       w(ixI^S,rho_)=rhoj
       w(ixI^S,mom(1))=zero
       w(ixI^S,mom(2))=rhoj*vj
       w(ixI^S,e_)=one/(mhd_gamma-one)+0.5d0*rhoj*vj**2.0d0+0.5*B_strength**2.0d0
       w(ixI^S,mag(1))=0.0d0
       w(ixI^S,mag(2))=B_strength
    else where
       ! Reflective:
       !   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,rho_) 
       !   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,e_) 
       !   w(ixI^S,mom(1)) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,mom(1))
       !   w(ixI^S,mom(2)) =-w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,mom(2))
    end where

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! e.g. refine for negative first coordinate x < 0 as
    !
    if (minval(dabs(x(ix^S,1))) < 0.1.and.minval(dabs(x(ix^S,2))) < 0.1) refine=1

  end subroutine specialrefine_grid

end module mod_usr
