module mod_usr
  use mod_mhd
  implicit none
  double precision :: amp, tilt_pc
  double precision :: rhoj, eta, vj, x0, delta_x, driv_trw, PoI, smoothness, j_h, BB2
  double precision :: jet_speed, jet_time, tilt, tilt_deg, jet_width, plasma_beta
  double precision :: xc, yc
  logical :: driv_gaussian, driv_tanh, orig_setup, cir
contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2D")

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid
    usr_aux_output    => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call mhd_activate()
    call params_read(par_files)

  end subroutine usr_init

!> Read parameters from a file
  subroutine params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /my_parameters/ jet_speed,plasma_beta,jet_time,jet_width,tilt_deg,driv_gaussian, driv_tanh, orig_setup, cir
  do n = 1, size(files)
    open(unitpar, file=trim(files(n)), status="old")
    read(unitpar, my_parameters, end=113)
113 close(unitpar)
      end do
  end subroutine params_read


  subroutine initglobaldata_usr

    mhd_gamma=1.4d0
    rhoj=mhd_gamma
    eta=3.d0
    vj=jet_speed ! jet Mach speed
    delta_x = jet_width*2
    x0 = 0.0d0 ! centre pt of guassian
    driv_trw = jet_width ! transition width for driver
    PoI = 0.5d0*(-jet_width+(-jet_width+driv_trw)) ! pnt of inflection
    smoothness = 0.1d0*driv_trw ! smoothness of tanh driver
    j_h = jet_width 
    tilt = (dpi/180.0d0)*tilt_deg ! convert deg to radians
    xc = 0.0d0 
    yc = 0.0d0 
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: magtotal(ix^S,1:ndir), B_sq(ix^S)
    integer :: idir
    double precision :: trw_wall, mid_pt_wall, smoothness_wall, rho_bg, rho_base
    double precision:: rcloud(ixG^T)

    {^IFONED call mpistop("This is a multi-D MHD problem") }

    BB2 = dsqrt(two/Plasma_beta)
    w(ix^S,mag(1))=zero
    w(ix^S,mag(2))=BB2
    do idir=1,ndir
      if(B0field) then
        magtotal(ix^S,idir)=w(ix^S,mag(idir))+block%B0(ix^S,idir,0)
      else
        magtotal(ix^S,idir)=w(ix^S,mag(idir))
      endif
    end do
    ! B^2
    B_sq(ix^S)=sum((magtotal(ix^S,:))**2,dim=ndim+1)
    if(orig_setup) then
      where(dabs(x(ix^S,1))<0.05d0.and.x(ix^S,2)<0.00d0)
         w(ix^S,rho_)=rhoj
         w(ix^S,mom(1))=0.0d0
         w(ix^S,mom(2))=rhoj*vj
         w(ix^S,e_)=one/(mhd_gamma-one)+0.5d0*rhoj*vj**2.0d0+0.5d0*B_sq(ix^S)
      else where
         w(ix^S,rho_) = rhoj/eta
         w(ix^S,e_) = one/(mhd_gamma-one)+0.5d0*B_sq(ix^S)
         w(ix^S,mom(1)) = 0.0d0
         w(ix^S,mom(2)) = 0.0d0
      end where
    endif

    if(driv_gaussian) then
      where(dabs(x(ix^S,1))<jet_width.and.x(ix^S,2)<0.50d0)
         w(ix^S,rho_)=rhoj
         w(ix^S,mom(1))=0.0d0
         w(ix^S,mom(2))=w(ix^S,rho_)*vj*dexp(-((x(ix^S,1)-x0)/delta_x)**2.0d0)
         w(ix^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ix^S,rho_)*((w(ix^S,mom(1))/w(ix^S,rho_))**2.0d0+(w(ix^S,mom(2))/w(ix^S,rho_))**2.0d0)+0.5d0*B_sq(ix^S)
      else where
         w(ix^S,rho_) = rhoj/eta
         w(ix^S,e_) = one/(mhd_gamma-one)+0.5d0*B_sq(ix^S)
         w(ix^S,mom(1)) = 0.0d0
         w(ix^S,mom(2)) = 0.0d0
      end where
    endif

    if(driv_tanh) then
      ! outside of jet
      w(ix^S,rho_) = rhoj/eta
      w(ix^S,e_) = one/(mhd_gamma-one)
      w(ix^S,mom(1)) = 0.0d0
      w(ix^S,mom(2)) = 0.0d0

      ! LHS of jet
      where(x(ix^S,1)>-jet_width .and. x(ix^S,1) <= 0.0d0 .and. x(ix^S,2)<j_h)
        w(ix^S,rho_)=rhoj/eta+(rhoj-rhoj/eta)*0.5d0*(1.0d0+dtanh((x(ix^S,1)-PoI)/smoothness))
        w(ix^S,mom(1))=dsin(tilt)*w(ix^S,rho_)*(vj*0.5d0)*(1.0d0+dtanh((x(ix^S,1)-PoI)/smoothness))
        w(ix^S,mom(2))=dcos(tilt)*w(ix^S,rho_)*(vj*0.5d0)*(1.0d0+dtanh((x(ix^S,1)-PoI)/smoothness))
        w(ix^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ix^S,rho_)*((w(ix^S,mom(1))/w(ix^S,rho_))**2.0d0 +(w(ix^S,mom(2))/w(ix^S,rho_))**2.0d0)+0.5d0*B_sq(ix^S)
      end where
      ! RHS of jet
      where(x(ix^S,1)<jet_width .and. x(ix^S,1) >= 0.0d0 .and. x(ix^S,2)<j_h)
        w(ix^S,rho_)=rhoj/eta+(rhoj-rhoj/eta)*0.5d0*(1.0d0-dtanh((x(ix^S,1)+PoI)/smoothness))
        w(ix^S,mom(1))=dsin(tilt)*w(ix^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((x(ix^S,1)+PoI)/smoothness))
        w(ix^S,mom(2))=dcos(tilt)*w(ix^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((x(ix^S,1)+PoI)/smoothness))
        w(ix^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ix^S,rho_)*((w(ix^S,mom(1))/w(ix^S,rho_))**2.0d0+(w(ix^S,mom(2))/w(ix^S,rho_))**2.0d0)+0.5d0*B_sq(ix^S)
      end where
    endif

    if(cir) then
      rcloud(ix^S)=(x(ix^S,1)-xc)**2+(x(ix^S,2)-yc)**2
      w(ix^S,rho_)=rhoj/eta+(rhoj-rhoj/eta)*0.5d0*(1.0d0-dtanh((dsqrt(rcloud(ix^S))+PoI)/smoothness))
      w(ix^S,mom(1))=dsin(tilt)*w(ix^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((dsqrt(rcloud(ix^S))+PoI)/smoothness))
      w(ix^S,mom(2))=dcos(tilt)*w(ix^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((dsqrt(rcloud(ix^S))+PoI)/smoothness))
      w(ix^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ix^S,rho_)*((w(ix^S,mom(1))/w(ix^S,rho_))**2.0d0+(w(ix^S,mom(2))/w(ix^S,rho_))**2.0d0)+0.5d0*B_sq(ix^S)
     endif

    ! set wall
    trw_wall = 0.2d0 !transition width
    mid_pt_wall = 0.5d0*(xprobmax2-xprobmin2) ! pos of transtion reg
    smoothness_wall = trw_wall*0.001d0! how rapid is the change
    rho_bg = rhoj/eta
    rho_base = rho_bg*1e-2
!    where(x(ix^S,2)>(0.5d0*(xprobmax2-xprobmin2)-trw_wall))
!      w(ix^S,rho_) = rho_base+(rho_bg-rho_base)*0.5d0*(1.0d0-dtanh((x(ix^S,2)-mid_pt_wall)/smoothness_wall))
!    end where
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)

    ! special boundary types, user defined
    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer :: ixI^L, ix2
    double precision :: magtotal(ixG^S,1:ndir), B_sq(ixG^S)
    integer :: idir
    double precision:: rcloud(ixG^T)

    w(ixG^S,mag(1))=zero
    w(ixG^S,mag(2))=BB2
    do idir=1,ndir
      if(B0field) then
        magtotal(ixG^S,idir)=w(ixG^S,mag(idir))+block%B0(ixG^S,idir,0)
      else
        magtotal(ixG^S,idir)=w(ixG^S,mag(idir))
      endif
    end do
    ! B^2
    B_sq(ixG^S)=sum((magtotal(ixG^S,:))**2,dim=ndim+1)

    ixImin^DD=ixBmin^DD;
    ixImax^DD=ixBmin^D-1+nghostcells^D%ixImax^DD=ixBmax^DD;

    ! Outflow:
    if(orig_setup) then
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
         w(ixI^S,e_)=one/(mhd_gamma-one)+0.5d0*rhoj*vj**2.0d0+0.5*B_sq(ixI^S)
      else where
         ! Reflective:
         !   w(ixI^S,rho_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,rho_) 
         !   w(ixI^S,e_) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,e_) 
         !   w(ixI^S,mom(1)) = w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,mom(1))
         !   w(ixI^S,mom(2)) =-w(ixImin1:ixImax1,ixImax2+nghostcells:ixImax2+1:-1,mom(2))
      end where
    endif

    if(driv_gaussian) then
      where(dabs(x(ixI^S,1))<jet_width.and.x(ixI^S,2)<0.50d0)
         w(ixI^S,rho_)=rhoj
         w(ixI^S,mom(1))=0.0d0
         w(ixI^S,mom(2))=w(ixI^S,rho_)*vj*dexp(-((x(ixI^S,1)-x0)/delta_x)**2.0d0)
         w(ixI^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ixI^S,rho_)*((w(ixI^S,mom(1))/w(ixI^S,rho_))**2.0d0+(w(ixI^S,mom(2))/w(ixI^S,rho_))**2.0d0)+0.5d0*B_sq(ixI^S)
      else where
         w(ixI^S,rho_) = rhoj/eta
         w(ixI^S,e_) = one/(mhd_gamma-one)+0.5d0*B_sq(ixI^S)
         w(ixI^S,mom(1)) = 0.0d0
         w(ixI^S,mom(2)) = 0.0d0
      end where
    endif

    if(driv_tanh) then
      ! outside of jet
      w(ixI^S,rho_) = rhoj/eta
      w(ixI^S,e_) = one/(mhd_gamma-one)
      w(ixI^S,mom(1)) = 0.0d0
      w(ixI^S,mom(2)) = 0.0d0

      ! LHS of jet
      where(x(ixI^S,1)>-jet_width .and. x(ixI^S,1) <= 0.0d0 .and. x(ixI^S,2)<j_h)
        w(ixI^S,rho_)=rhoj/eta+(rhoj-rhoj/eta)*0.5d0*(1.0d0+dtanh((x(ixI^S,1)-PoI)/smoothness))
        w(ixI^S,mom(1))=dsin(tilt)*w(ixI^S,rho_)*(vj*0.5d0)*(1.0d0+dtanh((x(ixI^S,1)-PoI)/smoothness))
        w(ixI^S,mom(2))=dcos(tilt)*w(ixI^S,rho_)*(vj*0.5d0)*(1.0d0+dtanh((x(ixI^S,1)-PoI)/smoothness))
        w(ixI^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ixI^S,rho_)*((w(ixI^S,mom(1))/w(ixI^S,rho_))**2.0d0 +(w(ixI^S,mom(2))/w(ixI^S,rho_))**2.0d0)+0.5d0*B_sq(ixI^S)
      end where
      ! RHS of jet
      where(x(ixI^S,1)<jet_width .and. x(ixI^S,1) >= 0.0d0 .and. x(ixI^S,2)<j_h)
        w(ixI^S,rho_)=rhoj/eta+(rhoj-rhoj/eta)*0.5d0*(1.0d0-dtanh((x(ixI^S,1)+PoI)/smoothness))
        w(ixI^S,mom(1))=dsin(tilt)*w(ixI^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((x(ixI^S,1)+PoI)/smoothness))
        w(ixI^S,mom(2))=dcos(tilt)*w(ixI^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((x(ixI^S,1)+PoI)/smoothness))
        w(ixI^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ixI^S,rho_)*((w(ixI^S,mom(1))/w(ixI^S,rho_))**2.0d0+(w(ixI^S,mom(2))/w(ixI^S,rho_))**2.0d0)+0.5d0*B_sq(ixI^S)
      end where
    endif

    if(cir) then
      rcloud(ixI^S)=(x(ixI^S,1)-xc)**2+(x(ixI^S,2)-yc)**2
      w(ixI^S,rho_)=rhoj/eta+(rhoj-rhoj/eta)*0.5d0*(1.0d0-dtanh((dsqrt(rcloud(ixI^S))+PoI)/smoothness))
      w(ixI^S,mom(1))=dsin(tilt)*w(ixI^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((dsqrt(rcloud(ixI^S))+PoI)/smoothness))
      w(ixI^S,mom(2))=dcos(tilt)*w(ixI^S,rho_)*(vj*0.5d0)*(1.0d0-dtanh((dsqrt(rcloud(ixI^S))+PoI)/smoothness))
      w(ixI^S,e_)=one/(mhd_gamma-one)+0.5d0*w(ixI^S,rho_)*((w(ixI^S,mom(1))/w(ixI^S,rho_))**2.0d0+(w(ixI^S,mom(2))/w(ixI^S,rho_))**2.0d0)+0.5d0*B_sq(ixI^S)
     endif

  end subroutine specialbound_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    integer :: idirmin,idir,ix^D

    double precision:: gradrho(ixG^T),rho(ixG^T),drho(ixG^T)
    double precision:: gradp(ixG^T),dp(ixG^T), ggrid(ixI^S), p(ixG^T)
    double precision:: kk,kk0,grhomax,kk1
    integer         :: idims
    logical, save   :: firstrun=.true.

    w(ixO^S,nw+1) = w(ixO^S,e_)
    w(ixO^S,nw+2) = w(ixO^S,mom(1))
    w(ixO^S,nw+3) = w(ixO^S,mom(2))
    ! output temperature
    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+4)=pth(ixO^S)/w(ixO^S,rho_)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
      else
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
      endif
    end do
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+5)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    ! output divB1
    call get_normalized_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+6)=divb(ixO^S)

    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+7)=pth(ixO^S)*two/B2(ixO^S)
!    if(mype==0) print*, 'pressure', pth(ixO^S) 
!    if (mype==0) print*, 'btotal', B2(ixO^S) 

     rho(ixI^S)=w(ixI^S,rho_)
     gradrho(ixO^S)=zero
     do idims=1,ndim
       select case(typegrad)
       case("central")
         call gradient(rho,ixI^L,ixO^L,idims,drho)
       case("limited")
         call gradientS(rho,ixI^L,ixO^L,idims,drho)
       end select
       gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
     enddo

     gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
     kk=5.0d0
     kk0=0.01d0
     kk1=1.0d0
!     grhomax=1.d-8

     grhomax=100.0d0!MAXVAL(gradrho(ixO^S))

  ! putting the schlierplot of density in nwauxio=1
     w(ixO^S,nw+8)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
   
     w(ixO^S,nw+9)=dsqrt(mhd_gamma*pth(ixO^S)/w(ixO^S,rho_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='e m1 m2 Te Alfv divB beta schrho cs'

  end subroutine specialvarnames_output

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
