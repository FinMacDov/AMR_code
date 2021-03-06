module mod_usr
  use mod_mhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya,BB1,BB2,BB3
  double precision, allocatable :: pa(:),ra(:),ya(:)
  integer, parameter :: jmax=8000

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2.5D")

    mhd_gamma=1.66666667d0
    mhd_eta=zero ! This gives idea MHD

    unit_length        = 1.d8                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters
    integer :: iv

    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade
    kx=dpi/(xprobmax1-xprobmin1)
    ly=kx*dcos(theta)
    SRadius=6.96d10/unit_length ! Solar radius
    !SRadius=69.61d0 ! Solar radius
    
    BB1=0.0d0/unit_magneticfield
    BB2=60.0d0/unit_magneticfield
    BB3=0.0d0/unit_magneticfield

    !=> To allow output to be in physical uints
    length_convert_factor = unit_length
    time_convert_factor = unit_time
    w_convert_factor(1) = unit_density
    w_convert_factor(5) = unit_pressure

    do iv = 2,4
      w_convert_factor(iv) = unit_velocity
      w_convert_factor(iv+4) = unit_magneticfield
    enddo 

    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

  end subroutine initglobaldata_usr

  subroutine inithdstatic
  !! initialize the table in a vertical line through the global domain
    use mod_global_parameters

    integer :: j,na,nb,ibc
    double precision, allocatable :: Ta(:),gg(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa

    rpho=1.151d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=2.0!0.2d0 ! height of initial transition region
    wtra=0.02d0 ! width of initial transition region
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
    
    !> creates temperture profile
    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))
    do j=1,jmax
       ya(j)=(dble(j)-0.5d0)*dya-gzone
       if(ya(j)>htra) then
         Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
       else
         Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
       endif
       gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo
    !! solution of hydrostatic equation
    nb=int(gzone/dya)
    ra(1)=rpho
    pa(1)=rpho*Tpho
    invT=gg(1)/Ta(1) !<1/H(z)
    invT=0.d0
    !>using scale height for HS equation
    do j=2,jmax
       invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
       pa(j)=pa(1)*dexp(invT*dya)
       ra(j)=pa(j)/Ta(j)
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
        write(*,*)'dimensionless vars:'
        write(*,*)'T =', unit_time, 's'
        write(*,*)'L =', unit_length, 'cm'
        write(*,*)'V =', unit_velocity, 'cm/s'
        write(*,*)'rho =', unit_numberdensity, 'g/cm3'
        write(*,*)'rho =', unit_density, 'cm-3'
        write(*,*)'B =', unit_magneticfield, 'G'
        write(*,*)'T =', unit_temperature, 'K'
        write(*,*)'p =', unit_pressure, 'dyn cm-2'
      endif
      first=.false.
    endif
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)= ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}
    w(ixO^S,mom(:))=zero
    if(B0field) then
      w(ixO^S,mag(:))=zero
    else
    select case(iprob)
      case(1,11)
      w(ixO^S,mag(1))=BB1
      w(ixO^S,mag(2))=BB2
      w(ixO^S,mag(3))=BB3
      case(2,12)
      !origninal setup
      w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
      w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
      w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
     case default
        call mpistop('iprob to implement')
    endselect
    endif
 
   ! Set tracer
    if(mhd_n_tracer>0) then
        w(ixO^S,tracer(1))=0.0d0
    end if

    if(mhd_glm) w(ixO^S,psi_)=0.d0

    call phys_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: delydelx, j_origx, rho_j, deltax, deltay
    double precision ::jet_w, jet_h, A, switch, delta_x, delta_y, phase, switch_off, endtime
    integer :: ind^D, test
    integer :: ix^D,idir,ixInt^L
    logical, save :: first=.true.

    select case(iB)
    case(3)
      !! fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        select case(iprob)
          case(1,11)
          w(ixO^S,mag(1))=BB1
          w(ixO^S,mag(2))=BB2
          w(ixO^S,mag(3))=BB3
          case(2,12)
          !origninal setup
          w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
          w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
          w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
         case default
            call mpistop('iprob to implement')
        endselect
      endif
      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      !to drive jet
      jet_w = 3.5d7/2.0d0/unit_length !< (350 km)/unit_length
      A = 5.0d6/unit_velocity
      delta_x = (jet_w)/3.d0 !< This defines the width of guass dist.
                                   !< divided by 3 as guass dist = 0 
                                   !< after 3 sigma. 
      j_origx = (abs(xprobmax1)-abs(xprobmin1))/2.0d0
      rho_j = 1.0d-9/unit_density
      phase = 2.0d0*dpi/10.0d0
!      phase = 2.0d0*dpi/(60.0d0/unit_time)
      endtime = 240.0d0/unit_time
      switch_off = endtime/2
!     note qt has no dim
     test = 2
     if(test == 0) then
!      call mhd_to_primitive(ixI^L,ixO^L,w,x)
      do ind1=ixOmin1,ixOmax1
        do ind2=ixOmin2,ixOmax2
          if( x(ind^D,1)<=jet_w .and. x(ind^D,1)>=-jet_w) then
            w(ind^D,rho_) = rho_j 
            if(mhd_n_tracer>0) then
              w(ind^D,tracer(1))=100.0d0
            endif
              w(ind^D,mom(1))=zero
              w(ind^D,mom(2))= A
            else
              w(ind^D,mom(2))=zero
            endif
        end do
      end do
      endif

!     do ind1 = ixOmin1,ixOmax1
!     print*, w(ind1,ixOmax2,mom(2))*unit_velocity, x(ind1,ixOmax2,1), qt*unit_time
!     enddo

     if(test==1) then
      do ind1=ixOmin1,ixOmax1
        do ind2=ixOmin2,ixOmax2
          if( x(ind^D,1)<=jet_w .and. x(ind^D,1)>=-jet_w .and. x(ind^D,2)<=jet_h) then
            w(ind^D,rho_) = rho_j 
            if(mhd_n_tracer>0) then
              w(ind^D,tracer(1))=100.0d0
            endif
            if (qt >= switch_off) then
              if (qt-endtime > 0.0d0) then
              w(ind^D,mom(2))= 0.0d0
              else  
              w(ind^D,mom(2))= -A*dtanh(phase*(qt-endtime))*& 
                               dexp(-((x(ind1,ind2,1)-j_origx)/delta_x)**2)  
              endif
              else
                w(ind^D,mom(2))= A*dtanh(phase*qt)*&
                               dexp(-((x(ind1,ind2,1)-j_origx)/delta_x)**2)
            endif
          endif
        end do
      end do
      endif

     if(test==2) then
      !to drive jet
      jet_w = (xprobmax1-xprobmin1)/domain_nx1 !This makes the jet 1 cell radius 
      jet_h = (xprobmax2-xprobmin2)/domain_nx2
      A = 5.0d6/unit_velocity
      deltax = (xprobmax1-xprobmin1)/domain_nx1
      deltay = (xprobmax2-xprobmin2)/domain_nx2
      do ind1=ixOmin1,ixOmax1
        do ind2=ixOmin2,ixOmax2
          if( x(ind^D,1)<=jet_w .and. x(ind^D,1)>=-jet_w .and. x(ind^D,2)<=jet_h) then
            switch = 2.0d0*dpi*qt/10.0d0 ! will peak around 4.25 s
           print*, x(ixOmin1,ind2,1), ind1, ind2, qt*unit_time, jet_w, x(ixOmax1,ind2,1), x(ind1,ind2,1)
            w(ind^D,mom(2))= A*dtanh(switch)*&
                             dexp(-(dabs(((x(ind1,ind2,1)-x(ixOmin1,ind2,1)/deltax)))**2&
                             +(dabs((x(ind1,ind2,2)-x(ind1,ixOmin2,2))/deltay))**2))  
            w(ind^D,rho_) = 1.0d-9/unit_density 
            if(mhd_n_tracer>0) then
              w(ind^D,tracer(1))=100.0d0
            endif
          end if
        end do
      end do
      endif


      if(mhd_glm) w(ixO^S,psi_)=0.d0
      call mhd_to_conserved(ixI^L,ixO^L,w,x)

!     do ind1 = ixOmin1,ixOmax1
!     print*, (w(ind1,ixOmax2,mom(2))*unit_velocity)/(rho_j), qt*unit_time
!     enddo 
   case(4)
      ixInt^L=ixO^L;
      ixIntmin2=ixOmin2-1;ixIntmax2=ixOmin2-1;
      call phys_get_pthermal(w,x,ixI^L,ixInt^L,pth)
      ixIntmin2=ixOmin2-1;ixIntmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixInt^L,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      enddo
      !> fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      !> zero normal gradient extrapolation
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
                    (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
               +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
      enddo
      if(mhd_glm) w(ixO^S,psi_)=0.d0
      call phys_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,2)=ggrid(ixO^S)

  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
  ! calculate background heating bQ
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S)

    bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/5.d0)

  end subroutine getbQ

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif

  end subroutine special_refine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_radiative_cooling
    use mod_geometry

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S), wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    integer :: idirmin,idir,ix^D

    double precision:: gradrho(ixG^T),rho(ixG^T),drho(ixG^T)
    double precision:: kk,kk0,grhomax,kk1
    integer                            :: idims

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=(pth(ixO^S)/w(ixO^S,rho_))*unit_temperature

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
    w(ixO^S,nw+2)=unit_velocity*dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    ! output divB1
    call get_normalized_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)

    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)

     rho(ixI^S)=w(ixI^S,rho_)
     gradrho(ixO^S)=zero
     do idims=1,ndim
       select case(typegrad)
       case("central")
         call gradient(rho*unit_density,ixI^L,ixO^L,idims,drho)
       case("limited")
         call gradientS(rho*unit_density,ixI^L,ixO^L,idims,drho)
       end select
       gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
     enddo

     gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
     kk=5.0d0
     kk0=0.01d0
     kk1=1.0d0
     grhomax=1.d-8

!     grhomax=MAXVAL(gradrho(ixO^S))

  ! putting the schlierplot of density in nwauxio=1
     w(ixO^S,nw+5)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
   
     w(ixO^S,nw+6)=unit_velocity*dsqrt(mhd_gamma*pth(ixO^S)/w(ixO^S,rho_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te Alfv divB beta schrho cs'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or
  ! linear force-free background field
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)


    select case(iprob)
      case(1,11)
      wB0(ixO^S,1)=BB1
      wB0(ixO^S,2)=BB2
      wB0(ixO^S,3)=BB3
      case(2,12)
      !origninal setup
      wB0(ixO^S,1)=-B0*dcos(kx*x(ixO^S,1)-dpi/2.0d0)*dexp(-ly*x(ixO^S,2))*dcos(theta)
      wB0(ixO^S,2)=+B0*dsin(kx*x(ixO^S,1)-dpi/2.0d0)*dexp(-ly*x(ixO^S,2))
      wB0(ixO^S,3)=-B0*dcos(kx*x(ixO^S,1)-dpi/2.0d0)*dexp(-ly*x(ixO^S,2))*dsin(theta)
     case default
        call mpistop('iprob to implement')
    endselect


  end subroutine specialset_B0

end module mod_usr

