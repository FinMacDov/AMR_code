module mod_usr
  use mod_mhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya
  double precision, allocatable :: pa(:),ra(:),ya(:),Ha(:)
  integer, parameter :: jmax=8000

  double precision :: s0,s1,Bv,B_y,y_r,BB1,BB2,BB3
  logical :: driver, tanh_profile, c7_profile, driver_kuz, driver_random,&
      integrate, derivative, derivative_2
  logical :: driver_injetion, driver_injetion_1, jet_cont, jet_switch
  double precision :: randphase(10), randA(10), randP(10) 
  integer :: nxmodes

contains
  subroutine usr_init()
    call set_coordinate_system("Cartesian_2.5D")

    unit_length= 1.d9  ! cm
    unit_temperature = 1.d6  ! K
    unit_numberdensity = 1.d9  ! cm-3,cm-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    !>aaded
    usr_add_aux_names   => specialvarnames_output
    usr_var_for_errest  => Te_for_errest
 
    call mhd_activate()
    call params_read(par_files)
  end subroutine usr_init

!> Read parameters from a file
  subroutine params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /my_switches/ driver, driver_kuz, driver_random, tanh_profile,&
      c7_profile, integrate,derivative, derivative_2, driver_injetion,&
      driver_injetion_1, jet_cont, jet_switch
  do n = 1, size(files)
   open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_switches, end=111)
111    close(unitpar)
    end do
  end subroutine params_read

!> Can base AMR on Te
  subroutine Te_for_errest(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iflag,w,x,var)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,iflag
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: var(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2)

    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,pth)
    var(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
  end subroutine Te_for_errest

  subroutine initglobaldata_usr()
    integer :: iv
    integer,dimension(:),allocatable:: seed
    integer::  seed_size,ix
    real:: randphaseP1(1:10), randA1(1:10), randP1(1:10)

    heatunit=unit_pressure/unit_time !3.697693390805347E-003 erg*cm-3/s,erg*cm-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) !cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 !the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
!    kx=dpi/(xprobmax1-xprobmin1)
    kx=dpi/((xprobmax1-xprobmin1)/2.0d0)-0.5d0
    ly=kx*dcos(theta)
!    SRadius=6.96d10/unit_length ! Solar radius
    SRadius=69.61d0 ! Solar radius

    !=> in Guass
    BB1=0.0d0/unit_magneticfield
    BB2=60.0d0/unit_magneticfield
    BB3=0.0d0/unit_magneticfield

    !=>Konkol et al. 2012 and Kuzma et al. 2017
    !=>s1 = s and s0 = a ("a" is not sound speed)
    B_y = 8.0d0/unit_magneticfield
    y_r = 1.0d9/unit_length
    Bv  = 6.0d0/unit_magneticfield
    s0  = -1.5d8/unit_length 
    s1  =  (Bv-B_y)*(y_r-s0)**2!<= see own notes
    !<= s1 fixes the value of |B| = 8G @ (0,y_r)     

    !=> input for random drivers
    if(driver_random)then
    nxmodes=10
    randphase(1:10) = zero
    randA(1:10)     = zero
    randP(1:10)     = zero
    if(nxmodes>10) call mpistop('too many modes, edit nxmodes')

    if(mype==0)then
      call random_seed(SIZE=seed_size)
      allocate(seed(seed_size))
      call random_seed(GET=seed(1:seed_size))
      call random_number(randphaseP1(1:nxmodes))
      randphase(1:nxmodes)=-dpi+2.0d0*dpi*dble(randphaseP1(1:nxmodes))
      call random_number(randA1(1:nxmodes))
      randA(1:nxmodes)=(1.0d3+floor(randA1(1:nxmodes)*1.01d5))/unit_velocity !1m/s-1km/s
      call random_number(randP1(1:nxmodes))
      randP(1:nxmodes)=(10.0d0+floor(randP1(1:nxmodes)*101.0d0))/unit_time
    endif
    call MPI_BARRIER(icomm,ierrmpi)
    if(npe>1)then
         call MPI_BCAST(randphase,10,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
         call MPI_BCAST(randA,10,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
         call MPI_BCAST(randP,10,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif
    if(mype==0)then
!      print *,'number of modes=',nxmodes
      open(123,file='phaseinfo',form='formatted')
      write(123,*) nxmodes
      do ix=1,nxmodes
          write(123,*) ix,randphase(ix),randA(ix)*unit_velocity,&
             randP(ix)*unit_time
      enddo
      close(123)
    endif
    endif

! Should leave this alone if posssible
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
!    use mod_global_parameters
  !! initialize the table in a vertical line through the global domain
    integer :: j,na,nb,ibc
    double precision, allocatable :: Ta(:),gg(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
  
    rpho=1.151d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=0.2d0 ! height of initial transition region
    wtra=0.02d0 ! width of initial transition region 
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**&
       3
 
    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax),Ha(jmax))
   !Setting up atmospheric profiles
   !=> set up of tanh profile
   if(tanh_profile) then    
    !=> creates temperture profile
    do j=1,jmax
       ya(j)=(dble(j)-0.5d0)*dya-gzone
      !<=remove for steeper T profile 
      if(ya(j)>htra) then
         Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
       else
         Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
       endif
       gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo
    !!=> solution of hydrostatic equation
    ra(1)=rpho
    pa(1)=rpho*Tpho
    invT=gg(1)/Ta(1) !<1/H(y)
    invT=0.d0
   endif

   ! methods for constructing 1D profiles
   if(integrate) then
   !=>scale height for HS equation
    do j=2,jmax
       invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
       pa(j)=pa(1)*dexp(invT*dya)
       ra(j)=pa(j)/Ta(j)
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya !<= Residual
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do
    endif
   if(derivative)then
    do j=2,jmax
      pa(j)=(pa(j-1)+dya*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dya*(gg(j)+&
         gg(j-1))/Ta(j)/4.d0)
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
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do
   endif

   if(derivative_2)then
    do j=2,jmax
      pa(j)=ra(j)*Ta(j)
    end do
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do
   endif

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    ! initialize one grid
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: res
    integer :: ix1,ix2,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
        write(*,*)'dimensionless vars:'
        write(*,*)'T =', unit_time, 's'
        write(*,*)'L =', unit_length, 'cm'
        write(*,*)'V =', unit_velocity, 'cm s-1'
        write(*,*)'rho =', unit_numberdensity, 'cm-3'
        write(*,*)'rho =', unit_density, 'g cm-3'
        write(*,*)'B =', unit_magneticfield, 'G'
        write(*,*)'T =', unit_temperature, 'K'
        write(*,*)'p =', unit_pressure, 'dyn cm-2'

        if(tanh_profile) write(*,*)'tanh atmos =', tanh_profile
        if(c7_profile) write(*,*)'C7 atmos =', c7_profile
        if(driver) write(*,*)'driver Fedun=', driver
        if(driver_random) write(*,*)'driver random =', driver_random
        if(driver_kuz) write(*,*)'driver Kuzma =', driver_kuz
        if(driver_injetion) write(*,*)'jet driver 1  =', driver_injetion
        if(driver_injetion_1) write(*,*)'jet driver 2  =', driver_injetion_1
        if(jet_cont) write(*,*)'jet_cont =',jet_cont 
        if(jet_switch) write(*,*) 'jet_switch =', jet_switch
      endif
      first=.false.
    endif

    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
        na=floor((x(ix1,ix2,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix1,ix2,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix1,ix2,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix1,ix2,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    end do
    end do
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(:))=zero
    if(B0field) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=zero
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2))*dcos(theta)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))= B0*dsin(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2))*dsin(theta)
    endif

    if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)=0.d0

    call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, iB, ixImin1,&
       ixImin2,ixImax1,ixImax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2),ggrid(ixImin1:ixImax1,&
       ixImin2:ixImax2),invT(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: delydelx
    integer :: ix1,ix2,idir,ixIntmin1,ixIntmin2,ixIntmax1,ixIntmax2

    select case(iB)
    case(3)
      !! fixed zero velocity
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) =-w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))/w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=0.d0
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(1))=-B0*dcos(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*dcos(theta)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(2))= B0*dsin(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(3))=-B0*dcos(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*dsin(theta)
      endif
      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)=0.d0
      call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case(4)
      ixIntmin1=ixOmin1;ixIntmin2=ixOmin2;ixIntmax1=ixOmax1;ixIntmax2=ixOmax2;
      ixIntmin2=ixOmin2-1;ixIntmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixIntmin1,&
         ixIntmin2,ixIntmax1,ixIntmax2,pth)
      ixIntmin2=ixOmin2-1;ixIntmax2=ixOmax2;
      call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixIntmin1,ixIntmin2,&
         ixIntmax1,ixIntmax2,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin1:ixOmax1,ixOmin2-1)=w(ixOmin1:ixOmax1,ixOmin2-1,&
         rho_)/pth(ixOmin1:ixOmax1,ixOmin2-1)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin1:ixOmax1,ixOmin2-1)=tmp(ixOmin1:ixOmax1,&
           ixOmin2-1)+0.5d0*(ggrid(ixOmin1:ixOmax1,ix2)+ggrid(ixOmin1:ixOmax1,&
           ix2-1))*invT(ixOmin1:ixOmax1,ixOmin2-1)
        w(ixOmin1:ixOmax1,ix2,p_)=pth(ixOmin1:ixOmax1,&
           ixOmin2-1)*dexp(tmp(ixOmin1:ixOmax1,ixOmin2-1)*dxlevel(2))
        w(ixOmin1:ixOmax1,ix2,rho_)=w(ixOmin1:ixOmax1,ix2,&
           p_)*invT(ixOmin1:ixOmax1,ixOmin2-1)
      enddo
      !> fixed zero velocity
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) =-w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))/w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      !> zero normal gradient extrapolation
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* (-w(ixOmin1:ixOmax1,ix2-2,&
           mag(:))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
      enddo
      if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)=0.d0
      call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr

  subroutine gravity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,wCT,x,gravity_field)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ndim)

    double precision                :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    gravity_field=0.d0
    call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,x)
    gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=ggrid(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine gravity

  subroutine getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=usr_grav*(SRadius/(SRadius+&
       x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))**2
  end subroutine

  subroutine special_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2),&
       bQgrid(ixImin1:ixImax1,ixImin2:ixImax2)

    ! add global background heating bQ
    call getbQ(bQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,qtC,wCT,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+qdt*bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

  end subroutine special_source

  subroutine getbQ(bQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
  ! calculate background heating bQ
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: bQgrid(ixImin1:ixImax1,ixImin2:ixImax2)

    bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=bQ0*dexp(-x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)/5.d0)

  end subroutine getbQ

  subroutine special_refine_grid(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif

  end subroutine special_refine_grid

  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       B2(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2),dRdT(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: ens(ixImin1:ixImax1,ixImin2:ixImax2),&
       divb(ixImin1:ixImax1,ixImin2:ixImax2),wlocal(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    double precision :: Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
       curlvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    integer :: idirmin,idir,ix1,ix2

    wlocal(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,pth)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,mag(idir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
           idir,0)
      else
        Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,mag(idir))
      endif
    end do
    ! B^2
    B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum((Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=dsqrt(B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

    ! output divB1
    call get_normalized_divb(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=divb(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! output the plasma beta p*2/B**2
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*two/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    ! output heating rate
    call getbQ(ens,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,global_time,wlocal,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5)=ens(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! store the cooling rate 
    if(mhd_radiative_cooling)call getvar_cooling(ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wlocal,x,ens)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6)=ens(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    ! store current
    call get_current(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idirmin,curlvec)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6+idir)=curlvec(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
    end do
  
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Alfv divB beta bQ rad j1 j2 j3'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)

    wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2))*dcos(theta)
    wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=+B0*dsin(kx*x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
    wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2))*dsin(theta)

  end subroutine specialset_B0

end module mod_usr
