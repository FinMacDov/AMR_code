
module mod_usr
  use mod_hd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:), pubc(:),rubc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya,BB1,BB2,&
     BB3
  double precision, allocatable :: pa(:),ra(:),ya(:),Ha(:)
  integer, parameter :: jmax=8000

  double precision, allocatable :: p_profile(:), rho_profile(:),&
      Te_profile(:) 
  integer, parameter :: kmax=8000
  double precision :: s0,s1,Bv,B_y,y_r, dyk
  logical :: driver, tanh_profile, c7_profile, driver_kuz, driver_random 
  double precision :: randphase(10), randA(10), randP(10)
  integer :: nxmodes

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    hd_gamma=1.66666667d0

    unit_length        = 1.d8 !cm = 1 Mm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm-3

!    unit_pressure = unit_temperature*unit_density

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_var_for_errest  => Te_for_errest

    call set_coordinate_system("Cartesian")
    call hd_activate()
    call params_read(par_files)
  end subroutine usr_init

!> Read parameters from a file
  subroutine params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /my_switches/ driver, driver_kuz, driver_random, tanh_profile,&
      c7_profile
  do n = 1, size(files)
   open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_switches, end=111)
111    close(unitpar)
    end do
  end subroutine params_read

  subroutine initglobaldata_usr()
    use mod_global_parameters
    integer :: iv
    integer,dimension(:),allocatable:: seed
    integer::  seed_size,ix
    real:: randphaseP1(1:10), randA1(1:10), randP1(1:10)
 
   ! unit_pressure = unit_temperature*unit_density
    heatunit=unit_pressure/unit_time !3.697693390805347E-003 erg*cm-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax1-xprobmin1)/dble(jmax) !cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 !the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade
    kx=dpi/((xprobmax1-xprobmin1)/2.0d0)
    ly=kx*dcos(theta)
    SRadius=6.96d10/unit_length ! Solar radius
    !SRadius=69.61d0 ! Solar radius
    
    !=> in Guass
    BB1=0.0d0/unit_magneticfield
    BB2=10.0d0/unit_magneticfield
    BB3=0.0d0/unit_magneticfield
    
    !=>Konkol et al. 2012 and Kuzma et al. 2017
    !=>s1 = s and s0 = a ("a" is not sound speed)
    B_y = 8.0d0/unit_magneticfield
    y_r = 1.0d9/unit_length
    Bv  = 6.0d0/unit_magneticfield
    s0  = -1.5d8/unit_length 
    s1  =  (Bv-B_y)*(y_r-s0)**2!<= see own notes
    !<= s1 fixes the value of |B| = 8G @ (0,y_r)     
    
    !=> To allow output to be in physical uints
    length_convert_factor = 1.0d0!unit_length
    time_convert_factor = unit_time
    w_convert_factor(1) = unit_density
    w_convert_factor(5) = unit_pressure

    do iv = 2,4
      w_convert_factor(iv) = unit_velocity
      w_convert_factor(iv+4) = unit_magneticfield
    enddo 

    !=> hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

  end subroutine initglobaldata_usr

  subroutine inithdstatic
  !! initialize the table in a vertical line through the global domain
    use mod_global_parameters

    integer :: j,na,nb,ibc, i, k
    double precision, allocatable :: Ta(:),gg(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
    integer :: step2,nb_pts 
    double precision:: step1
    rpho=1.151d15/unit_numberdensity ! number density at the bottom
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.8d6/unit_temperature ! estimated temperature in the top
    htra=1.95d0!0.2d0 ! height of initial transition region
    wtra=0.01d0!0.02d0 ! width of initial transition region
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**&
       3
   
   !=> set up of tanh profile
   if(tanh_profile) then    
    !=> creates temperture profile
    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))
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
    nb=int(gzone/dya)
    ra(1)=rpho
    pa(1)=rpho*Tpho
    invT=gg(1)/Ta(1) !<1/H(y)
    invT=0.d0
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

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif
   endif

 !=> set up of c7 profile
   if(c7_profile) then   
    allocate(ya(kmax),Ta(kmax),gg(kmax),pa(kmax),ra(kmax),Ha(kmax))
   !=> Data imported here is in SI units
!   open (unit = 11, file ="atmos_data/c7/fort/rho.dat", status='old')
!   open (unit = 12, file ="atmos_data/c7/fort/Temp.dat", status='old')
!   open (unit = 13, file ="atmos_data/c7/fort/y.dat", status='old')

   open (unit = 11, file ="atmos_data/c7/1dinterp/c7_rho.dat", status='old')
   open (unit = 12, file ="atmos_data/c7/1dinterp/c7_Te.dat", status='old')
   open (unit = 13, file ="atmos_data/c7/1dinterp/c7_y.dat", status='old')

  do i=1,kmax  
   read(11,*) ra(i) !kg m-3
   read(12,*) Ta(i) !K
   read(13,*) ya(i) !0-10Mm
   ra(i) = ra(i)*0.001d0/unit_density !SI to cgs to dimensionless 
   Ta(i) = Ta(i)/unit_temperature
   gg(i)=usr_grav*(SRadius/(SRadius+ya(i)))**2
  end do 
   close(11)
   close(12)
   close(13)

    dyk = ya(2) ! Cells size for C7 data
    pa(1)=ra(1)*Ta(1)
    invT=gg(1)/Ta(1) !<1/H(z)
    invT=0.d0
    !=>used scale height for HS equation
    do i=2,kmax
       invT=invT+(gg(i)/Ta(i)+gg(i-1)/Ta(i-1))*0.5d0
       Ha(i)=-1.0d0/((gg(i)/Ta(i)+gg(i-1)/Ta(i-1))*0.5d0)
       pa(i)=pa(1)*dexp(invT*dyk)
       ra(i)=pa(i)/Ta(i)
       if(mype==0)then
       endif
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dyk+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dyk !<= Residual
    rhob=ra(na)+res/dyk*(ra(na+1)-ra(na))
    pb=pa(na)+res/dyk*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(1,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0))/dyk+0.5d0)
      res=gzone-dx(1,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dyk
      rbc(ibc)=ra(na)+res/dyk*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dyk*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif
 endif
  end subroutine inithdstatic

  subroutine initonegrid_usr(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    ! initialize one grid
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: res,step1
    integer :: ix1,na,i
    integer :: nb_pts
    logical, save :: first=.true.
    double precision :: width, A, y0,x0

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
        write(*,*)'dimensionless vars:'
        write(*,*)'T =', unit_time, 's'
        write(*,*)'L =', unit_length, 'cm'
        write(*,*)'V =', unit_velocity, 'cm s-1'
        write(*,*)'rho =', unit_numberdensity, 'g cm-3'
        write(*,*)'rho =', unit_density, 'g cm-3'
        write(*,*)'B =', unit_magneticfield, 'G'
        write(*,*)'T =', unit_temperature, 'K'
        write(*,*)'p =', unit_pressure, 'dyn cm-2'
        write(*,*)'tanh atmos =', tanh_profile
        write(*,*)'C7 atmos =', c7_profile
        write(*,*)'driver =', driver
        write(*,*)'driver Fedun=', driver
        write(*,*)'driver random =', driver_random
        write(*,*)'driver Kuzma=', driver_kuz
      endif
      first=.false.
    endif
      
    if(.NOT.firstprocess)then
    if(tanh_profile) then
    do ix1=ixOmin1,ixOmax1
        na=floor((x(ix1,1)-xprobmin1+gzone)/dya+0.5d0)
        res=x(ix1,1)-xprobmin1+gzone-(dble(na)-0.5d0)*dya
        w(ix1,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix1,p_)=pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    end do
    endif

    if(c7_profile) then
!    dyk = ya(2) ! Cells size for C7 data
!    nb_pts = domain_nx2+2*nghostcells 
!=> interpolation to obtain rho & p of HSE
    do ix1=ixOmin1,ixOmax1
        na=floor((x(ix1,1)-xprobmin1+gzone)/dyk+0.5d0)
        res=x(ix1,1)-xprobmin1+gzone-(dble(na)-0.5d0)*dyk
        w(ix1,rho_)=ra(na)+(one-cos(dpi*res/dyk))/two*(ra(na+1)-ra(na))
        w(ix1,p_)=pa(na)+(one-cos(dpi*res/dyk))/two*(pa(na+1)-pa(na))
    enddo
    endif
    endif
 
!    if(.NOT.firstprocess)then
     call phys_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
!    endif
  end subroutine initonegrid_usr

  subroutine Te_for_errest(ixImin1,ixImax1,ixOmin1,ixOmax1,iflag,w,x,var)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1,iflag
    double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw),x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: var(ixImin1:ixImax1)
    double precision :: pth(ixImin1:ixImax1)

    call hd_get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
    var(ixOmin1:ixOmax1) = pth(ixOmin1:ixOmax1)/w(ixOmin1:ixOmax1,rho_)
  end subroutine Te_for_errest

  subroutine specialbound_usr(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixOmin1,ixOmax1, iB, ixImin1,ixImax1
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: v_sum
    double precision :: pth(ixImin1:ixImax1),tmp(ixImin1:ixImax1),&
       ggrid(ixImin1:ixImax1),invT(ixImin1:ixImax1)
    double precision :: delydelx, x0, y0, width
    double precision ::jet_w, jet_h, A, period, deltax, deltay, step1
    integer :: ind1, na, i
    integer :: ix1,idir,ixIntmin1,ixIntmax1
    integer :: nb_pts
    double precision :: rand_driv(10)


    select case(iB)
    case(1)
      !! fixed zero velocity
      do idir=1,ndir
        w(ixOmin1:ixOmax1,mom(idir)) =-w(ixOmax1+nghostcells:ixOmax1+1:-1,&
           mom(idir))/w(ixOmax1+nghostcells:ixOmax1+1:-1,rho_)
      end do

      do ix1=ixOmin1,ixOmax1
        w(ix1,rho_)=rbc(ix1)
        w(ix1,p_)=pbc(ix1)
        if(mype==1)then
         write(*,*) rbc(ix1)*unit_density,pbc(ix1)*unit_pressure
        endif
      enddo      
      call phys_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    case(2)
      ixIntmin1=ixOmin1-1;ixIntmax1=ixOmax1;
      call hd_get_pthermal(w,x,ixImin1,ixImax1,ixIntmin1,ixIntmax1,pth)
      ixIntmin1=ixOmin1-1;ixIntmax1=ixOmax1;
      call getggrav(ggrid,ixImin1,ixImax1,ixIntmin1,ixIntmax1,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,rho_)/pth(ixOmin1:ixOmax1)
      tmp=0.d0
      do ix1=ixOmin1,ixOmax1
        tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+&
           0.5d0*(ggrid(ixOmin1:ixOmax1)+&
           ggrid(ixOmin1:ixOmax1))*invT(ixOmin1:ixOmax1)
        w(ixOmin1:ixOmax1,p_)=pth(ixOmin1:ixOmax1)*dexp(tmp(&
           ixOmin1:ixOmax1)*dxlevel(1))
        w(ixOmin1:ixOmax1,rho_)=w(ixOmin1:ixOmax1,p_)*invT(ixOmin1:ixOmax1)
      enddo

      !> fixed zero velocity
      do idir=1,ndir
        w(ixOmin1:ixOmax1,mom(idir)) =-w(ixOmin1-1:ixOmin1-nghostcells:-1,&
           mom(idir))/w(ixOmin1-1:ixOmin1-nghostcells:-1,rho_)
      end do
      call phys_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr

  subroutine gravity(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,ndim)

    double precision                :: ggrid(ixImin1:ixImax1)

    gravity_field=0.d0
    call getggrav(ggrid,ixImin1,ixImax1,ixOmin1,ixOmax1,x)
    gravity_field(ixOmin1:ixOmax1,1)=ggrid(ixOmin1:ixOmax1)

  end subroutine gravity

  subroutine getggrav(ggrid,ixImin1,ixImax1,ixOmin1,ixOmax1,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out)   :: ggrid(ixImin1:ixImax1)

    ggrid(ixOmin1:ixOmax1)=usr_grav*(SRadius/(SRadius+x(ixOmin1:ixOmax1,&
       1)))**2
  end subroutine

  subroutine specialvar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixImin1:ixImax1),B2(ixImin1:ixImax1),&
       tmp2(ixImin1:ixImax1),dRdT(ixImin1:ixImax1)
    double precision :: ens(ixImin1:ixImax1),divb(ixImin1:ixImax1)
    double precision :: Btotal(ixImin1:ixImax1,1:ndir),curlvec(ixImin1:ixImax1,&
       1:ndir)
    integer :: idirmin,idir,ix1

    double precision:: gradrho(ixGlo1:ixGhi1),rho(ixGlo1:ixGhi1),&
       drho(ixGlo1:ixGhi1)
    double precision:: kk,kk0,grhomax,kk1
    integer                            :: idims

    ! output temperature
    call hd_get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
    w(ixOmin1:ixOmax1,nw+1)=unit_temperature*pth(ixOmin1:ixOmax1)/w(&
       ixOmin1:ixOmax1,rho_)

     rho(ixImin1:ixImax1)=w(ixImin1:ixImax1,rho_)
     gradrho(ixOmin1:ixOmax1)=zero
     do idims=1,ndim
       select case(typegrad)
       case("central")
         call gradient(rho*unit_density,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,&
            drho)
       case("limited")
         call gradientS(rho*unit_density,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,&
            drho)
       end select
       gradrho(ixOmin1:ixOmax1)=gradrho(ixOmin1:ixOmax1)+&
          drho(ixOmin1:ixOmax1)**2.0d0
     enddo

     gradrho(ixOmin1:ixOmax1)=dsqrt(gradrho(ixOmin1:ixOmax1))
     kk=5.0d0
     kk0=0.01d0
     kk1=1.0d0
     grhomax=1000!10.0d0

  ! putting the schlierplot of density in nwauxio=1
     w(ixOmin1:ixOmax1,nw+2)=dexp(-kk*(gradrho(ixOmin1:ixOmax1)-&
        kk0*grhomax)/(kk1*grhomax-kk0*grhomax))
   
    w(ixOmin1:ixOmax1,nw+3)=unit_velocity*dsqrt(hd_gamma*pth(&
       ixOmin1:ixOmax1)/w(ixOmin1:ixOmax1,rho_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te schrho cs'

  end subroutine specialvarnames_output


end module mod_usr
