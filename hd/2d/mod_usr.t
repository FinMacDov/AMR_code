module mod_usr
  use mod_hd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:), pubc(:),rubc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya,BB1,BB2,BB3
  double precision, allocatable :: dpa(:),pa(:),ra(:),ya(:),Ha(:)
  integer, parameter :: jmax=8000

  double precision, allocatable :: p_profile(:), rho_profile(:), Te_profile(:) 
  double precision :: s0,s1,Bv,B_y,y_r, dyk
  logical :: driver, tanh_profile, c7_profile, driver_kuz, driver_random, testing
  logical :: derivative, integrate
  double precision :: randphase(10), randA(10), randP(10)
  integer :: nxmodes

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2.5D")

    hd_gamma=1.66666667d0

    unit_length        = 1.d8                                         ! cm = 1 Mm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

!    unit_pressure = unit_temperature*unit_density

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_var_for_errest  => Te_for_errest

    call hd_activate()
    call params_read(par_files)
  end subroutine usr_init

!> Read parameters from a file
  subroutine params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /my_switches/ driver, driver_kuz, driver_random, tanh_profile, & 
                         c7_profile, derivative, integrate
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
    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade
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
    enddo 

    !=> hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic
    
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
      randphase(1:nxmodes)=-dpi+two*dpi*dble(randphaseP1(1:nxmodes))
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
      print *,'number of modes=',nxmodes
      open(123,file='phaseinfo',form='formatted')
      write(123,*) nxmodes
      do ix=1,nxmodes
          write(123,*) ix,randphase(ix),randA(ix)*unit_velocity,randP(ix)*unit_time
      enddo
      close(123)
    endif
    endif

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
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
   
   !=> set up of tanh profile
    !=> creates temperture profile
    allocate(ya(jmax),Ta(jmax),gg(jmax),dpa(jmax),pa(jmax),ra(jmax),Ha(jmax))
  if(tanh_profile) then   
        dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
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
    ra(1)=rpho
    pa(1)=rpho*Tpho
   endif

 !=> set up of c7 profile
   if(c7_profile) then   
   !=> Data imported here is in SI units
   open (unit = 11, file ="atmos_data/c7/1dinterp/c7_rho.dat", status='old')
   open (unit = 12, file ="atmos_data/c7/1dinterp/c7_Te.dat", status='old')
   open (unit = 13, file ="atmos_data/c7/1dinterp/c7_y.dat", status='old')

  do i=1,jmax  
   read(11,*) ra(i) !kg m-3
   read(12,*) Ta(i) !K
   read(13,*) ya(i) !0-10Mm
   ra(i) = ra(i)*1d-3/unit_density !SI to cgs to dimensionless 
   Ta(i) = Ta(i)/unit_temperature
   gg(i)=usr_grav*(SRadius/(SRadius+ya(i)))**2
  end do 
   close(11)
   close(12)
   close(13)
   
   dya = ya(2) ! Cells size for C7 data
   pa(1)=ra(1)*Ta(1)
   endif
 
   !!=> solution of hydrostatic equation

    if(integrate) then
    invT=gg(1)/Ta(1) !<1/H(y)
    invT=0.d0
    !=>scale height for HS equation
    do j=2,jmax
       invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
       pa(j)=pa(1)*dexp(invT*dya)
       ra(j)=pa(j)/Ta(j)
    end do
    endif

   if(derivative)then
    pa(1)=ra(1)*Ta(1)
    do j=2,jmax
      pa(j)=(pa(j-1)+dya*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dya*(gg(j)+gg(j-1))/Ta(j)/4.d0)
      ra(j)=pa(j)/Ta(j)
    end do
!    do j=2,jmax-1
!      dpa(j)=pa(j+1)-2.0d0*pa(j)+pa(j-1)
!      ra(j)=-dpa(j)/gg(j)
!    end do
   endif

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
     print*,'dy',dya     
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res,step1
    integer :: ix^D,na,i,idims
    integer :: nb_pts
    logical, save :: first=.true.
    double precision :: width, A, y0,x0

    double precision:: gradp(ixG^T),dp(ixG^T), ggrid(ixI^S), p(ixG^T)

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
        write(*,*)'dimensionless vars:'
        write(*,*)'T =', unit_time, 's'
        write(*,*)'L =', unit_length, 'cm'
        write(*,*)'V =', unit_velocity, 'cm s-1'
        write(*,*)'rho =', unit_numberdensity, 'g cm-3'
        write(*,*)'rho =', unit_density, 'g cm-3'
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
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix^D,p_)=pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}

!     p(ixI^S)=w(ixI^S,p_)
!     gradp(ixO^S)=zero
!    do idims=1,ndim
!       select case(typegrad)
!       case("central")
!         call gradient(p,ixI^L,ixO^L,idims,dp)
!       case("limited")
!         call gradientS(p,ixI^L,ixO^L,idims,dp)
!       end select
!       gradp(ixO^S)=gradp(ixO^S)+dp(ixO^S)**2.0d0
!     enddo
!     gradp(ixO^S)=dsqrt(gradp(ixO^S))
!     if(mype==0) then
!     write(*,*) gradp(ixO^S) 
!     endif
!     call getggrav(ggrid,ixI^L,ixO^L,x)
!     w(ixO^S,rho_) = gradp(ixO^S)/ggrid(ixO^S) 


    if(testing) then
!    dyk = ya(2) ! Cells size for C7 data
!    nb_pts = domain_nx2+2*nghostcells 
!=> interpolation to obtain rho & p of HSE
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix^D,p_)=pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    enddo
!        write(*,*)na,floor(x(ix^D,2)+0.5d0),pa(na),ra(na),dya,floor(x(ix^D,2)-xprobmin2+gzone+0.5d0) 
    enddo

    endif
    endif

      if(driver_kuz)then
      width = 0.1 !Mm
      A = 2e5/unit_velocity !<=2 km/s
      x0 = (xprobmax1-abs(xprobmin1))/2.0d0+0.0d0!<=x origin
      y0 =0.7d0 !<=y origin
      w(ixI^S,mom(2))= A*dexp(-((x(ixI^S,1)-x0)**2+(x(ixI^S,2)-y0)**2)/&
      width**2)
      endif

    if(.NOT.firstprocess)then
     call phys_to_conserved(ixI^L,ixO^L,w,x)
    endif
  end subroutine initonegrid_usr

  subroutine Te_for_errest(ixI^L,ixO^L,iflag,w,x,var)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)
    double precision :: pth(ixI^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    var(ixO^S) = pth(ixO^S)/w(ixO^S,rho_)
  end subroutine Te_for_errest

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: v_sum
    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: delydelx, x0, y0, width
    double precision ::jet_w, jet_h, A, period, deltax, deltay, step1
    integer :: ind^D, na, i
    integer :: ix^D,idir,ixInt^L
    integer :: nb_pts
    double precision :: rand_driv(10)


    select case(iB)
    case(3)
      !! fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do

      !! fixed gravity stratification of density and pressure pre-determined in initial condition

      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
        if(mype==1)then
         write(*,*) rbc(ix2)*unit_density,pbc(ix2)*unit_pressure
        endif
      enddo

      !=> Driver
      if(driver) then 
      jet_w = (xprobmax1-xprobmin1)/domain_nx1 !<= 1 cell radius 
      jet_h = (xprobmax2-xprobmin2)/domain_nx2
      A = 5.0d4/unit_velocity !500 m/s !5.0d6/unit_velocity!
      deltax = (xprobmax1-xprobmin1)/domain_nx1
      deltay = (xprobmax2-xprobmin2)/domain_nx2
      x0 =(xprobmax1-abs(xprobmin1))/2.0d0+0.0d0 !<=x origin
      y0 =-gzone !<=y origin
      period = 30.0d0/unit_time
        w(ixO^S,mom(2))= A*dsin(2.0d0*dpi*qt/period)*&
                         dexp(-(((x(ixO^S,1)-x0)/deltax)**2&
                         +((x(ixO^S,2)-y0)/deltay)**2))
      endif

      if(driver_random)then
      deltax = (xprobmax1-xprobmin1)/domain_nx1
      deltay = (xprobmax2-xprobmin2)/domain_nx2
      x0 =(xprobmax1-abs(xprobmin1))/2.0d0+0.0d0 !<=x origin
      y0 =-gzone !<=y origin
      do i = 1,10
         rand_driv(i) = randA(i)*sin(2*dpi*qt/randP(i)+randphase(i)) 
!         write(*,*)i,mype,rand_driv(i),randA(i),randP(i),randphase(i) 
      enddo

        w(ixO^S,mom(2))= sum(rand_driv)*&
                         dexp(-(((x(ixO^S,1)-x0)/deltax)**2&
                         +((x(ixO^S,2)-y0)/deltay)**2))
      endif
      
      call phys_to_conserved(ixI^L,ixO^L,w,x)
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

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    integer :: idirmin,idir,ix^D

    double precision:: gradrho(ixG^T),rho(ixG^T),drho(ixG^T)
    double precision:: gradp(ixG^T),dp(ixG^T), ggrid(ixI^S), p(ixG^T)
    double precision:: kk,kk0,grhomax,kk1
    integer         :: idims
    logical, save   :: firstrun=.true.

    ! output temperature
    call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=unit_temperature*pth(ixO^S)/w(ixO^S,rho_)

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
     grhomax=1000!10.0d0

  ! putting the schlierplot of density in nwauxio=1
     w(ixO^S,nw+2)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))

     w(ixO^S,nw+3)=gradrho(ixO^S)

   
    w(ixO^S,nw+4)=unit_velocity*dsqrt(hd_gamma*pth(ixO^S)/w(ixO^S,rho_))

     p(ixI^S)=pth(ixI^S)
     gradp(ixO^S)=zero
     do idims=1,ndim
       select case(typegrad)
       case("central")
         call gradient(p,ixI^L,ixO^L,idims,dp)
       case("limited")
         call gradientS(p,ixI^L,ixO^L,idims,dp)
       end select
       gradp(ixO^S)=gradp(ixO^S)+dp(ixO^S)**2.0d0
     enddo
     gradp(ixO^S)=dsqrt(gradp(ixO^S))

   
   call getggrav(ggrid,ixI^L,ixO^L,x)
!gradp is fine
   w(ixO^S,nw+5)=-gradp(ixO^S)+w(ixO^S,rho_)*ggrid(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te schrho delrho cs fb'

  end subroutine specialvarnames_output


end module mod_usr
