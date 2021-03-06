module mod_usr
  use mod_mhd
  implicit none
  integer :: npts
  double precision, allocatable :: pbc(:),rbc(:)
  double precision, allocatable :: pa(:),ra(:),ya(:)
  double precision :: gzone,dya
  double precision, allocatable :: p_profile(:), rho_profile(:)
  logical :: tanh_profile, c7_profile
  !> Name of temperture profile chosen 
  character(len=std_len), private  :: Te_profile

contains

  subroutine usr_init()
    integer :: iv    
    call set_coordinate_system("Cartesian_2D")


   SI_unit = .False.
   if(SI_unit)then
    unit_length        = 1.d7 ! m
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d3 ! m-3,m-3
   else
    unit_length        = 1.d8 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm-3,cm-3
   endif

    usr_set_parameters      => initglobaldata_usr
    usr_init_one_grid       => initonegrid_usr
    usr_special_bc          => specialbound_usr
    usr_aux_output          => specialvar_output
    usr_add_aux_names       => specialvarnames_output 
    usr_var_for_errest      => p_for_errest

    call mhd_activate()
    call params_read(par_files)
  end subroutine usr_init

!> Read parameters from a file
  subroutine params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /my_switches/ tanh_profile,c7_profile, /atmos_list/ npts,Te_profile
  do n = 1, size(files)
   open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_switches, end=111)
111    close(unitpar)
        open(unitpar, file=trim(files(n)), status="old")
        read(unitpar, atmos_list, end=112)
112     close(unitpar)
  end do
  end subroutine params_read

  subroutine initglobaldata_usr()
    use mod_global_parameters
 
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    if(c7_profile) call inithdstatic

  end subroutine initglobaldata_usr

  subroutine inithdstatic
  !! initialize the table in a vertical line through the global domain
    use mod_global_parameters
    use mod_solar_atmosphere

    integer :: j,na,ibc, i, k
    double precision, allocatable :: Ta(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa

    allocate(ya(npts),Ta(npts),pa(npts),ra(npts))

 !=> set up of c7 profile
   if(c7_profile) then   
    call solar_atmosphere_init(mhd_gamma,He_abundance)
    do i=1,npts  
     Ta(i) = Te_intrpl(i)/unit_temperature
     ra(i) = 1.0d0/Ta(i)
     pa(i)=ra(i)*Ta(i)
    end do 
  !number density in cgs units
    dya = abs(z_pos(2)-z_pos(1)) ! Cells size for C7 data
   endif

   !=> methods
   !=>scale height for HS equation
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

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
!    use mod_global_parameters
  ! initialize one grid
    integer :: iv
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: htra, wtra, rpho
    logical, save:: first=.true., units=.true.
    double precision :: res
    integer :: ix1,ix2,na,i
    
    if(units)then
      !=> To allow output to be in physical uints
      length_convert_factor = unit_length
      time_convert_factor = unit_time
      w_convert_factor(1) = unit_density
      w_convert_factor(4) = unit_pressure
      do iv = 2,3
        w_convert_factor(iv) = unit_velocity
        w_convert_factor(iv+3) = unit_magneticfield
      enddo 
    endif
    if (first) then
       if (mype==0) then
          print *,'YOKOYAMA and SHIBATA 2001 ApJ'
          print *, 'unit_length = ', unit_length
          print *, 'unit_time = ', unit_time
          print *, 'unit_density = ', unit_density
          print *, 'unit_pressure = ', unit_pressure
          print *, 'unit_velocity = ', unit_velocity
          print *, 'unit_magneticfield = ', unit_magneticfield
       end if
       first=.false.
    end if

    if(tanh_profile) then
      rpho=1.d5 ! number density at the bottom relaxla
      htra=2.0d0 ! height of initial transition region
      wtra=0.06d0 ! width of initial transition region 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=1.d0+&
         (rpho-1.d0)*(1.d0-tanh((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)-htra)/wtra))/2.d0
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)=1.d0
    endif
    
     if(c7_profile) then
      do ix2=ixImin2,ixImax2
      do ix1=ixImin1,ixImax1
        na=floor((x(ix1,ix2,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix1,ix2,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix1,ix2,rho_)=ra(na)+(one-cos(dpi*res/dya))/2.0d0*(ra(na+1)-ra(na))
        w(ix1,ix2,p_)=pa(na)+(one-cos(dpi*res/dya))/2.0d0*(pa(na+1)-pa(na))
      end do
      end do
    endif

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(:))=zero
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=zero
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))=10.0d0/unit_magneticfield

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
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: ix1,ix2, ixAmin1,ixAmin2,ixAmax1,ixAmax2


  end subroutine specialbound_usr

  subroutine p_for_errest(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iflag,w,x,var)
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,iflag
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: var(ixImin1:ixImax1,ixImin2:ixImax2)

    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,var)
    
  end subroutine p_for_errest

  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       B2(ixImin1:ixImax1,ixImin2:ixImax2),divb(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
       current_o(ixImin1:ixImax1,ixImin2:ixImax2,3)
    integer :: idir,idirmin

    ! output temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,e_)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(2))
    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,pth)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)
    if(B0field) then
      Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndir,0)
    else
      Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))
    endif
    ! B^2
    B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum((Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5)=dsqrt(B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))*unit_velocity
    ! output divB1
    call divvector(Btotal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6)=0.5d0*divb(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/dsqrt(B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))/(1.0d0/dxlevel(1)+1.0d0/dxlevel(2))
    ! output the plasma beta p*2/B**2
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+7)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*two/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='e m1 m2 Te Alfv divB beta'
  end subroutine specialvarnames_output


end module mod_usr
