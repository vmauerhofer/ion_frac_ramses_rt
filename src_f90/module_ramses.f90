module module_ramses

  use module_constants, only : kB, mp, XH, mSi, mMg, planck, clight,cmtoA

  implicit none

  private 


  ! Stop pretending this would work in 2D 
  integer(kind=4),parameter :: ndim = 3
  integer(kind=4),parameter :: twondim = 6
  integer(kind=4),parameter :: twotondim= 8 

  ! conversion factors (units)
  logical                        :: conversion_scales_are_known = .False. 
  real(kind=8)                   :: dp_scale_l,dp_scale_d,dp_scale_t
  real(kind=8)                   :: dp_scale_T2,dp_scale_zsun,dp_scale_nH
  real(kind=8)                   :: dp_scale_nHe,dp_scale_v,dp_scale_m,dp_scale_pf

  ! cooling-related stuff -------------------------------------------------------------
  type cooling_table
     integer(kind=4)          :: n11
     integer(kind=4)          :: n22
     real(kind=8),allocatable :: nH(:)    
     real(kind=8),allocatable :: T2(:)    
     real(kind=8),allocatable :: metal(:,:)  
     real(kind=8),allocatable :: cool(:,:)  
     real(kind=8),allocatable :: heat(:,:)  
     real(kind=8),allocatable :: cool_com(:,:)  
     real(kind=8),allocatable :: heat_com(:,:)  
     real(kind=8),allocatable :: cool_com_prime(:,:)  
     real(kind=8),allocatable :: heat_com_prime(:,:)  
     real(kind=8),allocatable :: metal_prime(:,:)  
     real(kind=8),allocatable :: cool_prime(:,:)  
     real(kind=8),allocatable :: heat_prime(:,:)  
     real(kind=8),allocatable :: mu(:,:)  
     real(kind=8),allocatable :: spec(:,:,:)  ! last dimension (6) is n_e, n_HI, n_HII, n_HeI, n_HeII, n_HeIII
  end type cooling_table
  type(cooling_table) :: cooling
  type cool_interp
     integer(kind=4)  :: n_nH
     real(kind=8)     :: nH_start,nH_step
     integer(kind=4)  :: n_T2
     real(kind=8)     :: T2_start,T2_step
  end type cool_interp
  type(cool_interp)  :: cool_int
  logical            :: cooling_is_read = .False. 
  ! ----------------------------------------------------------------------------------

  ! particle-related stuff -----------------------------------------------------------
  character(30) :: ParticleFields(20)  ! array of particle fields (e.g. (/'pos','vel','mass','iord','level'/) for a DM-only run)
  ! conformal time things
  integer(kind=4),parameter             :: n_frw = 1000
  real(KIND=8),dimension(:),allocatable :: aexp_frw,hexp_frw,tau_frw,t_frw
  ! ----------------------------------------------------------------------------------

   ! A few things from module_domain, to avoid having additional files ---------------
  type shell
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: r_inbound,r_outbound
  end type shell

  type cube
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: size            ! convention: size is the full size of the cube, corners are x+-size/2
  end type cube

  type sphere
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: radius
  end type sphere

  type slab                                   ! infinite slab in xy direction (make use of periodic boundaries)
     real(kind=8)              :: zc          ! position of the slab in the z direction
     real(kind=8)              :: thickness   ! thickness of the slab in the z direction 
  end type slab

  type domain
     character(10) :: type  ! one of the possible shapes ('shell', 'cube', 'sphere', 'slab')
     type(shell)   :: sh
     type(cube)    :: cu
     type(sphere)  :: sp
     type(slab)    :: sl
  end type domain
  ! ----------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [ramses] of the parameter file
  ! --------------------------------------------------------------------------
  ! ramses options (not guessable from outputs)
  character(20)            :: star_reading_method  = 'default'  ! "default" : like originally in rascas,  "teyssier" : like in Romain's github Ramses version, and changes minit,  "fromlist", from an external list
  character(200)           :: star_info_file       = '.'        ! Where the file with star lists is, in case star_reading_method = 'fromlist'
  real(kind=8)             :: t_sne_Myr            = 1d1        ! Time before supernovae in Romain's version of Ramses
  real(kind=8)             :: eta_sn               = 2d-1       ! In Teyssier version :  if(t > t_sne_myr), minit = m/(1 - eta_sn)
  logical                  :: self_shielding       = .true.     ! if true, reproduce self-shielding approx made in ramses to compute nHI. 
  logical                  :: ramses_rt            = .false.    ! if true, read ramses-RT output and compute nHI and T accordingly.
  logical                  :: read_rt_variables    = .false.    ! if true, read RT variables (e.g. to compute heating terms)
  logical                  :: use_initial_mass     = .false.    ! if true, use initial masses of star particles instead of mass at output time
  logical                  :: cosmo                = .true.     ! if false, assume idealised simulation
  logical                  :: use_proper_time      = .false.    ! if true, use proper time instead of conformal time for cosmo runs. 
  logical                  :: QuadHilbert       = .false.  ! if true, do not use hilbert indexes for now ... 

  ! miscelaneous
  logical                  :: verbose        = .false. ! display some run-time info on this module
  ! RT variable indices
  integer(kind=4) :: itemp  = 5 ! index of thermal pressure
  integer(kind=4) :: imetal = 6 ! index of metallicity 
  integer(kind=4) :: ihii   = 7 ! index of HII fraction 
  integer(kind=4) :: iheii  = 8 ! index of HeII fraction 
  integer(kind=4) :: iheiii = 9 ! index of HeIII fraction
  integer(kind=4) :: iGamma1=10 ! index of first photon bin
  ! --------------------------------------------------------------------------

  
  public :: compute_csn_in_box, compute_csn_in_domain, read_ramses_params, print_ramses_params, get_ncpu, get_nSEDgroups
  public :: ramses_get_T_nhi_cgs, ramses_get_metallicity, ramses_get_nh_cgs, ramses_get_fractions, ramses_get_flux, ramses_read_stars_in_domain
  public :: domain_constructor_from_scratch, get_param_real
  
  !==================================================================================
contains

  ! ----------------
  ! public functions 
  ! ----------------


  subroutine ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nhi)

    implicit none 

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nhi(nleaf), temp(nleaf)
    real(kind=8),allocatable    :: boost(:)
    real(kind=8)                :: xhi
    integer(kind=4) :: ihx,ihy,i
    real(kind=8)    :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4) :: if1,if2,jf1,jf2

    real(kind=8),allocatable,dimension(:)    :: mu

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if(ramses_rt)then
       ! ramses RT
       allocate(mu(1:nleaf))
       nhi  = ramses_var(1,:) * dp_scale_nh  * (1.d0 - ramses_var(ihii,:))   ! nb of H atoms per cm^3
       mu   = XH * (1.d0+ramses_var(ihii,:)) + 0.25d0*(1.d0-XH)*(1.d0 + ramses_var(iheii,:) + 2.d0*ramses_var(iheiii,:)) ! assumes no metals
       mu   = 1.0d0 / mu   
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2            ! T/mu [ K ]
       temp = temp * mu                                                      ! This is now T (in K) with no bloody mu ... 
       deallocate(mu)
    else
       ! ramses standard

       if (.not. cooling_is_read) then
          call read_cooling(repository,snapnum)
          cooling_is_read = .True.
       end if

       nhi  = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
       temp = ramses_var(itemp,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]

       allocate(boost(nleaf))
       if (self_shielding) then
          do i=1,nleaf
             boost(i)=MAX(exp(-nhi(i)/0.01),1.0D-20) ! same as hard-coded in RAMSES. 
          end do
       else
          boost = 1.0d0
       end if

    
       ! compute the ionization state and temperature using the 'cooling' tables
       do i = 1, nleaf
          xx  = min(max(log10(nhi(i)/boost(i)),cooling%nh(1)),cooling%nh(cooling%n11))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          ! neutral H density 
          f = dxx1 * dyy1 * cooling%spec(if2,jf2,2) + dxx2 * dyy1 * cooling%spec(if1,jf2,2) &
               & + dxx1 * dyy2 * cooling%spec(if2,jf1,2) + dxx2 * dyy2 * cooling%spec(if1,jf1,2)
          xhi = 10.0d0**(f-xx)  ! this is xHI = nHI/nH where the n's include the boost. 
          nhi(i) = xhi * nhi(i) ! nHI (cm^-3) (boost-free)
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
       end do

       deallocate(boost)

    endif
    
    return

  end subroutine ramses_get_T_nhi_cgs

  

  subroutine ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nh)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nh(nleaf)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    nh = ramses_var(1,:) * dp_scale_nh ! [ H / cm^3 ]

    return

  end subroutine ramses_get_nh_cgs


  subroutine ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)

    implicit none

    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: metallicity(nleaf)

    if (nvar < 6) then
       print*,'No metals !!! '
       stop
    end if
    metallicity = ramses_var(imetal,:) 

    return

  end subroutine ramses_get_metallicity


  subroutine ramses_get_fractions(nleaf,nvar,ramses_var,x)

    implicit none

    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: x(3,nleaf)

    x(1,:) = ramses_var(ihii,:)
    x(2,:) = ramses_var(iheii,:)
    x(3,:) = ramses_var(iheiii,:)

    return

  end subroutine ramses_get_fractions
  

  subroutine ramses_get_flux(repository,snapnum,nleaf,nvar,nGroups,ramses_var,flux)

    implicit none
    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar, nGroups
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: flux(nGroups,nleaf)
    integer(kind=4)             :: i
    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if
    do i=1,nGroups
       flux(i,:) = ramses_var(iGamma1+4*(i-1),:)*dp_scale_pf
    end do

    return

  end subroutine ramses_get_flux


  function ramses_get_box_size_cm(repository,snapnum)

    implicit none

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    real(kind=8)                :: ramses_get_box_size_cm 

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    ramses_get_box_size_cm = get_param_real(repository,snapnum,'boxlen') * dp_scale_l  ! [ cm ] 

    return

  end function ramses_get_box_size_cm


  function get_ncpu(repository,snapnum)

    implicit none

    integer(kind=4)            :: get_ncpu
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    get_ncpu = nint(get_param_real(repository,snapnum,'ncpu'))

    return
  end function get_ncpu


  function get_nSEDgroups(repository,snapnum)

    implicit none

    integer(kind=4)            :: get_nSEDgroups
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    get_nSEDgroups = nint(get_param_rt_real(repository,snapnum,'nGroups'))

    return
  end function get_nSEDgroups


  function get_param_real(repository,snapnum,param,default_value,rt_info)

    implicit none 

    real(kind=8)                    :: get_param_real
    character(512),intent(in)       :: repository
    integer(kind=4),intent(in)      :: snapnum
    character(*),intent(in)         :: param
    real(kind=8),optional,intent(in):: default_value
    logical,optional,intent(in)     :: rt_info
    logical(kind=4)                 :: not_ok
    character(512)                  :: filename, infofile
    character(512)                  :: line,name,value
    integer(kind=4)                 :: i
    integer(kind=4),parameter       :: param_unit = 13

    not_ok = .true.
    infofile='/info_'
    if(present(rt_info)) then ! read info_rt file
       if(rt_info)  infofile='/info_rt_'
    endif
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,trim(infofile),snapnum,'.txt'
    open(unit=param_unit,file=filename,status='old',form='formatted')
    do 
       read(param_unit,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))

       if (trim(name) .eq. trim(param)) then 
          read(value,*) get_param_real
          not_ok = .false.
          exit
       end if

    end do
2   close (param_unit)  

    if (not_ok) then 
       write(6,*) '--> parameter not found in infoxxx.txt : ',trim(param)
       if(present(default_value)) then
          get_param_real = default_value
       else
          stop
       endif
    end if

    return

  end function get_param_real


  function get_param_rt_real(repository,snapnum,param)

    implicit none

    real(kind=8)               :: get_param_rt_real
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum
    character(*),intent(in)    :: param
    logical(kind=4)            :: not_ok
    character(1000)            :: nomfich
    character(512)             :: line,name,value
    integer(kind=4)            :: i
    integer(kind=4),parameter  :: param_unit = 13

    not_ok = .true.
    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,'/info_rt_',snapnum,'.txt'
    open(unit=param_unit,file=nomfich,status='old',form='formatted')
    do 
       read(param_unit,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))

       if (trim(name) .eq. trim(param)) then 
          read(value,*) get_param_rt_real
          not_ok = .false.
       end if

    end do
2   close (param_unit)  

    if (not_ok) then 
       write(6,*) '> parameter not found in infoxxx.txt :',trim(param)
       stop
    end if

    return

  end function get_param_rt_real


  subroutine read_conversion_scales(repository,snapnum)

    implicit none 

    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    ! set global variables 
    dp_scale_l    = get_param_real(repository,snapnum,'unit_l')
    dp_scale_d    = get_param_real(repository,snapnum,'unit_d')
    dp_scale_t    = get_param_real(repository,snapnum,'unit_t')
    dp_scale_pf   = get_param_rt_real(repository,snapnum,'unit_pf')
    dp_scale_nH   = XH/mp * dp_scale_d      ! convert mass density (code units) to numerical density of H atoms [/cm3]
    dp_scale_nHe  = (1d0-XH)/mp * dp_scale_d ! mass dens to He/cm3
    dp_scale_v    = dp_scale_l/dp_scale_t   ! -> converts velocities into cm/s
    dp_scale_T2   = mp/kB * dp_scale_v**2   ! -> converts P/rho to T/mu, in K
    dp_scale_zsun = 1.d0/0.0127     
    dp_scale_m    = dp_scale_d * dp_scale_l**3 ! convert mass in code units to cgs. 

    return

  end subroutine read_conversion_scales

  subroutine read_cosmo_params(repository,snapnum,omega_0,lambda_0,little_h)

    implicit none 

    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum
    real(kind=8),intent(out)   :: omega_0,lambda_0,little_h

    omega_0  = get_param_real(repository,snapnum,'omega_m')
    lambda_0 = get_param_real(repository,snapnum,'omega_l')
    little_h = get_param_real(repository,snapnum,'H0') / 100.0d0

    return

  end subroutine read_cosmo_params

  subroutine get_fields_from_header(dir,ts,nfields)

    implicit none

    character(1000),intent(in)  :: dir
    integer(kind=4),intent(in)  :: ts
    integer(kind=4),intent(out) :: nfields
    character(2000)             :: nomfich,line
    integer(kind=4) :: i

    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(dir),'/output_',ts,'/header_',ts,'.txt'
    open(unit=50,file=nomfich,status='old',action='read',form='formatted')
    read(50,*) ! total nb of particles
    read(50,*)
    read(50,*) ! nb of DM particles
    read(50,*)
    read(50,*) ! nb of star particles
    read(50,*)
    read(50,*) ! nb of sinks
    read(50,*)
    read(50,*) ! Field list
    read(50,'(a)') line
    close(50)

    ! parse the Field list ...
    nfields = 0
    do
       i    = scan(line,' ') ! find a blank
       nfields = nfields + 1
       ParticleFields(nfields) = trim(adjustl(line(:i)))
       line = trim(adjustl(line(i:)))
       if (len_trim(line) == 0) exit
    end do

    return

  end subroutine get_fields_from_header

  !*****************************************************************************************************************

  subroutine read_cooling(repository,snapnum)

    implicit none

    character(1000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    character(1024)            :: filename
    integer(kind=4)            :: n1,n2

    ! initialize cooling variables
    call clear_cooling

    ! read cooling variables from cooling.out file
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,"/cooling_",snapnum,".out"
    open(unit=44,file=filename,form='unformatted')
    read(44) n1,n2
    cooling%n11 = n1
    cooling%n22 = n2
    allocate(cooling%nH(n1),cooling%T2(n2))
    allocate(cooling%cool_com(n1,n2),cooling%heat_com(n1,n2),cooling%cool_com_prime(n1,n2),cooling%heat_com_prime(n1,n2))
    allocate(cooling%cool(n1,n2),cooling%heat(n1,n2),cooling%mu(n1,n2))
    allocate(cooling%cool_prime(n1,n2),cooling%heat_prime(n1,n2),cooling%metal_prime(n1,n2))
    allocate(cooling%metal(n1,n2),cooling%spec(n1,n2,6))
    read(44)cooling%nH
    read(44)cooling%T2
    read(44)cooling%cool
    read(44)cooling%heat
    read(44)cooling%cool_com
    read(44)cooling%heat_com
    read(44)cooling%metal
    read(44)cooling%cool_prime
    read(44)cooling%heat_prime
    read(44)cooling%cool_com_prime
    read(44)cooling%heat_com_prime
    read(44)cooling%metal_prime
    read(44)cooling%mu
    read(44)cooling%spec
    close(44)

    ! define useful quantities for interpolation 
    cool_int%n_nh     = n1
    cool_int%nh_start = minval(cooling%nh)
    cool_int%nh_step  = cooling%nh(2) - cooling%nh(1)
    cool_int%n_t2     = n2
    cool_int%t2_start = minval(cooling%t2)
    cool_int%t2_step  = cooling%t2(2) - cooling%t2(1)

    return

  end subroutine read_cooling

  !*****************************************************************************************************************

  subroutine clear_cooling

    implicit none

    if (cooling%n11 > 0 .or. cooling%n22 > 0) then 
       cooling%n11 = 0
       cooling%n22 = 0
       deallocate(cooling%nH,cooling%T2,cooling%metal)
       deallocate(cooling%heat_com,cooling%cool_com,cooling%heat_com_prime,cooling%cool_com_prime)
       deallocate(cooling%cool,cooling%heat,cooling%metal_prime,cooling%cool_prime)
       deallocate(cooling%heat_prime,cooling%mu,cooling%spec)
    end if

    cool_int%n_nh = 0
    cool_int%nh_start = 0.0d0
    cool_int%nh_step  = 0.0d0
    cool_int%n_t2 = 0
    cool_int%t2_start = 0.0d0
    cool_int%t2_step  = 0.0d0

    return

  end subroutine clear_cooling



  !==================================================================================
  ! STARS utilities 

  subroutine ramses_read_stars_in_domain(repository,snapnum,selection_domain,star_pos,star_age,star_mass,star_vel,star_met)

    !$ use OMP_LIB

    implicit none

    character(1000),intent(in)             :: repository
    integer(kind=4),intent(in)             :: snapnum
    type(domain),intent(in)                :: selection_domain
    real(kind=8),allocatable,intent(inout) :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:)
    real(kind=8),allocatable               :: star_pos_all(:,:),star_age_all(:),star_mass_all(:),star_vel_all(:,:),star_met_all(:)
    integer(kind=4)                        :: nstars
    real(kind=8)                           :: omega_0,lambda_0,little_h,omega_k,H0
    real(kind=8)                           :: aexp,stime,time_cu,boxsize
    integer(kind=4)                        :: ncpu,ilast,icpu,npart,i,ifield,nfields
    character(1000)                        :: filename
    integer(kind=4),allocatable            :: id(:)
    real(kind=8),allocatable               :: age(:),m(:),x(:,:),v(:,:),mets(:),imass(:),skipy(:)
    real(kind=8)                           :: temp(3)
    integer(kind=4)                        :: rank, iunit, ilast_all

    if(star_reading_method == 'fromlist') then

       open(unit=20, file=star_info_file, status='old', form='unformatted')

       read(20) nstars
       allocate(star_pos(3,nstars), star_age(nstars), star_mass(nstars), star_vel(3,nstars), star_met(nstars))

       read(20) star_pos(1,:) !code units
       read(20) star_pos(2,:)
       read(20) star_pos(3,:)
       read(20) star_age(:)   !Myr
       if(use_initial_mass) then
          read(20) star_mass(:)  !g
          read(20)
       else
          read(20)
          read(20) star_mass(:)  !g
       end if
       read(20) star_vel(1,:) ! cm/s
       read(20) star_vel(2,:)
       read(20) star_vel(3,:)
       read(20) star_met(:)   !M fraction of metals
       close(20)

       ! select particles in domain
       allocate(x(3,nstars), age(nstars), m(nstars), v(3,nstars), mets(nstars))
       ilast = 1
       do i=1,nstars
          temp = star_pos(:,i)
          if (domain_contains_point(temp,selection_domain)) then ! it is inside the domain
             x(:,ilast) = star_pos(:,i)
             age(ilast) = star_age(i)
             m(ilast) = star_mass(i)
             v(:,ilast) = star_vel(:,i)
             mets(ilast) = star_met(i)

             ilast = ilast+1
          end if
       end do

       !Change size of star_xxx
       deallocate(star_pos, star_age, star_mass, star_vel, star_met)
       nstars = ilast - 1
       allocate(star_pos(3,nstars), star_age(nstars), star_mass(nstars), star_vel(3,nstars), star_met(nstars))

       do i=1,3
          star_pos(i,:) = x(i,1:nstars)
          star_vel(i,:) = v(i,1:nstars)
       end do
       star_age(:) = age(1:nstars)
       star_mass(:) = m(1:nstars)
       star_met(:) = mets(1:nstars)

       deallocate(x, age, m, v, mets)

    else

       ! get cosmological parameters to convert conformal time into ages
       call read_cosmo_params(repository,snapnum,omega_0,lambda_0,little_h)
       omega_k = 0.0d0
       h0      = little_h * 100.0d0
       call ct_init_cosmo(omega_0,lambda_0,omega_k,h0)
       ! compute cosmic time of simulation output (Myr)
       aexp  = get_param_real(repository,snapnum,'aexp') ! exp. factor of output
       stime = ct_aexp2time(aexp) ! cosmic time

       ! read units
       if (.not. conversion_scales_are_known) then
          call read_conversion_scales(repository,snapnum)
          conversion_scales_are_known = .True.
       end if

       if(.not.cosmo)then
          ! read time
          time_cu = get_param_real(repository,snapnum,'time') ! code unit
          write(*,*)'Time simu [Myr] =',time_cu, time_cu*dp_scale_t/(365.*24.*3600.*1d6)
          boxsize = get_param_real(repository,snapnum,'boxlen') !!!* dp_scale_l  ! [ cm ]
          write(*,*)'boxlen =',boxsize
       endif

       ! read stars
       nstars = get_tot_nstars(repository,snapnum)
       if (nstars == 0) then
          write(*,*) 'ERROR : no star particles in output '
          stop
       end if
       allocate(star_pos(3,nstars),star_age(nstars),star_mass(nstars),star_vel(3,nstars),star_met(nstars))
       allocate(star_pos_all(3,nstars),star_age_all(nstars),star_mass_all(nstars),star_vel_all(3,nstars),star_met_all(nstars))
       ! get list of particle fields in outputs
       call get_fields_from_header(repository,snapnum,nfields)
       ncpu  = get_ncpu(repository,snapnum)
       ilast_all = 1

       !$OMP PARALLEL &
       !$OMP DEFAULT(private) &
       !$OMP SHARED(ncpu, repository, snapnum, ParticleFields, nfields, selection_domain) &
       !$OMP SHARED(h0, stime, dp_scale_t, dp_scale_m, dp_scale_v, boxsize, time_cu, aexp, cosmo, use_initial_mass, use_proper_time) &
       !$OMP SHARED(ilast_all, star_pos_all, star_age_all, star_vel_all, star_mass_all, star_met_all, star_reading_method)
       !$OMP DO
       do icpu = 1, ncpu
          rank = 1
          !$ rank = OMP_GET_THREAD_NUM()
          iunit=10+rank*2
          write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository), '/output_', snapnum, '/part_', snapnum, '.out', icpu
          open(unit=iunit,file=filename,status='old',form='unformatted')
          
          if(star_reading_method == 'Teyssier') then

             read(iunit) !ncpu2
             read(iunit) !ndim2
             read(iunit) npart
             read(iunit) !localseed
             read(iunit) !nstar_tot
             read(iunit) !mstar_tot
             read(iunit) !mstar_lost
             read(iunit) !nsink

             allocate(age(1:npart))
             allocate(x(1:npart,1:ndim),m(npart),imass(npart))
             allocate(id(1:npart))
             allocate(mets(1:npart))
             allocate(v(1:npart,1:ndim))
             allocate(skipy(1:npart))

             ! Read position
             do i=1,ndim
                read(iunit) x(1:npart,i)
             end do
             ! Read velocity
             do i=1,ndim
                read(iunit) v(1:npart,i)
             end do
             ! Read mass
             read(iunit) m(1:npart)
             ! Read identity
             read(iunit) id(1:npart)
             ! Read level
             read(iunit)
             ! Read family
             read(iunit)
             ! Read tag
             read(iunit)
             ! Read birth epoch
             read(iunit) age(1:npart)
             ! Read metallicity
             read(iunit) mets(1:npart)

          else  !default,  like sphinx

             read(iunit)
             read(iunit)
             read(iunit)npart
             read(iunit)
             read(iunit)
             read(iunit)
             read(iunit)
             read(iunit)
             allocate(age(1:npart))
             allocate(x(1:npart,1:ndim),m(npart),imass(npart))
             allocate(id(1:npart))
             allocate(mets(1:npart))
             allocate(v(1:npart,1:ndim))
             allocate(skipy(1:npart))
             do ifield = 1,nfields
                select case(trim(ParticleFields(ifield)))
                case('pos')
                   do i = 1,ndim
                      read(iunit) x(1:npart,i)
                   end do
                case('vel')
                   do i = 1,ndim
                      read(iunit) v(1:npart,i)
                   end do
                case('mass')
                   read(iunit) m(1:npart)
                case('iord')
                   read(iunit) id(1:npart)
                case('level')
                   read(iunit)
                case('family')
                   read(iunit)
                case('tag')
                   read(iunit)
                case('tform')
                   read(iunit) age(1:npart)
                case('metal')
                   read(iunit) mets(1:npart)
                case('metalII')
                   read(iunit)
                case('imass')
                   read(iunit) imass(1:npart)
                case default
                   ! Note: we presume here that the unknown field is an 1d array of size 1:npart
                   read(iunit) skipy(1:npart)
                   print*,'Error, Field unknown: ',trim(ParticleFields(ifield))
                end select
             end do

          end if
          close(iunit)

          if(.not.cosmo)then
             x=x/boxsize
          endif

          ! save star particles within selection region
          ilast = 0
          do i = 1,npart
             if (age(i).ne.0.0d0) then ! This is a star
                temp(:) = x(i,:)
                if (domain_contains_point(temp,selection_domain)) then ! it is inside the domain
                   ilast = ilast + 1
                   if(cosmo)then
                      if (use_proper_time) then
                         star_age(ilast) = (stime - ct_proptime2time(age(i),h0))*1.d-6 ! Myr
                      else
                         ! Convert from conformal time to age in Myr
                         star_age(ilast) = (stime - ct_conftime2time(age(i)))*1.d-6 ! Myr
                      end if
                   else
                      ! convert from tborn to age in Myr
                      star_age(ilast)   = max(0.d0, (time_cu - age(i)) * dp_scale_t / (365.d0*24.d0*3600.d0*1.d6))
                   endif
                   if (use_initial_mass) then
                      !From Romain Teyssier version : have to divide mass by 1-eta_sn
                      if(star_reading_method == 'Teyssier') then
                         star_mass(ilast) = m(i) * dp_scale_m ! [g]
                         if(star_age(ilast) > t_sne_Myr) star_mass(ilast) = m(i) * dp_scale_m / (1d0 - eta_sn)
                      !Normal way, like initially in Rascas
                      else
                         star_mass(ilast) = imass(i) * dp_scale_m ! [g]
                      end if

                   else
                      star_mass(ilast) = m(i)     * dp_scale_m ! [g]
                   end if
                   star_pos(:,ilast) = x(i,:)              ! [code units]
                   star_vel(:,ilast) = v(i,:) * dp_scale_v ! [cm/s]
                   star_met(ilast) = mets(i)
                end if
             end if
          end do

          deallocate(age,m,x,id,mets,v,skipy,imass)

          !$OMP CRITICAL
          if(ilast .gt. 0) then
             star_age_all(ilast_all:ilast_all+ilast-1) = star_age(1:ilast)
             star_mass_all(ilast_all:ilast_all+ilast-1) = star_mass(1:ilast)
             star_pos_all(1:3,ilast_all:ilast_all+ilast-1) = star_pos(1:3,1:ilast)
             star_vel_all(1:3,ilast_all:ilast_all+ilast-1) = star_vel(1:3,1:ilast)
             star_met_all(ilast_all:ilast_all+ilast-1) = star_met(1:ilast)
          endif
          ilast_all = ilast_all + ilast
          !$OMP END CRITICAL
       end do
       !$OMP END PARALLEL

       deallocate(star_age,star_pos,star_vel,star_met,star_mass)

       ! resize star arrays
       nstars = ilast_all-1
       print*,'Nstars in domain =',nstars
       ! ages
       allocate(star_age(nstars))
       star_age = star_age_all(1:nstars)
       deallocate(star_age_all)
       ! masses
       allocate(star_mass(nstars))
       star_mass = star_mass_all(1:nstars)
       deallocate(star_mass_all)
       ! positions
       allocate(star_pos(3,nstars))
       do i = 1,nstars
          star_pos(:,i) = star_pos_all(:,i)
       end do
       deallocate(star_pos_all)
       ! velocities
       allocate(star_vel(3,nstars))
       do i = 1,nstars
          star_vel(:,i) = star_vel_all(:,i)
       end do
       deallocate(star_vel_all)
       ! metals
       allocate(star_met(nstars))
       star_met = star_met_all(1:nstars)
       deallocate(star_met_all)

    end if

    return
  end subroutine ramses_read_stars_in_domain
  

  function get_tot_nstars(dir,ts)

    implicit none 

    integer(kind=4),intent(in) :: ts
    character(1000),intent(in) :: dir
    character(2000)            :: nomfich, useless
    integer(kind=4)            :: get_tot_nstars

    get_tot_nstars = 0
    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(dir),'/output_',ts,'/header_',ts,'.txt'
    open(unit=50,file=nomfich,status='old',action='read',form='formatted')
    if(star_reading_method == 'Teyssier') then
       read(50,*) ! total nb of particles
       read(50,*)
       read(50,*) ! nb of DM particles
       read(50,*)
       read(50,*) ! nb of star particles 
       read(50,*) ! nb of DM particles
       read(50,*)
       read(50,*) ! nb of star particles 
       read(50,*) useless, get_tot_nstars
    else
       read(50,*) ! total nb of particles
       read(50,*)
       read(50,*) ! nb of DM particles
       read(50,*)
       read(50,*) ! nb of star particles 
       read(50,*) get_tot_nstars
    end if

    close(50)

    return

  end function get_tot_nstars


  !*************************************************************************
  subroutine compute_csn_in_domain(repository, snapnum, dom, n_elements, elements, nIons, csn)

    use module_spectra

    character(2000),intent(in)           :: repository 
    type(domain),intent(in)              :: dom
    integer(kind=4),intent(in)           :: n_elements, elements(n_elements), nIons(n_elements), snapnum
    real(kind=8),allocatable             :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:), star_L(:), star_L_tot(:), csn_help(:,:)
    integer(kind=4)                      :: nstars, nSEDgroups, i,j,k, nProp
    real(kind=8),allocatable,intent(out) :: csn(:,:)
    character(10)                        :: type
    integer(kind=4),allocatable          :: prop(:)

    prop = get_csn_indices(n_elements,elements,nIons)
    nProp = size(prop)

    nSEDgroups = get_nOptBins()
    allocate(csn(nSEDgroups,nProp)) ; csn = 0d0

    call ramses_read_stars_in_domain(repository,snapnum,dom,star_pos,star_age,star_mass,star_vel,star_met) ; deallocate(star_vel)
    nstars = size(star_age)
    print*, 'number of stars : ', nstars

    !initialize SED properties
    call init_SED_table()

    allocate(star_L(nSEDgroups), star_L_tot(nSEDgroups), csn_help(nSEDgroups,nProp)) ; star_L_tot = 0d0

    do i=1,nstars
       call inp_sed_table(star_age(i)/1d3, star_met(i), 1, .false., star_L(:))    !Third variable :  1 for L[#photons/s],  3 for mean energy in bin,  2+2*Iion for mean cross-section,  Iion: 1:HI,2:HeI,3:HeII,4:CI, etc, see module_spectra.f90
       star_L(:) = star_mass(i)*star_L(:)
       do j=1,nProp
          call inp_sed_table(star_age(i)/1d3, star_met(i), 2+prop(j)*2, .false., csn_help(:,j))
       end do

       do j=1,nSEDgroups
          csn(j,:) = csn(j,:) + star_L(j)*csn_help(j,:)
          star_L_tot(j) = star_L_tot(j) + star_L(j)
       end do
    end do

    call deallocate_table()

    print*, 'csn (nSEDgroups * nIons)'
    do j=1,nSEDgroups
       csn(j,:) = csn(j,:)/star_L_tot(j)
       print*, csn(j,:)
    end do
    deallocate(star_age, star_met, star_mass, csn_help)


  end subroutine compute_csn_in_domain
  !*************************************************************************

  !*************************************************************************
  subroutine compute_csn_in_box(repository, snapnum, n_elements, elements, nIons, csn)

    use module_spectra

    character(2000),intent(in)           :: repository 
    type(domain)                         :: dom
    integer(kind=4),intent(in)           :: n_elements, elements(n_elements), nIons(n_elements), snapnum
    real(kind=8),allocatable             :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:), star_L(:), star_L_tot(:), csn_help(:,:)
    integer(kind=4)                      :: nstars, nSEDgroups, i,j,k, nProp
    real(kind=8),allocatable,intent(out) :: csn(:,:)
    character(10)                        :: type
    integer(kind=4),allocatable          :: prop(:)


    prop = get_csn_indices(n_elements,elements,nIons)
    nProp = size(prop)

    nSEDgroups = get_nOptBins()                              !The groups of which the routine computes csn are defined in the section [spectra] of the parameters,  not in the simulation
    allocate(csn(nSEDgroups,nProp)) ; csn = 0d0

    type = 'sphere'
    call domain_constructor_from_scratch(dom, type, 5d-1, 5d-1, 5d-1, 1d3)

    call ramses_read_stars_in_domain(repository,snapnum,dom,star_pos,star_age,star_mass,star_vel,star_met) ; deallocate(star_vel)
    if(star_mass(1) <= 0d0) then
       print*, 'Problem with star mass, try using use_initial_mass = False  in the [ramses] section of the parameter file'
       stop
    end if
    if(star_age(1) <= 0d0) then
       print*, 'Problem with star age, try using use_proper_time = True  in the [ramses] section of the parameter file'
       stop
    end if
    
    nstars = size(star_age)
    print*, 'number of stars : ', nstars

    !initialize SED properties
    call init_SED_table()

    allocate(star_L(nSEDgroups), star_L_tot(nSEDgroups), csn_help(nSEDgroups,nProp)) ; star_L_tot = 0d0

    do i=1,nstars
       call inp_sed_table(star_age(i)/1d3, star_met(i), 1, .false., star_L(:))    !Third variable :  1 for L[#photons/s],  3 for mean energy in bin,  2+2*Iion for mean cross-section,  Iion: 1:HI,2:HeI,3:HeII,4:CI, etc, see module_spectra.f90
       star_L(:) = star_mass(i)*star_L(:)
       do j=1,nProp
          call inp_sed_table(star_age(i)/1d3, star_met(i), 2+prop(j)*2, .false., csn_help(:,j))
       end do

       do j=1,nSEDgroups
          csn(j,:) = csn(j,:) + star_L(j)*csn_help(j,:)
          star_L_tot(j) = star_L_tot(j) + star_L(j)
       end do
    end do

    call deallocate_table()

    print*, 'csn (nSEDgroups * nIons)'
    !print*, 'star_L_tot, ', star_L_tot
    do j=1,nSEDgroups
       csn(j,:) = csn(j,:)/star_L_tot(j)
       print*, csn(j,:)
    end do
    deallocate(star_age, star_met, star_mass, csn_help)


  end subroutine compute_csn_in_box
  !*************************************************************************

  

  ! conformal time utils :
    function ct_conftime2time(tau)

    ! return look-back time in yr
    
    implicit none 
    
    real(kind=8),intent(in) :: tau
    real(kind=8)            :: ct_conftime2time
    integer(kind=4)         :: i
    
    
    ! locate bracketing conf. times
    i = 1
    do while(tau_frw(i) > tau .and. i < n_frw)
       i = i + 1
    end do
    ! Interploate time
    ct_conftime2time = t_frw(i) * (tau-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
         & t_frw(i-1)        * (tau-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
    
    return

  end function ct_conftime2time
  
  function ct_proptime2time(tau,h0)
    ! return look-back time in yr
    
    implicit none 
    real(kind=8),intent(in) :: tau,h0
    real(kind=8)            :: ct_proptime2time
    
    ct_proptime2time = tau / (h0 / 3.08d19) / (365.25*24.*3600.)
    return
  end function ct_proptime2time
  
  
  function ct_aexp2time(aexp)

    ! return look-back time in yr

    implicit none

    real(kind=8),intent(in) :: aexp
    real(kind=8)            :: ct_aexp2time
    integer(kind=4)         :: i

    ! find bracketting aexp's 
     i = 1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i = i + 1
     end do
     ! Interploate time
     ct_aexp2time = t_frw(i) * (aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)    * (aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))

    return
    
  end function ct_aexp2time

  
  subroutine ct_init_cosmo(omega_m,omega_l,omega_k,h0)
    
    ! h0 is in km/s/Mpc

    implicit none 
    real(kind=8),intent(in) :: omega_m,omega_l,omega_k,h0
    real(kind=8)            :: time_tot

    allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
    allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
    call ct_friedman(omega_m,omega_l,omega_k,1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
    ! convert time to yr
    t_frw = t_frw / (h0 / 3.08d19) / (365.25*24.*3600.)

    return
    
  end subroutine ct_init_cosmo



  subroutine ct_friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
       & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

    implicit none
    integer::ntable
    real(kind=8)::O_mat_0, O_vac_0, O_k_0
    real(kind=8)::alpha,axp_min,age_tot
    real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
    ! ######################################################!
    ! This subroutine assumes that axp = 1 at z = 0 (today) !
    ! and that t and tau = 0 at z = 0 (today).              !
    ! axp is the expansion factor, hexp the Hubble constant !
    ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
    ! time, and t the look-back time, both in unit of 1/H0. !
    ! alpha is the required accuracy and axp_min is the     !
    ! starting expansion factor of the look-up table.       !
    ! ntable is the required size of the look-up table.     !
    ! ######################################################!
    real(kind=8)::axp_tau, axp_t
    real(kind=8)::axp_tau_pre, axp_t_pre
    real(kind=8)::dtau,dt
    real(kind=8)::tau,t
    integer::nstep,nout,nskip

    !  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
    !     write(*,*)'Error: non-physical cosmological constants'
    !     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
    !     write(*,*)'The sum must be equal to 1.0, but '
    !     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
    !     stop
    !  end if

    axp_tau = 1.0D0
    axp_t = 1.0D0
    tau = 0.0D0
    t = 0.0D0
    nstep = 0

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

    end do

    age_tot=-t
!!$    write(*,666)-t
!!$666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

    nskip=nstep/ntable

    axp_t = 1.d0
    t = 0.d0
    axp_tau = 1.d0
    tau = 0.d0
    nstep = 0
    nout=0
    t_out(nout)=t
    tau_out(nout)=tau
    axp_out(nout)=axp_tau
    hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

       if(mod(nstep,nskip)==0)then
          nout=nout+1
          t_out(nout)=t
          tau_out(nout)=tau
          axp_out(nout)=axp_tau
          hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
       end if

    end do
    t_out(ntable)=t
    tau_out(ntable)=tau
    axp_out(ntable)=axp_tau
    hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  contains
    function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
      implicit none 
      real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
      dadtau = axp_tau*axp_tau*axp_tau *  &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
           &     O_k_0   * axp_tau )
      dadtau = sqrt(dadtau)
      return
    end function dadtau
    
    function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
      implicit none
      real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
      dadt   = (1.0D0/axp_t)* &
           &   ( O_mat_0 + &
           &     O_vac_0 * axp_t*axp_t*axp_t + &
           &     O_k_0   * axp_t )
      dadt = sqrt(dadt)
      return
    end function dadt
    
  end subroutine ct_friedman


  subroutine read_ramses_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:8) == '[ramses]') then
          section_present = .true.
          exit
       end if
    end do
    ! read section if present
    if (section_present) then 
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:1) == '[') exit ! next section starting... -> leave
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case('star_reading_method')
             read(value,*) star_reading_method
          case('star_info_file')
            write(star_info_file,'(a)') trim(value)
          case('t_sne_Myr')
             read(value,*) t_sne_Myr
          case('eta_sn')
             read(value,*) eta_sn
          case ('self_shielding')
             read(value,*) self_shielding
          case ('ramses_rt')
             read(value,*) ramses_rt
          case ('read_rt_variables')
             read(value,*) read_rt_variables
          case ('verbose')
             read(value,*) verbose
          case ('use_initial_mass')
             read(value,*) use_initial_mass
          case ('cosmo')
             read(value,*) cosmo
          case ('use_proper_time')
             read(value,*) use_proper_time
          case('itemp') ! index of thermal pressure
             read(value,*) itemp
          case('imetal')! index of metallicity  
             read(value,*) imetal
          case('ihii') ! index of HII fraction 
             read(value,*) ihii
          case ('iheii') ! index of HeII fraction 
             read(value,*) iheii
          case('iheiii') ! index of HeIII fraction 
             read(value,*) iheiii
          case('iGamma1') ! index of first bin of photons  
             read(value,*) iGamma1
          case('QuadHilbert') ! True if simulation was run with -DQUADHILBERT option  
             read(value,*) QuadHilbert
          end select
       end do
    end if
    close(10)
    return

  end subroutine read_ramses_params


  
  subroutine print_ramses_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then
       write(*,'(a)')            '--------------------------------------------------------------------------------'
       write(unit,'(a,a,a)') '[ramses]'
       write(unit,'(a,a)')      '  star_reading_method  = ', star_reading_method
       write(unit,'(a,a)')      '  star_info_file       = ', trim(star_info_file)
       write(unit,'(a,ES10.3)') '  t_sne_Myr       = ', t_sne_Myr
       write(unit,'(a,ES10.3)') '  eta_sn          = ', eta_sn
       write(unit,'(a,L1)')     '  self_shielding       = ', self_shielding
       write(unit,'(a,L1)')     '  ramses_rt            = ', ramses_rt
       write(unit,'(a,L1)')     '  read_rt_variables    = ', read_rt_variables
       write(unit,'(a,L1)')     '  use_initial_mass     = ', use_initial_mass
       write(unit,'(a,L1)')     '  cosmo                = ', cosmo
       write(unit,'(a,L1)')     '  use_proper_time      = ', use_proper_time
       write(unit,'(a,L1)')     '  QuadHilbert       = ',QuadHilbert
       write(unit,'(a,L1)')     '  verbose              = ', verbose
       write(unit,'(a,i2)')     '  itemp                = ', itemp
       write(unit,'(a,i2)')     '  imetal               = ', imetal
       write(unit,'(a,i2)')     '  ihii                 = ', ihii
       write(unit,'(a,i2)')     '  iheii                = ', iheii
       write(unit,'(a,i2)')     '  iheiii               = ', iheiii
       write(unit,'(a,i2)')     '  iGamma1              = ', iGamma1
       write(*,'(a)')            '--------------------------------------------------------------------------------'
    else
       write(*,'(a)')            '--------------------------------------------------------------------------------'
       write(*,'(a,a,a)') '[ramses]'
       write(*,'(a,a)')      '  star_reading_method  = ', star_reading_method
       write(*,'(a,a)')      '  star_info_file       = ', trim(star_info_file)
       write(*,'(a,ES10.3)') '  t_sne_Myr       = ', t_sne_Myr
       write(*,'(a,ES10.3)') '  eta_sn          = ', eta_sn
       write(*,'(a,L1)')     '  self_shielding       = ', self_shielding
       write(*,'(a,L1)')     '  ramses_rt            = ', ramses_rt
       write(*,'(a,L1)')     '  read_rt_variables    = ', read_rt_variables
       write(*,'(a,L1)')     '  use_initial_mass     = ', use_initial_mass
       write(*,'(a,L1)')     '  cosmo                = ', cosmo
       write(*,'(a,L1)')     '  use_proper_time      = ', use_proper_time
       write(*,'(a,L1)')     '  QuadHilbert          = ',QuadHilbert
       write(*,'(a,L1)')     '  verbose              = ', verbose
       write(*,'(a,i2)')     '  itemp                = ', itemp
       write(*,'(a,i2)')     '  imetal               = ', imetal
       write(*,'(a,i2)')     '  ihii                 = ', ihii
       write(*,'(a,i2)')     '  iheii                = ', iheii
       write(*,'(a,i2)')     '  iheiii               = ', iheiii
       write(*,'(a,i2)')     '  iGamma1              = ', iGamma1
       write(*,'(a)')            '--------------------------------------------------------------------------------'

    end if
    
    return
  end subroutine print_ramses_params


  !--------------------------------------------------------------------------------------------------
  ! domain constructors
  !--------------------------------------------------------------------------------------------------

  subroutine domain_constructor_from_scratch(dom,type,xc,yc,zc,r,r_inbound,r_outbound,size,thickness)

    implicit none
    character(10),intent(in)         :: type
    real(kind=8),intent(in),optional :: xc,yc,zc
    real(kind=8),intent(in),optional :: r                     ! parameters for sphere
    real(kind=8),intent(in),optional :: r_inbound,r_outbound  ! parameters for shell
    real(kind=8),intent(in),optional :: size                  ! parameters for cube
    real(kind=8),intent(in),optional :: thickness             ! parameters for slab
    type(domain),intent(out)         :: dom
    logical                          :: ok

    select case(trim(type))

    case('sphere')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a sphere domain are missing...'
          stop
       endif
       dom%type=type
       dom%sp%center(1)=xc
       dom%sp%center(2)=yc
       dom%sp%center(3)=zc
       dom%sp%radius=r

    case('shell')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r_inbound).and.present(r_outbound)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a shell domain are missing...'
          stop
       endif
       dom%type=type
       dom%sh%center(1)=xc
       dom%sh%center(2)=yc
       dom%sh%center(3)=zc
       dom%sh%r_inbound=r_inbound
       dom%sh%r_outbound=r_outbound

    case('cube')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(size)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a cube domain are missing...'
          stop
       endif
       dom%type=type
       dom%cu%center(1)=xc
       dom%cu%center(2)=yc
       dom%cu%center(3)=zc
       dom%cu%size=size


    case('slab')
       ! check if optional argument required for sphere are present
       ok = present(zc).and.present(thickness)
       if (.not.ok) then
          print *,'ERROR: arguments to construct a slab domain are missing...'
          stop
       endif
       dom%type=type
       dom%sl%zc=zc
       dom%sl%thickness=thickness

    case default
       print *,'ERROR: type not defined ',trim(type)
       stop
    end select

    return

  end subroutine domain_constructor_from_scratch

  function domain_contains_point(x,dom)
    ! -> returns T/F if point xyz is in domain dom.
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    logical                              :: domain_contains_point
    real(kind=8)                         :: rr,dx,dy,dz
    domain_contains_point=.false.
    select case(trim(dom%type))

    case('sphere')
       ! correct cell's position for periodic boundaries 
       dx = x(1)-dom%sp%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if

       dy = x(2)-dom%sp%center(2)
       if (dy > 0.5d0) then 
          dy = dy -1.0d0 
       else if (dy < -0.5d0) then 
          dy = dy + 1.0d0
       end if

       dz = x(3)-dom%sp%center(3)
       if (dz > 0.5d0) then 
          dz = dz -1.0d0 
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if

       !rr = (x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2
       rr = dx**2 + dy**2 + dz**2
       if(rr<dom%sp%radius*dom%sp%radius)domain_contains_point=.true.

    case('shell')
       ! correct positions for periodic boundaries 
       dx = x(1)-dom%sh%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if
       dy = x(2)-dom%sh%center(2)
       if (dy > 0.5d0) then 
          dy = dy -1.0d0 
       else if (dy < -0.5d0) then 
          dy = dy + 1.0d0
       end if
       dz = x(3)-dom%sh%center(3)
       if (dz > 0.5d0) then 
          dz = dz -1.0d0 
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if
       rr = dx*dx + dy*dy + dz*dz
       !!rr = (x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2
       if((rr>dom%sh%r_inbound*dom%sh%r_inbound).and.(rr<dom%sh%r_outbound*dom%sh%r_outbound))domain_contains_point=.true.

    case('cube')

       ! correct positions for periodic boundaries 

       dx = x(1)-dom%cu%center(1)
       if (dx > 0.5d0) then 
          dx = dx -1.0d0 
       else if (dx < -0.5d0) then 
          dx = dx + 1.0d0
       end if

       if (abs(dx) < dom%cu%size*0.5d0) then 

          dx = x(2)-dom%cu%center(2)
          if (dx > 0.5d0) then 
             dx = dx -1.0d0 
          else if (dx < -0.5d0) then 
             dx = dx + 1.0d0
          end if

          if (abs(dx) < dom%cu%size*0.5d0) then 

             dx = x(3)-dom%cu%center(3)
             if (dx > 0.5d0) then 
                dx = dx -1.0d0 
             else if (dx < -0.5d0) then 
                dx = dx + 1.0d0
             end if

             if (abs(dx) < dom%cu%size*0.5d0) then

                domain_contains_point=.true.
             end if
          end if
       end if
       
    case('slab')
       dz = x(3) - dom%sl%zc
       if (dz > 0.5d0) then 
          dz = dz - 1.0d0
       else if (dz < -0.5d0) then 
          dz = dz + 1.0d0
       end if
       if(abs(dz) < dom%sl%thickness*0.5d0)domain_contains_point=.true.
    end select
    return
  end function domain_contains_point


end module module_ramses
!==================================================================================
!==================================================================================
    
