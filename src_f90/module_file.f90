module module_file

  use krome_commons
  use krome_main
  use krome_user
  use module_parallel_mpi
  use module_krome
  use module_ramses
  use module_cell
  use module_spectra

  implicit none

  type file
     integer(kind=4)               :: ID
     integer(kind=4)               :: status  ! 0 if not done, 1 if done
  end type file

  integer(kind=4)                  :: ncpu, ncell, ncoarse, ngridmax
  real(kind=8),allocatable         :: xg(:,:)      ! grids position
  real(kind=8),allocatable         :: var(:,:)
  integer,allocatable              :: headl(:,:),taill(:,:),numbl(:,:),numbtot(:,:)
  integer,allocatable              :: headb(:,:),tailb(:,:),numbb(:,:)
  integer,allocatable              :: nbor(:,:)    ! neighboring father cells
  integer,allocatable              :: son(:)       ! sons grids
  integer,allocatable              :: cpu_map(:)  ! domain decomposition
  integer,allocatable              :: next(:)      ! next grid in list
  real(KIND=8),dimension(1:3)      :: xbound=(/0d0,0d0,0d0/) 
  integer,parameter                :: twotondim = 8, ndim = 3, twondim = 6
  logical                          :: first=.true., first2=.true., first3=.true.
  integer(kind=4)                  :: U_precision=8, RT_precision=8 ! RT-precision in RAMSES output

  public :: init_files, compute_file

contains  


  subroutine init_files(repository, snapnum, filegrid)


    implicit none

    character(2000),intent(in)              :: repository
    integer(kind=4),intent(in)              :: snapnum
    type(file),allocatable,intent(inout)    :: filegrid(:)
    integer(kind=4)                         :: i

    ncpu = get_ncpu(repository, snapnum)

    allocate(filegrid(ncpu))
    do i=1,ncpu
       filegrid(i)%ID = i
       filegrid(i)%status = 0
    end do

  end subroutine init_files



  subroutine compute_file(repository, snapnum, icpu, output_path, repository_restart)

    implicit none

    character(2000),intent(in)              :: repository, output_path, repository_restart
    integer(kind=4),intent(in)              :: snapnum, icpu
    real(kind=8),allocatable                :: ramses_var(:,:)
    real(kind=8),allocatable,dimension(:,:) :: cells_rt, fractions
    real(kind=8),allocatable,dimension(:)   :: nH, Tgas, mets, nHI
    integer(kind=4)                         :: nSEDgroups, nvar, nleaf
    character(1000)                         :: nomfich, nomfich2
    logical                                 :: exist


    write(nomfich,'(a,a,a,a,a,i5.5,a,i5.5)') trim(output_path),'/',trim(element_names(elements(n_elements))),trim(roman_num(n_ions(n_elements))),'_',snapnum,'.out',icpu
    inquire(file=nomfich, exist=exist)

    if(exist) then
       print*, 'file ', icpu, ' already computed, passing to the next'
    else

       call read_hydro(repository, snapnum, icpu, nvar, nleaf, ramses_var, repository_restart)
       call init_cells(repository, snapnum, nvar, nleaf, ramses_var)
       deallocate(ramses_var)

       call compute_cells(icpu)

       call write_ion_T_ne_files(snapnum, icpu, output_path)

    end if

  end subroutine compute_file


  subroutine read_hydro(repository, snapnum, icpu, nvar, nleaf, ramses_var, repository_restart)

    ! Reads the hydro part of the Ramses-RT output

    implicit none 
    character(2000),intent(in)                :: repository, repository_restart
    integer(kind=4),intent(in)                :: snapnum, icpu 
    integer(kind=4),intent(inout)             :: nleaf, nvar
    real(kind=8),allocatable, intent(inout)   :: ramses_var(:,:)
    integer(kind=4)                           :: ileaf,nleaf_test,k,ivar,iloop
    real(kind=8)                              :: time1,time2,time3,rate
    integer(kind=8)                           :: c1,c2,c3,cr
    character(1000)                           :: filename 
    logical                                   :: ok_cell
    integer(kind=4)                           :: i,j,ilevel,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)                           :: idim,ind,iu1,iu2,iu3,iu4,rank
    ! stuff read from AMR files
    integer(kind=4)                           :: ngridmax,ngrid_current
    real(kind=8),allocatable                  :: xg(:,:)        ! grids position
    integer(kind=4),allocatable               :: son(:,:)       ! sons grids
    real(KIND=8),dimension(1:3)               :: xbound=(/0d0,0d0,0d0/)  
    integer(kind=4),allocatable               :: ngridfile(:,:),ngridlevel(:,:),ngridbound(:,:)
    integer(kind=4)                           :: ngrida,ncpused
    logical,allocatable                       :: ref(:,:)
    real(kind=8)                              :: dx,boxlen
    integer(kind=4)                           :: ix,iy,iz,nvarH,nvarRT,nvarRT_restart
    real(kind=8),allocatable                  :: xc(:,:),xp(:,:,:)
    ! stuff read from the HYDRO files
    real(kind=8),allocatable                  :: var(:,:,:)
    real(kind=4),allocatable                  :: var_sp(:)
    logical                                   :: cellInDomain
    real(kind=8),dimension(3)                 :: xx
    logical,allocatable                       :: cpu_is_useful(:)
    integer(kind=4)                           :: U_precision=8 ! hydro-precision in RAMSES output
    integer(kind=4)                           :: RT_precision=8 ! RT-precision in RAMSES output


    call get_nleaf_per_cpu(repository,snapnum,icpu,nleaf)
    nvar     = get_nvar(repository,snapnum,repository_restart)

    allocate(ramses_var(nvar,nleaf))

    ! Check whether the ramses output is in single or double precision
    U_precision = nint(get_param_real(repository,snapnum,'U_precision',default_value=8d0))
    RT_precision = nint(get_param_real(repository,snapnum,'rtprecision',default_value=8d0,rt_info=.true.))

    rank = 1
    iu1 = 10+rank*4
    iu2 = 10+rank*4+1
    iu3 = 10+rank*4+2
    iu4 = 10+rank*4+3

    ileaf=1
    iloop=0

    ! verify AMR input file -> already done above in get_nleaf_new
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    ! Open AMR file and skip header
    open(unit=iu1,file=filename,form='unformatted',status='old',action='read')
    read(iu1)ncpu
    read(iu1)      !ndim
    read(iu1)nx,ny,nz
    read(iu1)nlevelmax
    read(iu1)ngridmax
    read(iu1)nboundary
    read(iu1)ngrid_current
    read(iu1)boxlen
    do i=1,13
       read(iu1)
    end do
    !twotondim=2**ndim
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
    if(allocated(ngridfile)) deallocate(ngridfile,ngridlevel)
    allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
    allocate(ngridlevel(1:ncpu,1:nlevelmax))
    if(nboundary>0)then
       if(allocated(ngridbound)) deallocate(ngridbound)
       allocate(ngridbound(1:nboundary,1:nlevelmax))
    endif
    ! Read grid numbers
    read(iu1)ngridlevel
    ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
    read(iu1)
    if(nboundary>0)then
       do i=1,2
          read(iu1)
       end do
       read(iu1)ngridbound
       ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
    endif
    read(iu1)
    ! ROM: comment the single follwing line for old stuff
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)

    if(allocated(xc)) deallocate(xc)
    allocate(xc(1:twotondim,1:ndim))


    ! open hydro file and get nvarH
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=iu2,file=filename,form='unformatted',status='old',action='read')
    read(iu2)
    read(iu2)nvarH
    read(iu2)
    read(iu2)
    read(iu2)
    read(iu2)

    ! Open RT file and get nvarRT
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
    open(unit=iu3,file=filename,status='old',form='unformatted')
    read(iu3)
    read(iu3)nvarRT
    read(iu3)
    read(iu3)
    read(iu3)
    read(iu3)

    ! Same for restart
    if(use_rt_restart) then
       ! Open RT file and get nvarRT
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository_restart),'/output_',snapnum+1,'/rt_',snapnum+1,'.out',icpu
       open(unit=iu4,file=filename,status='old',form='unformatted')
       read(iu4)
       read(iu4)nvarRT_restart
       read(iu4)
       read(iu4)
       read(iu4)
       read(iu4)
    else
       nvarRT_restart = 0
    end if

    !ncoarse = nx*ny*nz
    !ncell   = ncoarse+twotondim*ngridmax

    ! Loop over levels
    do ilevel=1,nlevelmax

       ! Geometry
       dx=0.5**ilevel
       do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc(ind,1)=(dble(ix)-0.5D0)*dx
          xc(ind,2)=(dble(iy)-0.5D0)*dx
          xc(ind,3)=(dble(iz)-0.5D0)*dx
       end do

       ! Allocate work arrays
       if(allocated(xg)) then 
          deallocate(xg,son,var,xp,ref)
       endif
       if(allocated(var_sp)) deallocate(var_sp)
       ngrida=ngridfile(icpu,ilevel)
       if(ngrida>0)then
          allocate(xg(1:ngrida,1:ndim))
          allocate(son(1:ngrida,1:twotondim))
          allocate(var(1:ngrida,1:twotondim,1:nvarh+nvarRT+nvarRT_restart))
          if(rt_Precision.eq.4 .or. U_Precision.eq.4) allocate(var_sp(1:ngrida))
          allocate(xp(1:ngrida,1:twotondim,1:ndim))
          allocate(ref(1:ngrida,1:twotondim))
          ref=.false.
       endif


       ! Loop over domains
       do j=1,nboundary+ncpu

          ! Read AMR data
          if(ngridfile(j,ilevel)>0)then
             read(iu1) ! Skip grid index
             read(iu1) ! Skip next index
             read(iu1) ! Skip prev index
             ! Read grid center
             do idim=1,ndim
                if(j.eq.icpu)then
                   read(iu1)xg(:,idim)
                else
                   read(iu1)
                endif
             end do
             read(iu1) ! Skip father index
             do ind=1,2*ndim
                read(iu1) ! Skip nbor index
             end do
             ! Read son index
             do ind=1,twotondim
                if(j.eq.icpu)then
                   read(iu1)son(:,ind)
                else
                   read(iu1)
                end if
             end do
             ! Skip cpu map
             do ind=1,twotondim
                read(iu1)
             end do
             ! Skip refinement map
             do ind=1,twotondim
                read(iu1)
             end do
          endif

          ! Read HYDRO data
          read(iu2)
          read(iu2)
          read(iu3)
          read(iu3)
          if(use_rt_restart) then
             read(iu4)
             read(iu4)
          end if
          if(ngridfile(j,ilevel)>0)then
             ! Read hydro variables
             do ind=1,twotondim
                do ivar=1,nvarh
                   if(j.eq.icpu)then
                      if(U_precision.eq.4) then
                         read(iu2) var_sp(:)
                         var(:,ind,ivar) = var_sp(:)
                      else
                         read(iu2)var(:,ind,ivar)
                      end if
                   else
                      read(iu2)
                   end if
                end do
                ! Restart before normal RT, since it's lower energies
                if(use_rt_restart) then
                   do ivar=1,nvarRT_restart
                      if(j.eq.icpu)then
                         if(rt_Precision.eq.4) then
                            read(iu4) var_sp(:)
                            var(:,ind,nvarH+ivar) = var_sp(:)
                         else
                            read(iu4)var(:,ind,nvarH+ivar)
                         end if
                      else
                         read(iu4)
                      end if
                   end do
                end if
                !
                do ivar=1,nvarRT
                   if(j.eq.icpu)then
                      if(rt_Precision.eq.4) then
                         read(iu3) var_sp(:)
                         var(:,ind,nvarH+nvarRT_restart+ivar) = var_sp(:)
                      else
                         read(iu3)var(:,ind,nvarh+nvarRT_restart+ivar)
                      end if
                   else
                      read(iu3)
                   end if
                end do
             end do
          end if

       enddo ! end loop over domains

       ! Get leaf cells and store data
       if(ngrida>0)then
          ! Loop over cells
          do ind=1,twotondim
             ! Compute cell center
             do i=1,ngrida
                xp(i,ind,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                xp(i,ind,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                xp(i,ind,3)=(xg(i,3)+xc(ind,3)-xbound(3))
             end do
             ! Check if cell is refined
             do i=1,ngrida
                ref(i,ind)=son(i,ind)>0.and.ilevel<nlevelmax
             end do
             ! Store leaf cells
             do i=1,ngrida
                ok_cell= .not.ref(i,ind)
                if(ok_cell)then
                   !if(.not.ref(i,ind))then
                   cellInDomain=.true.
                   if(cellInDomain)then
                      do ivar = 1,nvar
                         ramses_var(ivar,ileaf) = var(i,ind,ivar)
                      end do
                      ileaf=ileaf+1
                   endif
                endif
             enddo
          end do
       endif


    enddo ! end loop over levels


    close(iu1)
    close(iu2)
    close(iu3)
    close(iu4)

    nleaf_test = ileaf-1
    if(nleaf /= nleaf_test) then
       print*,'Problem with nleaf in module_file. nleaf from ramses,  ileaf-1 : ', nleaf, nleaf_test
       stop
    end if

    return

  end subroutine read_hydro



  function get_nvar(repository,snapnum,repository_restart)
    implicit none 
    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository, repository_restart
    character(1000)             :: filename
    integer(kind=4)             :: get_nvar,icpu,nvarRT,nvarRT_restart

    icpu = 1
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=11,file=filename,form='unformatted',status='old',action='read')
    read(11)
    read(11)get_nvar
    close(11)
    ! Open RT file and get nvarRT
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/rt_',snapnum,'.out',icpu
    open(unit=12,file=filename,status='old',form='unformatted')
    read(12)
    read(12)nvarRT
    close(12)

    if(use_rt_restart) then
       ! Open RT file and get nvarRT
       write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository_restart),'/output_',snapnum+1,'/rt_',snapnum+1,'.out',icpu
       open(unit=13,file=filename,status='old',form='unformatted')
       read(13)
       read(13)nvarRT_restart
       close(13)
    else
       nvarRT_restart = 0
    end if

    get_nvar = get_nvar + nvarRT + nvarRT_restart

    return
  end function get_nvar


  subroutine get_nleaf_per_cpu(repository,snapnum,icpu,nleaf)
    ! purpose: get the number of leaf cells in this icpu file

    implicit none 

    integer(kind=4),intent(in)  :: snapnum,icpu
    character(1000),intent(in)  :: repository
    integer(kind=4),intent(out) :: nleaf

    character(1000)             :: filename 
    logical                     :: ok
    integer(kind=4)             :: i,j,ilevel,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)             :: idim,ind,iu1,rank
    integer(kind=4)             :: ngridmax,ngrid_current,ngrida
    integer,allocatable         :: son(:,:)       ! sons grids
    real(KIND=8),dimension(1:3) :: xbound=(/0d0,0d0,0d0/)  
    integer, allocatable        :: ngridfile(:,:),ngridlevel(:,:),ngridbound(:,:)
    logical, allocatable        :: ref(:,:)
    real(kind=8)                :: boxlen

    rank = 1
    iu1 = 10+rank*3

    ! verify AMR input file
    write(filename,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    inquire(file=filename, exist=ok)
    if(.not. ok)then
       write(*,*)'File '//TRIM(filename)//' not found'    
       stop
    end if

    ! Open AMR file and skip header
    open(unit=iu1,file=filename,form='unformatted',status='old',action='read')
    read(iu1)ncpu
    read(iu1)      !ndim
    read(iu1)nx,ny,nz
    read(iu1)nlevelmax
    read(iu1)ngridmax
    read(iu1)nboundary
    read(iu1)ngrid_current
    read(iu1)boxlen
    do i=1,13
       read(iu1)
    end do
    !twotondim=2**ndim
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
    allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
    allocate(ngridlevel(1:ncpu,1:nlevelmax))
    if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

    ! Read grid numbers
    read(iu1)ngridlevel
    ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
    read(iu1)
    if(nboundary>0)then
       do i=1,2
          read(iu1)
       end do
       read(iu1)ngridbound
       ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
    endif
    read(iu1)
    ! ROM: comment the single follwing line for old stuff
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)
    read(iu1)

    nleaf=0

    ! Loop over levels
    do ilevel=1,nlevelmax

       ! Allocate work arrays
       if(allocated(son)) then 
          deallocate(son,ref)
       endif
       ngrida=ngridfile(icpu,ilevel)
       if(ngrida>0)then
          allocate(son(1:ngrida,1:twotondim))
          allocate(ref(1:ngrida,1:twotondim))
          ref=.false.
       endif

       ! Loop over domains
       do j=1,nboundary+ncpu
          ! Read AMR data
          if(ngridfile(j,ilevel)>0)then
             read(iu1) ! Skip grid index
             read(iu1) ! Skip next index
             read(iu1) ! Skip prev index
             do idim=1,ndim ! Skip grid center
                read(iu1)
             end do
             read(iu1) ! Skip father index
             do ind=1,2*ndim
                read(iu1) ! Skip nbor index
             end do
             ! Read son index
             do ind=1,twotondim
                if(j.eq.icpu)then
                   read(iu1)son(:,ind)
                else
                   read(iu1)
                end if
             end do
             ! Skip cpu map
             do ind=1,twotondim
                read(iu1)
             end do
             ! Skip refinement map
             do ind=1,twotondim
                read(iu1)
             end do
          endif
       enddo

       ! Count leaf cells
       if(ngrida>0)then
          ! Loop over cells
          do ind=1,twotondim
             ! Check if cell is refined
             do i=1,ngrida
                ref(i,ind)=son(i,ind)>0.and.ilevel<nlevelmax
                if (.not.ref(i,ind))then
                   nleaf=nleaf+1
                endif
             end do
          end do
       endif

    enddo

    close(iu1) 
    return
  end subroutine get_nleaf_per_cpu


end module module_file

