module module_cell

  use krome_commons
  use krome_main
  use krome_user
  use module_parallel_mpi
  use module_krome
  use module_ramses


  implicit none

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [UVB] of the parameter file
  ! If you use restart files (when Ramses was restarted to propagate a bin0) (see parameters in ion_fractions.f90), this UVB is ignored.
  ! --------------------------------------------------------------------------
  real(kind=8)              :: threshold    = 1d2     ! Threshold density above which we don't include the UVB for bin0. If negative, no threshold
  real(kind=8)              :: bin1_factor  = 1d1     ! The intensity of bin0 is always at least bin1_factor * bin1. I advise againt putting it smaller than 1 (especially 0).
  real(kind=8)              :: UVB          = 7.8d7   ! Integral of the UVB at the z of the simulation, between the two lowest energies in [spectra] parameters. 7.8d7 is for z=3, HM12, from 6-13.6 eV
  ! --------------------------------------------------------------------------

  type cell
     real(kind=8)                                :: nHI
     real(kind=8)                                :: nHII
     real(kind=8)                                :: nHeII
     real(kind=8)                                :: nHeIII
     real(kind=8)                                :: T
     real(kind=8)                                :: Z
     real(kind=8),dimension(nPhotoRea)           :: rates      !nPhotoRea is the number of photoionization reactions in the chemical network of Krome
     real(kind=8),dimension(nmols-natoms-3)      :: den_ions   !Density of the ionization stages of interest, e.g. Si0, Si+, O0, O+, O++, O+++
     real(kind=8)                                :: ne
  end type cell

  real(kind=8),allocatable,dimension(:,:),public :: csn        ! Effective cross-sections of photoionization, weighted by number of photons
  type(cell),allocatable                         :: cellgrid(:)

  public :: init_cells, init_csn, compute_cells, write_ion_T_ne_files, read_UVB_params, print_UVB_params

contains

  subroutine init_csn(repository, snapnum, csn_file, nGroups, number_ions)
    !----------------------------------------------------------------------------
    ! Compute the effective cross-sections, and save them in a file. If file already exists, it simply reads the file
    !----------------------------------------------------------------------------

    implicit none

    character(2000),intent(in)              :: repository, csn_file
    integer(kind=4),intent(in)              :: snapnum, nGroups, number_ions
    integer(kind=4)                         :: i,j, test_nGroups, test_number_ions
    logical                                 :: exist

    ! nGroups: number of energy bins for the radiation. 4 for SPHINX10 (with UVB), 3 for SPHINX20 (with UVB)

    inquire(file=csn_file, exist=exist)

    ! Reading the existing file, and checking for consistency
    if(exist) then
       open(unit=10, file=csn_file, form='unformatted', action='read')
       read(10) test_nGroups
       if(test_nGroups /= nGroups) then
          print*, 'Problem, the number of photon groups in the csn_file is not correct. Check your parameters or your csn_file'
          call stop_mpi
       end if
       read(10) test_number_ions
       if(test_number_ions /= number_ions) then
          print*, 'Problem, the number of ions in the csn_file is not correct. Check your parameters or your csn_file'
          call stop_mpi
       end if

       allocate(csn(nGroups, number_ions))
       
       do i=1,number_ions
          read(10) (csn(j,i), j=1,nGroups)
       end do
       close(10)

       print*, 'csn (nSEDgroups * nIons)'
       do j=1,nGroups
          print*, csn(j,:)
       end do

    else
       print*, 'Beginning computation of csn, takes a few minutes, up to 1 hour for large simulation boxes.'
       call compute_csn_in_box(repository, snapnum, n_elements, elements, n_ions, csn) ! n_elements and elements defined in module_krome.f90

       open(unit=10, file=csn_file, form='unformatted', action='write')
       write(10) nGroups
       write(10) number_ions
       do i=1,number_ions
          write(10) (csn(j,i), j=1,nGroups)
       end do
       close(10)
       
    end if

  end subroutine init_csn


  subroutine init_cells(repository, snapnum, nvar, ncell, ramses_var)
    ! Sets the properties of all gas cells. Densities, temperature, photoionization rates, metallicity.

    use module_ramses
    use module_krome
    use module_spectra

    implicit none

    character(2000),intent(in)                    :: repository
    integer(kind=4),intent(in)                    :: snapnum, nvar, ncell
    real(kind=8),intent(in),dimension(ncell,nvar) :: ramses_var    ! Cell properties read directly from the Ramses-RT output
    real(kind=8),allocatable                      :: cells_rt(:,:)
    real(kind=8),dimension(ncell)                 :: nH, Tgas, mets, nHI
    real(kind=8),dimension(3,ncell)               :: fractions
    integer(kind=4)                               :: nSEDgroups, i, j, k
    real(kind=8)                                  :: bin0

    if(use_rt_restart) then
       nSEDgroups = get_nOptBins()
    else
       nSEDgroups = get_nOptBins() - 1     ! Remove 1 because in this case the non-ionising group is set by a UVB model
    end if
    allocate(cells_rt(nSEDgroups,ncell))

    ! From ramses_var, those methods provide desired quantities in the correct (cgs) units
    call ramses_get_nh_cgs(repository,snapnum,ncell,nvar,ramses_var,nH)
    call ramses_get_T_nhi_cgs(repository,snapnum,ncell,nvar,ramses_var,Tgas,nHI)
    call ramses_get_metallicity(ncell,nvar,ramses_var,mets)
    call ramses_get_fractions(ncell,nvar,ramses_var,fractions) ! Ionization fractions of H and He
    call ramses_get_flux(repository,snapnum,ncell,nvar,nSEDgroups,ramses_var,cells_rt)

    allocate(cellgrid(ncell))
    do i=1,ncell
       cellgrid(i)%nHI = nH(i)*(1d0-fractions(1,i))
       cellgrid(i)%nHII = nH(i)*fractions(1,i)
       cellgrid(i)%nHeII = 7.895d-2*nH(i)*fractions(2,i)
       cellgrid(i)%nHeIII = 7.895d-2*nH(i)*fractions(3,i)
       cellgrid(i)%T = Tgas(i)
       cellgrid(i)%Z = mets(i)
       

       if(use_rt_restart) then
          do j=1,nPhotoRea
             cellgrid(i)%rates(j) = sum(cells_rt(:,i)*csn(:,j)) ! Photoionization rates are products and photon fluxes and effective cross-sections
          end do
       else

          ! In this case, we use our model of UV background
          if(cellgrid(i)%nHI < threshold) then
             bin0 = max(UVB, bin1_factor*cells_rt(1,i))
          else
             bin0 = bin1_factor*cells_rt(1,i)
          end if

          do j=1,nPhotoRea
             cellgrid(i)%rates(j) = bin0*csn(1,j) + sum(cells_rt(:,i)*csn(2:nSEDgroups+1,j))
          end do
       end if

       cellgrid(i)%den_ions(:) = 0d0 ! Initializing the ionization fractions of metals (our end desired quantities) to 0
    end do

    deallocate(cells_rt)

  end subroutine init_cells


  subroutine compute_cells(icpu)
    ! Main part of the code, that computes the ionization fractions with KROME. icpu corresponds to the part of the output written by one core used in the Ramses-RT run.

    implicit none

    integer(kind=4),intent(in)         :: icpu
    integer(kind=4)                    :: i, j, k, l, ion_state(n_elements), ncell, non_zero_index(n_elements)
    real(kind=8)                       :: densities(nmols), n_ion_save(n_elements)
    logical                            :: krome_converged

    ncell = size(cellgrid)

    do i=1,ncell

       if(cellgrid(i)%T > 1d10) then
          cellgrid(i)%den_ions(:) = 1d-18
       else

          call krome_set_photoBin_rates(cellgrid(i)%rates) ! Uses my custom Krome routine to set the rates by hand

          densities(:)               = 1d-18
          densities(krome_idx_H)     = max(cellgrid(i)%nHI,1d-18)
          densities(krome_idx_Hj)    = max(cellgrid(i)%nHII,1d-18)
          densities(krome_idx_He)    = max(7.895d-2*(cellgrid(i)%nHI + cellgrid(i)%nHII) - cellgrid(i)%nHeII - cellgrid(i)%nHeIII,1d-18)
          densities(krome_idx_Hej)   = max(cellgrid(i)%nHeII,1d-18)
          densities(krome_idx_Hejj)  = max(cellgrid(i)%nHeIII,1d-18)
          densities(krome_idx_E)     = densities(krome_idx_Hj) + densities(krome_idx_Hej) + 2*densities(krome_idx_Hejj)


          do j=1,n_elements
             n_ion_save(j) = cellgrid(i)%Z/0.0134*abundances(j)*(densities(krome_idx_H) + densities(krome_idx_Hj))
             non_zero_index(j) = get_non_zero_index(j,cellgrid(i)%T,ion_state(j))
             densities(non_zero_index(j)) = max(n_ion_save(j), 1d-18)
             densities(krome_idx_E) = densities(krome_idx_E) + (ion_state(j)-1)*densities(non_zero_index(j))
          end do
          densities(krome_idx_E) = max(densities(krome_idx_E),1d-18)

          call krome_equilibrium(densities, cellgrid(i)%T, icpu, i, krome_converged)

          cellgrid(i)%ne = max(densities(krome_idx_E), 1d-18)

          l=0
          do j=1,n_elements
             !Rare cases of bugs :
             if(abs( sum(densities(indices(j,1:n_ions(j)+1))) - n_ion_save(j)) / n_ion_save(j) > 0.01) then
                print*, 'small bug in file ', icpu
                print*, 'temperature, nHI, nHII, photorates, metallicity'
                print*, cellgrid(i)%T, cellgrid(i)%nHI, cellgrid(i)%nHII, cellgrid(i)%rates(:), cellgrid(i)%Z
                cellgrid(i)%den_ions(l+1:l+n_ions(j)) = 1d-18
                if(ion_state(j) <= n_ions(j)) cellgrid(i)%den_ions(l+ion_state(j)) = max(n_ion_save(j),1d-18)
             !No bug
             else
                do k=1,n_ions(j)
                   if(indices(j,k) < 1) then
                      print*, 'Problem ! Trying to update a ionization state that does not exist'
                      call stop_mpi
                   end if
                   cellgrid(i)%den_ions(k+l) = max(densities(indices(j,k)), 1d-18)
                end do
             end if
             l = l + n_ions(j)
          end do
       end if
    end do
  end subroutine compute_cells
  

  subroutine write_ion_T_ne_files(snapnum, icpu, output_path)
    
    implicit none

    character(2000),intent(in)              :: output_path
    integer(kind=4),intent(in)              :: snapnum, icpu
    character(1000)                         :: nomfich
    integer(kind=4)                         :: i,j,k,l,start, ncell

    ncell = size(cellgrid)

    l=0
    do j=1,n_elements
       do k=min(n_ions(j),2),n_ions(j)
          write(nomfich,'(a,a,a,a,a,i5.5,a,i5.5)') trim(output_path),'/',trim(element_names(elements(j))),trim(roman_num(k)),'_',snapnum,'.out',icpu
          open(unit=10, file=nomfich, form='unformatted', action='write')
          write(10) ncell
          write(10) (cellgrid(i)%den_ions(k+l), i=1,ncell)
          close(10)
       end do
       l = l + n_ions(j)
    end do

    write(nomfich,'(a,a,i5.5,a,i5.5)') trim(output_path),'/T_ne_',snapnum,'.out',icpu
    open(unit=10, file=nomfich, form='unformatted', action='write')
    write(10) ncell
    write(10) (cellgrid(i)%T, i=1,ncell)
    write(10) (cellgrid(i)%ne, i=1,ncell)
    close(10)

    deallocate(cellgrid)

  end subroutine write_ion_T_ne_files
  

  subroutine read_UVB_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (mesh)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    logical         :: ndomain_present 

    section_present = .false.
    ndomain_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:5) == '[UVB]') then
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
          case('threshold')
             read(value,*) threshold
             if(threshold<0) threshold = 1d18
          case('bin1_factor')
             read(value,*) bin1_factor
          case ('UVB')
             read(value,*) UVB
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_UVB_params


  subroutine print_UVB_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    character(100) :: fmt

    if (present(unit)) then
       write(*,'(a)')            '--------------------------------------------------------------------------------'
       write(unit,'(a,a,a)')     '[UVB]'
       write(unit,'(a,ES10.3)')  '  UVB          = ',UVB
       write(unit,'(a,ES10.3)')  '  threshold    = ',threshold
       write(unit,'(a,ES10.3)')  '  bin1_factor  = ',bin1_factor
       write(*,'(a)')            '--------------------------------------------------------------------------------'
    else
       write(*,'(a)')            '--------------------------------------------------------------------------------'
       write(*,'(a,a,a)')     '[UVB]'
       write(*,'(a,ES10.3)')  '  UVB          = ',UVB
       write(*,'(a,ES10.3)')  '  threshold    = ',threshold
       write(*,'(a,ES10.3)')  '  bin1_factor  = ',bin1_factor
       write(*,'(a)')            '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_UVB_params
    
end module module_cell


