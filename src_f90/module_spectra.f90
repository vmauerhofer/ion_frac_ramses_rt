MODULE spectrum_integrator_module
  !_________________________________________________________________________
  implicit none

  PUBLIC integrateSpectrum, f1, fLambda, fdivLambda, fSig, fSigLambda,   &
       fSigdivLambda, trapz1, get_prop

  PRIVATE   ! default

CONTAINS

  !*************************************************************************
  FUNCTION integrateSpectrum(X, Y, N, e0, e1, species, func, doPrint)

    ! Integrate spectral weighted function in energy interval [e0,e1]
    ! X      => Wavelengths [angstrom]
    ! Y      => Spectral luminosity per angstrom at wavelenghts [XX A-1] 
    ! N      => Length of X and Y
    ! e0,e1  => Integrated interval [ev]
    ! species=> ion species, used as an argument in fx
    ! func   => Function which is integrated (of X, Y, species)
    !-------------------------------------------------------------------------
    use module_constants,only:clight,evtoerg, planck

    real(kind=8):: integrateSpectrum, X(N), Y(N), e0, e1
    integer :: N, species
    interface
       real(kind=8) function func(wavelength,intensity,species)
         real(kind=8)::wavelength,intensity
         integer::species
       end function func
    end interface!----------------------------------------------------------
    real(kind=8),dimension(:),allocatable:: xx, yy, f
    real(kind=8):: la0, la1
    integer :: i
    logical,optional::doPrint
    !-------------------------------------------------------------------------
    integrateSpectrum=0.
    if(N .le. 2) RETURN
    ! Convert energy interval to wavelength interval
    la0 = X(1) ; la1 = X(N)
    if(e1.gt.0) la0 = max(la0, 1.d8 * planck * clight / e1 / evtoerg)                
    if(e0.gt.0) la1 = min(la1, 1.d8 * planck * clight / e0 / evtoerg)
    if(la0 .ge. la1) RETURN         
    ! If we get here, the [la0, la1] inverval is completely within X
    allocate(xx(N)) ; allocate(yy(N)) ; allocate(f(N))
    xx =  la0   ;   yy =  0.   ;   f = 0.
    i=2
    do while ( i.lt.N .and. X(i).le.la0 )
       i = i+1                      !              Below wavelength interval
    enddo                           !   X(i) is now the first entry .gt. la0 
    ! Interpolate to value at la0
    yy(i-1) = Y(i-1) + (xx(i-1)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
    f(i-1)  = func(xx(i-1), yy(i-1), species)
    do while ( i.lt.N .and. X(i).le.la1 )              ! Now within interval
       xx(i) = X(i) ; yy(i) = Y(i) ; f(i) = func(xx(i),yy(i),species)
       i = i+1
    enddo                          ! i=N or X(i) is the first entry .gt. la1
    xx(i:) = la1                   !             Interpolate to value at la1
    yy(i) = Y(i-1) + (xx(i)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
    f(i)  = func(xx(i),yy(i),species) 

    integrateSpectrum = trapz1(xx,f,i)
    deallocate(xx) ; deallocate(yy) ; deallocate(f)

  END FUNCTION integrateSpectrum


  !*************************************************************************
  ! FUNCTIONS FOR USE WITH integrateSpectrum:
  ! lambda  => wavelengths in Angstrom
  ! f       => function of wavelength (a spectrum in some units)
  ! species => 1=HI, 2=HeI or 3=HeII
  !_________________________________________________________________________
  FUNCTION f1(lambda, f, species)
    real(kind=8):: f1, lambda, f
    integer :: species
    f1 = f
  END FUNCTION f1

  FUNCTION fLambda(lambda, f, species)
    real(kind=8):: fLambda, lambda, f
    integer :: species
    fLambda = f * lambda
  END FUNCTION fLambda

  FUNCTION fdivLambda(lambda, f, species)
    real(kind=8):: fdivlambda, lambda, f
    integer :: species
    fdivLambda = f / lambda
  END FUNCTION fdivLambda

  FUNCTION fSig(lambda, f, species)
    real(kind=8):: fSig, lambda, f
    integer :: species
    fSig = f * getCrosssection_Hui(lambda,species)
  END FUNCTION fSig

  FUNCTION fSigLambda(lambda, f, species)
    real(kind=8):: fSigLambda, lambda, f
    integer :: species
    fSigLambda = f * lambda * getCrosssection_Hui(lambda,species)
  END FUNCTION fSigLambda

  FUNCTION fSigdivLambda(lambda, f, species)
    real(kind=8):: fSigdivLambda, lambda, f
    integer :: species
    fSigdivLambda = f / lambda * getCrosssection_Hui(lambda,species)
  END FUNCTION fSigdivLambda
  !_________________________________________________________________________


  !*************************************************************************
  FUNCTION trapz1(X,Y,N,cum)

    ! Integrates function Y(X) along the whole interval 1..N, using a very 
    ! simple staircase method and returns the result.
    ! Optionally, the culumative integral is returned in the cum argument.
    !-------------------------------------------------------------------------
    integer :: N,i
    real(kind=8):: trapz1
    real(kind=8):: X(N),Y(N)
    real(kind=8),optional::cum(N)
    real(kind=8),allocatable::cumInt(:)
    !-------------------------------------------------------------------------
    allocate(cumInt(N))
    cumInt(:)=0.d0
    if (N.le.1) RETURN
    do i=2,N
       cumInt(i)= cumInt(i-1) + abs(X(i)-X(i-1)) * (Y(i)+Y(i-1)) / 2.d0
    end do
    trapz1 = cumInt(N)
    if(present(cum)) cum=cumInt
    deallocate(cumInt)
  END FUNCTION trapz1


  !*************************************************************************
  FUNCTION getCrosssection_Hui(lambda, species)

    ! Gives an atom-photon cross-section of given species at given wavelength, 
    ! as given by Hui and Gnedin (1997).
    ! lambda  => Wavelength in angstrom
    ! species => 1=HI, 2=HeI or 3=HeII
    ! returns :  photoionization cross-section in cm^2
    !------------------------------------------------------------------------
    use module_constants,only:clight, evtoerg, ionEvs, planck
    real(kind=8)      :: lambda, getCrosssection_Hui
    integer           :: species
    real(kind=8)      :: E, E0, cs0, P, ya, yw, y0, y1, x, y
    !------------------------------------------------------------------------
    E = planck * clight/(lambda*1.d-8) / evtoerg         ! photon energy in ev
    if ( E .lt. ionEvs(species) ) then            ! below ionization energy
       getCrosssection_Hui=0.
       RETURN
    endif
    select case (species) ! Fitting parameters from Verner et al. 1998
    case(1)  !HI
       E0 = 4.298e-01 ; cs0 = 5.475e-14 ; P = 2.963e+00
       ya = 3.288e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(2)  !HeI
       E0 = 1.361e+01 ; cs0 = 9.492e-16 ; P = 3.188e+00
       ya = 1.469e+00 ; yw = 2.039e+00 ; y0 = 4.434e-01 ; y1 = 2.136e+00
    case(3)  !HeII
       E0 = 1.720e+00 ; cs0 = 1.369e-14 ; P = 2.963e+00
       ya = 3.288e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(4)  !CI
       E0 = 2.144e+00 ; cs0 = 5.027e-16 ; P = 5.101e+00
       ya = 6.216e+01 ; yw = 9.157e-02 ; y0 = 1.133e+00 ; y1 = 1.607e+00
    case(5)  !CII
       E0 = 4.058e-01 ; cs0 = 8.709e-18 ; P = 8.578e+00
       ya = 1.261e+02 ; yw = 2.093e+00 ; y0 = 4.929e+01 ; y1 = 3.234e+00
    case(6)  !CIII
       E0 = 4.614e+00 ; cs0 = 1.539e-14 ; P = 1.593e+01
       ya = 1.737e+00 ; yw = 5.922e+00 ; y0 = 4.378e-03 ; y1 = 2.528e-02
    case(7)  !CIV
       E0 = 3.506e+00 ; cs0 = 1.068e-16 ; P = 7.457e+00
       ya = 1.436e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(8)  !CV
       E0 = 4.624e+01 ; cs0 = 2.344e-16 ; P = 2.581e+00
       ya = 2.183e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(9)  !CVI
       E0 = 1.548e+01 ; cs0 = 1.521e-15 ; P = 2.963e+00
       ya = 3.288e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(10)  !NI
       E0 = 4.034e+00 ; cs0 = 8.235e-16 ; P = 3.928e+00
       ya = 8.033e+01 ; yw = 9.097e-02 ; y0 = 8.598e-01 ; y1 = 2.325e+00
    case(11)  !NII
       E0 = 6.128e-02 ; cs0 = 1.944e-18 ; P = 8.773e+00
       ya = 8.163e+02 ; yw = 1.043e+01 ; y0 = 4.280e+02 ; y1 = 2.030e+01
    case(12)  !NIII
       E0 = 2.420e-01 ; cs0 = 9.375e-19 ; P = 9.156e+00
       ya = 2.788e+02 ; yw = 1.850e+00 ; y0 = 1.877e+02 ; y1 = 3.999e+00
    case(13)  !NIV
       E0 = 5.494e+00 ; cs0 = 1.690e-14 ; P = 1.706e+01
       ya = 1.714e+00 ; yw = 7.904e+00 ; y0 = 6.415e-03 ; y1 = 1.937e-02
    case(14)  !NV
       E0 = 4.471e+00 ; cs0 = 8.376e-17 ; P = 6.003e+00
       ya = 3.297e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(15)  !NVI
       E0 = 6.943e+01 ; cs0 = 1.519e-16 ; P = 2.315e+00
       ya = 2.627e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(16)  !OI
       E0 = 1.240e+00 ; cs0 = 1.745e-15 ; P = 1.764e+01
       ya = 3.784e+00 ; yw = 7.589e-02 ; y0 = 8.698e+00 ; y1 = 1.271e-01
    case(17)  !OII
       E0 = 1.386e+00 ; cs0 = 5.967e-17 ; P = 8.943e+00
       ya = 3.175e+01 ; yw = 1.934e-02 ; y0 = 2.131e+01 ; y1 = 1.503e-02
    case(18)  !OIII
       E0 = 1.723e-01 ; cs0 = 6.753e-16 ; P = 6.822e+00
       ya = 3.852e+02 ; yw = 1.191e-01 ; y0 = 3.839e-03 ; y1 = 4.569e-01
    case(19)  !OIV
       E0 = 2.044e-01 ; cs0 = 8.659e-19 ; P = 8.785e+00
       ya = 4.931e+02 ; yw = 3.143e+00 ; y0 = 3.328e+02 ; y1 = 4.285e+01
    case(20)  !OV
       E0 = 2.854e+00 ; cs0 = 1.642e-14 ; P = 2.647e+01
       ya = 1.792e+00 ; yw = 2.836e+01 ; y0 = 3.036e-02 ; y1 = 5.554e-02
    case(21)  !OVI
       E0 = 7.824e+00 ; cs0 = 6.864e-17 ; P = 5.495e+00
       ya = 3.210e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(22)  !OVII
       E0 = 8.709e+01 ; cs0 = 1.329e-16 ; P = 2.336e+00
       ya = 2.535e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(23)  !OVIII
       E0 = 2.754e+01 ; cs0 = 8.554e-16 ; P = 2.963e+00
       ya = 3.288e+01 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(24)  !NeI
       E0 = 4.870e+00 ; cs0 = 4.287e-15 ; P = 8.355e+00
       ya = 5.798e+00 ; yw = 2.434e-01 ; y0 = 4.236e-02 ; y1 = 5.873e+00
    case(25)  !NeII
       E0 = 1.247e+01 ; cs0 = 1.583e-15 ; P = 7.810e+00
       ya = 3.935e+00 ; yw = 6.558e-02 ; y0 = 1.520e+00 ; y1 = 1.084e-01
    case(26)  !NeIII
       E0 = 7.753e-01 ; cs0 = 5.708e-18 ; P = 1.005e+01
       ya = 6.725e+01 ; yw = 4.633e-01 ; y0 = 7.654e+01 ; y1 = 2.023e+00
    case(27)  !NeIV
       E0 = 5.566e+00 ; cs0 = 1.685e-15 ; P = 3.056e+00
       ya = 6.409e+02 ; yw = 8.290e-03 ; y0 = 5.149e+00 ; y1 = 6.687e+00
    case(28)  !NeV
       E0 = 1.248e+00 ; cs0 = 2.430e-18 ; P = 8.999e+00
       ya = 1.066e+02 ; yw = 6.855e-01 ; y0 = 9.169e+01 ; y1 = 3.702e-01
    case(29)  !NeVI
       E0 = 1.499e+00 ; cs0 = 9.854e-19 ; P = 8.836e+00
       ya = 1.350e+02 ; yw = 1.656e+00 ; y0 = 1.042e+02 ; y1 = 1.435e+00
    case(30)  !MgI
       E0 = 1.197e+01 ; cs0 = 1.372e-10 ; P = 1.574e+01
       ya = 2.228e-01 ; yw = 2.805e-01 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(31)  !MgII
       E0 = 8.139e+00 ; cs0 = 3.278e-18 ; P = 3.610e+00
       ya = 4.341e+07 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(32)  !MgIII
       E0 = 1.086e+01 ; cs0 = 5.377e-16 ; P = 7.117e+00
       ya = 9.779e+00 ; yw = 2.604e+00 ; y0 = 4.860e+00 ; y1 = 3.722e+00
    case(33)  !MgIV
       E0 = 2.912e+01 ; cs0 = 1.394e-15 ; P = 6.487e+00
       ya = 2.895e+00 ; yw = 4.326e-02 ; y0 = 9.402e-01 ; y1 = 1.135e-01
    case(34)  !MgV
       E0 = 9.762e-01 ; cs0 = 1.728e-18 ; P = 1.006e+01
       ya = 9.184e+01 ; yw = 8.090e-01 ; y0 = 1.276e+02 ; y1 = 3.979e+00
    case(35)  !MgVI
       E0 = 1.711e+00 ; cs0 = 2.185e-18 ; P = 9.202e+00
       ya = 9.350e+01 ; yw = 6.325e-01 ; y0 = 1.007e+02 ; y1 = 1.729e+00
    case(36)  !AlI
       E0 = 1.381e+01 ; cs0 = 7.195e-18 ; P = 3.642e+00
       ya = 1.621e+03 ; yw = 3.166e-01 ; y0 = 2.041e-01 ; y1 = 4.753e-01
    case(37)  !AlII
       E0 = 2.048e-01 ; cs0 = 6.948e-20 ; P = 9.049e+00
       ya = 5.675e+02 ; yw = 4.615e-01 ; y0 = 9.149e+01 ; y1 = 6.565e-01
    case(38)  !AlIII
       E0 = 1.027e+01 ; cs0 = 4.915e-18 ; P = 3.477e+00
       ya = 1.990e+06 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(39)  !AlIV
       E0 = 3.130e+00 ; cs0 = 1.513e-17 ; P = 1.180e+01
       ya = 1.674e+01 ; yw = 5.342e+00 ; y0 = 3.994e+01 ; y1 = 4.803e+00
    case(40)  !AlV
       E0 = 2.414e+01 ; cs0 = 2.925e-16 ; P = 6.724e+00
       ya = 6.973e+00 ; yw = 1.000e-01 ; y0 = 3.495e+00 ; y1 = 2.701e-01
    case(41)  !AlVI
       E0 = 3.483e-01 ; cs0 = 1.962e-20 ; P = 2.084e+01
       ya = 1.856e+01 ; yw = 8.839e+00 ; y0 = 5.675e-02 ; y1 = 2.768e-01
    case(42)  !SiI
       E0 = 2.317e+01 ; cs0 = 2.506e-17 ; P = 3.546e+00
       ya = 2.057e+01 ; yw = 2.837e-01 ; y0 = 1.672e-05 ; y1 = 4.207e-01
    case(43)  !SiII
       E0 = 2.556e+00 ; cs0 = 4.140e-18 ; P = 1.191e+01
       ya = 1.337e+01 ; yw = 1.570e+00 ; y0 = 6.634e+00 ; y1 = 1.272e-01
    case(44)  !SiIII
       E0 = 1.659e-01 ; cs0 = 5.790e-22 ; P = 1.336e+01
       ya = 1.474e+02 ; yw = 8.626e-01 ; y0 = 9.613e+01 ; y1 = 6.442e-01
    case(45)  !SiIV
       E0 = 1.288e+01 ; cs0 = 6.083e-18 ; P = 3.353e+00
       ya = 1.356e+06 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(46)  !SiV
       E0 = 7.761e-01 ; cs0 = 8.863e-19 ; P = 9.980e+00
       ya = 1.541e+02 ; yw = 1.303e+00 ; y0 = 2.009e+02 ; y1 = 4.537e+00
    case(47)  !SiVI
       E0 = 6.305e+01 ; cs0 = 7.293e-17 ; P = 2.400e+00
       ya = 1.558e+02 ; yw = 2.989e-03 ; y0 = 1.115e+00 ; y1 = 8.051e-02
    case(48)  !S I
       E0 = 1.808e+01 ; cs0 = 4.564e-14 ; P = 1.361e+01
       ya = 1.000e+00 ; yw = 6.385e-01 ; y0 = 9.935e-01 ; y1 = 2.486e-01
    case(49)  !S II
       E0 = 8.787e+00 ; cs0 = 3.136e-16 ; P = 1.281e+01
       ya = 3.442e+00 ; yw = 7.354e-01 ; y0 = 2.782e+00 ; y1 = 1.788e-01
    case(50)  !S III
       E0 = 2.027e+00 ; cs0 = 6.666e-18 ; P = 8.611e+00
       ya = 5.454e+01 ; yw = 4.109e+00 ; y0 = 1.568e+01 ; y1 = 9.421e+00
    case(51)  !S IV
       E0 = 2.173e+00 ; cs0 = 2.606e-18 ; P = 8.655e+00
       ya = 6.641e+01 ; yw = 1.863e+00 ; y0 = 1.975e+01 ; y1 = 3.361e+00
    case(52)  !S V
       E0 = 1.713e-01 ; cs0 = 5.072e-22 ; P = 1.307e+01
       ya = 1.986e+02 ; yw = 7.880e-01 ; y0 = 9.424e+01 ; y1 = 6.265e-01
    case(53)  !S VI
       E0 = 1.413e+01 ; cs0 = 9.139e-18 ; P = 3.626e+00
       ya = 1.656e+03 ; yw = 0.000e+00 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(54)  !FeI
       E0 = 5.461e-02 ; cs0 = 3.062e-19 ; P = 7.923e+00
       ya = 2.671e+07 ; yw = 2.069e+01 ; y0 = 1.382e+02 ; y1 = 2.481e-01
    case(55)  !FeII
       E0 = 1.761e-01 ; cs0 = 4.365e-15 ; P = 5.204e+00
       ya = 6.298e+03 ; yw = 1.141e+01 ; y0 = 9.272e+01 ; y1 = 1.075e+02
    case(56)  !FeIII
       E0 = 1.698e-01 ; cs0 = 6.107e-18 ; P = 8.055e+00
       ya = 1.555e+03 ; yw = 8.698e+00 ; y0 = 1.760e+02 ; y1 = 1.847e+01
    case(57)  !FeIV
       E0 = 2.544e+01 ; cs0 = 3.653e-16 ; P = 6.538e+00
       ya = 8.913e+00 ; yw = 5.602e-01 ; y0 = 0.000e+00 ; y1 = 0.000e+00
    case(58)  !FeV
       E0 = 7.256e-01 ; cs0 = 1.523e-21 ; P = 1.767e+01
       ya = 3.736e+01 ; yw = 5.064e+01 ; y0 = 8.871e+01 ; y1 = 5.280e-02
    case(59)  !FeVI
       E0 = 2.656e+00 ; cs0 = 5.259e-19 ; P = 1.632e+01
       ya = 1.450e+01 ; yw = 1.558e+01 ; y0 = 3.361e+01 ; y1 = 3.743e-03

    end select

    x = E/E0 - y0
    y = sqrt(x**2+y1**2)

    getCrosssection_Hui = &
         cs0 * ((x-1.)**2 + yw**2) * y**(0.5*P-5.5)/(1.+sqrt(y/ya))**P

  END FUNCTION getCrosssection_Hui

  function get_prop(n_elements,elements,nIons)

    implicit none

    integer(kind=4),intent(in)  :: n_elements, elements(n_elements), nIons(n_elements)
    integer(kind=4),allocatable :: get_prop(:)
    integer(kind=4)             :: nProp, i, j, k

    do i=1,n_elements
       if(elements(i) /= 6 .and. elements(i) /= 7 .and. elements(i) /= 8 .and. elements(i) /= 10 .and. elements(i) /= 12 .and. elements(i) /= 13 .and. elements(i) /= 14 .and. elements(i) /= 16 .and. elements(i) /= 26) then
          print*, 'The photoionization is not yet implemented for elements other than Carbon(6), Nitrogen(7), Oxygen(8), Neon(10, Magnesium(12), Aluminium(13) and Silicone(14), Sulfur(16) and Iron(26), please use those ones, sorry'
          stop
       end if
       if(nIons(i) > 6 .and. elements(i) /= 8) then
          print*, 'The photoionization is not yet implemented for nIons>6, sorry (except for oxygen)'
          stop
       end if
    end do
    
    nProp = sum(nIons)
    allocate(get_prop(nProp))

    j=1
    do i=1,n_elements
       if(elements(i) == 6) then
          get_prop(j) = 4
       else if(elements(i) == 7) then
          get_prop(j) = 10
       else if(elements(i) == 8) then
          get_prop(j) = 16
       else if(elements(i) == 10) then
          get_prop(j) = 24
       else if(elements(i) == 12) then
          get_prop(j) = 30
       else if(elements(i) == 13) then
          get_prop(j) = 36
       else if(elements(i) == 14) then
          get_prop(j) = 42
       else if(elements(i) == 16) then
          get_prop(j) = 48
       else
          get_prop(j) = 54
       end if

       do k=1,Nions(i)-1
          get_prop(j+k) = get_prop(j)+k
       end do

       j = Nions(i)+j
    end do
  end function get_prop
    


END MODULE spectrum_integrator_module





module module_spectra

  use module_constants

  implicit none

  public init_SED_table, inp_SED_table, get_nOptBins, get_csn_indices, deallocate_table, read_spectra_params, print_spectra_params, use_rt_restart

  private

  ! Light properties for different spectral energy distributions----------
  integer::SED_nA, SED_nZ=8           ! Number of age bins and Z bins
  ! Age and z logarithmic intervals and lowest values:
  real(kind=8),parameter::SED_dlgA=0.02d0
  real(kind=8)::SED_dlgZ
  real(kind=8)::SED_lgA0, SED_lgZ0
  real(kind=8),allocatable,dimension(:,:,:,:)::SED_table
  real(kind=8),allocatable,dimension(:)::SED_ages, SED_zeds

  ! SED_table: iAges, imetallicities, igroups, properties 
  !                                         (Lum, Lum-acc, egy, csn, cse).
  ! Lum is photons per sec per solar mass (eV per sec per solar mass in 
  ! the case of SED_isEgy=true). Lum-acc is accumulated lum.

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [spectra] of the parameter file
  ! --------------------------------------------------------------------------
  character(200)             ::   sed_dir            = './sed'          ! Where the file with star lists is, in case star_reading_method = 'fromlist'
  logical                    ::   use_rt_restart     = .false.          ! whether to use Harley's restarted Sphinx20
  integer(kind=4)            ::   nSEDgroups_restart = 2                ! number of photon groups
  logical                    ::   verbose            = .false.          ! display some run-time info on this module
  integer(kind=4)            ::   nSEDgroups         = 1                ! number of photon groups
  real(kind=8),allocatable   ::   group_decomp(:)                       ! energy intervals of the photon groups, in eV
  ! --------------------------------------------------------------------------
  integer(kind=4)            :: nSEDgroups_tot



contains

  !#################################################################################################
  SUBROUTINE init_SED_table()

    ! Initiate SED properties table, which gives photon luminosities, 
    ! integrated luminosities, average photon cross sections and energies of 
    ! each photon group as a function of stellar population age and
    ! metallicity.  The SED is read from a directory specified by sed_dir.
    !-------------------------------------------------------------------------

    use spectrum_integrator_module

    integer:: nAges, nzs, nLs              ! # of bins of age, z, wavelength
    integer,parameter::nIons=57             !BAD !!! This is hardcoded,   should think about that
    real(kind=8),allocatable::ages(:), Zs(:), Ls(:), rebAges(:)
    real(kind=8),allocatable::SEDs(:,:,:)           ! SEDs f(lambda,age,met)
    real(kind=8),allocatable::tbl(:,:,:), tbl2(:,:,:), reb_tbl(:,:,:)
    integer::i,ia,iz,ip,ii,dum
    character(len=128)::fZs, fAges, fSEDs                        ! Filenames
    logical::ok,okAge,okZ
    real(kind=8)::dlgA, pL0, pL1, tmp
    integer::locid,ncpu2,ierr
    integer::nv=3+2*nIons  ! # vars in SED table: L,Lacc,egy,nions*(csn,cse)
    integer,parameter::tag=1132
    integer::dummy_io,info2
    real(kind=8), parameter:: Gyr2sec = 3.15569d+16       !       Gyr to sec conversion constant

    !-------------------------------------------------------------------------
    inquire(FILE=TRIM(sed_dir)//'/all_seds.dat', exist=ok)
    if(.not. ok)then 
       write(*,*)'Cannot access SED directory ',TRIM(sed_dir)
       write(*,*)'Directory '//TRIM(sed_dir)//' not found'
       write(*,*)'You need to set the RAMSES_SED_DIR envvar' // &
            ' to the correct path, or use the namelist.'
       stop
    end if
    write(fZs,'(a,a)')   trim(sed_dir),"/metallicity_bins.dat"
    write(fAges,'(a,a)') trim(sed_dir),"/age_bins.dat"
    write(fSEDs,'(a,a)') trim(sed_dir),"/all_seds.dat"
    inquire(file=fZs, exist=okZ)
    inquire(file=fAges, exist=okAge)
    inquire(file=fSEDs, exist=ok)
    if(.not. ok .or. .not. okAge .or. .not. okZ) then
       write(*,*) 'Cannot read SED files...'
       write(*,*) 'Check if SED-directory contains the files ',  &
            'metallicity_bins.dat, age_bins.dat, all_seds.dat'
       stop
    end if


    ! READ METALLICITY BINS-------------------------------------------------
    open(unit=10,file=fZs,status='old',form='formatted')
    read(10,'(i8)') nzs
    allocate(zs(nzs))
    do i = 1, nzs
       read(10,'(e14.6)') zs(i)
    end do
    close(10)
    ! READ AGE BINS---------------------------------------------------------
    open(unit=10,file=fAges,status='old',form='formatted')
    read(10,'(i8)') nAges
    allocate(ages(nAges))
    do i = 1, nAges
       read(10,'(e14.6)') ages(i)
    end do
    close(10)
    ages = ages*1.e-9                       !         Convert from yr to Gyr
    if(ages(1) .ne. 0.) ages(1) = 0.
    ! READ SEDS-------------------------------------------------------------
    open(unit=10,file=fSEDs,status='old',form='unformatted')
    read(10) nLs, dum
    allocate(Ls(nLs))
    read(10) Ls(:)
    allocate(SEDs(nLs,nAges,nzs))
    do iz = 1, nzs
       do ia = 1, nAges
          read(10) SEDs(:,ia,iz)
       end do
    end do
    close(10)

    ! Perform SED integration of luminosity, csn and egy per (age,Z) bin----
    allocate(tbl(nAges,nZs,nv))
    do ip = 1,nSEDgroups_tot                                ! Loop photon groups
       tbl=0.
       pL0 = group_decomp(ip) ; pL1 = group_decomp(ip+1)! eV interval of photon group ip
       do iz = 1, nzs                                     ! Loop metallicity
          do ia = 1,nAges                                ! Loop age
             tbl(ia,iz,1) = getSEDLuminosity(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
             tbl(ia,iz,3) = getSEDEgy(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
             do ii = 1,nIons
                tbl(ia,iz,2+ii*2) = getSEDcsn(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
                tbl(ia,iz,3+ii*2) = getSEDcse(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
             end do
          end do ! End age loop
       end do ! End Z loop

       dlgA = SED_dlgA ; SED_dlgZ = -SED_nz
       call rebin_log(dlgA, SED_dlgZ                                       &
            , tbl(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, nv        &
            , reb_tbl, SED_nA, SED_nZ, rebAges, SED_Zeds)
       SED_nA=SED_nA+1                              ! Make room for zero age
       if(ip .eq. 1 ) allocate(SED_table(SED_nA, SED_nZ, nSEDgroups_tot, nv))
       SED_table(1, :,ip,:) = reb_tbl(1,:,:)            ! Zero age properties
       SED_table(1, :,ip,2) = 0.                        !  Lacc=0 at zero age
       SED_table(2:,:,ip,:) = reb_tbl
       deallocate(reb_tbl)

       if(ip==1)then
          SED_lgZ0 = log10(SED_Zeds(1))                  ! Interpolation intervals
          SED_lgA0 = log10(rebAges(1))
          allocate(SED_ages(SED_nA))
          SED_ages(1)=0.d0 ; SED_ages(2:)=rebAges ;    ! Must have zero initial age
       endif

       ! Integrate the cumulative luminosities:
       SED_table(:,:,ip,2)=0.d0
       do iz = 1, SED_nZ ! Loop metallicity
          tmp = trapz1( SED_ages, SED_table(:,iz,ip,1), SED_nA, SED_table(:,iz,ip,2) )
          SED_table(:,iz,ip,2) = SED_table(:,iz,ip,2) * Gyr2sec
       end do

    end do ! End photon group loop

    if(verbose) print*, 'Initialization of SED_table : done !'

    deallocate(SEDs) ; deallocate(tbl)
    deallocate(ages) ; deallocate(rebAges) 
    deallocate(zs)
    deallocate(Ls)

  end subroutine init_SED_table
  !#################################################################################################

  !*************************************************************************
  function get_csn_indices(n_elements, elements, nIons)

    use spectrum_integrator_module

    implicit none

    integer(kind=4),intent(in)  :: n_elements, elements(n_elements), nIons(n_elements)
    integer(kind=4),allocatable :: get_csn_indices(:)

    get_csn_indices = get_prop(n_elements, elements, nIons)   !In spectrum_integrator_module

  end function get_csn_indices
  !*************************************************************************


  !*************************************************************************
  SUBROUTINE inp_SED_table(age, Z, nProp, same, ret)!, na, nz, SED_ages, SED_zeds)

    ! Compute SED property by interpolation from table.
    ! input/output:
    ! age   => Star population age [Gyrs]
    ! Z     => Star population metallicity [m_metals/m_tot]
    ! nprop => Number of property to fetch
    !          1=log(photon # intensity [# Msun-1 s-1]),
    !          2=log(cumulative photon # intensity [# Msun-1]),
    !          3=avg_egy, 2+2*iIon=avg_csn, 3+2*iIon=avg_cse
    ! same  => If true then assume same age and Z as used in last call.
    !          In this case the interpolation indexes can be recycled.
    ! ret   => The interpolated values of the sed property for every photon
    !          group
    !-------------------------------------------------------------------------
    use module_constants
    real(kind=8), intent(in):: age, Z
    !integer(kind=4), intent(in):: na, nz
    !real(kind=8), intent(in):: SED_ages(na), SED_zeds(nz)
    real(kind=8):: lgAge, lgZ
    integer:: nProp
    logical:: same
    real(kind=8),dimension(:):: ret
    integer,save:: ia, iz
    real(kind=8),save:: da, da0, da1, dz, dz0, dz1
    !-------------------------------------------------------------------------
    ! ia, iz: lower indexes: 0<ia<sed_nA etc.
    ! da0, da1, dz0, dz1: proportional distances from edges:
    ! 0<=da0<=1, 0<=da1<=1 etc.

    if(.not. same) then
       if(age.le.0d0) then
          lgAge=-4d0
       else
          lgAge = log10(age)
       endif

       lgZ=log10(Z)
       ia = min(max(floor((lgAge-SED_lgA0)/SED_dlgA ) + 2, 1  ),  SED_nA-1 )
       da = SED_ages(ia+1)-SED_ages(ia)
       da0= min( max(   (age-SED_ages(ia)) /da,       0. ), 1.          )
       da1= min( max(  (SED_ages(ia+1)-age)/da,       0. ), 1.          )

       iz = min(max(floor((lgZ-SED_lgZ0)/SED_dlgZ ) + 1,   1  ),  SED_nZ-1 )
       dz = sed_Zeds(iz+1)-SED_Zeds(iz)
       dz0= min( max(   (Z-SED_zeds(iz)) /dz,         0. ),  1.         )
       dz1= min( max(  (SED_Zeds(iz+1)-Z)/dz,         0. ),  1.         )

       if (abs(da0+da1-1.0d0) > 1.0d-5 .or. abs(dz0+dz1-1.0d0) > 1.0d-5) then
          write(*,*) 'Screwed up the sed interpolation ... '
          write(*,*) da0+da1,dz0+dz1
          stop
       end if
    endif

    ret = da0 * dz0 * SED_table(ia+1, iz+1, :, nProp) + &
         da1 * dz0 * SED_table(ia,   iz+1, :, nProp) + &
         da0 * dz1 * SED_table(ia+1, iz,   :, nProp) + &
         da1 * dz1 * SED_table(ia,   iz,   :, nProp)

  END SUBROUTINE inp_SED_table


  !#################################################################################################
  function get_nOptBins()

    implicit none

    integer(kind=4) :: get_nOptBins
 
    get_nOptBins = nSEDgroups_tot

    return

  end function get_nOptBins
  !#################################################################################################
  

  !#################################################################################################
  subroutine deallocate_table()

    implicit none

    deallocate(SED_table, SED_ages, SED_zeds)

  end subroutine deallocate_table
  !#################################################################################################


  !#################################################################################################
  subroutine read_spectra_params(pfile)

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
       if (line(1:9) == '[spectra]') then
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
          case('sed_dir')
             write(sed_dir,'(a)') trim(value)
          case ('use_rt_restart')
             read(value,*) use_rt_restart
          case('nSEDgroups')
             read(value,*) nSEDgroups
          case('nSEDgroups_restart')
             read(value,*) nSEDgroups_restart
          case('group_decomp')
             if(use_rt_restart) then
                nSEDgroups_tot = nSEDgroups+nSEDgroups_restart
             else
                nSEDgroups_tot = nSEDgroups
             end if
             allocate(group_decomp(nSEDgroups_tot+1))
             read(value,*) group_decomp
          case ('verbose')
             read(value,*) verbose
          end select
       end do
    end if
    close(10)
    return
  end subroutine read_spectra_params
  !#################################################################################################


  !#################################################################################################
  subroutine print_spectra_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then
       write(*,'(a)')            '--------------------------------------------------------------------------------'
       write(unit,'(a,a,a)') '[spectra]'
       write(unit,'(a,a)')      '  sed_dir            = ', trim(sed_dir)
       write(unit,'(a,L1)')     '  use_rt_restart     = ', use_rt_restart
       write(unit,'(a,i5)')     '  nSEDgroups         = ', nSEDgroups
       if(use_rt_restart) then
          write(unit,'(a,i5)')     '  nSEDgroups_restart = ', nSEDgroups_restart
       end if
       write(unit,2000)         '  group_decomp       = ', group_decomp
       write(unit,'(a,L1)')     '  verbose            = ', verbose
       write(*,'(a)')            '--------------------------------------------------------------------------------'
    else
       write(*,'(a)')            '--------------------------------------------------------------------------------'
       write(*,'(a,a,a)') '[spectra]'
       write(*,'(a,a)')      '  sed_dir               = ', trim(sed_dir)
       write(*,'(a,L1)')     '  use_rt_restart        = ', use_rt_restart
       write(*,'(a,i5)')     '  nSEDgroups            = ', nSEDgroups
       if(use_rt_restart) then
          write(*,'(a,i5)')     '  nSEDgroups_restart    = ', nSEDgroups_restart
       end if
       write(*,2000)         '  group_decomp          = ', group_decomp
       write(*,'(a,L1)')     '  verbose               = ', verbose
       write(*,'(a)')            '--------------------------------------------------------------------------------'

    end if

    2000 format (a,1000(ES10.3,1x))

    return
  end subroutine print_spectra_params
  !#################################################################################################




  !*************************************************************************
  ! START PRIVATE SUBROUTINES AND FUNCTIONS*********************************

  !*************************************************************************
  FUNCTION getSEDLuminosity(X, Y, N, e0, e1)

    ! Compute and return luminosity in energy interval (e0,e1) [eV]
    ! in SED Y(X). Assumes X is in Angstroms and Y in Lo/Angstroms/Msun.
    ! (Lo=[Lo_sun], Lo_sun=[erg s-1]. total solar luminosity is
    ! Lo_sun=10^33.58 erg/s)
    ! returns: Photon luminosity in, [# s-1 Msun-1],
    !                             or [eV s-1 Msun-1] if SED_isEgy=true
    !-------------------------------------------------------------------------
    use module_constants,only:clight, planck, evtoerg
    use spectrum_integrator_module
    real(kind=8):: getSEDLuminosity, X(n), Y(n), e0, e1
    integer :: N, species
    real(kind=8),parameter :: const=1.0e-8/planck/clight
    real(kind=8),parameter :: Lsun=3.8256d33 ! Solar luminosity [erg/s]
    ! const is a div by ph energy => ph count.  1e-8 is a conversion into
    ! cgs, since wly=[angstrom] h=[erg s-1], c=[cm s-1]
    !-------------------------------------------------------------------------
    species          = 1                   ! irrelevant but must be included
    !  Photon number per sec per Msun
    getSEDLuminosity = const &
         * integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
    getSEDLuminosity = getSEDLuminosity*Lsun  ! Scale by solar luminosity

  END FUNCTION getSEDLuminosity

  !*************************************************************************
  FUNCTION getSEDEgy(X, Y, N, e0, e1)

    ! Compute average energy, in eV, in energy interval (e0,e1) [eV] in SED
    ! Y(X). Assumes X is in Angstroms and Y is energy weight per angstrom
    ! (not photon count).
    !-------------------------------------------------------------------------
    use module_constants,only:clight, evtoerg, planck
    use spectrum_integrator_module
    real(kind=8):: getSEDEgy, X(N), Y(N), e0, e1, norm
    integer :: N,species
    real(kind=8),parameter :: const=1.d8*planck*clight/evtoerg! energy conversion
    !-------------------------------------------------------------------------
    species      = 1                       ! irrelevant but must be included
    norm         = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
    getSEDEgy    = const * &
         integrateSpectrum(X, Y, N, e0, e1, species, f1) / norm
  END FUNCTION getSEDEgy

  !*************************************************************************
  FUNCTION getSEDcsn(X, Y, N, e0, e1, species)

    ! Compute and return average photoionization
    ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
    ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
    ! angstrom (not photon #).
    ! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
    !-------------------------------------------------------------------------
    use spectrum_integrator_module
    use module_constants,only:ionEVs
    real(kind=8):: getSEDcsn, X(N), Y(N), e0, e1, norm
    integer :: N, species
    !-------------------------------------------------------------------------
    if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
       getSEDcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
    endif
    norm     = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
    getSEDcsn= integrateSpectrum(X, Y, N, e0, e1, species, fSigLambda)/norm
  END FUNCTION getSEDcsn

  !************************************************************************
  FUNCTION getSEDcse(X, Y, N, e0, e1, species)

    ! Compute and return average energy weighted photoionization
    ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
    ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
    ! angstrom (not photon #).
    ! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
    !-------------------------------------------------------------------------
    use spectrum_integrator_module
    use module_constants,only:ionEVs
    real(kind=8):: getSEDcse, X(N), Y(N), e0, e1, norm
    integer :: N, species
    !-------------------------------------------------------------------------
    if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
       getSEDcse=0. ; RETURN    ! [e0,e1] below ionization energy of species
    endif
    norm      = integrateSpectrum(X, Y, N, e0, e1, species, f1)
    getSEDcse = integrateSpectrum(X, Y, N, e0, e1, species, fSig) / norm
  END FUNCTION getSEDcse

  !*************************************************************************
  SUBROUTINE rebin_log(xint_log, yint_log,                                 &
       data,       nx,       ny,     x,     y,     nz,           &
       new_data, new_nx, new_ny, new_x, new_y      )

    ! Rebin the given 2d data into constant logarithmic intervals, using
    ! linear interpolation.
    ! xint_log,  => x and y intervals in the rebinned data. If negative,
    ! yint_log      these values represent the new number of bins.
    ! data       => The 2d data to be rebinned
    ! nx,ny      => Number of points in x and y in the original data
    ! nz         => Number of values in the data
    ! x,y        => x and y values for the data
    ! new_data  <=  The rebinned 2d data
    ! new_nx,ny <=  Number of points in x and y in the rebinned data
    ! new_x,y   <=  x and y point values for the rebinned data
    !-------------------------------------------------------------------------
    real(kind=8):: xint_log, yint_log
    integer,intent(in):: nx, ny, nz
    integer::new_nx, new_ny
    real(kind=8),intent(in):: x(nx),y(ny)
    real(kind=8),intent(in):: data(nx,ny,nz)
    real(kind=8),dimension(:,:,:),allocatable:: new_data
    real(kind=8),dimension(:),allocatable:: new_x, new_lgx, new_y, new_lgy
    real(kind=8):: dx0, dx1, dy0, dy1, x_step, y_step
    real(kind=8):: x0lg, x1lg, y0lg, y1lg
    integer :: i, j, ix, iy, ix1, iy1
    !-------------------------------------------------------------------------
    if(allocated(new_x)) deallocate(new_x)
    if(allocated(new_y)) deallocate(new_y)
    if(allocated(new_data)) deallocate(new_data)

    ! Find dimensions of the new_data and initialize it
    x0lg = log10(x(1));   x1lg = log10(x(nx))
    y0lg = log10(y(1));   y1lg = log10(y(ny))

    if(xint_log .lt. 0 .and. nx .gt. 1) then
       new_nx=int(-xint_log)                        ! xint represents wanted
       xint_log = (x1lg-x0lg)/(new_nx-1)            !     number of new bins
    else
       new_nx =int((x1lg-x0lg)/xint_log) + 1
    endif
    allocate(new_x(new_nx)) ; allocate(new_lgx(new_nx))
    do i = 0, new_nx-1                              !  initialize the x-axis
       new_lgx(i+1) = x0lg + i*xint_log
    end do
    new_x=10.d0**new_lgx

    if(yint_log .lt. 0 .and. ny .gt. 1) then        ! yint represents wanted
       new_ny=int(-yint_log)                        !     number of new bins
       yint_log = (y1lg-y0lg)/(new_ny-1)
    else
       new_ny = int((y1lg-y0lg)/yint_log) + 1
    endif
    allocate(new_y(new_ny)) ; allocate(new_lgy(new_ny))
    do j = 0, new_ny-1                              !      ...and the y-axis
       new_lgy(j+1) = y0lg + j*yint_log
    end do
    new_y=10.d0**new_lgy

    ! Initialize new_data and find values for each point in it
    allocate(new_data(new_nx, new_ny, nz))
    do j = 1, new_ny
       call locate(y, ny, new_y(j), iy)
       ! y(iy) <= new_y(j) <= y(iy+1)
       ! iy is lower bound, so it can be zero but not larger than ny
       if(iy < 1) iy=1
       if (iy < ny) then
          iy1  = iy + 1
          y_step = y(iy1) - y(iy)
          dy0  = max(new_y(j) - y(iy),    0.0d0)  / y_step
          dy1  = min(y(iy1)   - new_y(j), y_step) / y_step
       else
          iy1  = iy
          dy0  = 0.0d0 ;  dy1  = 1.0d0
       end if

       do i = 1, new_nx
          call locate(x, nx, new_x(i), ix)
          if(ix < 1) ix=1
          if (ix < nx) then
             ix1  = ix+1
             x_step = x(ix1)-x(ix)
             dx0  = max(new_x(i) - x(ix),    0.0d0)  / x_step
             dx1  = min(x(ix1)   - new_x(i), x_step) / x_step
          else
             ix1  = ix
             dx0  = 0.0d0 ;  dx1  = 1.0d0
          end if

          if (abs(dx0+dx1-1.0d0) .gt. 1.0d-6 .or.                          &
               abs(dy0+dy1-1.0d0) > 1.0d-6) then
             write(*,*) 'Screwed up the rebin interpolation ... '
             write(*,*) dx0+dx1,dy0+dy1
             stop
          end if

          new_data(i,j,:) =                                                &
               dx0 * dy0 * data(ix1,iy1,:) + dx1 * dy0 * data(ix, iy1,:) + &
               dx0 * dy1 * data(ix1,iy, :) + dx1 * dy1 * data(ix, iy, :)
       end do
    end do

    deallocate(new_lgx)
    deallocate(new_lgy)

  END SUBROUTINE rebin_log



  SUBROUTINE locate(xx,n,x,j)
    ! Locates position j of a value x in an ordered array xx of n elements
    ! After: xx(j) <= x <= xx(j+1) (assuming increasing order)
    ! j is lower bound, so it can be zero but not larger than n
    !-------------------------------------------------------------------------
    !use amr_commons,only:kind=8
    integer ::  n,j,jl,ju,jm
    real(kind=8)::  xx(n),x
    !-------------------------------------------------------------------------
    jl = 0
    ju = n+1
    do while (ju-jl > 1) 
       jm = (ju+jl)/2
       if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    j = jl
  END SUBROUTINE locate


end module module_spectra
