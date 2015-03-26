!---------------------------------------------------------------
program crossection
!---------------------------------------------------------------
  !
  !       Cross section for scattering of a particle by a potential
  !       Forward integration only, Numerov algorithm
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi = 3.141592653589792_dp
  !
  ! Units: energies in meV, distances in A
  ! hbar2m = hbar^2/2m (meV*A^2), with m = mass of the particle (proton)
  !
  real(dp), parameter :: hbar2m = 2.08
  !
  ! xmin = starting point for integration
  ! xmax = point at which the potential is negligible
  !
  real(dp) :: xmin, xmax
  ! energy range
  real(dp) :: elw, eup, de
  !
  integer :: lmax, mesh, i, l, i1, i2, n, ne
  real(dp) :: dx, ddx12, norm, e, q, kappa, lambda2, &
              rmin, rmax, r1, r2, tandelta, delta
  real(dp), external :: nl, jl
  real(dp), allocatable :: r(:), u(:), p(:), vpot(:), f(:), cross(:)
  character (len=80) :: fileout
  !
  integer :: ios
  ! read input data
  ! default values for LJ potential, sigma=3.57
  !xmin=1.8
  !xmax=18.0
  !
  !mesh = 2000
  !elw=0.1_dp
  !eup=3.5_dp
  !de=0.01_dp
  !fileout = 'mah'
  !lmax = 6
  !
  print '(a,$)', 'Number of grid points for radial integration > '
  read (*,*) mesh
  print '(a,$)','Radial integration starts at (A) > '
  read (*,*) xmin
  print '(a,$)','First matching point (A) > '
  read (*,*) xmax
  print '(a,$)', 'Max l for which phase shifts are calculated > '
  read (*,*) lmax
  print '(a,$)', 'Min and Max energy, Delta E (meV) >'
  read (*,*) elw, eup, de
  print '(a,$)', 'Output data written to file > '
  read (*,'(a)') fileout
  !
  open (unit=8,file=fileout,status='unknown',form='formatted')
  write (8,"('# mesh=',i4,', xmin=',f7.2,', lmax=',i1)") mesh, xmin, lmax
  write (8,"('#     E          sigma(E)     sigma_l(E), l=0, lmax')")

  open(unit=9, file='function.out', iostat=ios, status="unknown")
  if ( ios /= 0 ) stop "Error opening file energy.out"
  
  !
  ! allocate arrays
  !
  allocate ( r(0:mesh), u(0:mesh), p(0:mesh), vpot(0:mesh), f(0:mesh) )
  allocate ( cross(0:lmax) )
  !
  ! average half-wavelength on the energy range
  ! recalculating it at each energy produces some numerical noise
  !
  lambda2 = pi/sqrt( (eup+elw)/2.0_dp/hbar2m )
  !
  ne = nint((eup-elw)/de) + 1
  do n=1,ne
     !
     ! initialization
     !
     e = elw + (n-1)*de
     q = sqrt(e/hbar2m)
     !
     rmin  = xmin
     ! also add half a wavelength to allow matching to spherical waves
     rmax  = xmax + lambda2
     dx = (rmax-rmin)/mesh 
     ddx12=dx*dx/12.0_dp
     r2 = rmax
     i2 = mesh
     ! r1 is a second point, half a wavelength from r2=rmax, needed to
     ! calculate matching to spherical waves
     r1 = rmax - dx*nint(lambda2/dx)
     i1 = mesh - nint(lambda2/dx)
     !
     do i = 0, mesh
        r(i) = rmin + float(i) * dx
     end do
     !
     ! fill vpot with the potential
     !
     call potential ( mesh, r, vpot )
     !
     cross=0.0
     !
     do l=0,lmax
        !
        !
        ! set up the f-function used by the Numerov algorithm
        !
        do i=0,mesh
           f(i)=ddx12/hbar2m*(hbar2m*l*(l+1)/r(i)**2+vpot(i)-e)
        end do
        f = 1.0_dp - f
        !
        u = 0.0_dp
        !
        ! wave-function in the first two points 
        !
        call starting_points ( hbar2m, r(0), r(1), u(0), u(1) )
        !
        ! outward integration 
        !
        do i =1,mesh-1
           u(i+1)=((12.0_dp-10.0_dp*f(i))*u(i)-f(i-1)*u(i-1))/f(i+1)
        end do
        !
        ! normalization (simple integral)
        !
        norm = dot_product (u, u) * dx 
        u = u / sqrt(norm)
        !
        ! matching
        !
        kappa = r1*u(i2)/(r2*u(i1))
        tandelta = (kappa*jl(l,q*r1) - jl(l,q*r2) ) / &
                   (kappa*nl(l,q*r1) - nl(l,q*r2) ) 
        delta = atan(tandelta)
	      cross(l) = cross(l)+  4*pi/q**2 * (2*l+1)*sin(delta)**2
        !
        ! calculate the asymptotic wavefunction p
        !
        do i = 0, mesh
           p(i)= q*r(i) * ( jl(l,q*r(i))*cos(delta) - nl(l,q*r(i))*sin(delta) )
        end do
        ! normalize p so that it is equal to u for r=r2
        p = p / (p(i2)/u(i2) )
        !
        ! uncomment to write wavefunction, asymptotic wvf, eff.potential
        !
        if ( l==1 ) then
          do i=mesh/2,mesh
            !write (10+l,'(f10.4,3e16.8)') r(i), u(i), p(i), vpot(i)+hbar2m*l*(l+1)/r(i)**2
            write(unit=9, fmt=*, iostat=ios) r(i), u(i), p(i)
            if ( ios /= 0 ) stop "Write error in file unit 9"     
          enddo
        end if

        !
     enddo
     ! max l: 13, or else the lines will be wrapped
     write (8,'(2f12.6,13e11.3)') e, SUM(cross), cross(:)
     !
   end do
   close(unit=8)
   close(unit=9)
   !
end program crossection
!
! jl and nl functions - beware! unstable recurrence for large l
!--------------------------------------------------------------------
function jl (l, x )
  !--------------------------------------------------------------------
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: l
  real(dp) :: x, jl
  !
  real(dp) :: jm1, jp1
  integer :: lm
 
  jm1= cos (x) / x
  jl = sin (x) / x
  
  do lm = 0, l-1 
    jp1= (2*lm+1)/x*jl-jm1
    jm1= jl
    jl = jp1
  end do
  !
  return
end function jl
!
!--------------------------------------------------------------------
function nl (l, x )
  !--------------------------------------------------------------------
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: l
  real(dp) :: x, nl
  !
  real(dp) :: nm1, np1
  integer :: lm
 
  nm1= sin (x) / x
  nl =-cos (x) / x
  
  do lm = 0, l-1 
    np1= (2*lm+1)/x*nl-nm1
    nm1= nl
    nl = np1
  end do
  !
  return
end function nl
 !
subroutine potential ( mesh, r, v ) 
  ! 
  ! Lennard-Jones potential
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  !
  integer, intent(in) :: mesh
  real(dp), intent (in) :: r(0:mesh)
  real(dp), intent (out):: v(0:mesh)
  !
  ! epsilon and sigma are the parameters of Lennard-Jones potential
  ! epsilon in energy units (here:meV), sigma in length units (here: A)
  !
  real(dp), parameter :: epsilon = 5.9, sigma=3.57
  integer :: i

  do i = 0, mesh
     v(i) = epsilon * ( -2.0_dp*(sigma/r(i))**6 + (sigma/r(i))**12 )
  end do
  !
end subroutine potential
!
subroutine starting_points ( hbar2m, r0, r1, u0, u1 ) 
  ! 
  ! Starting wavefunction for Lennard-Jones potential
  ! assuming dominant r^12 term - not accurate, see Thjissen notes
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: epsilon = 5.9, sigma=3.57 ! as defined above
  !
  real(dp), intent (in) :: hbar2m, r0, r1
  real(dp), intent (out):: u0, u1
  !
  u0 = exp ( -sqrt( epsilon/hbar2m*sigma**12/25.0_dp ) * r0**(-5) )
  u1 = exp ( -sqrt( epsilon/hbar2m*sigma**12/25.0_dp ) * r1**(-5) )

end subroutine starting_points