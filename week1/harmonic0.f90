!---------------------------------------------------------------
program harmonic0
!---------------------------------------------------------------
  !
  !       Solution of the quantum harmonic oscillator
  !       Forward integration only, Numerov algorithm
  !       Eigenvalue search using the shooting method.
  !       Adimensional units: x = (mK/hbar^2)^(1/4) X
  !                           e = E/(hbar omega)
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  !
  integer :: mesh, i, icl
  integer :: nodes, hnodes, ncross, kkk, iterate
  real(dp) :: xmax, dx, ddx12, norm, arg
  real(dp) :: eup, elw, e
  real(dp), allocatable :: x(:), y(:), p(:), vpot(:), f(:)
  character (len=80) :: fileout
  !
  ! read input data
  !
  print '(a,$)','Max value for x (typical value: 10) ? '
  read (*,*) xmax
  print '(a,$)', 'Number of grid points (typically a few hundreds) ? '
  read (*,*) mesh
  !
  ! allocate arrays, initialize grid
  !
  allocate ( x(0:mesh), y(0:mesh), p(0:mesh), vpot(0:mesh), f(0:mesh) )
  dx =  xmax/mesh 
  ddx12=dx*dx/12.0_dp
  !
  ! set up the potential (must be even w.r.t. x=0)
  !
  do i = 0, mesh
     x(i) = float(i) * dx
     vpot(i) = 0.5_dp * x(i)*x(i)
  end do
  !
  print '(a,$)','output file name = '
  read (*,'(a)') fileout
  open (7,file=fileout,status='unknown',form='formatted')

999 continue ! this is the entry point for a new eigenvalue search
  !
  ! read number of nodes (stop if < 0)
  !         
  print '(a,$)', ' nodes (type -1 to stop) >> '
  read (*,*) nodes
  if (nodes < 0) then
     close(7)
     deallocate ( f, vpot, p, y, x )
     stop 
  end if
  !
  ! set initial lower and upper bounds to the eigenvalue
  !
  eup=maxval (vpot(:))
  elw=minval (vpot(:))
  !
  ! Set trial energy
  !
  print '(a,$)','Trial energy (0=search with bisection) ? '
  read (*,*) e
  if ( e == 0.0_dp ) then
     ! search eigenvalues with bisection
     e = 0.5_dp * (elw + eup)
     iterate = 1
  else
     ! test a single energy value
     iterate = 0
  endif
  
  kkk = 0
1 continue ! this is the entry point for the solution at fixed energy
  kkk = kkk + 1
  !
  ! set up the f-function used by the Numerov algorithm
  ! and determine the position of its last crossing, i.e. change of sign
  ! f < 0 means classically allowed   region
  ! f > 0 means classically forbidden region
  !
  f(0)=ddx12*(2.0_dp*(vpot(0)-e))
  icl=-1
  do i=1,mesh
     f(i)=ddx12*2.0_dp*(vpot(i)-e)
     ! beware: if f(i) is exactly zero the change of sign is not observed
     ! the following line is a trick to prevent missing a change of sign 
     ! in this unlikely but not impossible case:
     if ( f(i) == 0.0_dp) f(i)=1.d-20
     ! store the index 'icl' where the last change of sign has been found
     if ( f(i) /= sign(f(i),f(i-1)) ) icl=i
  end do
  
  if (icl >= mesh-2) then
     deallocate ( f, vpot, p, y, x )
     stop 'last change of sign too far'
  else if (icl < 1) then
     deallocate ( f, vpot, p, y, x )
     stop 'no classical turning point?'
  end if
  !
  ! f(x) as required by the Numerov algorithm
  !
  f = 1.0_dp - f
  y = 0.0_dp
  !
  ! determination of the wave-function in the first two points 
  !
  hnodes = nodes/2
  !
  ! beware the integer division: 1/2 = 0 !
  ! if nodes is even, there are 2*hnodes nodes
  ! if nodes is odd,  there are 2*hnodes+1 nodes (one is in x=0)
  ! hnodes is thus the number of nodes in the x>0 semi-axis (x=0 excepted)
  !
  if (2*hnodes == nodes) then
     ! even number of nodes: wavefunction is even
     y(0) = 1.0_dp
     ! assume f(-1) = f(1)
     y(1) = 0.5_dp*(12.0_dp-10.0_dp*f(0))*y(0)/f(1)
  else
     ! odd number of nodes: wavefunction is odd
     y(0) = 0.0_dp
     y(1) = dx
  end if
  !
  ! outward integration and count number of crossings
  !
  ncross=0
  do i =1,mesh-1
     y(i+1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
     if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
  end do
  !
  ! check number of crossings
  !
  print *, kkk, e, ncross, hnodes
  if ( iterate /= 0 ) then
     if (ncross > hnodes) then
        ! Too many crossings: current energy is too high, lower the upper bound
        eup = e
     else 
        ! Too few or correct number of crossings: current energy is too low,
        ! raise the lower bound
        elw = e
     end if
     ! New trial value:
     e = 0.5_dp * (eup+elw)
     ! Convergence criterion:
     if (eup-elw > 1.d-10) go to 1
  end if
  !
  ! ---- convergence has been achieved (or it wasn't required) -----
  ! Note that the wavefunction is not normalized: 
  ! the problem is the divergence at large |x| 
  !
  ! Calculation of the classical probability density for energy e:
  !
  norm = 0.0_dp
  p(icl:) = 0.0_dp
  do i=0,icl
     arg = (e - x(i)**2/2.0_dp)
     if ( arg > 0.0_dp) then
        p(i) = 1.0_dp/sqrt(arg)
     else
        p(i) = 0.0_dp
     end if
     norm = norm + 2.0_dp*dx*p(i)
  enddo
! The point at x=0 must be counted once:
  norm = norm - dx*p(0)
! Normalize p(x) so that  Int p(x)dx = 1
  p(:icl-1) = p(:icl-1)/norm
! x<0 region:
  do i=mesh,1,-1
     write (7,'(f7.3,3e16.8,f12.6)') & 
          -x(i), (-1)**nodes*y(i), y(i)*y(i), p(i), vpot(i)
  enddo
  ! x>0 region:
  do i=0,mesh
     write (7,'(f7.3,3e16.8,f12.6)') &
          x(i), y(i), y(i)*y(i), p(i), vpot(i)
  enddo

  go to 999

end program harmonic0
