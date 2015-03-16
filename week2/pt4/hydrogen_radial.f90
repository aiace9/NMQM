!---------------------------------------------------------------
program hydrogen
  !---------------------------------------------------------------
  !
  ! Find a solution with given N, L for an hydrogenic atom
  ! solving the radial schroedinger equation by Numerov method
  ! Atomic (Ry) units
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  integer mesh
  integer n, l, i
  real (dp) ::  zeta, zmesh, rmax, xmin, dx, e
  real (dp), allocatable :: r(:), sqr(:), r2(:), y(:), vpot(:)
  !
  ! initialize atomic charge (Z)
  !
  write (*,'(" Atomic Charge > ",$)')
  read (*,*) zeta
  !write (*,'(" Atomic charge = ",f10.6)') zeta
  if ( zeta < 1.0_dp) stop 'zeta should be >= 1'
  !
  ! initialize logarithmic mesh
  !
  zmesh=  zeta 
  rmax = 100.0_dp
  xmin = -8.0_dp
  dx   =  0.01_dp
  !
  mesh = (log(zmesh*rmax) - xmin)/dx
  !
  allocate ( r(0:mesh), sqr(0:mesh), r2(0:mesh), vpot(0:mesh), y(0:mesh) )
  !
  call do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
  !
  ! initialize the potential
  !
  call init_pot ( zeta, r, mesh, vpot)
  !
  ! open output file that will contain the wavefunctions
  !
  open (7,file='wfc.out',status='unknown',form='formatted')
  !
  ! read number of nodes (stop if nodes < 0)
  !
1 continue
  write (*,'(a,$)') ' n,l >> '
  read (*,*) n,l
  if (n < 1) then
     deallocate (y, vpot, r2, sqr, r)
     stop
  else if (n < l+1) then
     write(*,*) 'error in main: n < l+1 -> wrong number of nodes '
     deallocate (y, vpot, r2, sqr, r)
     stop
  else if (l < 0) then
     write(*,*) 'error in main: l < 0 unphysical '
     deallocate (y, vpot, r2, sqr, r)
     stop
  end if
  !
  ! solve the schroedinger equation in radial coordinates by Numerov method
  !
  call solve_sheq (n, l, e, mesh, dx, r, sqr, r2, vpot, zeta,y )
  !
  ! write out the eigenvalue energy to be compared with the external potential
  !
  write (*,'(" eigenvalue =",f15.8,",  eig*(n/zeta)^2 =",f15.8)') &
     e, e*(n/zeta)**2
  write(7,"('#       r             R(r)          X(r)=rR(r)          e')")
  do i=0,mesh
     write (7,*) r(i),y(i)/sqr(i), y(i)*sqr(i), e
  end do
  write (7,'(/)')
  
  go to 1
  
end program hydrogen
!
!--------------------------------------------------------------------
subroutine do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
!--------------------------------------------------------------------
  !
  ! initialize radial grid
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent (in) :: mesh
  real (dp), intent (in) :: zmesh, xmin, dx
  real (dp), intent (out) :: r(0:mesh), sqr(0:mesh), r2(0:mesh)
  !
  integer :: i
  real(dp) :: x
  !
  do i=0,mesh
     x = xmin + dx * i
     r(i)  = exp(x)/zmesh
     sqr(i)= sqrt(r(i))
     r2(i) = r(i) * r(i)
  end do
  write(*,'(/" radial grid information:")')
  write(*,'(" dx =",f10.6,", xmin =",f10.6,", zmesh =",f10.6)') &
       dx,xmin,zmesh
  write(*,'(" mesh =",i6,", r(0) =",f10.6,", r(mesh) =",f10.6)') &
       mesh, r(0), r(mesh)
  write(*,*) 
  !
  return
end subroutine do_mesh
!---------------------------------------------------------------------
subroutine solve_sheq ( n, l, e, mesh, dx, r, sqr, r2, vpot, zeta, y )
  !---------------------------------------------------------------------
  !
  ! solve the schroedinger equation in radial coordinates on a 
  ! logarithmic grid by Numerov method - atomic (Ry) units
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200), maxiter=10000
  real(dp), parameter :: eps=1.0D-10
  integer, intent(in) :: mesh, n,l
  real(dp), intent(in) :: dx, r(0:mesh), sqr(0:mesh), r2(0:mesh), & 
       vpot(0:mesh), zeta
  real(dp), intent(out) :: e, y(0:mesh)
  
  integer :: i, j, icl, nodes, ncross, kkk
  real (dp) :: ddx12, sqlhf, x2l2, ycusp, dfcusp, fac, norm, eup, elw, de
  real (dp), allocatable :: f(:)
  
  allocate ( f(0:mesh) )
  ddx12=dx*dx/12.0_dp
  sqlhf = (l+0.5_dp)**2
  x2l2 = 2*l+2
  !
  ! set (very rough) initial lower and upper bounds to the eigenvalue
  !
  eup = vpot(mesh)
  elw=minval (sqlhf/r2(:)+vpot(:))
  if (eup-elw < eps) then
     write (*,*) elw, eup
     write (*,*) 'solve_sheq: lower and upper bounds are equal'
     stop
  end if
  e = 0.5_dp * (elw + eup)

  do kkk = 1, maxiter
     !
     ! set up the f-function and determine the position of its last
     ! change of sign
     ! f < 0 (approximately) means classically allowed   region
     ! f > 0         "         "        "      forbidden   "
     !
     icl = -1
     f(0) = ddx12 *( sqlhf + r2(0) * (vpot(0)-e) )
     do i=1,mesh
        f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
        !
        ! beware: if f(i) is exactly zero the change of sign is not observed
        ! the following line is a trick to prevent missing a change of sign 
        ! in this unlikely but not impossible case:
        !
        if( f(i) == 0.0_dp ) f(i)=1.d-20
        if( f(i) /= sign(f(i),f(i-1)) ) icl=i
     end do

     if ( icl < 0 .or. icl >= mesh-2 ) then
        !
        ! classical turning point not found or too far away
        ! no panic: it may follow from a bad choice of eup in
        ! the first iterations. Update e and eup and re-try
        eup = e
        e = 0.5_dp * (eup+elw)
        cycle
     end if
     !
     ! f function as required by numerov method
     !
     f = 1.0_dp - f
     y = 0
     !
     ! determination of the wave-function in the first two points 
     ! (asymptotic behaviour - second term depends upon the potential)
     !
     nodes = n - l - 1
     y(0) = r(0)**(l+1) *(1.0_dp - 2.0_dp*zeta*r(0)/x2l2) / sqr(0)
     y(1) = r(1)**(l+1) *(1.0_dp - 2.0_dp*zeta*r(1)/x2l2) / sqr(1)
     !
     ! outward integration, count number of crossings
     !
     ncross=0
     do i =1, icl-1
        y(i+1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
        if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
     end do
     fac = y(icl) 
     !
     ! check number of crossings
     !
     if ( ncross /= nodes ) then
        if ( ncross > nodes ) then
           eup = e
        else
           elw = e
        end if
        e = 0.5_dp * (eup+elw)
        cycle
     end if
     !
     ! determination of the wave-function in the last two points 
     ! assuming y(mesh+1) = 0 and y(mesh) = dx
     !
     y(mesh) = dx
     y(mesh-1) = (12.0_dp-10.0_dp*f(mesh))*y(mesh)/f(mesh-1) 
     !
     ! inward integration 
     !
     do i = mesh-1,icl+1,-1
        y(i-1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i+1)*y(i+1))/f(i-1)
        if (y(i-1) > 1.0d10) then
           do j=mesh,i-1,-1
              y(j) = y(j)/y(i-1)
           end do
        end if
     end do
     !
     ! rescale function to match at the classical turning point (icl)
     !
     fac = fac/y(icl)
     y (icl:) = y(icl:)*fac
     !
     ! normalization - note the change of variable:
     !  \int f(r)dr => \sum_i f_i r_i Delta x
     !
     norm = 0.d0
     do i=1,mesh
        norm = norm + y(i)*y(i) * r2(i) * dx
     end do
     norm = sqrt(norm)
     y = y / norm
     !
     ! find the value of the cusp at the matching point (icl)
     !
     i = icl
     ycusp = (y(i-1)*f(i-1)+f(i+1)*y(i+1)+10.0_dp*f(i)*y(i)) / 12.0_dp
     dfcusp = f(i)*(y(i)/ycusp - 1.0_dp)
     !
     ! eigenvalue update using perturbation theory
     !
     de = dfcusp/ddx12 * ycusp*ycusp * dx 
     if (de > 0.0_dp) elw = e
     if (de < 0.0_dp) eup = e
     !
     ! prevent e to go out of bounds, i.e. e > eup or e < elw 
     ! (might happen far from convergence)
     !
     e = max( min (e+de,eup),elw)
     !
     ! convergence check
     !
     if (abs(de) < eps) exit
     !
  end do
  !
  ! was convergence achieved ?
  !
  if ( abs(de) > 1.d-10 ) then
     if ( ncross /= nodes ) then
        write (*,*) e, elw, eup, ncross, nodes, icl
     else
        write (*,*) e, de
     end if
     write (*,*) ' error in solve_sheq: too many iterations'
     stop
  else 
  ! ---- convergence has been achieved -----
     write (*,'(" convergence achieved at iter #",i3," de = ", e10.4)') kkk,de
  end if
  return
end subroutine solve_sheq
!
!--------------------------------------------------------------------
subroutine init_pot ( zeta, r, mesh, vpot )
!--------------------------------------------------------------------
  !
  ! initialize potential
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent (in) :: mesh
  real (dp), intent(in) :: zeta, r(0:mesh)
  real (dp), intent(out):: vpot(0:mesh)
  real (dp), parameter :: qe = 1
  real (dp), parameter :: Q = 0.05
  integer :: i

  open (7,file='pot.out',status='unknown',form='formatted')
  do i =0,mesh
     vpot (i) = - 2.0_dp*zeta/r(i) - zeta*(qe**2)*exp(-Q*r(i))/r(i)
     write (7,*) r(i),vpot(i)
  end do
  close(7)
  
  return
end subroutine init_pot
