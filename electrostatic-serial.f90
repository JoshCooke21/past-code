program electrostatic
  !use mpi
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)

  real(dp), dimension(:,:), allocatable :: phi

  real(dp) :: r_outer, r_inner, l_outer, l_inner, phi_0, phi_old, phi_change, w_parameter, k_step, tol
  integer :: i, k, Nr, Nz, istat=0, rinner_coord, linner_coord, loop_counter
  integer :: A_rcoord, A_zcoord, B_rcoord, B_zcoord, C_rcoord, C_zcoord
  logical :: exit_condition

  Nr=200
  Nz=400

  w_parameter = 1.5_dp
  tol=0.000001_dp

  r_inner=1._dp
  r_outer=10._dp

  l_inner=5._dp
  l_outer=20._dp

  phi_0=1000._dp

  !create the grid
  allocate(phi(0:Nz+1,0:Nr+1), stat=istat)
  if (istat.ne.0) stop 'Error allocating phi'

  !set the inner conductor to phi_0
  rinner_coord = floor(real(Nr,dp)*r_inner/r_outer)
  k_step = r_outer / real(Nr,dp)

  linner_coord = floor(real(Nz,dp)*l_inner/l_outer)

  phi(:,:) = 0._dp
  phi(:linner_coord,:rinner_coord)=phi_0

  do k=rinner_coord+1,Nr !Along z=0
    phi(0,k) = rAxisPot( phi_0, r_outer, r_inner, k*k_step )
  end do

  do loop_counter = 1, 1000000
    exit_condition = .true.

    !loop through the grid, updating each point (ignoring i=0, which doesn't change)
    ! loop over r=0
    do i=linner_coord+1,Nz

      phi_old    = phi(i,0)
      phi_change = zAxisPot( phi(i+1,0), phi(i-1,0), phi(i,1) )
      phi(i,0)   = phi_old + w_parameter * (phi_change - phi_old)

      !If the % change in phi is OVER the chosen tolerance, then don't exit
      if (abs((phi(i,0)-phi_old)/phi_old).gt.tol) exit_condition = .false.
    end do

    !loop over region A
    do i=1,linner_coord
      do k=rinner_coord+1,Nr

        phi_old    = phi(i,k)
        phi_change = normalPot( phi(i+1,k), phi(i-1,k), phi(i,k+1), phi(i,k-1), real(k,dp) )
        phi(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

        if (abs((phi(i,k)-phi_old)/phi_old).gt.tol) exit_condition = .false.
      end do
    end do

    !loop over region B
    do i=linner_coord+1,Nz
      do k=1,Nr

        phi_old    = phi(i,k)
        phi_change = normalPot( phi(i+1,k), phi(i-1,k), phi(i,k+1), phi(i,k-1), real(k,dp) )
        phi(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

        if (abs((phi(i,k)-phi_old)/phi_old).gt.tol) exit_condition = .false.
      end do
    end do

    !This only triggers if NO changes in phi have been over the tolerance
    if (exit_condition) exit

  end do

  print*, 'num loops:', loop_counter

  A_rcoord = 0
  B_rcoord = floor(real(Nr,dp)*6._dp/r_outer)
  C_rcoord = floor(real(Nr,dp)*5._dp/r_outer)

  A_zcoord = floor(real(Nz,dp)*7.5_dp/l_outer)
  B_zcoord = floor(real(Nz,dp)*5._dp/l_outer)
  C_zcoord = floor(real(Nz,dp)*12.5_dp/l_outer)

  print*, 'Potential at A = ', phi(A_zcoord,A_rcoord)
  print*, 'Potential at B = ', phi(B_zcoord,B_rcoord)
  print*, 'Potential at C = ', phi(C_zcoord,C_rcoord)

  phi = transpose(phi)
  call write_pgm(phi)

  deallocate(phi, stat=istat)
  if (istat.ne.0) stop 'Error deallocating phi'

contains
  function normalPot(phi_up, phi_down, phi_right, phi_left, k)
    implicit none
    real(dp), intent(in) :: phi_up, phi_down, phi_right, phi_left, k
    real(dp) :: normalPot

    normalPot = 0.25_dp * (phi_up + phi_down + phi_right + phi_left) + (phi_right - phi_left)/(8._dp*k)

  end function normalPot

  function zAxisPot(phi_up, phi_down, phi_right)
    implicit none
    real(dp), intent(in) :: phi_up, phi_down, phi_right
    real(dp) :: zAxisPot

    zAxisPot = 2._dp * phi_right/3._dp + (phi_up + phi_down)/6._dp

  end function zAxisPot

  function rAxisPot(phi_0, r_outer, r_inner, r_current)
    implicit none
    real(dp), intent(in) :: phi_0, r_outer, r_inner, r_current
    real(dp) :: rAxisPot

    rAxisPot = phi_0 * (log(r_outer) - log(r_current))/(log(r_outer) - log(r_inner))

  end function rAxisPot

  subroutine write_pgm(grid)
    implicit none

    real(dp), dimension(:,:), intent(in) :: grid
    integer, dimension(:,:), allocatable :: pixels

    real(dp) :: Pmax, Pmin
    integer  :: i, j, max_greys, Nx, Ny, out_unit, ierr=0

    print*, 'Writing to file'

    Nx=size(grid,1)
    Ny=size(grid,2)

    !set up the pgmfile
    allocate(pixels(Nx,Ny),stat=ierr)
    if (ierr/=0) stop 'Error in allocating pixels'

    max_greys=255
    Pmax=maxval(grid)
    Pmin=minval(grid)

    do i=1,Nx
      do j=1,Ny
        pixels(i,j)=int((grid(i,j)-Pmin)*max_greys/(Pmax-Pmin)) !Tmin<T<Tmax
      end do
    end do

    out_unit=10
    open(file='potential.pgm',unit=out_unit,status='unknown')
    write (out_unit,11) 'P2'                 !pgm magic number
    write (out_unit,12) Nx,Ny                !width, height
    write (out_unit,13) max_greys            !max gray value
    do j=1,Ny
      do i=1,Nx
        write (out_unit,14) pixels(i,j)
      end do
    end do
    close (unit=out_unit)

    11 format(a2)
    12 format(i10,1x,i10)
    13 format(i10)
    14 format (15(1x,i3))
  end subroutine

end program
