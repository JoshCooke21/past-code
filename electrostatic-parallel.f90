program electrostatic
  use mpi
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)


  real(dp) :: r_outer, r_inner, l_outer, l_inner, phi_0, phi_old, phi_change, w_parameter, k_step, tol
  integer :: i, k, Nr, Nz, istat=0, rinner_coord, linner_coord, loop_counter
  integer :: A_rcoord, A_zcoord, B_rcoord, B_zcoord, C_rcoord, C_zcoord
  logical :: local_exitcond

  real(dp), dimension(:,:), allocatable :: phi_local, phi_gather, phi_temp, phi_final
  real(dp), dimension(:), allocatable :: send1, send2, recv1, recv2
  real(dp) :: time1, time2
  integer :: numranks, myrank, bloc_width, lbound, ubound, linner_rank, local_linner_coord
  integer :: send1_req, send2_req, rec1_req, rec2_req
  logical :: linner_check, global_exitcond

  integer, dimension(MPI_STATUS_SIZE)   :: status

  call MPI_INIT(istat)
  if (istat.ne.0) stop 'Error initializing MPI'

  call MPI_COMM_SIZE(MPI_COMM_WORLD, numranks, istat)
  if (istat.ne.0) stop 'Error getting numranks'

  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, istat)
  if (istat.ne.0) stop 'Error getting myrank'

  time1 = MPI_wtime()

  w_parameter = 1.5_dp
  tol=0.0000001_dp
  phi_0=1000._dp

  r_inner=1._dp
  r_outer=10._dp

  l_inner=5._dp
  l_outer=20._dp

  Nr=400
  Nz=800
  !Split along z (i)
  bloc_width = Nz / numranks

  rinner_coord = floor(real(Nr,dp)*r_inner/r_outer)
  k_step = r_outer / real(Nr,dp)

  linner_coord = floor(real(Nz,dp)*l_inner/l_outer)

  lbound = myrank * bloc_width + 1
  ubound = (myrank + 1) * bloc_width


  !Finds in which thread the z=5mm line is (top of the inner conductor and region A)
  linner_check = (lbound.le.linner_coord).and.(ubound.ge.linner_coord)
  if (linner_check) linner_rank = myrank

  !create the grid
  allocate(phi_local(0:bloc_width+1,0:Nr+1), stat=istat)
  if (istat.ne.0) stop 'Error allocating phi_local'

  !set the inner conductor to phi_0
  phi_local(:,:) = 0._dp
  if (myrank.eq.linner_rank) then
    local_linner_coord = mod(linner_coord-1,bloc_width) + 1
    phi_local(0:local_linner_coord, 0:rinner_coord) = phi_0
  else if (myrank.lt.linner_rank) then
    phi_local(:,0:rinner_coord) = phi_0
  end if

  !Set the perfect conductor case, which will only be on rank=0
  if (myrank.eq.0) then
    do k=rinner_coord+1,Nr !Along z=0
      phi_local(0,k) = rAxisPot( phi_0, r_outer, r_inner, k*k_step )
    end do
  end if

  !Allocate the arrays for sending/receiving
  allocate(send1(0:Nr+1), stat=istat)
  if (istat.ne.0) stop 'Error allocating send1'
  allocate(send2(0:Nr+1), stat=istat)
  if (istat.ne.0) stop 'Error allocating send2'
  allocate(recv1(0:Nr+1), stat=istat)
  if (istat.ne.0) stop 'Error allocating recv1'
  allocate(recv2(0:Nr+1), stat=istat)
  if (istat.ne.0) stop 'Error allocating recv2'

  do loop_counter = 1, 1000000
    !local_exitcond will later be reduced into a single variable
    local_exitcond = .true.

    send1(0:Nr+1) = phi_local(1         ,0:Nr+1)
    send2(0:Nr+1) = phi_local(bloc_width,0:Nr+1)


    !Send between cores
    if (myrank.ne.0) then
      call MPI_ISEND(send1, Nr+2, MPI_DOUBLE_PRECISION, myrank-1, 21, MPI_COMM_WORLD, send1_req, istat)
      if (istat.ne.0) stop 'Error sending down'
      call MPI_IRECV(recv1, Nr+2, MPI_DOUBLE_PRECISION, myrank-1, 22, MPI_COMM_WORLD, rec1_req, istat)
      if (istat.ne.0) stop 'Error receiving from down'
    end if
    if ( myrank.ne.(numranks-1) ) then
      call MPI_ISEND(send2, Nr+2, MPI_DOUBLE_PRECISION, myrank+1, 22, MPI_COMM_WORLD, send2_req, istat)
      if (istat.ne.0) stop 'Error sending up'
      call MPI_IRECV(recv2, Nr+2, MPI_DOUBLE_PRECISION, myrank+1, 21, MPI_COMM_WORLD, rec2_req, istat)
      if (istat.ne.0) stop 'Error receiving from up'
    end if

    !Latency hiding

    if (myrank.eq.linner_rank) then !Rank with z=5mm (can contain both regions A and B)

      !loop over region A, excluding bottom row that needs data from other threads
      do i=2,local_linner_coord
        do k=rinner_coord+1,Nr

          phi_old    = phi_local(i,k)
          phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
          phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

          !If the % change in phi is OVER the chosen tolerance, then don't exit
          if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
        end do
      end do


      do i=local_linner_coord+1,bloc_width-1

        ! loop over z=0
        phi_old    = phi_local(i,0)
        phi_change = zAxisPot( phi_local(i+1,0), phi_local(i-1,0), phi_local(i,1) )
        phi_local(i,0)   = phi_old + w_parameter * (phi_change - phi_old)

        if (abs((phi_local(i,0)-phi_old)/phi_old).gt.tol) local_exitcond = .false.

        !loop over region B, excluding top row that needs data from other threads
        do k=1,Nr

          phi_old    = phi_local(i,k)
          phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
          phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

          if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
        end do

      end do

    else if (myrank.lt.linner_rank) then !Ranks below z=5mm (region A)

      !loop over region A, excluding rows that need data from other threads
      do i=2,bloc_width-1
        do k=rinner_coord+1,Nr

          phi_old    = phi_local(i,k)
          phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
          phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

          !If the % change in phi is OVER the chosen tolerance, then don't exit
          if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
        end do
      end do

    else if (myrank.gt.linner_rank) then !Ranks above z=5mm (region B)

      do i=2,bloc_width-1

        ! loop over r=0
        phi_old    = phi_local(i,0)
        phi_change = zAxisPot( phi_local(i+1,0), phi_local(i-1,0), phi_local(i,1) )
        phi_local(i,0)   = phi_old + w_parameter * (phi_change - phi_old)

        if (abs((phi_local(i,0)-phi_old)/phi_old).gt.tol) local_exitcond = .false.

        !loop over region B
        do k=1,Nr

          phi_old    = phi_local(i,k)
          phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
          phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

          if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
        end do

      end do

    end if


    !Wait for send/receive
    if (myrank.ne.0) then
      call MPI_WAIT(send1_req, status, istat)
      if(istat.ne.0) stop 'Error waiting on send down'
      call MPI_WAIT(rec1_req, status, istat)
      if(istat.ne.0) stop 'Error waiting on receive from down'

      phi_local(0, 0:Nr+1) = recv1
    end if
    !On rank 0, the bottom row (z=0) doesn't change (as set before the looping)

    if ( myrank.ne.(numranks-1) ) then
      call MPI_WAIT(send2_req, status, istat)
      if(istat.ne.0) stop 'Error waiting on send up'
      call MPI_WAIT(rec2_req, status, istat)
      if(istat.ne.0) stop 'Error waiting on receive from up'

      phi_local(bloc_width+1, 0:Nr+1) = recv2
    end if
    !On maxrank, the top row is already set to 0.

    !do the bits that needed send/receive
    if (myrank.eq.linner_rank) then !Rank with z=5mm (can contain both regions A and B)

      !region A i=1, which was excluded from above
      i=1
      do k=rinner_coord+1,Nr

        phi_old    = phi_local(i,k)
        phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
        phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

        !If the % change in phi is OVER the chosen tolerance, then don't exit
        if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
      end do

      i=bloc_width
      ! r=0
      phi_old    = phi_local(i,0)
      phi_change = zAxisPot( phi_local(i+1,0), phi_local(i-1,0), phi_local(i,1) )
      phi_local(i,0)   = phi_old + w_parameter * (phi_change - phi_old)

      if (abs((phi_local(i,0)-phi_old)/phi_old).gt.tol) local_exitcond = .false.

      !loop over region B, but only top row that needed data from other threads
      do k=1,Nr

        phi_old    = phi_local(i,k)
        phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
        phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

        if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
      end do

    else if (myrank.lt.linner_rank) then !Ranks below z=5mm (region A)

      !loop over region A, but only the rows that need data from other threads (i=1 and bloc_width)
      do i=1,bloc_width, bloc_width-1 !Jump from i=1 to i=bloc_width
        do k=rinner_coord+1,Nr

          phi_old    = phi_local(i,k)
          phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
          phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

          !If the % change in phi is OVER the chosen tolerance, then don't exit
          if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
        end do
      end do

    else if (myrank.gt.linner_rank) then !Ranks above z=5mm (region B)

      do i=1,bloc_width, bloc_width-1 !Jump from i=1 to i=bloc_width
        ! loop over k=0
        phi_old    = phi_local(i,0)
        phi_change = zAxisPot( phi_local(i+1,0), phi_local(i-1,0), phi_local(i,1) )
        phi_local(i,0)   = phi_old + w_parameter * (phi_change - phi_old)

        if (abs((phi_local(i,0)-phi_old)/phi_old).gt.tol) local_exitcond = .false.

        !loop over region B, but only the rows that need data from other threads (i=1 and bloc_width)
        do k=1,Nr

          phi_old    = phi_local(i,k)
          phi_change = normalPot( phi_local(i+1,k), phi_local(i-1,k), phi_local(i,k+1), phi_local(i,k-1), real(k,dp) )
          phi_local(i,k)   = phi_old + w_parameter * (phi_change - phi_old)

          if (abs((phi_local(i,k)-phi_old)/phi_old).gt.tol) local_exitcond = .false.
        end do

      end do

    end if

    !combine the exit checks together
    call MPI_ALLREDUCE(local_exitcond, global_exitcond, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, istat)
    if (istat.ne.0) stop 'Error reducing exit conditions'

    !This only triggers if NO changes in any phi have been over the tolerance
    if (global_exitcond) exit

    !This print statement can be used to check if the program is running or not
    ! if ((myrank.eq.0).and.(mod(loop_counter,1000).eq.0))  print*, loop_counter

  end do

  allocate(phi_gather(Nr+1,Nz), stat=istat)
  if (istat.ne.0) stop 'Error allocating phi_gather'
  allocate(phi_temp(Nr+1,bloc_width), stat=istat)
  if(istat.ne.0) stop 'Error allocating phi_temp'
  allocate(phi_final(0:Nr, 0:Nz), stat=istat)
  if(istat.ne.0) stop 'Error allocating phi_final'

  !The transpose is necessary due to the contiguous nature of MPI_GATHER
  phi_temp = transpose(phi_local(1:bloc_width, 0:Nr))

  !Gather all grids into a single array
  call MPI_GATHER(phi_temp, bloc_width*(Nr+1), MPI_DOUBLE_PRECISION, phi_gather, &
  & bloc_width*(Nr+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, istat)
  if (istat.ne.0) stop 'Error gathering temperatures'

  if (myrank.eq.0) then
    phi_final = 0._dp
    phi_final(0:Nr, 1:Nz) = phi_gather
    phi_final(0:rinner_coord,0) = phi_0
    do k=rinner_coord+1,Nr !z=0 is lost in the gather
      phi_final(k,0) = rAxisPot( phi_0, r_outer, r_inner, k*k_step )
    end do

    print*, 'num loops:', loop_counter

    A_rcoord = 0
    B_rcoord = floor(real(Nr,dp)*6._dp/r_outer)
    C_rcoord = floor(real(Nr,dp)*5._dp/r_outer)

    A_zcoord = floor(real(Nz,dp)*7.5_dp/l_outer)
    B_zcoord = floor(real(Nz,dp)*5._dp/l_outer)
    C_zcoord = floor(real(Nz,dp)*12.5_dp/l_outer)
    print*, 'Potential at A = ', phi_final(A_rcoord,A_zcoord)
    print*, 'Potential at B = ', phi_final(B_rcoord,B_zcoord)
    print*, 'Potential at C = ', phi_final(C_rcoord,C_zcoord)

    call write_pgm(phi_final)

    time2 = MPI_wtime()
    print*, 'Time taken: ', time2-time1
  end if

  deallocate(phi_local, stat=istat)
  if (istat.ne.0) stop 'Error deallocating phi_local'
  deallocate(send1, stat=istat)
  if (istat.ne.0) stop 'Error allocating send1'
  deallocate(send2, stat=istat)
  if (istat.ne.0) stop 'Error allocating send2'
  deallocate(recv1, stat=istat)
  if (istat.ne.0) stop 'Error allocating recv1'
  deallocate(recv2, stat=istat)
  if (istat.ne.0) stop 'Error allocating recv2'
  deallocate(phi_gather, stat=istat)
  if (istat.ne.0) stop 'Error deallocating phi_gather'
  deallocate(phi_temp, stat=istat)
  if (istat.ne.0) stop 'Error deallocating phi_temp'
  deallocate(phi_final, stat=istat)
  if (istat.ne.0) stop 'Error deallocating phi_final'

  call MPI_FINALIZE(istat)
  if (istat.ne.0) stop 'Error finalizing MPI'

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
