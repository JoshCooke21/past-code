module grapheneFunctions
    implicit none

    integer, parameter :: dp = selected_real_kind(15,300)
    real(dp), parameter :: PI=4._dp*atan(1._dp)
    complex*16, parameter :: imag = (0._dp,1._dp)

    integer :: istat=0,eigenvalues_unit=21
    real(dp) :: self_interaction=0._dp, hopping_1=2.7_dp, hopping_2=0.2_dp

contains
    subroutine runprogram(make_hamiltonian, num_atoms, num_neigh)
        implicit none

        complex*16 :: make_hamiltonian
        integer, intent(in) :: num_atoms, num_neigh

        integer :: loop_counter, num_points, i, j
        real(dp) :: ka, d_ka

        ! zheev variables
        complex*16, allocatable, dimension(:,:) :: hamiltonian
        complex*16, allocatable, dimension(:) :: WORK
        real(dp), allocatable, dimension(:) :: RWORK, eigenvalues
        integer :: LWORK, INFO


        ! Allocate the zheev arrays
        LWORK=2*num_atoms-1
        allocate( WORK(LWORK), stat=istat)
        if (istat.ne.0) stop 'Error allocating WORK'
        allocate( RWORK(3*num_atoms-2), stat=istat)
        if (istat.ne.0) stop 'Error allocating RWORK'
        allocate( eigenvalues(num_atoms), stat=istat)
        if (istat.ne.0) stop 'Error allocating eigenvalues'

        ! Allocate the Hamiltonian matrix
        allocate( hamiltonian(num_atoms, num_atoms), stat=istat)
        if (istat.ne.0) stop 'Error allocating Hamiltonian'

        ! Initialize some variables
        ka=0._dp
        num_points=1000
        ! d_ka will increment ka along its range 0 -> pi
        d_ka = PI / real(num_points,dp)

        ! This loop runs along the range of ka
        do loop_counter = 1, num_points+1

            ! Initialize the Hamiltonian
            hamiltonian(:,:)= (0._dp , 0._dp)

            ! Calculate each element depending on the type of lattice
            do i=1,num_atoms
                do j=i,num_atoms
                    hamiltonian(i,j) = make_hamiltonian(ka, i, j, num_neigh)
                end do
            end do

            ! Find the eigenvalues of the matrix
            call zheev('N','U',num_atoms, hamiltonian, num_atoms, eigenvalues, WORK, LWORK, RWORK, INFO)
            if (INFO.ne.0) stop 'zheev not successful'

            ! Write the value of ka and the whole eigenvalue array
            write(unit=eigenvalues_unit, fmt=*, iostat=istat) ka, eigenvalues
            if (istat.ne.0) stop 'Error writing to file'

            ka = ka + d_ka
        end do

    end subroutine

    function square_hamiltonian( x_pos, index_1, index_2, dummy_nn )
        implicit none

        real(dp), intent(in) :: x_pos
        integer, intent(in) :: index_1, index_2, dummy_nn ! Dummy is not used for the square lattice

        complex*16 :: square_hamiltonian
        integer :: dummy_dash

        ! This line is entirely unnecessary, aside from avoiding a warning when compiling
        dummy_dash = dummy_nn

        if ( index_2.eq.index_1 ) then
            ! If on the diagonal
            square_hamiltonian = self_interaction
        else if ( index_2.eq.(index_1+1) ) then
            ! Left/right interactions
            square_hamiltonian = - hopping_1 * (1._dp + exp(-imag*x_pos))
        else if ( index_2.eq.(index_1+2) ) then
            ! Up/down interactions
            square_hamiltonian = - hopping_1
        else
            ! For when there are no valid interactions
            square_hamiltonian=0._dp
        end if

    end function square_hamiltonian

    function agnr_hamiltonian( x_pos, index_1, index_2, neighbours)
        implicit none

        real(dp), intent(in) :: x_pos
        integer, intent(in) :: index_1, index_2, neighbours

        complex*16 :: agnr_hamiltonian

        if ( index_2.eq.index_1 ) then
            ! On the diagonal
            agnr_hamiltonian = self_interaction
        else if ( index_2.eq.(index_1+2) ) then
            ! Interactions along the row (a1 to a2 etc)
            agnr_hamiltonian = - hopping_1
        else if ( index_2.eq.(index_1+1) ) then
            ! Interactions with atoms above/below
            ! The first if statement checks for j is even.
            !   If so, it then checks whether it is a multiple of 4.
            !       If so, then the interaction is with the next cell.
            !       Else the interaction is within the unit cell.
            ! If j is odd, there is no valid interaction
            if ( mod(index_2,2).eq.0 ) then
                if ( mod(index_2,4).eq.0 ) then
                    agnr_hamiltonian = - hopping_1 * exp(imag * x_pos)
                else
                    agnr_hamiltonian = - hopping_1
                end if
            else
                agnr_hamiltonian = 0._dp
            end if
        else
            agnr_hamiltonian = 0._dp
        end if

        if (neighbours.ge.2) then
            if (index_2.eq.(index_1+1)) then
                if ( mod(index_2,2).eq.1) then
                    agnr_hamiltonian = agnr_hamiltonian - hopping_2 * (1._dp + exp(-imag*x_pos))
                end if
            else if ( index_2.eq.(index_1+3) ) then
                if ( mod(index_2,2).eq.0 ) then
                    agnr_hamiltonian = agnr_hamiltonian - hopping_2 * (1._dp + exp(imag*x_pos))
                end if
            else if ( index_2.eq.(index_1+4) ) then
                agnr_hamiltonian = agnr_hamiltonian - hopping_2
            end if
        end if

    end function agnr_hamiltonian

    function zgnr_hamiltonian( x_pos, index_1, index_2, neighbours )
        implicit none

        real(dp), intent(in) :: x_pos
        integer, intent(in) :: index_1, index_2, neighbours

        complex*16 :: zgnr_hamiltonian

        if ( index_2.eq.index_1 ) then
            ! On the diagonal
            zgnr_hamiltonian = self_interaction
        else if ( index_2.eq.(index_1+1) ) then
            ! There isn't a consistent way of describing visually the relationship with reference to
            !   the unit cell.
            ! However, in the Hamiltonian, this represents the value along the 1st offset diagonal.
            !   If j is even, then a value is set. Else, the element is 0.
            if ( mod(index_2,2).eq.0 ) then
                zgnr_hamiltonian = - hopping_1 * (1 + exp(imag*x_pos) )
            else
                zgnr_hamiltonian = 0._dp
            end if
        else if ( index_2.eq.(index_1+2) ) then
            ! Along the 2nd offset diagonal, the pattern is more obscure.
            ! Wherever j%4 is 2 or 3, a hopping_1 term is set.
            !   E.g. j=3,6,7,10,11,etc
            ! Else element=0
            if ((mod(index_2,4).eq.2).or.(mod(index_2,4).eq.3) ) then
                zgnr_hamiltonian = -hopping_1
            else
                zgnr_hamiltonian = 0._dp
            end if
        else
            zgnr_hamiltonian = 0._dp
        end if

        if (neighbours.ge.2) then
            if ( index_2.eq.index_1 ) then
                zgnr_hamiltonian = zgnr_hamiltonian - 2._dp * hopping_2 * cos(x_pos)
            else if ( index_2.eq.(index_1+1) ) then
                if ( mod(index_2,2).eq.1 ) then
                    zgnr_hamiltonian = zgnr_hamiltonian - hopping_2 * (1._dp + exp(-imag*x_pos))
                end if
            else if ( index_2.eq.(index_1+3) ) then
                if ( mod(index_2,2).eq.0 ) then
                    zgnr_hamiltonian = zgnr_hamiltonian - hopping_2 * (1._dp + exp(imag*x_pos))
                end if
            end if

        end if

    end function zgnr_hamiltonian

end module grapheneFunctions
