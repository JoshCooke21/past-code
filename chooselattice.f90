program chooselattice
    use grapheneFunctions
    implicit none

    integer :: num_atoms, lattice_choice, nn_choice

    open(unit=eigenvalues_unit, file='data.txt', iostat=istat)
    if (istat.ne.0) stop 'Error opening data.txt'

    write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Graphene simulation'

    ! Request a lattice shape from the user
    print *, 'What lattice shape would you like?'
    print *, ' [1] Square lattice'
    print *, ' [2] Armchair lattice'
    print *, ' [3] Zig-zag lattice'
    read(*,*) lattice_choice
    print *, ' '
    print *, 'How many atoms would you like'
    read(*,*) num_atoms
    print *, ' '
    print *, 'How many nearest neighbours? (Enter as an integer, 1 or 2)'
    print *, 'If the square lattice was chosen, this will be defaulted to 1st NN'
    read(*,*) nn_choice

    ! Write the information to the data file
    write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Total atoms = ', num_atoms
    write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Self-interaction = ', self_interaction
    write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Hopping 1 = ', hopping_1
    write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Hopping 2 = ', hopping_2

    ! Depending on the lattice type chosen, the program uses a different Hamiltonian
    !   maker function.
    ! The lattice type is also written to the data file.
    select case(lattice_choice)
    case(1)
        print *, ' '
        print *, 'Square lattice shape chosen.'
        write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Lattice type = square'
        call runprogram(square_hamiltonian, num_atoms, 1)

    case(2)
        print *, ' '
        print *, 'Armchair (AGNR) lattice chosen.'
        write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Lattice type = AGNR'
        call runprogram(agnr_hamiltonian, num_atoms, nn_choice)
    case(3)
        print *, ' '
        print *, 'Zig-zag (ZGNR) lattice chosen.'
        write(unit=eigenvalues_unit, fmt=*, iostat=istat) '# Lattice type = ZGNR'
        call runprogram(zgnr_hamiltonian, num_atoms, nn_choice)

    case default
        ! In case the user inputs a non-valid choice, the default is set to a
        !   square lattice.
        print *, ' '
        print *, 'Invalid choice.'
        print *, 'Defaulting to square lattice.'
        lattice_choice=1

        call runprogram(square_hamiltonian, num_atoms, 1)
    end select

    close(unit=eigenvalues_unit)

end program chooselattice
