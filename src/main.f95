program main
    use precision_mod, only : wp
    use periodic_lattice_mod
    use lattice_mod, only : Lattice
    use math_mod, only : pi
    use io_mod, only : savetxt, str

    integer, parameter :: n = 101
    real(wp) :: omega = 1.0_wp
    real(wp) :: dx = 0.01_wp
    integer :: steps = 10000
    integer :: outstep = 1000

    class(Lattice), allocatable :: latt
    type(PeriodicLattice) :: platt
    real(wp), allocatable :: density(:)

    integer :: t

    ! density = [(sin(3*pi*dx*real(i-1,wp))+1, i = 1, n)]
    density = [(0.5*dx*real(i-1,wp), i = 1, n)]

    call PeriodicLattice_(latt,omega,density)
    deallocate(latt)
    allocate(latt, source=Lattice(omega,density))
    platt = PeriodicLattice(omega,density)

    do t = 0, steps
        call latt%collide_and_stream(density)

        if (mod(t,outstep)==0) then
            call savetxt("data_"//str(t),density)
        end if
    end do

end program