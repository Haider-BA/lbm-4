module bc_module
    use precision_mod, only : wp
    use boundary_condition_mod

    public UserFunction, convectiveBC

    type, extends(UserFunction), public :: convectiveBC
        real(wp) :: k
        real(wp) :: cinfty
    contains
        procedure, pass(params) :: eval => evaluate
    end type

contains

    function evaluate(x,params) result(val)
        real(wp), intent(in) :: x
        class(convectiveBC), intent(in) :: params
        real(wp) :: val

        val = params%k*(x-params%cinfty)
    end function
end module


program main
    use precision_mod, only : wp
    ! use lattice_mod, only : Lattice
    use sman_lattice_mod, only : SmanLattice, new_SmanLatticeWithBC
    ! use periodic_lattice_mod, only : PeriodicLattice
    use bc_module
    use math_mod, only : pi
    use io_mod, only : savetxt, str

    integer, parameter :: n = 100
    real(wp) :: omega = 1.0_wp
    real(wp) :: dx = 0.01_wp
    integer :: steps = 10000
    integer :: outstep = 1000

    real(wp), parameter :: area = 1.0_wp, dt = 1.0_wp
    ! class(Lattice), allocatable :: latt
    ! type(PeriodicLattice) :: platt
    type(SmanLattice) :: slatt
    real(wp), allocatable :: density(:), positions(:)

    integer :: t

    type(convectiveBC), target :: bc
    class(UserFunction), pointer :: pbc => null()


    bc%k = 0.01_wp
    bc%cinfty = 0.1_wp

    pbc => bc

    ! density = [(sin(3*pi*dx*real(i-1,wp))+1, i = 1, n)]
    ! density = [(0.5_wp*(dx*real(i-1,wp)+0.5_wp*dx), i = 1, n)]

    density = [(1.0, i=1,n)]

    ! call PeriodicLattice_(latt,omega,density)
    ! deallocate(latt)
    ! allocate(latt, source=Lattice(omega,density))
    ! platt = PeriodicLattice(omega,density)


    slatt = new_SmanLatticeWithBC(omega,density,area,dt,pbc)

    positions = slatt%positions()
    do t = 0, steps
        ! call latt%collide_and_stream(density)
        ! if (mod(t,outstep)==0) then
            ! call savetxt("data_"//str(t),density)
        ! end if
        call slatt%collide_and_stream(density)
        ! call slatt%setLeftBC(1.0_wp)
        if (mod(t,outstep)==0) then
            call savetxt("sdata_"//str(t),reshape([positions,density],[n,2]))
        end if
    end do
end program