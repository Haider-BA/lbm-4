module alattice_mod
    use precision_mod, only : wp

    implicit none
    private

    type, abstract, public :: aLattice
        integer :: n
    contains
    private
        procedure(int1), deferred :: collide_and_stream_simple
        procedure(int2), deferred :: collide_and_stream_density
    end type


    abstract interface
        subroutine int1(this)
            import :: aLattice
            class(aLattice), intent(inout) :: this
        end subroutine
        subroutine int2(this,dens)
            import :: aLattice, wp
            class(aLattice), intent(inout) :: this
            real(wp) :: dens(this%n)
        end subroutine
    end interface
end module

module lattice_mod
    use precision_mod, only : wp
    use cell_mod, only : LatticeCell, new_LatticeCell
    use alattice_mod

    implicit none
    private

    public Lattice

    type, extends(aLattice) :: Lattice
        class(LatticeCell), allocatable :: cell(:)
    contains
        private
        procedure, public :: print
        procedure, public :: density
        procedure :: collide_and_stream_simple
        procedure :: collide_and_stream_density
        generic, public :: collide_and_stream => collide_and_stream_simple, collide_and_stream_density
        procedure :: assign_
        generic :: assignment(=) => assign_
        ! procedure :: density
    end type

    interface Lattice
        module procedure new_Lattice
        module procedure new_LatticeWithDensity
    end interface

contains

    subroutine assign_(this,rhs)
        class(Lattice), intent(inout) :: this
        class(Lattice), intent(in) :: rhs
        this%n = rhs%n
        this%cell = rhs%cell
    end subroutine

    function new_Lattice(n,omega,dens) result(latt)
        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        type(Lattice) :: latt

        latt%n = n
        allocate(latt%cell(n))

        ! do i = 1, n
           ! latt%cell(i) = LatticeCell(omega,1.0_wp,dens)
        ! end do

        latt%cell(:) = LatticeCell(omega,1.0_wp,dens)

    end function

    function new_LatticeWithDensity(omega,dens) result(latt)
        real(wp), intent(in) :: omega
        real(wp) :: dens(:)
        type(Lattice) :: latt
        integer :: i
        latt%n = size(dens)
        allocate(latt%cell(latt%n))
        do i = 1, latt%n
            latt%cell(i) = LatticeCell(omega,1.0_wp,dens(i))
        end do
    end function

    subroutine print(this)
        class(Lattice), intent(in) :: this

        integer :: i
        do i = 1, this%n
            call this%cell(i)%print
        end do
    end subroutine

    function density(this) result(dens)
        class(Lattice), intent(in) :: this
        real(wp), allocatable :: dens(:)
        dens = this%cell%density()
    end function


    subroutine collide_and_stream_simple(this)
        class(Lattice), intent(inout) :: this
        integer :: i

        call this%cell(1)%simple_collision()
        do i = 2, this%n
            call this%cell(i)%simple_collision()
            call this%cell(i)%cell_swap(this%cell(i-1))
        end do
    end subroutine

    subroutine collide_and_stream_density(this,dens)
        class(Lattice), intent(inout) :: this
        real(wp) :: dens(this%n)
        integer :: i

        call this%cell(1)%simple_collision(dens(1))
        do i = 2, this%n
            call this%cell(i)%simple_collision(dens(i))
            call this%cell(i)%cell_swap(this%cell(i-1))
        end do
    end subroutine
end module

module periodic_mod
    use precision_mod
    use lattice_mod, only : Latticeg => Lattice
    use cell_mod, only : LatticeCell
    ! use alattice_mod
    ! use cell_mod, only: swap
    implicit none
    private
    public PeriodicLattice, assignment(=), Lattice

    type, public, extends(Latticeg) :: PeriodicLattice
    contains
        private
        procedure :: collide_and_stream_simple
        procedure :: collide_and_stream_density
    end type

    interface Lattice
        module procedure new_Lattice
        module procedure new_LatticeWithDensity
    end interface

    interface assignment(=)
        module procedure periodic_from_normal
    end interface

contains

    function new_Lattice(n,omega,dens) result(latt)
        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        type(PeriodicLattice) :: latt

        latt%n = n
        allocate(latt%cell(n))

        ! do i = 1, n
           ! latt%cell(i) = LatticeCell(omega,1.0_wp,dens)
        ! end do

        latt%cell(:) = LatticeCell(omega,1.0_wp,dens)

    end function

    function new_LatticeWithDensity(omega,dens) result(latt)
        real(wp), intent(in) :: omega
        real(wp) :: dens(:)
        type(PeriodicLattice) :: latt
        integer :: i
        latt%n = size(dens)
        allocate(latt%cell(latt%n))
        do i = 1, latt%n
            latt%cell(i) = LatticeCell(omega,1.0_wp,dens(i))
        end do
    end function

    subroutine collide_and_stream_simple(this)
        class(PeriodicLattice), intent(inout) :: this
        integer :: i

        call this%cell(1)%simple_collision()
        do i = 2, this%n
            call this%cell(i)%simple_collision()
            call this%cell(i)%cell_swap(this%cell(i-1))
        end do
        call this%cell(1)%cell_swap(this%cell(this%n)) ! periodic bc
    end subroutine

    subroutine collide_and_stream_density(this,dens)
        class(PeriodicLattice), intent(inout) :: this
        real(wp) :: dens(this%n)
        integer :: i

        call this%cell(1)%simple_collision(dens(1))
        do i = 2, this%n
            call this%cell(i)%simple_collision(dens(i))
            call this%cell(i)%cell_swap(this%cell(i-1))
        end do
        call this%cell(1)%cell_swap(this%cell(this%n)) ! periodic bc
    end subroutine

    subroutine periodic_from_normal(lhs,rhs)
        type(PeriodicLattice), intent(inout) :: lhs
        class(Latticeg), intent(in) :: rhs

        lhs%n = rhs%n
        lhs%cell = rhs%cell

    end subroutine

end module

program main
    use precision_mod, only : wp
    use periodic_mod, only : PeriodicLattice, Lattice
    ! use lattice_mod, only : Lattice
    use math_mod, only : pi
    use io_mod

    integer, parameter :: n = 101
    real(wp) :: omega = 1.0_wp
    real(wp) :: dx = 0.01_wp
    integer :: steps = 1000

    type(PeriodicLattice) :: latt
    real(wp), allocatable :: density(:)

    integer :: t


    density = [(sin(3*pi*dx*real(i-1,wp))+1, i = 1, n)]
    density = [(0.5*dx*real(i-1,wp), i = 1, n)]

    ! allocate(latt,source=PeriodicLattice(omega,density))
    latt = Lattice(omega,density)

    do t = 0, steps
        call latt%collide_and_stream(density)

        if (mod(t,100)==0) then
            call savetxt("data_"//str(t),density)
        end if
    end do

end program