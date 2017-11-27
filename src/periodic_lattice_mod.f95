module periodic_lattice_mod
    use precision_mod
    use lattice_mod, only : Lattice
    use cell_mod, only : LatticeCell

    implicit none
    private
    public PeriodicLattice, PeriodicLattice_

    type, public, extends(Lattice) :: PeriodicLattice
    contains
        private
        procedure :: collide_and_stream_simple
        procedure :: collide_and_stream_density
    end type

    interface PeriodicLattice
        module procedure new_Lattice
        module procedure new_LatticeWithDensity
    end interface

    interface PeriodicLattice_
        module procedure PeriodicLattice_simple
        module procedure PeriodicLattice_density
    end interface


contains

    subroutine PeriodicLattice_simple(latt,omega,dens)
        class(Lattice), intent(out), allocatable :: latt
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens(:)

        allocate(latt, source=PeriodicLattice(omega,dens))
    end subroutine

    subroutine PeriodicLattice_density(latt,n,omega,dens)
        class(Lattice), intent(out), allocatable :: latt
        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens

        allocate(latt, source=PeriodicLattice(n,omega,dens))
    end subroutine

    function new_Lattice(n,omega,dens) result(latt)
        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        type(PeriodicLattice) :: latt

        latt%n = n
        allocate(latt%cell(n))
        latt%cell(:) = LatticeCell(omega,dens)
    end function

    function new_LatticeWithDensity(omega,dens) result(latt)
        real(wp), intent(in) :: omega
        real(wp) :: dens(:)
        type(PeriodicLattice) :: latt
        integer :: i, n

        n = size(dens)
        latt%n = n
        allocate(latt%cell(n))
        do i = 1, n
            latt%cell(i) = LatticeCell(omega,dens(i))
        end do
    end function

    subroutine collide_and_stream_simple(this)
        class(PeriodicLattice), intent(inout) :: this
        integer :: i

        call this%cell(1)%collide()
        do i = 2, this%n
            call this%cell(i)%collide()
            call this%cell(i)%swap_with(this%cell(i-1))
        end do
        call this%cell(1)%swap_with(this%cell(this%n)) ! periodic bc

    end subroutine

    subroutine collide_and_stream_density(this,dens)
        class(PeriodicLattice), intent(inout) :: this
        real(wp), intent(out) :: dens(this%n)
        integer :: i

        call this%cell(1)%collide(dens(1))
        do i = 2, this%n
            call this%cell(i)%collide(dens(i))
            call this%cell(i)%swap_with(this%cell(i-1))
        end do
        call this%cell(1)%swap_with(this%cell(this%n)) ! periodic bc
    end subroutine

end module