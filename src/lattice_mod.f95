module lattice_mod
    use precision_mod, only : wp
    use abstract_lattice_mod, only : AbstractLattice
    use cell_mod, only : LatticeCell

    implicit none
    private

    public Lattice

    type, extends(AbstractLattice), public :: Lattice
        class(LatticeCell), allocatable :: cell(:)
    contains
        private
        procedure, public :: density
        procedure :: collide_and_stream_simple
        procedure :: collide_and_stream_density
        generic, public :: collide_and_stream => collide_and_stream_simple, collide_and_stream_density
        procedure, public :: print => print_lattice
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
        class(AbstractLattice), intent(in) :: rhs
        select type(rhs)
            class is (Lattice)
                this%n = rhs%n
                this%cell = rhs%cell
            class default
                print *, "unsupported lattice class"
        end select
    end subroutine

    function new_Lattice(n,omega,dens) result(latt)
        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        type(Lattice) :: latt

        latt%n = n
        allocate(latt%cell(n))
        latt%cell = LatticeCell(omega,dens)
    end function

    function new_LatticeWithDensity(omega,dens) result(latt)
        real(wp), intent(in) :: omega
        real(wp) :: dens(:)
        type(Lattice) :: latt
        integer :: i, n
    
        n = size(dens)        
        latt%n = n
        allocate(latt%cell(n))
        do i = 1, n
            latt%cell(i) = LatticeCell(omega,dens(i))
        end do
    end function

    subroutine print_lattice(this)
        class(Lattice), intent(in) :: this
        integer :: i

        do i = 1, this%n
            call this%cell(i)%print
        end do
    end subroutine

    function density(this) result(dens)
        class(Lattice), intent(in) :: this
        real(wp), allocatable :: dens(:)
        dens = this%cell%zeroth_moment()
    end function


    subroutine collide_and_stream_simple(this)
        class(Lattice), intent(inout) :: this
        integer :: i

        call this%cell(1)%collide()
        do i = 2, this%n
            call this%cell(i)%collide()
            call this%cell(i)%swap_with(this%cell(i-1))
        end do
    end subroutine

    subroutine collide_and_stream_density(this,dens)
        class(Lattice), intent(inout) :: this
        real(wp) :: dens(this%n)
        integer :: i

        call this%cell(1)%collide(dens(1))
        do i = 2, this%n
            call this%cell(i)%collide(dens(i))
            call this%cell(i)%swap_with(this%cell(i-1))
        end do
    end subroutine
end module