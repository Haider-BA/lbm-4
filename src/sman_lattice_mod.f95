module sman_lattice_mod
    use precision_mod, only : wp
    use abstract_lattice_mod, only : AbstractLattice
    use sman_cell_mod, only : SmanCell

    implicit none
    private

    public SmanLattice

    type, extends(AbstractLattice), public :: SmanLattice
        real(wp) :: dt
        real(wp) :: area
        class(SmanCell), allocatable :: cell(:)
        real(wp), allocatable :: solid_velocity(:)
    contains
        private
        procedure, public :: density
        procedure :: collide_and_stream_simple
        procedure :: collide_and_stream_density
        ! procedure :: collide_and_stream_sman
        ! procedure :: collide_and_stream_sman_density
        generic, public :: collide_and_stream => collide_and_stream_simple, &
                                                 collide_and_stream_density
        procedure, public :: print => print_lattice
        ! procedure :: assign_
        ! generic :: assignment(=) => assign_
        ! procedure :: density
        procedure, public :: setLeftBC
        procedure, public :: setRightBC
    end type

    interface SmanLattice
        module procedure new_SmanLattice
        module procedure new_SmanLatticeWithDensity
    end interface

contains

    ! subroutine assign_(this,rhs)
    !     class(Lattice), intent(inout) :: this
    !     class(AbstractLattice), intent(in) :: rhs
    !     select type(rhs)
    !         class is (Lattice)
    !             this%n = rhs%n
    !             this%cell = rhs%cell
    !         class default
    !             print *, "unsupported lattice class"
    !     end select
    ! end subroutine

    function new_SmanLattice(n,omega,dens,dt,A) result(latt)
        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        real(wp), intent(in), optional :: dt, A
        type(SmanLattice) :: latt

        latt%n = n
        allocate(latt%cell(n))
        latt%cell = SmanCell(omega,dens)
        allocate(latt%solid_velocity(n+1))
        latt%solid_velocity = 0.0_wp
        if (present(dt)) then
            latt%dt = dt
        else
            latt%dt = 1.0_wp
        end if

        if (present(A)) then
            latt%area = A
        else
            latt%area = 1.0_wp
        end if
    end function

    function new_SmanLatticeWithDensity(omega,dens,dt,A) result(latt)
        real(wp), intent(in) :: omega
        real(wp) :: dens(:)
        real(wp), intent(in), optional :: dt, A
        type(SmanLattice) :: latt
        integer :: i, n

        n = size(dens)
        latt%n = n
        allocate(latt%cell(n))
        do i = 1, n
            latt%cell(i) = SmanCell(omega,dens(i))
        end do
        allocate(latt%solid_velocity(n+1))
        latt%solid_velocity = 0.0_wp
        if (present(dt)) then
            latt%dt = dt
        else
            latt%dt = 1.0_wp
        end if

        if (present(A)) then
            latt%area = A
        else
            latt%area = 1.0_wp
        end if
    end function



    subroutine print_lattice(this)
        class(SmanLattice), intent(in) :: this
        integer :: i

        do i = 1, this%n
            call this%cell(i)%print
        end do
    end subroutine

    function density(this) result(dens)
        class(SmanLattice), intent(in) :: this
        real(wp) :: dens(this%n)
        integer :: i
        dens = [(this%cell(i)%zeroth_moment(),i=1,this%n)]
    end function

    ! subroutine collide_and_stream_simple(this)
    !     class(SmanLattice), intent(inout) :: this
    !     integer :: i

    !     call this%cell(1)%collide()
    !     do i = 2, this%n
    !         call this%cell(i)%collide()
    !         call this%cell(i)%swap_with(this%cell(i-1))
    !     end do
    ! end subroutine

    ! subroutine collide_and_stream_density(this,dens)
    !     class(SmanLattice), intent(inout) :: this
    !     real(wp) :: dens(this%n)
    !     integer :: i

    !     call this%cell(1)%collide(dens(1))
    !     do i = 2, this%n
    !         call this%cell(i)%collide(dens(i))
    !         call this%cell(i)%swap_with(this%cell(i-1))
    !     end do
    ! end subroutine

    subroutine collide_and_stream_simple(this)
        class(SmanLattice), intent(inout) :: this
        real(wp) :: vs(2), dt
        integer :: i

        dt = this%dt
        vs = this%solid_velocity(1:2)
        call this%cell(1)%collide(vs,dt)
        do i = 2, this%n
            vs = this%solid_velocity(i:i+1)
            call this%cell(i)%collide(vs,dt)
            call this%cell(i)%swap_with(this%cell(i-1))
        end do
    end subroutine

    subroutine collide_and_stream_density(this,dens)
        class(SmanLattice), intent(inout) :: this
        real(wp) :: dens(this%n)
        integer :: i

        associate(dt => this%dt, vs => this%solid_velocity)
        ! vs = this%solid_velocity(1:2)
        call this%cell(1)%collide(vs(1:2),dt,dens(1))
        do i = 2, this%n
            ! vs = this%solid_velocity(i:i+1)
            call this%cell(i)%collide(vs(i:i+1),dt,dens(i))
            call this%cell(i)%swap_with(this%cell(i-1))
        end do
        end associate
    end subroutine

    subroutine setLeftBC(this,value)
        class(SmanLattice), intent(inout) :: this
        real(wp), intent(in) :: value

        call this%cell(1)%force_west(value)
    end subroutine

    subroutine setRightBC(this,value)
        class(SmanLattice), intent(inout) :: this
        real(wp), intent(in) :: value

        call this%cell(this%n)%force_east(value)
    end subroutine
end module