module sman_lattice_mod
    use precision_mod, only : wp
    use abstract_lattice_mod, only : AbstractLattice
    use sman_cell_mod, only : SmanCell
    use rootfinder_mod, only : secant
    use boundary_condition_mod


    implicit none
    private

    public SmanLattice, new_SmanLatticeWithBC

    type, extends(AbstractLattice), public :: SmanLattice
        real(wp) :: dt
        real(wp) :: area
        class(SmanCell), allocatable :: cell(:)
        real(wp), allocatable :: solid_velocity(:)
        class(DirichletCondition), allocatable :: dleft, dright
        class(FluxCondition), allocatable :: fleft, fright
    contains
        private
        procedure, public :: density
        procedure :: collide_and_stream_simple
        procedure :: collide_and_stream_density
        ! procedure :: collide_and_stream_sman
        ! procedure :: collide_and_stream_sman_density
        procedure :: collide_and_stream_density_with_bc
        ! procedure :: collide_and_stream_density_with_bc2
        generic, public :: collide_and_stream => collide_and_stream_simple, &
                                                 collide_and_stream_density, &
                                                 collide_and_stream_density_with_bc
        procedure, public :: print => print_lattice
        ! procedure :: assign_
        ! generic :: assignment(=) => assign_
        ! procedure :: density
        procedure, public :: setLeftBC
        procedure, public :: setRightBC
        procedure, public :: positions
    end type

    interface SmanLattice
        module procedure new_SmanLatticeWithBC
        ! module procedure new_SmanLattice
        ! module procedure new_SmanLatticeWithDensity
    end interface

contains


    pure function positions(this) result(x)
        class(SmanLattice), intent(in) :: this
        real(wp) :: x(this%n)
        integer :: i

        x(1) = this%cell(1)%length*0.5_wp
        do i = 2, this%n
            x(i) = x(i-1) + 0.5_wp*(this%cell(i-1)%length+this%cell(i)%length)
        end do
    end function
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

        latt%dleft = DirichletCondition(1.0_wp,latt%cell(1)%diffusivity(),dt,A)
    end function

    function new_SmanLatticeWithDensity(omega,dens,dt,A) result(latt)
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens(:)
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
        latt%dleft = DirichletCondition(1.0_wp,latt%cell(1)%diffusivity(),dt,A)
    end function

    function new_SmanLatticeWithBC(omega,dens,dt,A,bc1) result(latt)
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens(:)
        real(wp), intent(in) :: dt, A
        class(UserFunction), pointer :: bc1
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
        latt%dt = dt

        latt%area = A

        allocate(latt%fleft)
        ! print *, allocated(latt%fleft)
        latt%fleft = new_Flux(dt,A,bc1)
        ! else if (present(bc2)) then
        !     latt%fright = FluxCondition(dt,A,bc2)
        ! end if

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
        real(wp), intent(out) :: dens(this%n)
        integer :: i

        associate(dt => this%dt, vs => this%solid_velocity, cell => this%cell, &
            left => this%fleft)

        call this%cell(1)%collide(vs(1:2),dt,dens(1))
        do i = 2, this%n
            call this%cell(i)%collide(vs(i:i+1),dt,dens(i))
            call this%cell(i)%swap_with(this%cell(i-1))
        end do

        if (allocated(this%fleft)) then
            left%old = this%cell(1)%pdf(2)
            left%pdf = this%cell(1)%pdf([1,3])
            left%grid(1) = 0.0_wp
            left%grid(2) = this%cell(1)%length*0.5_wp
            do i = 2, 4
                left%grid(i+1) = left%grid(i) + 0.5_wp*(this%cell(i-1)%length+this%cell(i)%length)
            end do
            ! print *, left%grid
        ! left%grid = [0.0_wp,(0.5*this%length(i),i = 1, 4)]
            left%conc(:) = [0.0_wp,(this%cell(i)%zeroth_moment(), i = 1, 4)]
            this%cell(1)%pdf(2) = secant(left,left%old)
        print *, this%cell(1)%pdf(2), left%old
        end if
        end associate
    end subroutine



    subroutine collide_and_stream_density_with_bc(this,dens,cwall)
        class(SmanLattice), intent(inout) :: this
        real(wp), intent(out) :: dens(this%n), cwall
        integer :: i

        associate(dt => this%dt, vs => this%solid_velocity, &
            left => this%dleft)

        call this%cell(1)%collide(vs(1:2),dt,dens(1))
        do i = 2, this%n
            call this%cell(i)%collide(vs(i:i+1),dt,dens(i))
            call this%cell(i)%swap_with(this%cell(i-1))
        end do


        left%old = this%cell(1)%pdf(2)
        left%pdf = this%cell(1)%pdf([1,3])
        left%grid(1) = 0.0_wp
        left%grid(2) = this%cell(1)%length*0.5_wp
        do i = 2, 4
            left%grid(i+1) = left%grid(i) + 0.5_wp*(this%cell(i-1)%length+this%cell(i)%length)
        end do
        left%conc(:) = [cwall,(this%cell(i)%zeroth_moment(), i = 1, 4)]

        this%cell(1)%pdf(2) = secant(left,left%old) ! this is where bc actually happens
        end associate
    end subroutine

    subroutine collide_and_stream_density_with_bc2(this,dens)
        class(SmanLattice), intent(inout) :: this
        real(wp), intent(out) :: dens(this%n)
        integer :: i

        associate(dt => this%dt, vs => this%solid_velocity, &
            left => this%dleft)

        call this%cell(1)%collide(vs(1:2),dt,dens(1))
        do i = 2, this%n
            call this%cell(i)%collide(vs(i:i+1),dt,dens(i))
            call this%cell(i)%swap_with(this%cell(i-1))
        end do

        print *, "hi ivan from bc"
        ! left%flux = bc
        ! left%dt = this%dt
        ! left%old = this%cell(1)%pdf(2)
        ! left%pdf = this%cell(1)%pdf([1,3])
        ! left%grid(1) = 0.0_wp
        ! left%grid(2) = this%cell(1)%length*0.5_wp
        ! do i = 2, 4
        !     left%grid(i+1) = left%grid(i) + 0.5_wp*(this%cell(i-1)%length+this%cell(i)%length)
        ! end do
        ! print *, left%grid
        ! ! left%grid = [0.0_wp,(0.5*this%length(i),i = 1, 4)]
        ! left%conc(:) = [0.0_wp,(this%cell(i)%zeroth_moment(), i = 1, 4)]
        ! this%cell(1)%pdf(2) = secant(this%left,left%old)
        ! print *, this%cell(1)%pdf(2), left%old
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