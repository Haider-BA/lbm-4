module lattice_mod
    use precision_mod, only : wp
    use cell_mod, only : LatticeCell, new_LatticeCell

    implicit none   
    private

    public Lattice

    type :: Lattice
        integer :: n
        class(LatticeCell), allocatable :: cell(:)
    contains
        procedure :: print
        procedure :: density
        procedure, private :: collide_and_stream_simple
        procedure, private :: collide_and_stream_density
        generic :: collide_and_stream => collide_and_stream_simple, collide_and_stream_density
        ! procedure :: density
    end type

    interface Lattice
        module procedure new_Lattice
        module procedure new_LatticeWithDensity
    end interface

contains

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

program main
    use precision_mod, only : wp
    use lattice_mod, only : Lattice
    use math_mod, only : PI
    use io_mod

    integer, parameter :: n = 101
    real(wp) :: omega = 1.0_wp
    real(wp) :: dx = 0.01_wp
    integer :: steps = 1000

    type(Lattice) :: latt
    real(wp), allocatable :: density(:)

    integer :: t

    density = [(sin(2*PI*dx*real(i-1,wp))**2, i = 1, n)]

    latt = Lattice(omega,density)

    do t = 0, steps
        call latt%collide_and_stream(density)

        if (mod(t,100)==0) then
            call savetxt("data_"//str(t),density)
        end if
    end do

end program