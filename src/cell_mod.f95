module cell_mod
    use precision_mod, only : wp

    implicit none
    private
    public LatticeCell, new_LatticeCell

    real(wp), parameter :: W(3) = [2.0_wp/3.0_wp, 1.0_wp/6.0_wp, 1.0_wp/6.0_wp]
    real(wp), parameter :: C_SQR = 1.0_wp/3.0_wp
    real(wp), parameter :: CC(3) = [0.0_wp,1.0_wp,-1.0_wp]

    type :: LatticeCell
        real(wp) :: omega
        real(wp) :: length
        real(wp) :: pdf(3)
        real(wp) :: c(3)
    contains
        procedure, public :: density
        procedure, public :: first_moment
        procedure, public :: simple_collision
        procedure, public :: diffusivity
        procedure, public :: assign_cell
        procedure, public :: print
        procedure, public :: cell_swap
        generic :: assignment(=) => assign_cell
    end type

    interface LatticeCell
        module procedure new_LatticeCell
    end interface

contains

    elemental subroutine assign_cell(lhs,rhs)
        class(LatticeCell), intent(inout) :: lhs
        class(LatticeCell), intent(in) :: rhs
        select type(rhs)
            class is (LatticeCell)
                lhs%omega = rhs%omega
                lhs%length = rhs%length
                lhs%pdf = rhs%pdf
                lhs%c = rhs%c
            class default
                ! stop "assign_cell: unsupported class"
        end select
    end subroutine    

    pure elemental function new_LatticeCell(omega,length,dens) result(cell)
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: length
        real(wp), intent(in) :: dens
        type(LatticeCell) :: cell

        cell%omega = omega
        cell%length = length
        cell%pdf = W*dens
        cell%c = CC
    end function


    pure elemental function density(this) result(dens)
        class(LatticeCell), intent(in) :: this
        real(wp) :: dens
        dens = sum(this%pdf)
    end function


    pure elemental function first_moment(this) result(moment)
        class(LatticeCell), intent(in) :: this
        real(wp) :: moment
        moment = sum(this%c*this%pdf)
    end function


    ! subroutine perform_collision(this,dens)
    !     class(LatticeCell), intent(inout) :: this
    !     real(wp), intent(out), optional :: dens 

    !     ! calculate pre-collision moments from PDF, Eq. 13
    !     m0 = this%density()
    !     m1 = this%first_moment() - m0*c(1)

    !     ! calculate post-collision moments
    !     ! m0 = m0, first moment is conserved
    !     m1 = m1*(1.0_wp-omega)

    !     ! update cell size and lattice velocity

    !     ! output density if necessary
    !     if (present(dens)) then
    !         dens = m0
    !     end if
    ! end subroutine

    subroutine simple_collision(this,dens)
        class(LatticeCell), intent(inout) :: this
        real(wp), intent(out), optional :: dens

        real(wp) :: feq(3), m0, one_minus_omega

        m0 = this%density()
        feq = W*m0
        one_minus_omega = 1.0_wp - this%omega

        this%pdf(1) = one_minus_omega*this%pdf(1) + this%omega*feq(1)
        this%pdf(2) = one_minus_omega*this%pdf(3) + this%omega*feq(3)
        this%pdf(3) = one_minus_omega*this%pdf(2) + this%omega*feq(2)

        if (present(dens)) then
            dens = m0
        end if
    end subroutine

    subroutine cell_swap(this,previous)
        class(LatticeCell), intent(inout) :: this
        class(LatticeCell), intent(inout) :: previous

        call swap(this%pdf(2),previous%pdf(3))
    end subroutine

    subroutine swap(f1,f2)
        real(wp), intent(inout) :: f1, f2
        real(wp) :: temp
        temp = f1
        f1 = f2
        f2 = temp
    end subroutine swap

    pure elemental function diffusivity(this) result(d)
        class(LatticeCell), intent(in) :: this
        real(wp) :: d
        d = C_SQR*(1.0_wp/this%omega - 0.5_wp)
    end function

    subroutine print(this)
        class(LatticeCell), intent(in) :: this

        print *, this%omega, this%pdf
    end subroutine
end module