module abstract_cell_mod
    use precision_mod, only : wp

    implicit none
    private

    type, abstract, public :: AbstractCell
        real(wp) :: omega
        real(wp) :: one_minus_omega
        real(wp) :: pdf(3)
    contains
        procedure(symm), deferred :: assign
        generic :: assignment(=) => assign
        procedure(swap), deferred :: swap_with
    end type

    abstract interface
        elemental subroutine symm(lhs,rhs)
            import :: AbstractCell
            class(AbstractCell), intent(inout) :: lhs
            class(AbstractCell), intent(in) :: rhs
        end subroutine
        subroutine swap(this,previous)
            import :: AbstractCell
            class(AbstractCell), intent(inout) :: this
            class(AbstractCell), intent(inout) :: previous
        end subroutine
    end interface
end module

module cell_mod
    use precision_mod, only : wp
    use abstract_cell_mod, only : AbstractCell

    implicit none
    private
    public LatticeCell, AbstractCell

    real(wp), parameter :: W(3) = [2.0_wp/3.0_wp, 1.0_wp/6.0_wp, 1.0_wp/6.0_wp]
    real(wp), parameter :: C_SQR = 1.0_wp/3.0_wp

    type, extends(AbstractCell), public :: LatticeCell
    contains
        procedure, public :: zeroth_moment
        procedure, public :: first_moment
        procedure, public :: collide_simple
        generic :: collide => collide_simple
        procedure, public :: diffusivity
        procedure, public :: assign => assign_cell
        procedure, public :: print
        procedure, public :: swap_with => swap_with_simple
    end type

    interface LatticeCell
        module procedure new_LatticeCell
    end interface

contains

    elemental subroutine assign_cell(lhs,rhs)
        class(LatticeCell), intent(inout) :: lhs
        class(AbstractCell), intent(in) :: rhs
        select type(rhs)
            type is (LatticeCell)
                lhs%omega = rhs%omega
                lhs%pdf = rhs%pdf
        end select
    end subroutine    

    pure elemental function new_LatticeCell(omega,dens) result(cell)
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        type(LatticeCell) :: cell
        cell%omega = omega
        cell%one_minus_omega = 1.0_wp - omega
        cell%pdf = W*dens
    end function


    pure elemental function zeroth_moment(this) result(dens)
        class(LatticeCell), intent(in) :: this
        real(wp) :: dens
        dens = sum(this%pdf)
    end function


    pure elemental function first_moment(this) result(moment)
        class(LatticeCell), intent(in) :: this
        real(wp) :: moment
        moment = this%pdf(2)-this%pdf(3)
    end function


    subroutine collide_simple(this,dens)
        class(LatticeCell), intent(inout) :: this
        real(wp), intent(out), optional :: dens

        real(wp) :: feq(3), m0

        ! density
        m0 = this%zeroth_moment()
        if (present(dens)) dens = m0

        ! equilibrium distribution function
        feq = W*m0

        ! collide and swap in place
        this%pdf(1) = this%one_minus_omega*this%pdf(1) + this%omega*feq(1)
        this%pdf(2) = this%one_minus_omega*this%pdf(3) + this%omega*feq(3)
        this%pdf(3) = this%one_minus_omega*this%pdf(2) + this%omega*feq(2)
    end subroutine


    subroutine swap_with_simple(this,previous)
        class(LatticeCell), intent(inout) :: this
        class(AbstractCell), intent(inout) :: previous
        real(wp) :: temp
        temp = this%pdf(2)
        this%pdf(2) = previous%pdf(3)
        previous%pdf(3) = temp
    end subroutine


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