module sman_cell_mod
    use precision_mod, only : wp
    use cell_mod, only : LatticeCell, AbstractCell

    implicit none
    private

    real(wp) :: C_SQR = 1.0_wp/3.0_wp

    type, extends(LatticeCell), public :: SmanCell
        private
        real(wp) :: c(3)
        real(wp) :: length
    contains
        private
        procedure, public :: first_moment
        procedure :: collide_sman
        generic, public :: collide => collide_sman
        procedure, public :: assign => assign_sman_cell
        procedure :: swap_with_sman
        procedure, public :: swap_with => swap_with_sman
        procedure, public :: force_east
        procedure, public :: force_west
    end type

    interface SmanCell
        module procedure new_SmanCell
        module procedure new_SmanCellComplex
    end interface
contains

    elemental subroutine assign_sman_cell(lhs,rhs)
        class(SmanCell), intent(inout) :: lhs
        class(AbstractCell), intent(in) :: rhs
        select type(rhs)
            class is (LatticeCell)
                lhs%omega = rhs%omega
                lhs%one_minus_omega = rhs%one_minus_omega
                lhs%pdf = rhs%pdf
                lhs%c = [0.0_wp,1.0_wp,-1.0_wp]
                lhs%length = 1.0_wp
            class is (SmanCell)
                lhs%omega = rhs%omega
                lhs%one_minus_omega = rhs%one_minus_omega
                lhs%pdf = rhs%pdf
                lhs%c = rhs%c
                lhs%length = rhs%length
        end select
    end subroutine

    pure function weights_(c) result(w)
        real(wp), intent(in) :: c(3)
        real(wp) :: w(3)
        real(wp) :: c2, c3, csum

        c2 = abs(c(2))
        c3 = abs(c(3))
        csum = c2 + c3

        w(1) = 1.0_wp - C_SQR/(c2*c3)
        w(2) = C_SQR/(c2*csum)
        w(3) = C_SQR/(c3*csum)
    end function


    pure elemental function new_SmanCell(omega,dens) result(cell)
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        type(SmanCell) :: cell
        real(wp) :: t(3)

        cell%omega = omega
        cell%one_minus_omega = 1.0_wp - omega
        cell%length = 1.0_wp
        cell%c = [0.0_wp,1.0_wp,-1.0_wp]
        t = weights_(cell%c)
        cell%pdf = t*dens
    end function

    pure function new_SmanCellComplex(omega,dens,length,c) result(cell)
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: dens
        real(wp), intent(in) :: length
        real(wp), intent(in) :: c(3)
        type(SmanCell) :: cell
        real(wp) :: t(3)

        cell%omega = omega
        cell%one_minus_omega = 1.0_wp - omega
        cell%length = length
        cell%c = c
        t = weights_(cell%c)
        cell%pdf = t*dens
    end function

    pure elemental function first_moment(this) result(moment)
        class(SmanCell), intent(in) :: this
        real(wp) :: moment
        moment = sum(this%c*this%pdf)
    end function

    subroutine collide_sman(this,vel,dt,dens)
        class(SmanCell), intent(inout) :: this
        real(wp), intent(in) :: vel(2)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: dens
        real(wp) :: m0, m1, t(3), c2,c3, csum

        ! pre-collision moments (calculated from pdf's), Eq. 13
        m0 = this%zeroth_moment()
        m1 = this%first_moment() - m0*this%c(1)
        if (present(dens)) dens = m0

        ! update length and velocity vectors, Eq. 14
        this%length = this%length + (vel(1)-vel(2))*dt
        this%c(1) = 0.5*sum(vel)
        this%c(2) = this%c(1) + this%length/dt
        this%c(3) = this%c(1) - this%length/dt

        ! post-collision moments, ! Eq. 15
        ! m0 = m0, mass is conserved in this case
        m1 = m1*this%one_minus_omega

        ! new weight factors, Eq. 17, see above
        c2 = abs(this%c(2))
        c3 = abs(this%c(3))
        csum = c2 + c3

        t(1) = 1.0_wp - C_SQR/(c2*c3)
        t(2) = C_SQR/(c2*csum)
        t(3) = C_SQR/(c3*csum)

        ! post-collision pdf's, Eq. 19
        this%pdf(1) = t(1)*m0
        this%pdf(3) = t(2)*m0 + m1/csum ! NOTE: These are opposite on purpose.
        this%pdf(2) = t(3)*m0 - m1/csum
    end subroutine


    subroutine swap_with_sman(this,previous)
        class(SmanCell), intent(inout) :: this
        class(AbstractCell), intent(inout) :: previous
        real(wp) :: temp
        real(wp) :: c2, c3, csum, A, B2, B3

        select type(previous)
            class is (SmanCell)
                c2 = abs(previous%c(2))
                c3 = abs(this%c(3))
                csum = c2 + c3
                A = (c2 - c3)/csum
                B2 = 2.0_wp*c2/csum
                B3 = 2.0_wp*c3/csum

                temp = this%pdf(2)
                this%pdf(2) = B2*previous%pdf(3) + A*temp
                previous%pdf(3) = B3*temp - A*previous%pdf(3)
            class default
                print *, "unsupported swapping of classes"
        end select
    end subroutine

    subroutine force_west(this,value)
        class(SmanCell), intent(inout) :: this
        real(wp), intent(in) :: value

        this%pdf(2) = value - this%pdf(1) - this%pdf(3)
    end subroutine

    subroutine force_east(this,value)
        class(SmanCell), intent(inout) :: this
        real(wp), intent(in) :: value

        this%pdf(3) = value - this%pdf(1) - this%pdf(2)
    end subroutine
end module