module sman_cell_mod
    use precision_mod, only : wp
    use cell_mod, only : LatticeCell

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
        procedure, public :: collide
    end type

    interface SmanCell
        module procedure new_SmanCell
    end interface
contains

    elemental subroutine assign_cell(lhs,rhs)
        class(LatticeCell), intent(inout) :: lhs
        class(LatticeCell), intent(in) :: rhs

        lhs%omega = rhs%omega
        lhs%pdf = rhs%pdf
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

    pure elemental function first_moment(this) result(moment)
        class(SmanCell), intent(in) :: this
        real(wp) :: moment
        moment = sum(this%c*this%pdf)
    end function

    subroutine collide(this,vel,dt,dens)
        class(SmanCell), intent(inout) :: this
        real(wp), intent(in) :: vel(2)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: dens
        real(wp) :: m0, m1, c2, c3, csum, t1, t2, t3

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
        t = weights_()

        ! post-collision pdf's, Eq. 19
        this%pdf(1) = t(1)*m0
        this%pdf(2) = t(2)*m0 + m1/csum
        this%pdf(3) = t(3)*m0 - m1/csum
    end subroutine

end module