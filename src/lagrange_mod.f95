module lagrange_mod
    use precision_mod, only: wp

    implicit none
    private

    public LagrangePolynomial, linspace

    type, public :: LagrangePolynomial
        private
        integer :: n
        real(wp), allocatable :: x(:), y(:)
    contains
            procedure :: p => evaluate
            procedure :: d_dx
    end type

    interface LagrangePolynomial
        module procedure new_LagrangePolynomial
        module procedure new_LagrangePolynomial_safe
    end interface

contains

        pure function linspace(a,b,n) result(vals)
            real(wp), intent(in) :: a
            real(wp), intent(in) :: b
            integer, intent(in) :: n
            real(wp) :: vals(n)
            integer :: i
            real(wp) :: dx

            dx = (b-a)/real(n-1,wp)
            vals = [(a+real(i-1,wp)*dx, i = 1, n)]
        end function linspace

        pure elemental function evaluate(this,x) result(val)
            class(LagrangePolynomial), intent(in) :: this
            real(wp), intent(in) :: x
            real(wp) :: val, prod
            integer :: i, j

            val = 0.0_wp
            do i = 1, this%n
                prod = 1.0_wp
                do j = 1, this%n
                    if (i /= j) then
                        prod = prod*(x-this%x(j))/(this%x(i)-this%x(j))
                    end if
                end do
                val = val + this%y(i)*prod
            end do
        end function

        pure elemental function d_dx(this,x) result(dp)
            class(LagrangePolynomial), intent(in) :: this
            real(wp), intent(in) :: x
            real(wp) :: sum, prod
            real(wp) :: dp
            integer :: j, l, m

            dp = 0.0_wp
            do j = 1, this%n
                sum = 0.0_wp
                do l = 1, this%n
                    prod = 1.0_wp
                    do m = 1, this%n
                        if (m==l .or. m == j) then
                            continue
                        else
                            prod = prod*(x - this%x(m))/(this%x(j)-this%x(m))
                        end if
                    end do
                    if (l == j) then
                        continue
                    else
                        sum = sum + prod/(this%x(j)-this%x(l))
                    end if
                end do

                dp = dp + sum*this%y(j)
            end do
        end function

        function new_LagrangePolynomial(x,y) result(p)
            real(wp), intent(in) :: x(:),y(:)
            type(LagrangePolynomial) :: p
            integer :: n

            n = min(size(x),size(y))
            p%n = n

            allocate(p%x(n))
            allocate(p%y(n))
            p%x = x
            p%y = y
        end function

        function new_LagrangePolynomial_safe(n,x,y) result(p)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n),y(n)
            type(LagrangePolynomial) :: p

            p%n = n
            allocate(p%x(n))
            allocate(p%y(n))
            p%x = x
            p%y = y
        end function
end module


program interpolate
    use precision_mod, only : wp
    use lagrange_mod, only : LagrangePolynomial, linspace
    use math_mod, only : pi
    implicit none

    integer, parameter :: n = 5


    real(wp), allocatable :: x(:), y(:)
    type(LagrangePolynomial) :: p
    integer :: i


    x = linspace(0.0_wp,1.0_wp,n)
    y = x**2 + 3.0_wp*x + 1.0_wp

    p = LagrangePolynomial(n,x,y)

    do i = 1, n
        print *, x(i), y(i), 2.0_wp*x(i)+3.0_wp, p%d_dx(x(i))
    end do
end program
