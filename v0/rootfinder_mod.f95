module rootfinder_mod

    use precision_mod, only: wp

    implicit none

    private
    public user_function, simple_function, secant, newton

    type, abstract, public :: user_function
    contains
        procedure(function_evaluation), deferred, pass(params) :: eval
    end type

    abstract interface
        function function_evaluation(x,params) result(fx)
            import :: newton_function, wp
            real(wp) :: x
            class(user_function) :: params
            real(wp) :: fx
        end function
    end interface

    abstract interface
        pure function simple_function(x)
            use precision_mod, only: wp
            real(wp), intent(in) :: x
            real(wp) :: newton_func
        end function
    end interface

contains

    function secant(f,x0,tol,max_iter) result(p)
        class(user_function), intent(in) :: f
        real(wp), intent(in) :: x0
        real(wp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter

        real(wp) :: tol_
        integer :: max_iter_

        real(wp) :: p0, p1, q0, q1, p

        ! set tolerance and maximum number of iterations
        tol_ = 1.0e-7_wp
        if (present(tol)) tol_ = tol
        max_iter_ = 100
        if (present(max_iter)) max_iter_ = max_iter

        ! print warning
        if (tol_ <= tiny(1.0_wp)) print *, "[secant] Tolerance is too small (",tol," <= 0)."
        if (maxiter_ < 1) print *, "[secant] 'max_iter' must be greater than 0."

        ! secant method
        p0 = x0
        if (p0 >= 0) then
            p1 = x0*1.0001_wp + 0.0001_wp
        else
            p1 = x0*1.0001_wp - 0.0001_wp
        end if

        q0 = f%eval(p0)
        q1 = f%eval(p1)

        do i = 1, maxiter_
            if (q1 == q0) then
                if (p1 /= p0) then
                    print *, "[secant] Tolerance of ", p1 - p0," reached."
                end if
                p = 0.5_wp*(p1 + p0)
                return
            else
                p = p1 - q1*(p1 - p0)/(q1 - q0)
                if (abs(p - p1) < tol_) return ! found zero
            end if
            p0 = p1
            q0 = q1
            p1 = p
            q1 = f%eval(p1)
        end do
    end function

    function newton(f,x0,fprime,tol,maxiter,fprime2) result(p)
        procedure(simple_function), pointer, intent(in) :: f
        real(wp), intent(in) :: x0
        procedure(simple_function), optional :: fprime
        real(wp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        procedure(simple_function), optional :: fprime2
        real(wp) :: p

        real(wp) :: tol_ = 1.48e-7_wp
        integer :: maxiter_ = 100
        integer :: i

        real(wp) :: p0, p1, q0, q1
        real(wp) :: fval, fder, fder2, discr

        if (present(tol)) tol_ = tol
        if (present(maxiter)) maxiter_ = maxiter

        if (tol_ <= 0) print *, "[newton] Tolerance is too small (",tol," <= 0)"
        if (maxiter_ < 1) print *, "[newton] maxiter must be greater than 0"

        if (present(fprime)) then
            ! newton method
            p0 = x0
            do i = 1, maxiter_
                fder = fprime(p0)
                if (fder == 0) then
                    print *, "[newton] derivative was zero."
                    p = p0
                    return
                end if
                fval = f(p0)
                if (.not. present(fprime2)) then
                    ! Newton step
                    p = p0 - fval/fder
                else
                    ! Halley's method
                    fder2 = fprime2(p)
                    discr = fder**2 - 2.0_wp*fval*fder2
                    if (discr < 0 ) then
                        p = p0 - fder/fder2
                    else
                        p = p0 - 2._wp*fval/(fder + sign(sqrt(discr),fder))
                    end if
                end if
                if (abs(p - p0) < tol_) return ! found zero
                p0 = p
            end do
        else
            ! secant method
            p0 = x0
            if (p0 >= 0) then
                p1 = x0*1.0001_wp + 0.0001_wp
            else
                p1 = x0*1.0001_wp - 0.0001_wp
            end if

            q0 = f(p0)
            q1 = f(p1)

            do i = 1, maxiter_
                if (q1 == q0) then
                    if (p1 /= p0) then
                        print *, "[newton] Tolerance of ", p1 - p0," reached"
                    end if
                    p = 0.5_wp*(p1 + p0)
                    return
                else
                    p = p1 - q1*(p1 - p0)/(q1 - q0)
                    if (abs(p - p1) < tol_) return ! found zero
                end if
                p0 = p1
                q0 = q1
                p1 = p
                q1 = f(p1)
            end do
        end if

        print *, "[newton] Failed to converge after ", maxiter_," iterations, &
            & value is ", p

    end function newton


end module