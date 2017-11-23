module surrogate_mod

    implicit none
    private

    ! This stateless type serves only for purposes of extension by other types.
    ! In such a role, it can serve as a substitute for the child type when that
    ! type is inaccessible because of Fortran’s prohibition against circular references.

    type, abstract, public :: surrogate
    end type
end module


module integrators_mod
    
    use surrogate_mod, only : surrogate ! Substitute for integrable_model (avoid circular reference)

    implicit none
    private

    type, abstract, public :: integrator ! Abstract time integration strategy
    contains
        procedure(integrator_interface), nopass, deferred :: integrate ! Abstract integration procedure
    end type

    abstract interface
        subroutine integrator_interface(this,dt)
            import :: surrogate
            class(surrogate), intent(inout) :: this ! integrand
            real, intent(in) :: dt ! time step
        end subroutine
    end interface
end module


module integrable_model_mod
    use surrogate_mod, only : surrogate
    use integrators_mod, only: integrator
    implicit none

    private
    public :: integrate

    type, extends(surrogate), abstract, public:: integrable_model
        class(integrator), allocatable :: quadrature
    contains
        procedure :: integrate
        procedure(time_derivative), deferred :: d_dt
        procedure(symmetric_operator), deferred :: add
        procedure(symmetric_assignment), deferred :: assign
        procedure(asymmetric_operator), deferred :: multiply
        generic :: operator(+) => add
        generic :: operator(*) => multiply
        generic :: assignment(=) => assign
    end type integrable_model

    abstract interface
        function time_derivative(this) result(dthis_dt)
            import :: integrable_model
            class(integrable_model), intent(in) :: this
            class(integrable_model), allocatable :: dthis_dt
        end function
        function symmetric_operator(lhs,rhs) result(total)
            import :: integrable_model
            class(integrable_model), intent(in) :: lhs, rhs
            class(integrable_model), allocatable :: total
        end function
        function asymmetric_operator(lhs,rhs) result(total)
            import :: integrable_model
            class(integrable_model), intent(in) :: lhs
            real, intent(in) :: rhs
            class(integrable_model), allocatable :: total
        end function
        subroutine symmetric_assignment(lhs,rhs)
            import :: integrable_model
            class(integrable_model), intent(in) :: rhs
            class(integrable_model), intent(inout) :: lhs
        end subroutine
    end interface

contains
    
    subroutine integrate(model,dt)
        class(integrable_model), intent(inout) :: model
        real, intent(in) :: dt
        call model%quadrature%integrate(model,dt)
    end subroutine
end module

module explicit_euler_mod

    use surrogate_mod, only : surrogate ! integrable model parent
    use integrators_mod, only : integrator ! time integration strategy
    use integrable_model_mod, only : integrable_model ! abstract integrand

    implicit none
    private

    type, extends(integrator), public :: explicit_euler ! 1st-order explicit time integrator
    contains
        procedure, nopass :: integrate
    end type

 ! $5$
contains

    subroutine integrate(this,dt)
        class(surrogate), intent(inout) :: this
        real, intent(in) :: dt
        select type(this)
            class is (integrable_model)
                this = this + this%d_dt()*dt ! Explicit Euler formula
            class default
                stop "integrate: unsupported class"
        end select
    end subroutine
end module


module rk2_mod

    use surrogate_mod, only : surrogate ! integrable model parent
    use integrators_mod, only : integrator ! time integration strategy
    use integrable_model_mod, only : integrable_model ! abstract integrand

    implicit none
    private

    type, extends(integrator), public :: rk2 ! 2nd-order explicit time integrator
    contains
        procedure, nopass :: integrate
    end type

contains

    subroutine integrate(this,dt)
        class(surrogate), intent(inout) :: this
        real, intent(in) :: dt
        class(integrable_model), allocatable :: this_half
        select type(this)
            class is (integrable_model)
                allocate(this_half,source=this)
                this_half = this + this%d_dt()*(0.5*dt) ! predictor step
                this = this + this_half%d_dt()*dt ! corrector step
            class default
                stop "integrate: unsupported class"
        end select
    end subroutine
end module

module rk4_mod

    use surrogate_mod, only : surrogate ! integrable model parent
    use integrators_mod, only : integrator ! time integration strategy
    use integrable_model_mod, only : integrable_model ! abstract integrand

    implicit none
    private

    type, extends(integrator), public :: rk4 ! 2nd-order explicit time integrator
    contains
        procedure, nopass :: integrate
    end type

contains

    subroutine integrate(this,dt)
        class(surrogate), intent(inout) :: this
        real, intent(in) :: dt
        class(integrable_model), allocatable :: half, half2, full
        select type(this)
            class is (integrable_model)
                allocate(half,source=this)
                allocate(half2,source=this)
                allocate(full,source=this)
                half = this + this%d_dt()*(0.5*dt)
                half2 = this + half%d_dt()*(0.5*dt)
                full = this + half2%d_dt()*dt
                this = this + (this%d_dt()+half%d_dt()*2.0+half2%d_dt()*2.0+full%d_dt())*(dt/6.)
            class default
                stop "integrate: unsupported class"
        end select
    end subroutine
end module

module lorenz_mod

    use integrators_mod, only : integrator
    use integrable_model_mod, only : integrable_model

    implicit none
    private
    public :: integrable_model

    type, extends(integrable_model), public :: lorenz
        private
        real, allocatable :: state(:) ! solution vector
        real :: sigma, rho, beta ! Lorenz parameters
    contains
        procedure, public :: d_dt => dlorenz_dt
        procedure, public :: add => add_lorenz
        procedure, public :: multiply => multiply_lorenz
        procedure, public :: assign => assign_lorenz
        procedure, public :: construct
        procedure, public :: output
    end type

contains

    subroutine construct(this,initial_state,s,r,b,strategy)
        class(lorenz), intent(out) :: this
        real, intent(in) :: initial_state(:)
        real, intent(in) :: s, r, b
        class(integrator), intent(in) :: strategy

        this%state = initial_state
        this%sigma = s
        this%rho = r
        this%beta = b
        allocate(this%quadrature, source=strategy)
    end subroutine

    function dlorenz_dt(this) result(dv_dt)
        class(lorenz), intent(in) :: this
        class(integrable_model), allocatable :: dv_dt
        type(lorenz), allocatable :: local_dv_dt

        allocate(local_dv_dt)
        allocate(local_dv_dt%state(size(this%state)))

        local_dv_dt%state(1) = this%sigma*( this%state(2) -this%state(1))
        local_dv_dt%state(2) = this%state(1)*(this%rho-this%state(3))-this%state(2)
        local_dv_dt%state(3) = this%state(1)*this%state(2)-this%beta*this%state(3)

        call move_alloc(local_dv_dt,dv_dt)
    end function

    function add_lorenz(lhs,rhs) result(sum)
        class(lorenz), intent(in) :: lhs
        class(integrable_model), intent(in) :: rhs
        class(integrable_model), allocatable :: sum
        type(lorenz), allocatable :: local_sum

        select type(rhs)
            class is (lorenz)
                allocate(local_sum)
                local_sum%state = lhs%state + rhs%state
                local_sum%sigma = lhs%sigma
                local_sum%rho = lhs%rho
                local_sum%beta = lhs%beta
            class default
                stop "add_lorenz: unsupported class"
        end select
        call move_alloc(local_sum,sum)
    end function

    function multiply_lorenz(lhs,rhs) result(product)
        class(lorenz), intent(in) :: lhs
        real, intent(in) :: rhs
        class(integrable_model), allocatable :: product
        type(lorenz), allocatable :: local_product

        allocate(local_product)
        local_product%state = lhs%state*rhs
        call move_alloc(local_product,product)
    end function

    subroutine assign_lorenz(lhs,rhs)
        class(lorenz), intent(inout) :: lhs
        class(integrable_model), intent(in) :: rhs
        select type(rhs)
            class is (lorenz)
                lhs%state = rhs%state
                lhs%sigma = rhs%sigma
                lhs%rho = rhs%rho
                lhs%beta = rhs%beta
            class default
                stop "assign_lorenz: unsupported class"
        end select
    end subroutine


    function output(this) result(coordinates)
        class(lorenz), intent(in) :: this
        real, allocatable :: coordinates(:)
        coordinates = this%state
    end function

end module

module ode_solver_mod
    use explicit_euler_mod, only: explicit_euler
    use rk2_mod, only: rk2
    use rk4_mod, only: rk4
end module

program main
    use lorenz_mod, only : lorenz
    use ode_solver_mod
    implicit none

    type(rk4) :: integrator
    type(lorenz) :: attractor
    integer :: step

    real, parameter :: sigma = 10., rho = 28., beta = 8./3., dt = 0.0001
    integer, parameter :: num_steps = int((80)/dt)

    real, parameter :: y0(3) = [100.,60.,30.]

    call attractor%construct(y0,sigma,rho,beta,integrator)


    print *, attractor%output()
    do step = 1, num_steps
        call attractor%integrate(dt)
        print *, attractor%output()
    end do
end program
