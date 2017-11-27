module boundary_condition_mod
    use precision_mod, only : wp
    use rootfinder_mod, only : UserFunction, secant
    use lagrange_mod, only : LagrangePolynomial

    implicit none
    private

    public DirichletCondition, FluxCondition, UserFunction,new_Flux

    type, extends(UserFunction) :: DirichletCondition
        real(wp) :: cwall
        real(wp) :: diff, dt, area
        real(wp) :: pdf(2)
        real(wp) :: old
        real(wp), allocatable :: grid(:)
        real(wp), allocatable :: conc(:)
    contains
        procedure, pass(params) :: eval => evaluate_dirichlet_bc
        procedure :: assign => assign_dbc
        generic :: assignment(=) => assign
    end type

    type, extends(UserFunction) :: FluxCondition
        class(UserFunction), pointer :: flux => null()
        real(wp) :: dt, area
        real(wp) :: pdf(2)
        real(wp) :: old
        real(wp), allocatable :: grid(:)
        real(wp), allocatable :: conc(:)
    contains
        procedure, pass(params) :: eval => evaluate_flux_bc
        procedure :: assign => assign_fbc
        generic :: assignment(=) => assign
    end type

    interface DirichletCondition
        module procedure new_Dirichlet
    end interface

    interface FluxCondition
        module procedure new_Flux
    end interface

contains

    function new_Dirichlet(cwall,diff,dt,area) result(bc)
        real(wp), intent(in) :: cwall
        real(wp), intent(in) :: diff, dt, area
        type(DirichletCondition) :: bc

        bc%cwall = cwall
        bc%diff = diff
        bc%dt = dt
        bc%area = area
        bc%pdf = 0.0_wp
        bc%old = 0.0_wp
        allocate(bc%grid(5))
        allocate(bc%conc(5))
    end function

    function new_Flux(dt,area,flux) result(bc)
        real(wp), intent(in) :: dt, area
        class(UserFunction), pointer :: flux
        type(FluxCondition) :: bc

        bc%dt = dt
        bc%area = area
        bc%flux => flux
        bc%old = 0.0_wp
        bc%pdf = 0.0_wp
        allocate(bc%grid(5))
        allocate(bc%conc(5))
    end function

    subroutine assign_dbc(lhs,rhs)
        class(DirichletCondition), intent(inout) :: lhs
        class(UserFunction), intent(in) :: rhs
        select type(rhs)
            class is (DirichletCondition)
                lhs%cwall = rhs%cwall
                lhs%dt = rhs%dt
                lhs%area = rhs%area
                lhs%diff = rhs%diff
                lhs%pdf = rhs%pdf
                lhs%old = rhs%old
                call move_alloc(rhs%grid,lhs%grid)
                call move_alloc(rhs%conc,lhs%conc)
            class default
                print *, "not supported"
        end select
        ! lhs%eval => rhs%eval
    end subroutine

    subroutine assign_fbc(lhs,rhs)
        class(FluxCondition), intent(inout) :: lhs
        class(UserFunction), intent(in) :: rhs
        select type(rhs)
            class is (FluxCondition)
                lhs%dt = rhs%dt
                lhs%area = rhs%area
                lhs%pdf = rhs%pdf
                lhs%old = rhs%old
                lhs%flux => rhs%flux
                call move_alloc(rhs%grid,lhs%grid)
                call move_alloc(rhs%conc,lhs%conc)
            class default
                print *, "not supported"
        end select
        ! lhs%eval => rhs%eval
    end subroutine

    function evaluate_dirichlet_bc(x,params) result(fx)
        real(wp), intent(in) :: x
        class(DirichletCondition), intent(in) :: params
        real(wp) :: fx

        type(LagrangePolynomial) :: p
        real(wp) :: flux, conc(size(params%grid)), c0

        c0 = x + sum(params%pdf) ! zeroth_moment
        conc(1) = params%cwall
        conc(2) = c0
        conc(3:) = params%conc(3:)
        p = LagrangePolynomial(params%grid,conc) ! TODO: use barycentric interpolant
        ! print *, conc
        flux = params%diff*params%area*p%d_dx(params%grid(1))
        ! print *, flux
        fx = x - params%old + flux*params%dt
        ! print *, fx
    end function

    function evaluate_flux_bc(x,params) result(fx)
        real(wp), intent(in) :: x
        class(FluxCondition), intent(in) :: params
        real(wp) :: fx
        type(LagrangePolynomial) :: p
        real(wp) :: flux, conc(size(params%grid)), cw

        conc = params%conc
        conc(2) = x + sum(params%pdf) ! zeroth_moment
        p = LagrangePolynomial(params%grid(2:),conc(2:)) ! TODO: use barycentric interpolant
        cw = p%p(params%grid(1)) ! extrapolate wall concentration
        flux = params%flux%eval(cw)
        fx = x - params%old + flux*params%dt
        print *, fx
    end function
end module