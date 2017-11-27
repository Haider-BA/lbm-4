program heat_equation
    
    use lbm1d_mod

    implicit none

    integer, parameter :: n = 101
    integer, parameter :: tmax = 1000
    integer, parameter :: tout = 100
    character(len=*), parameter :: myfile = "T"

    real(wp) :: temp(n), x(n)
    real(wp) :: omega, dx
    real(wp) :: diff, Biot

    integer :: t, i

    type(Lattice) :: temperature
    procedure(bc), pointer :: rbc => null()

    ! physical parameters
    Biot = 0.1_wp
    diff = 1.0_wp

    ! grid
    dx = 1.0_wp/(n-1)
    x = [(dx*real((i-1),wp), i=1, n)]
    
    ! initial density
    temp = 1.0_wp

    ! initialize lattice
    t = 0
    omega = diffToOmega(diff)
    temperature = Lattice(n,omega,temp)


    print *, "omega =", temperature%getOmega()
    print *, "diff =", temperature%getDiffusivity()

    ! assign right BC
    rbc => convection

    temp = temperature%getDensity()
    call savetxt(myfile,t,temp,dx)

    ! time loop
    do t = 1, tmax

        ! print output
        if (mod(t,tout) == 0) then
            temp = temperature%getDensity()
            call savetxt(myfile,t,temp,dx)
        end if
        call temperature%collide()
        call temperature%stream(periodic=.false.)
        !call moisture%setLeftBC(1.0_wp)
        !call moisture%setRightBC(0.4_wp)
        call temperature%setRightBC(rbc)
        stop
        
        call temperature%prepNextStep()
    
    end do

    temp = temperature%getDensity()
    call savetxt(myfile,t,temp,dx)

contains

    function convection(T) result(zero)
        ! Function to calculate the value of the boundary condition
        !
        ! The boundary functions is:
        !   $D \frac{\partial u}{\partial x} = 0$
        ! $$
        real(wp), intent(in) :: T
        real(wp) :: zero

        zero = temperature%getRightGrad(T) - Biot*T

    end function

end program

 