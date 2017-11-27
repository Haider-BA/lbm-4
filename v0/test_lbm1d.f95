program test_lbm1d
    
    use lbm1d_mod

    implicit none

    integer, parameter :: n = 101
    integer, parameter :: tmax = 100000
    integer, parameter :: tout = 1000
    character(len=*), parameter :: myfile = "density"

    real(wp) :: dens(n), x(n)
    real(wp) :: omega, dx

    integer :: t, i

    type(Lattice) :: moisture
    procedure(bc), pointer :: rbc => null()

    ! assign right BC
    rbc => reactionbc

    ! grid
    dx = 1.0_wp/(n-1)
    x = [(dx*real((i-1),wp), i=1, n)]
    
    ! initial density
    dens = sin(2.0_wp*PI*x)**2

    ! initialize lattice
    t = 0
    omega = 1.0_wp
    moisture = Lattice(n,omega,dens)

    dens = moisture%getDensity()
    call savetxt(myfile,t,dens,dx)

    ! time loop
    do t = 1, tmax

        ! print output
        if (mod(t,tout) == 0) then
            dens = moisture%getDensity()
            call savetxt(myfile,t,dens,dx)
        end if
        call moisture%collide()
        call moisture%addSource(-0.000001_wp*moisture%getDensity())
        call moisture%stream()
        call moisture%setLeftBC(1.0_wp)
        !call moisture%setRightBC(0.4_wp)
        call moisture%setRightBC(rbc)
        
        call moisture%prepNextStep()
    
    end do

    dens = moisture%getDensity()
    call savetxt(myfile,t,dens,dx)

contains

    function reactionbc(u) result(zero)
        ! Function to calculate the value of the boundary condition
        !
        ! The boundary functions is:
        !   $D \frac{\partial u}{\partial x} = 0$
        ! $$
        real(wp), intent(in) :: u
        real(wp) :: zero

        zero = moisture%getDiffusivity()*moisture%getRightGrad(u)

    end function

end program

 