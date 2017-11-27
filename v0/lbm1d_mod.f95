module lbm1d_mod

    use precision_mod, only: wp
    use lbm_mod, only: Lattice, diffToOmega, bc
    use io_mod, only: savetxt

    implicit none

    real(wp), parameter :: PI = 3.14159265358979323846_wp

end module