module math_mod
    use precision_mod, only : wp
    implicit none
    private
    public pi

    real(wp), parameter :: pi = 4._wp*atan(1.0_wp)

end module