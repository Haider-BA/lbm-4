module precision_mod
    
    implicit none
    private
    public sp,dp,wp

    integer, parameter :: sp = kind(1.0e0)
    integer, parameter :: dp = kind(1.0d0)

    integer, parameter :: wp = dp
end module