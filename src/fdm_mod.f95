module fdm_mod
    use precision_mod, only: wp

    implicit none
    private
    public three_point_difference

contains

    pure function three_point_difference(f,x,i) result(df)
        real(wp), intent(in) :: f(3)
        real(wp), intent(in) :: x(3)
        integer, intent(in) :: i
        real(wp) :: df

        real(wp) :: h1, h2, hsum, hprod

        h1 = x(2) - x(1)
        h2 = x(3) - x(1)
        hsum = h1 + h2
        hprod = h1*h2

        df = 0.0_wp
        select case (i)
            case (1)
                df = -(2*h1+h2)/(h1*hsum)*f(1) + hsum/hprod*f(2) - h1/(hsum*h2)*f(3)
            case (2)
                df = -h2/(h1*hsum)*f(1) - (h1-h2)/hprod*f(2) + h1/(hsum*h2)*f(3);
            case (3)
                df = h2/(h1*hsum)*f(1) - hsum/hprod*f(2) - (h1+2*h2)/(hsum*h2)*f(3)
            case default
                return
        end select 
    end function

    ! pure function four_point_difference(f,x,i) result(df)
    !     real(wp), intent(in) :: f(4)
    !     real(wp), intent(in) :: x(4)
    !     integer, intent(in) :: i
    !     real(wp) :: df

    !     df = 0.0_wp ! TO DO

    ! end function

end module