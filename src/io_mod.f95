module io_mod

    use precision_mod, only : wp
    use utils_mod, only : str

    implicit none
    private

    public savetxt, str

    interface savetxt
        module procedure savetxt_simple
        module procedure savetxt_with_step
    end interface

contains

    subroutine savetxt_simple(fname,data,numbered)
        character(len=*), intent(in) :: fname
        real(wp), intent(in) :: data(:)
        logical, intent(in) :: numbered
        integer :: i, unit

        open(newunit=unit,file=fname//".out")
        if (numbered) then
            do i = 1, size(data)
                write(unit,*) i, data(i)
            end do
        else
            do i = 1, size(data)
                write(unit,*) data(i)
            end do
        end if
        close(unit)
    end subroutine


    subroutine savetxt_with_step(fname,data,dx)
        character(len=*), intent(in) :: fname
        real(wp), intent(in) :: data(:)
        real(wp), intent(in), optional :: dx 
        real(wp) :: dx_
        integer :: i, unit

        dx_ = 1.0_wp
        if (present(dx)) dx_ = dx

        open(newunit=unit,file=fname//".out")
        do i = 1, size(data)
            write(unit,*) real(i-1,wp)*dx_, data(i)
        end do
        close(unit)
    end subroutine

end module