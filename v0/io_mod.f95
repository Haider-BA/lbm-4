module io_mod
    
    use precision_mod, only: wp
    use newunit_mod, only: newunit
    use utils_mod, only: str

    implicit none

    private
    public savetxt

    interface savetxt
        module procedure saveArrayToTxt
        module procedure saveArrayToTxtAtTime
        module procedure saveMultiArrayToTxt
        module procedure saveMultiArrayToTxtAtTime
    end interface

contains

    subroutine saveArrayToTxt(filename,data,dx)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: data(:)
        real(wp), intent(in), optional :: dx

        integer :: i, unit

        open(newunit(unit),file=filename//".txt")
        do i = 1, size(data)
            if (present(dx)) then
                write(unit,*) real(i-1,wp)*dx, data(i)
            else
                write(unit,*) data(i)
            end if
        end do
        close(unit)

    end subroutine saveArrayToTxt

    
    subroutine saveArrayToTxtAtTime(filename,t,data,dx)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: t
        real(wp), intent(in) :: data(:)
        real(wp), intent(in), optional :: dx

        integer :: i, unit

        open(newunit(unit),file=filename//"_"//str(t)//".txt")
        do i = 1, size(data)
            if (present(dx)) then
                write(unit,*) real(i-1,wp)*dx, data(i)
            else
                write(unit,*) data(i)
            end if
        end do
        close(unit)

    end subroutine saveArrayToTxtAtTime


    subroutine saveMultiArrayToTxt(filename,data,dx)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: data(:,:)
        real(wp), intent(in), optional :: dx

        integer :: i, unit

        open(newunit(unit),file=filename//".txt")
        do i = 1, size(data,dim=1)
            if (present(dx)) then
                write(unit,*) real(i-1,wp)*dx, data(i,:)
            else
                write(unit,*) data(i,:)
            end if
        end do
        close(unit)

    end subroutine saveMultiArrayToTxt


    subroutine saveMultiArrayToTxtAtTime(filename,t,data,dx)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: t
        real(wp), intent(in) :: data(:,:)
        real(wp), intent(in), optional :: dx

        integer :: i, unit

        open(newunit(unit),file=filename//"_"//str(t)//".txt")
        do i = 1, size(data)
            if (present(dx)) then
                write(unit,*) real(i-1,wp)*dx, data(i,:)
            else
                write(unit,*) data(i,:)
            end if
        end do
        close(unit)

    end subroutine saveMultiArrayToTxtAtTime

end module