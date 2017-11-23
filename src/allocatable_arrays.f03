module point_mod

    implicit none

    private

    public wp, Point_t

    ! kind parameter to define precision
    integer, parameter :: wp = kind(1.0e0)

    ! A simple datatype defining a point
    type :: Point_t
        real(wp) :: x, y
    contains
        procedure :: print
        procedure :: add
        procedure :: add_and_return_two
        generic :: operator(+) => add
    end type Point_t

    ! Constructors
    interface Point_t
        module procedure new_Point
        module procedure new_Point_vec    
    end interface

contains

    function new_Point(a_,b_) result(this)
        real(wp), intent(in) :: a_, b_
        type(Point_t) :: this

        this%x = a_
        this%y = b_
    end function

    function new_Point_vec(vec) result(this)
        real(wp), intent(in) :: vec(2)
        type(Point_t) :: this

        this%x = vec(1)
        this%y = vec(2)
    end function

    ! this pure elemental stuff is the shizzle ;)
    pure elemental function add(lhs,rhs) result(sum)
        class(Point_t), intent(in) :: lhs
        class(Point_t), intent(in) :: rhs
        type(Point_t) :: sum

        sum%x = lhs%x + rhs%x
        sum%y = lhs%y + rhs%y
    end function

    subroutine print(this)
        class(Point_t), intent(in) :: this

        print *, "x = ", this%x, "y = ", this%y
    end subroutine


    elemental function add_and_return_two(lhs,rhs,two) result(sum)
        class(Point_t), intent(in) :: lhs,rhs
        type(Point_t) :: sum
        real(wp), value :: two
        sum%x = lhs%x + rhs%x
        sum%y = lhs%y + rhs%y
        
            two = 2.0_wp
        
    end function
end module

program main
    
    use point_mod

    implicit none

    integer :: i
    real(wp), allocatable :: v(:),two(:)

    type(Point_t) :: p1, p2, p3
    type(Point_t), allocatable :: pointset(:), other_set(:)

    ! 6 elements are allocated
    v = [sin([3.,4.,5.]),7.,8.,9.]
    print *, v
    print *, size(v)

    ! shortened down to 3 elements
    v = v(1:3)
    print *, v
    print *, size(v)

    ! Test point type
    p1 = Point_t(2.,3.) ! 1st constructor
    p2 = Point_t([8.,9.]) ! 2nd constructor
    call p1%print
    call p2%print

    ! Test summation operator
    p3 = p1 + p2 ! operator(+) - add
    call p3%print

    ! Test a point array (implicitly allocated)
    pointset = [p1,p2,p3]
    print *, pointset
    print *, size(pointset)

    ! Test summation operator on array of points
    pointset = pointset + pointset 
    print *, pointset
    print *, size(pointset)

    ! Try enlarging the point set with the old valued 
    pointset = [pointset,p1,p2,p3]
    do i = 1, size(pointset)
        print *, "Point ", i,":"
        call pointset(i)%print
    end do

    allocate(two(size(pointset)))
    other_set = pointset%add_and_return_two(pointset,two)
    print *, two

end program