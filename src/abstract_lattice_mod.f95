module abstract_lattice_mod
    use precision_mod, only : wp

    implicit none
    private

    type, abstract, public :: AbstractLattice
        integer :: n
    contains
    private
        procedure(int1), deferred :: collide_and_stream_simple
        procedure(int2), deferred :: collide_and_stream_density
        procedure(dens_interface), deferred :: density
    end type

    abstract interface
        subroutine int1(this)
            import :: AbstractLattice
            class(AbstractLattice), intent(inout) :: this
        end subroutine
        subroutine int2(this,dens)
            import :: AbstractLattice, wp
            class(AbstractLattice), intent(inout) :: this
            real(wp) :: dens(this%n)
        end subroutine
        function dens_interface(this) result(dens)
            import :: AbstractLattice, wp
            class(AbstractLattice), intent(in) :: this
            real(wp) :: dens(this%n)
        end function
    end interface

end module