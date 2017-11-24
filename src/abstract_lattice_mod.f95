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
    end interface
end module