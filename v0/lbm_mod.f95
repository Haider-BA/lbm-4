module lbm_mod
    
    use precision_mod, only: wp
    use rootfinder_mod, only: newton, bc => newton_func

    implicit none

    private 
    public Lattice, diffToOmega, bc

    ! lattice weights
    real(wp), parameter :: W(0:2) = [2._wp/3._wp,1._wp/6._wp,1._wp/6._wp]

    ! (lattice speed of sound)**2
    real(wp), parameter :: CSQR = 1._wp/3._wp

    ! lattice directions
    ! <-----o----->
    ! 2     0     1

    type :: Lattice
        private
        integer :: n
        real(wp) :: constOmega
        real(wp), allocatable :: f(:,:,:)
        integer :: old, new
        !procedure(bc), pointer, nopass :: leftbc => null()
        !procedure(bc), pointer, nopass :: rightbc => null()
    contains
        procedure :: addSource
        procedure :: collide
        procedure :: getDensity
        procedure :: getDensityAt
        procedure :: getDiffusivity
        procedure :: getLeftGrad
        procedure :: getRightGrad
        procedure :: getOmega
        procedure :: getSize
        procedure :: getTau
        procedure :: setEquilibriumValues_array
        procedure :: setEquilibriumValues_scalar
          generic :: setEquilibriumValues => setEquilibriumValues_array, &
                                             setEquilibriumValues_scalar
        procedure :: setLeftBC_dirichlet
        procedure :: setRightBC_dirichlet
        procedure :: setRightBC_robyn
          generic :: setLeftBC => setLeftBC_dirichlet
          generic :: setRightBC => setRightBC_dirichlet, setRightBC_robyn                                     
        procedure :: setOmega
        procedure :: setOmegaFromDiffusivity
        procedure :: stream
        procedure :: prepNextStep
    end type


    interface Lattice
        module procedure newLattice
        module procedure newLatticeWithRho
        module procedure newLatticeWithRhoArray
    end interface


contains


    function newLattice(n,omega) result(latt)
        ! Returns a new lattice.
        !
        ! Initializes a new lattice object based on user provided
        ! lattice size, relaxation frequency.
        !
        ! Parameters
        ! ----------
        !   n : integer
        !       The lattice size
        !   omega : scalar
        !       The value of the relaxation frequency.

        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        type(Lattice) :: latt

        latt%n = n
        latt%constOmega = omega
        allocate(latt%f(0:n+1,0:2,2))

        latt%new = 1
        latt%old = 2

        latt%f(:,:,:) = 0.0_wp

    end function newLattice


    function newLatticeWithRho(n,omega,rho) result(latt)
        ! Returns a new lattice.
        !
        ! Initializes a new lattice object based on user provided
        ! lattice size, relaxation frequency and density (optional).
        ! If the density is provided the PDF is initialized with
        ! the equilibrium values.
        !
        ! Parameters
        ! ----------
        !   n : integer
        !       The lattice size
        !   omega : scalar
        !       The value of the relaxation frequency.
        !   rho : scalar, optional
        !       A single initial lattice density.

        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: rho
        type(Lattice) :: latt

        latt%n = n
        latt%constOmega = omega
        allocate(latt%f(0:n+1,0:2,2))

        latt%new = 1
        latt%old = 2

        latt%f(1:n,0,:) = W(0)*rho
        latt%f(1:n,1,:) = W(1)*rho
        latt%f(1:n,2,:) = W(2)*rho

    end function newLatticeWithRHo


    function newLatticeWithRhoArray(n,omega,rho) result(latt)
        ! Returns a new lattice.
        !
        ! Initializes a new lattice object based on user provided
        ! lattice size, relaxation frequency and density (optional).
        ! If the density is provided the PDF is initialized with
        ! the equilibrium values.
        !
        ! Parameters
        ! ----------
        !   n : integer
        !       The lattice size
        !   omega : scalar
        !       The value of the relaxation frequency.
        !   rho : array-like
        !       The initial densities at all positions.

        integer, intent(in) :: n
        real(wp), intent(in) :: omega
        real(wp), intent(in) :: rho(n)
        type(Lattice) :: latt

        latt%n = n
        latt%constOmega = omega
        allocate(latt%f(0:n+1,0:2,2))

        latt%new = 1
        latt%old = 2

        latt%f(1:n,0,latt%old) = W(0)*rho(:)
        latt%f(1:n,1,latt%old) = W(1)*rho(:)
        latt%f(1:n,2,latt%old) = W(2)*rho(:)

        latt%f(:,:,latt%new) = latt%f(:,:,latt%old)
    end function newLatticeWithRHoArray

    subroutine stream(self,periodic)
        ! Performs streaming step on lattice

        class(Lattice), intent(inout) :: self
        logical, intent(in), optional :: periodic
        
        logical :: periodic_ = .false. ! Bounce-back is the default!!!
        integer :: i
        
        if (present(periodic)) periodic_ = periodic

        do i = 1, self%n
            self%f(i+1,1,self%new) = self%f(i,1,self%old) 
            self%f(i-1,2,self%new) = self%f(i,2,self%old)
        end do 

        ! Periodic and bounce-back boundary conditions
        if (periodic_) then
            self%f(1,1,self%new) = self%f(self%n,1,self%old)
            self%f(self%n,2,self%new) = self%f(1,2,self%old)        
        else ! bounce-back
            self%f(1,1,self%new) = self%f(1,2,self%old)
            self%f(self%n,2,self%new) = self%f(self%n,1,self%old)
        end if
    end subroutine stream


    subroutine collide(self)
        ! Performs the collision step on the lattice

        class(Lattice), intent(inout) :: self

        integer :: i
        real(wp) :: omega, oneMinOmega, dens, feq0,feq1,feq2

        oneMinOmega = 1.0_wp - self%constOmega
        omega = self%constOmega

        do i = 1, self%n
            dens = sum(self%f(i,:,self%old))
            feq0 = W(0)*dens
            feq1 = W(1)*dens
            feq2 = W(2)*dens

            self%f(i,0,self%old) = oneMinOmega*self%f(i,0,self%old) + omega*feq0
            self%f(i,1,self%old) = oneMinOmega*self%f(i,1,self%old) + omega*feq1
            self%f(i,2,self%old) = oneMinOmega*self%f(i,2,self%old) + omega*feq2
        end do
    end subroutine collide


    function getDensity(self) result(dens)
        class(Lattice), intent(in) :: self
        real(wp) :: dens(self%n)

        integer :: i

        do i = 1, self%n
            dens(i) = sum(self%f(i,:,self%old))
        end do
    end function getDensity


    function getDensityAt(self,pos) result(dens)
        class(Lattice), intent(in) :: self
        integer, intent(in) :: pos
        real(wp) :: dens

        dens = sum(self%f(pos,:,self%new))
    end function getDensityAt


    subroutine addSource(self,source)
        class(Lattice), intent(inout) :: self
        real(wp) :: source(getSize(self))

        integer :: i

        do i = 1, self%getSize()
            self%f(i,0,self%old) = self%f(i,0,self%old) + W(0)*source(i)
            self%f(i,1,self%old) = self%f(i,1,self%old) + W(1)*source(i)
            self%f(i,2,self%old) = self%f(i,2,self%old) + W(2)*source(i)
        end do
    end subroutine addSource


    pure function getSize(self) result(size)
        ! Returns the lattice size.
        !
        ! Parameters
        ! ----------
        !   self : Lattice
        !       A lattice object
        !
        ! Returns
        ! -------
        !   size : integer
        !       The lattice size.
        
        class(Lattice), intent(in) :: self
        integer :: size
        size = self%n
    end function getSize


    pure function getOmega(self) result(omega)
        ! Returns the lattice relaxation frequency.
        !
        ! Parameters
        ! ----------
        !   self : Lattice
        !       Our lattice object.
        !
        ! Returns
        ! -------
        !   omega : scalar
        !       The value of the relaxation frequency set initially.
        class(Lattice), intent(in) :: self
        real(wp) :: omega
        omega = self%constOmega
    end function getOmega


    pure function getTau(self) result(tau)
        class(Lattice), intent(in) :: self
        real(wp) :: tau
        tau = 1._wp/self%constOmega
    end function getTau


    subroutine prepNextStep(self)
        ! Prepares the lattice for the next timestep
        !
        ! Swaps the source and destination indexes of the lattice
        ! for the next time step. This way we do not have to care about 
        ! PDF values being ever over-written.
        !
        ! Parameters
        ! ----------
        ! self : Lattice
        !   A Lattice object that stores the PDF

        class(Lattice), intent(inout) :: self
        integer :: temp

        temp = self%old
        self%old = self%new
        self%new = temp
    end subroutine prepNextStep


    pure function getDiffusivity(self) result(D)
        class(Lattice), intent(in) :: self
        real(wp) :: D
        D = CSQR*(self%getTau() - 0.5_wp)
    end function getDiffusivity


    elemental function diffToOmega(diff) result(omega)
        ! Calculates the value of the relaxation rate/frequency.
        !
        ! Parameters
        ! ----------
        !   diff : scalar or array-like
        !       The value(s) of the diffusivity.
        ! 
        ! Returns
        ! -------
        !   omega : scalar or array-like
        !       The relaxation frequencies for the given
        !       diffusivity values.
        
        real(wp), intent(in) :: diff
        real(wp) :: omega

        omega = 1._wp/(diff/CSQR + 0.5_wp)
    end function diffToOmega


    subroutine setOmega(self,omega)
        class(Lattice), intent(inout) :: self
        real(wp) :: omega
        self%constOmega = omega
    end subroutine setOmega

    subroutine setOmegaFromDiffusivity(self,diff)
        class(Lattice), intent(inout) :: self
        real(wp) :: diff

        self%constOmega = diffToOmega(diff)

    end subroutine setOmegaFromDiffusivity


    subroutine setEquilibriumValues_scalar(self,rho)
        class(Lattice), intent(inout) :: self
        real(wp), intent(in) :: rho

        self%new = 1
        self%old = 2

        self%f(1:self%n,0,self%old) = W(0)*rho
        self%f(1:self%n,1,self%old) = W(1)*rho
        self%f(1:self%n,2,self%old) = W(2)*rho
    end subroutine setEquilibriumValues_scalar


    subroutine setEquilibriumValues_array(self,rho)
        class(Lattice), intent(inout) :: self
        real(wp), intent(in) :: rho(getSize(self))

        self%new = 1
        self%old = 2

        self%f(1:self%n,0,self%old) = W(0)*rho(:)
        self%f(1:self%n,1,self%old) = W(1)*rho(:)
        self%f(1:self%n,2,self%old) = W(2)*rho(:)
    end subroutine setEquilibriumValues_array


    subroutine setLeftBC_dirichlet(self,p)
        class(Lattice), intent(inout) :: self
        real(wp), intent(in) :: p
        real(wp) :: fl(0:2)

        fl = self%f(1,:,self%new)
        self%f(1,1,self%new) = p - sum(fl([0,2]))

    end subroutine setLeftBC_dirichlet


    subroutine setRightBC_dirichlet(self,p)
        class(Lattice), intent(inout) :: self
        real(wp), intent(in) :: p

        real(wp) :: fr(0:2)

        fr = self%f(self%n,:,self%new)
        self%f(self%n,2,self%new) = p - sum(fr([0,1]))
    end subroutine setRightBC_dirichlet


    subroutine setRightBC_robyn(self,h)
        class(Lattice), intent(inout) :: self
        procedure(bc), pointer, intent(in) :: h

        real(wp) :: fr(0:2)
        real(wp) :: rho

        ! old value as initial guess to root finder
        rho = sum(self%f(self%n,:,self%old))

        ! find the root
        rho = newton(h,rho)

        ! set the missing boundary value with the required value of rho
        fr = self%f(self%n,:,self%new)
        self%f(self%n,2,self%new) = rho - sum(fr([0,1]))

    end subroutine setRightBC_robyn


    function rgrad(l,c,r)
        ! Calculates the second order backward difference approximation for
        ! a first derivative on a uniform grid.
        !
        ! Parameters
        ! ----------
        !   l : scalar
        !       The value at pos = -2.
        !   c : scalar
        !       The value at pos = -1.
        !   r : scalar
        !       The value at pos = 0, where the approximation of
        !       the derivative is valid.
        !
        ! Returns
        ! -------
        !   rgrad : scalar
        !       Value of the derivative at pos = 0.

        real(wp), intent(in) :: l, c, r
        real(wp) :: rgrad

        rgrad = 0.5_wp*l - 2.0_wp*c + 1.5_wp*r
    end function rgrad


    function lgrad(l,c,r)
        ! Calculates the second order forward difference approximation for
        ! a first derivative on a uniform grid.
        !
        ! Parameters
        ! ----------
        !   l : scalar
        !       The value at pos = 0, where the approximation of
        !       the derivative is valid.
        !   c : scalar
        !       The value at pos = 1.
        !   r : scalar
        !       The value at pos = 2. 
        ! Returns
        ! -------
        !   lgrad : scalar
        !       Value of the derivative at pos = 0.

        real(wp), intent(in) :: l, c, r
        real(wp) :: lgrad
        lgrad = -1.5_wp*l + 2.0_wp*c - 0.5_wp*r
    end function lgrad


    function getRightGrad(self,rho) result(grad)
        class(Lattice), intent(in) :: self
        real(wp), intent(in), optional :: rho
        real(wp) :: grad

        real(wp) :: rholl, rhol, rho_

        rholl = self%getDensityAt(self%n-2)
        rhol = self%getDensityAt(self%n-1)
        if (present(rho)) then
            rho_ = rho
        else
            rho_ = self%getDensityAt(self%n)
        end if
        
        grad = rgrad(rholl,rhol,rho_)
    end function


    function getLeftGrad(self,rho) result(grad)
        class(Lattice), intent(in) :: self
        real(wp), intent(in), optional :: rho
        real(wp) :: grad
        
        real(wp) :: rho_, rhor, rhorr

        rhorr = self%getDensityAt(3)
        rhor = self%getDensityAt(2)
        if (present(rho)) then
            rho_ = rho
        else
            rho_ = self%getDensityAt(1)
        end if
        
        grad = lgrad(rho_,rhor,rhorr)
    end function
end module lbm_mod