module KuraNetwork
    use iso_fortran_env, only: real64
    implicit none

    type KuramotoNetwork
       integer :: nnodes
       integer :: nedges
       real(kind=real64)            :: couple_constant
       real(kind=real64), dimension(:),allocatable :: p 
       real(kind=real64), dimension(:),allocatable :: swap 
       integer,dimension(:),allocatable            :: Bt_csr !IncidenceMatrix 
    end type 

    interface KuramotoNetwork
        module procedure KNConstructor
    end interface
    ! TODO: Constructor that allocates tmp(nnodes)
contains
    type(KuramotoNetwork) function KNConstructor(nnodes,nedges,couple_constant)
        implicit none
        integer::nnodes
        integer::nedges
        real(real64)::couple_constant
        KNConstructor%nnodes = nnodes
        KNConstructor%nedges = nedges
        KNConstructor%couple_constant = couple_constant
        allocate(KNConstructor%p(nnodes)) 
        allocate(KNConstructor%swap(nedges)) 
        allocate(KNConstructor%Bt_csr(nedges*2))

    end function KNConstructor 

    subroutine f(this,t,x,dx)
        ! dx .= p .-  coupling_constant .* (B * sin.(B_t * x))   
        implicit none
        
        class(KuramotoNetwork), intent(inout) :: this
            
        real(kind=real64), intent(in)      :: t
        real(kind=real64), intent(in), dimension(:)     :: x
        real(kind=real64), intent(out), dimension(:)  :: dx
        integer :: i
        associate( nedges => this%nedges,&
                   nnodes => this%nnodes,&
                   col_index => this%Bt_csr,&
                   omega => this%p,&
                   couple_constant => this%couple_constant,&
                   y => this%swap) 
        y = 0 
        do concurrent (i =1:nedges) 
             y(i) = sin(x(col_index((i-1)*2+2)) - x(col_index((i-1)*2+1)))
        end do
        dx = 0
        do i = 1,nedges ! treat the csr as a csc to transpose it 
            associate(idx_s => col_index((i-1)*2+1),&
                  idx_k => col_index((i-1)*2+2)&
                  ) 
            dx(idx_s) = dx(idx_s) - y(i)
            dx(idx_k) = dx(idx_k) + y(i)
            end associate
        end do
        do concurrent (i = 1:nnodes) ! now dx = B*sin(B_t*x)
            dx(i) = omega(i) - couple_constant * dx(i)
        end do
        ! Now dx .= p .- coupling_constant .* (B*sin.(Bt*x))

    end associate
    end subroutine f 
end module KuraNetwork
