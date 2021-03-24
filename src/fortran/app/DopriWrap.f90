module dopriWrap 
  !======= Inclusions ===========
  use KuraNetwork
  use LoadNetwork
  use iso_fortran_env, only: real64
  !======= Declarations =========
  implicit none

  !type(KuramotoNetwork(:,:)),allocatable :: net 
  type(KuramotoNetwork)  :: net 
  real(real64),parameter :: dt = 2.0_real64

contains
    subroutine dopri_loadNet(database, systemid) 
    implicit none
    
    character(len=*), intent(in)     :: database
    integer, intent(in)              :: systemid
    integer                          :: state_length
    !real(kind=real64),dimension(:), allocatable    :: y,dy 
    call get_system_info(database,systemid)
    !call load_network(net) 
    net =  load_network() 
    return
end subroutine dopri_loadNet

  ! ----------------------------------------------------------------
  ! RhsFn provides the right hand side function for the
  ! ODE: dy/dt = f(t,y)
  ! Here wrapped for use with the dopri5 solver 
  ! ----------------------------------------------------------------
  subroutine dopri_RhsFn(N,t,y,dy,rpar,ipar)

    implicit none

    integer,        intent(in)               :: N
    real(real64)   ,intent(in)               :: t
    real(real64)   ,dimension(N),intent(in)  :: y
    real(real64)   ,dimension(N),intent(out) :: dy
    real(real64), dimension(:), intent(in),optional  :: rpar ! real params
    integer, dimension(:), intent(in),optional   :: ipar  ! integer params

    ! fill RHS vector fvec
    call f(net,t,y,dy)

    ! return success
    end subroutine dopri_RhsFn
end module dopriWrap 
