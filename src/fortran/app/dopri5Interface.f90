module dopri5Interface
!include 'dopri5.f'
use hairerDopri5
use iso_fortran_env, only: wp=>real64
implicit none
real(wp), dimension(1) :: DEFAULT_RTOL = (/ 1.0E-9_wp /)
real(wp), dimension(1) :: DEFAULT_ATOL = (/ 1.0E-9_wp /)
integer, dimension(1) :: DEFAULT_IPAR = (/ 0 /)
real(wp), dimension(1) :: DEFAULT_RPAR = (/ 0.0_wp /)
real(wp)               :: xout


type, public :: dopri5_t
    integer :: N ! Dimension of system
    real(wp)                :: h!
    integer :: NRDENS ! Number of components for which dense output is
                      ! required. Default is 0; for 0 < nrdens < n
                      ! components must be specified in iwork(21) -
                      ! iwork(nrdens+20). For nrdens = N this is done
                      ! automatically
    integer :: IOUT ! Output routine (0: no solout, 1 solout output, 2 dense out
    real(wp), dimension(1) :: RTOL
    real(wp), dimension(1) :: ATOL
    integer  :: ITOL ! itol = 0, both rtol and atol are
                      ! scalars. itol = 1, rtol and atol are vectors

    real(wp), dimension(:), allocatable :: WORK ! must be at least 
                      ! 8 * N + 5 * NRDENS + 21
    integer                :: LWORK ! Declared length of "WORK" 
    integer, dimension(:), allocatable :: IWORK ! must be at least
                      ! NRDENS + 21 
    integer                :: LIWORK ! Declared length of "IWORK" 
    
    real(wp), dimension(:), allocatable :: RPAR  
    integer, dimension(:), allocatable  :: IPAR  
    
    procedure(deriv_func), pointer, nopass  :: fcn ! rhs function
    procedure(solout_func), pointer, nopass :: solout ! rhs function
   
    contains
        procedure, public :: integrate  => dopri5_int

end type dopri5_t
interface dopri5_t
    module procedure initialize
end interface dopri5_t 

 abstract interface 
     subroutine deriv_func(N,x,y,f,rpar,ipar)
        use iso_fortran_env, only: wp => real64
        implicit none
        real(wp), intent(in)                :: x ! independent variable
        integer, intent(in)                 :: N ! Numberr or states 
        real(wp), dimension(N), intent(in)  :: y ! state vector
        real(wp), dimension(N), intent(out) :: f ! deriv vector dy/dx
        real(wp), dimension(:), intent(in),optional  :: rpar ! real params
        integer, dimension(:), intent(in),optional   :: ipar  ! integer params
    end subroutine
    
    subroutine solout_func(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
        use iso_fortran_env, only: wp => real64
        implicit none
        integer, intent(in)             :: nr
        integer, intent(in)             :: n, nd
        real(wp), intent(inout)                  :: xold
        real(wp), intent(in)                  :: x ! independent variable
        real(wp), dimension(n), intent(in)    :: y ! state vector
        real(wp), dimension(:), intent(inout)   :: con ! dense output vector
        integer                         :: lrc ! length of dense out
        real(wp), dimension(:), intent(in)    :: rpar ! real params
        integer, dimension(:), intent(in)     :: ipar  ! integer params
        integer, intent(inout)     :: icomp  ! integer params
        integer                        :: irtrn ! If < 0, dopri5 returns to 
                    ! vcalling program. If irtrn=2, the numerical
                 ! solution is changed in iout
    end subroutine
 end interface
contains
 type(dopri5_t) function initialize(n,fcn,rtol,atol,h,ndense,solout,iout) 
    implicit none
    integer, intent(in) :: n
    procedure(deriv_func) :: fcn
    procedure(solout_func), optional:: solout
    real(wp), optional, intent(in) :: h,ndense
    real(wp), optional, dimension(:), intent(in) :: rtol 
    integer, optional, dimension(:), intent(in) :: atol 
    integer, optional, intent(in) :: iout 
    
    initialize%n = n
    initialize%fcn => fcn
    if (present(h)) then
        initialize%h = h
    else
        initialize%h = 0.0
    end if
    if (present(rtol)) then
        if (size(rtol,dim=1)>1) then
            initialize%ITOL = 1 ! THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
                          ! RTOL(I)*ABS(Y(I))+ATOL(I).
        else
            initialize%ITOL = 0 ! Code keeps Y(I) BELOW RTOL*ABS(Y(I))+ATOL
        initialize%RTOL = rtol
        end if
    else
        initialize%ITOL = 0
        initialize%RTOL = DEFAULT_RTOL!  (/ 0 /)
    end if 
    if (present(atol)) then
        if ((size(atol,dim=1)>1)) then
                initialize%ITOL = 1 ! THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
                              ! RTOL(I)*ABS(Y(I))+ATOL(I).
        else
            initialize%ITOL = 0 ! Code keeps Y(I) BELOW RTOL*ABS(Y(I))+ATOL
        end if
        initialize%ATOL = atol
    else
        initialize%ITOL = 0
        initialize%ATOL = 1.0*initialize%RTOL ! Absolute tolerance default to real 
    end if 
    if (present(iout)) then
        initialize%IOUT = iout
        if (iout>0) then
            if (present(solout)) then
                initialize%SOLOUT => solout
            else
                initialize%SOLOUT => null()
            end if
        end if
    end if
    initialize%IPAR = DEFAULT_IPAR
    initialize%RPAR = DEFAULT_RPAR
    if (present(ndense)) then
        initialize%NRDENS = ndense
    else
        initialize%NRDENS = n
    end if
    
    associate(N => initialize%N, NRDENS => initialize%NRDENS)
        allocate(initialize%WORK(8*N+5*NRDENS+21)) !  
        initialize%LWORK = ubound(initialize%WORK,dim=1)
        initialize%WORK = 0.0_wp
        
        allocate(initialize%IWORK(NRDENS+21)) ! Hairer: "LIWORK" MUST BE AT LEAST NRDENS+21 
        initialize%LIWORK = ubound(initialize%IWORK,dim=1)
        initialize%IWORK = 0
    end associate 
 end function initialize 
 subroutine dopri5_int(this,x,y,xend,rtol,atol,iout,idid)
    implicit none
    class(dopri5_t), intent(inout) :: this
    real(wp), intent(inout)                :: x
    real(wp), dimension(:),intent(inout)   :: y
    real(wp), intent(in)                   :: xend
    real(wp), intent(inout),optional :: rtol
    real(wp), intent(inout),optional :: atol
    integer, intent(in),optional               :: iout !switch for solout 
    integer, intent(inout)                       :: idid
     
    !if (.not.present(atol)) then
    !   atol = this%atol(1)
    !end if 
    !if (.not.present(rtol)) then
    !   rtol = this%rtol(1)
    !end if 
    call DOPRI5(this%N,this%FCN,x,y,xend, &
                       (/ rtol /) , (/ atol /), 0, &
                       this%solout,this%IOUT, &
                       this%WORK,this%LWORK,this%IWORK,this%LIWORK,&
                       this%RPAR,this%IPAR,idid)
    end subroutine dopri5_int
    subroutine dopri5_solout(nr,xold,x,y,n,cont,icomp,nd,rpar,ipar,irtrn) 
        use iso_fortran_env, only: wp => real64
        implicit none
        !TODO: Where is the output in this function?
        integer, intent(in)             :: nr, n
        integer, intent(in)             :: nd
        real(wp), intent(inout)                  :: xold
        real(wp), intent(in)                  :: x ! independent variable
        real(wp), dimension(nr), intent(in)   :: y ! state vector
        real(wp), dimension(:), intent(inout)   :: cont ! dense output vector
        integer, intent(inout)                  :: icomp ! length of dense out
        real(wp), dimension(:), intent(in)    :: rpar ! real params
        integer, dimension(:), intent(in)     :: ipar  ! integer params
        integer                               :: irtrn ! If < 0,
        real(wp)                              :: currentValue
        ! radau returns to the calling program.  
        if (NR.eq.1) then
           print*, x,y(1),y(2),y(3)
           xout = 0.25_wp
        else
            do while (x.ge.xout)
                !do i = 1,ubound(y,dim=1)
                    print*, xout, CONTD5(1,xout,cont,(/ icomp /),nd),NR-1
                !end do
                xout = xout+0.25
            end do
        end if
    end subroutine dopri5_solout
end module
