program main

  !======= Inclusions ===========
  use iso_fortran_env, only: wp => real64

  use KuraNetwork
  use LoadNetwork
  use iso_varying_string, only: VARYING_STRING, var_str
  use DopriWrap
  use dopri5Interface 
  !======= Declarations =========
  implicit none
  ! local variables
  real(wp) :: tstart                   ! initial time
  real(wp) :: tend                     ! final time
  real(wp) :: rtol, atol               ! relative and absolute tolerance
  real(wp) :: dtout                    ! output time interval
  real(wp) :: tout                     ! output time
  real(wp) :: tcur                  ! current time
  integer :: nout                     ! number of outputs
  integer :: outstep                  ! output loop counter
  integer :: neq                  ! output loop counter
  type(dopri5_t)    :: prop
  logical :: status_ok
  integer :: idid  
  integer,dimension(:), allocatable  :: icomp 
  character(len=:), allocatable :: filename,buffer
  integer   :: testcase
  integer, dimension(:),allocatable  :: n,c 
  character, dimension(:),allocatable:: nt 
  character, dimension(:),allocatable:: cd
  type(VARYING_STRING) :: database 
  real(wp)                          :: start_time, stop_time
  real(wp)                          :: rate,theta
  real(wp), dimension(:), allocatable :: y, time
  real(wp), dimension(:,:), allocatable :: output
  integer                           :: start_clock, stop_clock,i,j
  type(system_table)                :: stable
  real(wp) :: tload, tintegrate, ttotal
  integer                           :: refeid
  real(wp) :: evref  
    
    allocate(character(512) :: filename)
    CALL get_environment_variable("KURABENCH_DB", filename)
    filename = trim(filename)
    
    if(command_argument_count().gt.0) then
        allocate(character(5) :: buffer)
        call get_command_argument(1,buffer)
        buffer = trim(buffer)
        read(buffer,'(I4.1)') testcase
        deallocate(buffer)
    else
        testcase=1 
    end if
        tstart=0.0_wp
    if(command_argument_count().gt.1) then
        allocate(character(5) :: buffer)
        call get_command_argument(2,buffer)
        buffer = trim(buffer)
        read(buffer,*) tend 
        deallocate(buffer)
    else
        tend=1.0_wp
    end if
    if(command_argument_count().gt.2) then
        allocate(character(10) :: buffer)
        call get_command_argument(3,buffer)
        buffer = trim(buffer)
        read(buffer,*) atol 
        deallocate(buffer)
    else
        atol=1.0d-6
    end if
    if(command_argument_count().gt.3) then
        allocate(character(10) :: buffer)
        call get_command_argument(4,buffer)
        buffer = trim(buffer)
        read(buffer,*) rtol 
        deallocate(buffer)
    else
        rtol=1.0d-6
    end if
  
  tcur   = tstart
  tout   = tstart
  dtout  = 2.0d0
  nout   = ceiling(tend/dtout)
  call cpu_time(start_time) 
  call dopri_loadNet(filename,testcase)
  neq = ln_table%nnodes
  allocate(y(neq))
  allocate(icomp(neq))
  
  icomp = (/(i, i=1,neq, 1)/)
  call load_initial_conditions(y)  
  call cpu_time(stop_time) 
  tload = stop_time - start_time
  prop = initialize(fcn=dopri_RhsFn,n=neq,solout=solout,iout=0)
  
  call cpu_time(start_time) 
     !now, perform the integration:
  call prop%integrate(tcur,y,tend,rtol,atol,idid=idid) 
  
  call cpu_time(stop_time) 
  tintegrate=stop_time - start_time
  !call write_experiment_states()
  call get_errorvref(refeid,evref,y)
  call write_experiment_header(atol, rtol, tstart, tend, &
                              tload, tintegrate, tload+tintegrate,&
                              refeid, evref)
  call write_final_state(y,tcur)
  call close_database_conn() 
end program main
