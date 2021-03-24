module LoadNetwork
    !use sqliteff
    use sqlite
    use iso_varying_string, only: VARYING_STRING, var_str
    use KuraNetwork  
    use iso_fortran_env, only: real64
    implicit none
    real(real64), parameter :: dx = 2.0_real64   
    real(real64), parameter :: pi = 3.141592653589793238462643383279_real64   
    type system_table
        integer :: id
        integer :: nnodes
        integer :: nedges
        real(kind=real64) :: coupling_constant
        type(SQLITE_DATABASE) :: database 
    end type
    type(system_table) :: ln_table
    type(SQLITE_DATABASE) :: db
    type(SQLITE_COLUMN),dimension(:),pointer :: statecol
    real(real64), dimension(:,:), allocatable :: outstates
    integer :: idx_out 
    real(real64), dimension(:), allocatable :: outtime 
    integer :: err 
    integer :: experimentid 
    integer, parameter :: nsteps = 11 
contains
     subroutine check_db_err(expected)
        implicit none
        integer, intent(in) :: expected
        err = sqlite3_error(db)
        if (err .ne. expected) then
            print*, "Error opening database!\n"
            print*, "expected: ", expected, "\n"
            print*, "got: ", err, "\n"
            stop -1
        end if
     end subroutine check_db_err
     subroutine get_system_info(database,systemid) 
        implicit none
        
        character(len=*), intent(in)     :: database
        integer, intent(in)              :: systemid

        integer :: status
        logical                                    :: finished
        type(VARYING_STRING) :: remaining
        character(512) :: command
        type(SQLITE_STATEMENT) :: stmt
        type(SQLITE_COLUMN),dimension(:),pointer          :: column
        !type(SqliteDatabase_t) :: connection
        real(kind=real64) :: coupling_constant
        integer :: nedges, nnodes
        call sqlite3_open(database, db)
        call check_db_err(0)
        write( command, "(a,I3.1,a)") &
"select nodes, edges, coupling_constant from systems where id=",systemid,";"
        
        !call sqlite3_do(db,command)
        allocate(column(3))
        call sqlite3_prepare(db, command, stmt, column) 
        call sqlite3_column_query( column(1), 'nodes', SQLITE_INT)
        call sqlite3_column_query( column(2), 'edges', SQLITE_INT)
        call sqlite3_column_query( column(3), 'coupling_constant', SQLITE_DOUBLE)
        call sqlite3_next_row(stmt, column, finished)  
        call sqlite3_get_column( column(1), ln_table%nnodes )
        call sqlite3_get_column( column(2), ln_table%nedges )
        call sqlite3_get_column( column(3), ln_table%coupling_constant )
        ln_table%id = systemid
        ln_table%nedges = ln_table%nedges/2 ! we only need the half edges
        deallocate(column) 
        call sqlite3_finalize( stmt ) 
    ! Storage for output at dt points
    ! In the Python and Julia versions, we don't include the time
    ! to write to the database in the timing; so I should do the 
    ! same here.
    !     -           -             !
        allocate(outstates(ln_table%nnodes,nsteps))
        allocate(outtime(nsteps))
        allocate(statecol(5))
    end subroutine get_system_info  

    type(KuramotoNetwork) function load_network()
        implicit none
        !type(KuramotoNetwork(:,:)), allocatable, intent(inout) :: network
        !type(KuramotoNetwork), intent(inout) :: network
        !type(system_table),         intent(in) :: intable
        type(VARYING_STRING)        :: remaining
        character(512)              :: command
        type(SQLITE_STATEMENT)     :: stmt
        type(SQLITE_COLUMN),dimension(:),pointer   :: column
        integer :: i
        logical :: finished

        !allocate(KuramotoNetwork(ln_table%nnodes,ln_table%nedges) :: network)
        !allocate(KuramotoNetwork :: network)
        load_network = KuramotoNetwork(ln_table%nnodes,ln_table%nedges,&
            ln_table%coupling_constant)
        !load_network%nedges = ln_table%nedges 
        !load_network%couple_constant = ln_table%coupling_constant 
        !allocate(load_network%p(ln_table%nnodes))
        !allocate(load_network%Bt_csr(ln_table%nedges*2))
        !allocate(load_network%swap(ln_table%nedges))
        !allocate(knswap(ln_table%nedges)) 
        !associate( connection => intable%database ) 
        !connection = ln_table%database
        !status = sqliteff_open(intable%database, connection)
        !if (status .ne. SQLITE_OK) then
        !    print*, "error loading database"
        !    stop -1
        !end if
        
        ! load the Bt_csr matrix
        write( command, "(a,a,a,I3.1,a,a)") &
            " select distinct x,y from ", &
            " (select source as x, dest as y from edges ", &
            " where system=", ln_table%id, &
            " and source < dest) ", &
            " order by x;"
        allocate(column(2))
        call sqlite3_column_query( column(1), 'x', SQLITE_INT)
        call sqlite3_column_query( column(2), 'y', SQLITE_INT)
        call sqlite3_prepare(db,command,stmt,column)
        call sqlite3_next_row(stmt,column,finished)
        i = 1 
        do while (.not. finished)
            call sqlite3_get_column(column(1), load_network%Bt_csr(i))
            call sqlite3_get_column(column(2), load_network%Bt_csr(i+1))
            i = i+2
            call sqlite3_next_row(stmt,column,finished)
        end do
        !call sqlite3_reset( stmt )
        deallocate(column)
        ! Load node parameters 
        write( command, "(a,a,a,I3.1,a)") &
            " select omega from ", &
            " nodes", &
            " where system=", ln_table%id, &
            " order by id;"
        allocate(column(1))
        call sqlite3_column_query( column(1), 'omega', SQLITE_DOUBLE)
        call sqlite3_prepare(db,command,stmt,column)
        call sqlite3_next_row(stmt,column,finished)
        i = 1
        do while (.not. finished) 
            call sqlite3_get_column(column(1), load_network%p(i))
            i = i+1
            call sqlite3_next_row(stmt,column,finished)
        end do
        deallocate(column)
        call sqlite3_finalize( stmt )
        !call write_experiment_header()
        !err = sqlite3_error(db)
        !print*, err 
    end function load_network 

    subroutine load_initial_conditions(x)
        implicit none
        integer :: status
        !type(system_table),         intent(in)         :: intable
        real(kind=real64),dimension(ln_table%nnodes),intent(inout)   :: x 
        character(512)              :: command
        type(SQLITE_STATEMENT)     :: stmt
        type(SQLITE_COLUMN),dimension(:),pointer   :: column
        integer :: i
        logical :: finished
        
        write( command, "(a,a,a,I3.1,a)") &
            " select initialCondition from ", &
            " nodes", &
            " where system=", ln_table%id, &
            " order by id;"
        allocate(column(1))
        call sqlite3_column_query( column(1), 'initialCondition', SQLITE_DOUBLE)
        call sqlite3_prepare(db,command,stmt,column)
        call sqlite3_next_row(stmt,column,finished)
        i = 1
        do while (.not. finished) 
            call sqlite3_get_column(column(1), x(i))
            outstates(i,1) = x(i)
            call sqlite3_next_row(stmt,column,finished)
            i = i+1
        end do
        outtime(1) = 0
        call sqlite3_finalize( stmt )
    end subroutine load_initial_conditions 
    !subroutine solout(me,nr,xold,x,y,irtrn,xout)
    subroutine solout(nr,xout,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
    use hairerDopri5 
    implicit none
    integer, intent(in)             :: nr
    integer, intent(in)             :: n, nd
    real(real64), intent(inout)                  :: xout
    real(real64), intent(in)                  :: x ! independent variable
    real(real64), dimension(n), intent(in)    :: y ! state vector
    real(real64), dimension(:), intent(inout)   :: con ! dense output vector
    integer                         :: lrc ! length of dense out
    real(real64), dimension(:), intent(in)    :: rpar ! real params
    integer, dimension(:), intent(in)     :: ipar  ! integer params
    integer, intent(inout)     :: icomp  ! integer params
    integer                    :: irtrn ! If < 0, dopri5 returns to 

    integer :: j

    if ( nr==1 ) then
        idx_out = 1
        do j=1,ln_table%nnodes
            outstates(:,idx_out) = y
        end do 
        outtime(idx_out) = 0
        xout = dx
        idx_out = 2
    else
        do
            if ( x<xout ) exit
            print*, xout
            do j=1,ln_table%nnodes
                !outstates(j,idx_out) = modulo(me%contd8(j,xout),2*pi)
                outstates(j,idx_out) = modulo(&
                    CONTD5(j,xout,con,(/ icomp /),nd),&
                    2*pi)
            end do 
            outtime(idx_out) = xout
            idx_out = idx_out + 1
            xout = xout + dx
        end do
    end if

    end subroutine solout
    subroutine get_errorvref(refeid,evref,y)
        implicit none
        integer, intent(out) :: refeid
        real(real64), dimension(:), optional, intent(in) :: y
        real(real64), intent(out):: evref
        real(real64), allocatable, dimension(:) :: errorstate
        type(VARYING_STRING)        :: remaining
        character(1024)              :: command
        type(SQLITE_STATEMENT)     :: stmt
        type(SQLITE_COLUMN),dimension(:),pointer   :: column
        integer :: i
        logical :: finished

        if (present(y)) then
            outstates(:,nsteps) = modulo(y,2.0*pi) 
        end if
        
        ! First get the reference trajectory from the table 
        write( command, "(a, I8.1, a)") &
            "select experimentid from experiments where system=",&
            ln_table%id,&
            " and tend=20 and atol=1e-12 and rtol=1e-10 and solver is 'radau';"
        allocate(column(1)) 
        call sqlite3_column_query( column(1), 'experimentid', SQLITE_INT)
        call check_db_err(0)
        call sqlite3_prepare(db, command, stmt, column )
        call check_db_err(0)
        call sqlite3_next_row(stmt, column ,finished)
        call sqlite3_get_column(column(1), refeid)
        ! Now compare the values to create an average error
        allocate(errorstate(ln_table%nnodes))
        deallocate(column) 
        call sqlite3_reset( stmt )
        allocate(column(1)) 
        call sqlite3_column_query( column(1), 'state', SQLITE_DOUBLE)
            write( command, "(a, I8.1, a, I8.1, a)") &
                "select state from states where experimentid=",&
                refeid,&
                " and time=20 order by idx;"
            call sqlite3_prepare(db, command, stmt, column )
        do i=1,ln_table%nnodes
            call sqlite3_next_row(stmt, column ,finished)
            call sqlite3_get_column(column(1), errorstate(i))
        end do
        evref = sum(abs(outstates(:,nsteps) - errorstate)) / ln_table%nnodes  
        deallocate(errorstate) 
        deallocate(column) 
        !! Load the results  
    end subroutine get_errorvref
    subroutine write_experiment_header(atol, rtol,tstart,tend,tload, tintegration, ttotal,refeid,evref)
        !type(system_table),         intent(in) :: intable
        !integer,         intent(in) :: systemid
        real(real64), intent(in) :: atol, rtol, tstart, tend, tload, tintegration, ttotal, evref
        integer, intent(in) :: refeid
        type(VARYING_STRING)        :: remaining
        character(1024)              :: command
        type(SQLITE_STATEMENT)     :: stmt
        type(SQLITE_COLUMN),dimension(:),pointer   :: column
        integer :: i
        logical :: finished
       character(len=50) :: githash_
       character(len=2048) :: gitchanges_
       character(len=50) :: host_
       integer :: system_
       character(512) :: clcommand
       character(1024) :: exe_directory 
       character(80) :: clcommand_
       character(80) :: runtime = "Fortran" 
       character(80) :: solver = "dop754" 
       character(40) :: date_string
       character(40) :: time_string
       character(40) :: stamp
       !type(SqliteDatabase_t),intent(inout)      :: database
       
       call execute_command_line('git rev-parse HEAD > /tmp/githashtmp')
       open(unit=71,file="/tmp/githashtmp")
       read(71,*) githash_
       close(71)
       
       call execute_command_line('git status -s -uno > /tmp/gitchanges')
       open(unit=71,file="/tmp/gitchanges")
       read(71,*) gitchanges_
       close(71)
       
       call execute_command_line('hostname > /tmp/ln_hostname')
       open(unit=71,file="/tmp/ln_hostname")
       read(71,*) host_
       close(71)
        
       system_ = ln_table%id 
       call get_command(clcommand)
       clcommand = trim(clcommand)
       i = LEN_TRIM(clcommand)
       if(i > 79) then
            clcommand_ = clcommand((i-79):i)
        else
            clcommand_ = trim(clcommand)
        end if
       call date_and_time( date = date_string, time = time_string )
       write( stamp, '(11a)' ) &
         date_string(1:4), '-', date_string(5:6), '-', date_string(7:8), ' ',&
         time_string(1:2), ':', time_string(3:4), ':', time_string(5:6)
       
       write( command, "(17a, I8.1, a, E10.4, a, E10.4, a, F20.15, a,&
          &F20.15, a, F20.15, a, F20.15, a, F20.15, a, I8.1, a, E25.15, a)") &
        "BEGIN TRANSACTION; INSERT INTO experiments (runtime, solver, datetime, githash, gitchanges, commandline, host, system,",&
        " atol, rtol, tstart, tend, tload, tintegration, ttotal, refeid, err_v_ref)",&
        " VALUES( '",&
        trim(runtime), "','",trim(solver), "','",&
        trim(stamp), "','", trim(githash_), "','", trim(gitchanges_), "','",&
        trim(clcommand_), "','",trim(host_), "',", system_, &
        ",", atol, ",", rtol, ",", tstart, ",", tend, ",", tload, ",",&
        tintegration, ",", ttotal, ",", refeid, ",", evref,&
        " ); END TRANSACTION;"
        !call sqlite3_begin( db )
        !err = sqlite3_error(db)
        !print*, err 
        call sqlite3_do(db,trim(command))
        call sqlite3_commit(db)
        !call sqlite3_commit( db )
        !err = sqlite3_error(db)
        ! Experimentid
        allocate(column(1))
        call sqlite3_prepare(db,'SELECT last_insert_rowid()', stmt, column )
        call check_db_err(0)
        call sqlite3_column_query( column(1), 'experimentid', SQLITE_INT)
        call sqlite3_next_row(stmt, column ,finished)
        call sqlite3_get_column( column(1), experimentid )

        deallocate(column)
    end subroutine write_experiment_header
     subroutine write_experiment_states()
        implicit none
        integer :: time
        integer :: idx
        type(VARYING_STRING)        :: remaining
        character(1024)              :: command
        type(SQLITE_STATEMENT)     :: stmt
        integer :: i,j
        logical :: finished
        type(SQLITE_COLUMN),dimension(:),pointer          :: column
        allocate(column(5))
        call sqlite3_begin( db )
        do i=1,nsteps
            do j=1,ln_table%nnodes
                call sqlite3_set_column(column(1),experimentid) 
                call sqlite3_set_column(column(2),outtime(i)) 
                call sqlite3_set_column(column(3),j) 
                call sqlite3_set_column(column(4),j) 
                call sqlite3_set_column(column(5),outstates(j,i))
                call sqlite3_insert(db, 'states', column)
            end do
        end do 
        call sqlite3_commit( db )
        deallocate(column) 
    end subroutine write_experiment_states
    subroutine write_final_state(y,t)
        implicit none
        real(real64), dimension(:), intent(in) :: y
        real(real64), intent(in) :: t
        type(VARYING_STRING)        :: remaining
        character(1024)              :: command
        type(SQLITE_STATEMENT)     :: stmt
        integer :: i,j
        logical :: finished
        type(SQLITE_COLUMN),dimension(:),pointer          :: column
        allocate(column(5))
        call sqlite3_begin( db )
        do j=1,ln_table%nnodes
            call sqlite3_set_column(column(1),experimentid) 
            call sqlite3_set_column(column(2),t) 
            call sqlite3_set_column(column(3),j) 
            call sqlite3_set_column(column(4),j) 
            call sqlite3_set_column(column(5),y(j))
            call sqlite3_insert(db, 'states', column)
        end do
        call sqlite3_commit( db )
        deallocate(column) 
    end subroutine write_final_state
    subroutine close_database_conn()
        call sqlite3_commit(db)
        call sqlite3_close(db)
    end subroutine close_database_conn
end module LoadNetwork

