!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                  ===================
!                     program NLSPMC
!                  ===================

!     MAIN PROGRAM FOR NON-LINEAR LEAST SQUARES FIT OF 2D-SPECTRA

!     This program reads in parameters and data for a nonlinear least-
!     squares fit using an EPRCGF-family slow-motional 2D spectral
!     calculation.

!     Modified to allow multiple components.  RC 7/98

!     ** Needs graphical interface **

!     Includes:
!        nlsdim.inc
!        nlsnam.inc
!        eprprm.inc
!        expdat.inc
!        parcom.inc
!        parmequ.inc
!        lmcomm.inc
!        stdio.inc
!        iterat.inc
! 	 miscel.inc	new variables controlling sites/spectra choice

!######################################################################

    program nlspmc
!c	include 'CXML_include.F90'

    implicit none

!      include 'nlsdim.inc'
!      include 'nlsnam.inc'
!      include 'eprprm.inc'
!      include 'expdat.inc'
!      include 'parcom.inc'
!      include 'parmequ.inc'
!      include 'lmcomm.inc'
!      include 'stdio.inc'
!      include 'iterat.inc'
    include 'stdio.inc'
    include 'limits.inc'
    include 'names.inc'
    include 'parms.inc'
    include 'parmequ.inc'
    include 'simparm.inc'
    include 'miscel.inc'

    integer :: i,iflg,ioerr,j,lth
    logical :: fexist
    character line*(linelg), token*(wordlg), scratch*(wordlg), &
    fileID*(wordlg), chr*2

    logical :: getlin
    external getlin

    write (luttyo,1000)
    write (luttyo,1010)

!    ----------------------
!    Initialize NLS system
!    ----------------------

    call nlsinit
    write (luttyo,1011)

!----------------------------------------------------------------------
!     Get next command from input stream (skip comments)
!----------------------------------------------------------------------

    25 if ((ihltcmd /= 0) .AND. (nfiles > 0)) then
    
    ! if ^c is struck, close all macro files and start over:
    
        if (luout /= luttyo) write(luttyo,1004)
        write(luout,1004)
        do 27 i=1,nfiles
            close (ludisk+i)
        27 END DO
        nfiles=0
        lucmd=luttyi
    !         call uncatchc(ihltcmd)
    !  reinitalize everything (later, maybe only do some things):
        call nlsinit
    end if
! input a line:
    if (getlin(line)) then
    
    !######################################################################
    ! get the first (next) token:
        call gettkn(line,token,lth)
        call touppr(token,lth)
    ! if done (empty line or comment), get next line:
        if (lth == 0 .OR. token == 'C' .OR. token == '/*') go to 25
    
    ! Now analyze the line received, first word is command, then extra
    ! information may follow or not.
    
    !----------------------------------------------------------------------
    !        ASSIGN command
    !  Assigns a given basis set to a given set of spectrum, site indices
    ! command removed since basis(i,j) selects site and spectrum info.
    !----------------------------------------------------------------------
    !         if (token.eq.'ASSIGN') then
    !            call assgnc(line)
    
    !----------------------------------------------------------------------
    !        AXIAL command
    !----------------------------------------------------------------------
        if (token == 'AXIAL') then
            call convtc(line,3)
        
        !----------------------------------------------------------------------
        !        BASIS command
        !----------------------------------------------------------------------
        else if (token == 'BASIS') then
        
            call basis_subroutine(line,luout)
        
        !----------------------------------------------------------------------
        !        CARTESIAN command
        !----------------------------------------------------------------------
        else if (token == 'CARTESIAN' .OR. token == 'CART') then
            call convtc(line,1)
        
        !----------------------------------------------------------------------
        !        COMPONENTS command (alias COMP)
        !----------------------------------------------------------------------
        else if (token == 'COMPONENTS' .OR. token == 'COMP') then
            call comps(line)
        
        !----------------------------------------------------------------------
        !        DATA command
        !----------------------------------------------------------------------
        else if (token == 'DATA') then
            call datac(line)
        
        !----------------------------------------------------------------------
        !        DEBUG command: open debug file identification
        !----------------------------------------------------------------------
        else if (token == 'DEBUG') then
            idebug=1
            open(ludeb,file=dbname(:lthfnm),status='unknown', &
            access='sequential',form='formatted',iostat=ioerr)
        
            ievec=0
            call gettkn(line,token,lth)
            call touppr(token,lth)
            if (token == 'VECTOR') ievec=1
        
        !----------------------------------------------------------------------
        !        ECHO command
        !----------------------------------------------------------------------
        else if (token == 'ECHO') then
            call gettkn(line,token,lth)
            if (lth > 0) then
                scratch=token
            else
                scratch=' '
            end if
            call touppr(scratch,3)
            if (scratch == 'ON') then
                luecho=luttyo
            else if (scratch == 'OFF') then
                luecho=0
            else
                call ungett(token,lth,line)
                write(luttyo,1060) line
            end if
        
        !----------------------------------------------------------------------
        !        FIT command
        !----------------------------------------------------------------------
        else if (token == 'FIT') then
            call fitc(line)
        
        !----------------------------------------------------------------------
        !        FIX command (alias REMOVE)
        !----------------------------------------------------------------------
        else if (token == 'FIX' .OR. token == 'REMOVE') then
            call fixc(line)
        
        !----------------------------------------------------------------------
        !        HELP command
        !----------------------------------------------------------------------
        else if (token == 'HELP') then
            call helpc(line)
        
        !----------------------------------------------------------------------
        !        LET command
        !----------------------------------------------------------------------
        else if (token == 'LET') then
            call letcmc(line)
        
        !----------------------------------------------------------------------
        !        LOG command: Set log file identification
        !----------------------------------------------------------------------
        else if (token == 'LOG') then
            call gettkn(line,fileID,lth)
            if (lth == 0) then
                write (luttyo,1020)
            
            else if (fileID == 'END' .OR. fileID == 'end') then
                if (luout == luttyo) then
                    write (luttyo,1021)
                else
                    close(lulog)
                    luout=luttyo
                end if
            
            else
                call setfil( fileID )
            
                open(lulog,file=lgname(:lthfnm),status='unknown', &
                access='sequential',form='formatted',iostat=ioerr)
            
                write(lulog,*)'log file of program nlspc.rc'
            
            !  remove appending to existing log file (RC 7/25/96)
            
            !               inquire(file=lgname(:lthfnm),exist=fexist)
            !               if (fexist) then
            !                                   * append the log file
            ! 30               read (lulog,'(a2)',end=32) chr
            !                  go to 30
            !               end if
                32 luout=lulog
            end if
        
        !----------------------------------------------------------------------
        !        PARMS command
        !----------------------------------------------------------------------
        else if (token == 'PARMS') then
            call parc(line)
        
        !----------------------------------------------------------------------
        !        QUIT command (alias EXIT)
        !----------------------------------------------------------------------
        else if (token == 'QUIT' .OR. token == 'EXIT') then
            goto 9999
        
        !----------------------------------------------------------------------
        !        READ command (alias CALL)
        !        open a new input file and set I/O units appropriately
        !----------------------------------------------------------------------
        
        else if (token == 'READ' .OR. token == 'CALL') then
        
        !       --- get filename
        
            call gettkn(line,fileID,lth)
        
            if (lth /= 0) then
            
            !       --- open input file if possible
            
                if (nfiles >= MXSPEC) then
                    write (luttyo,1050) fileID(:lth),MXSPEC
                else
                    nfiles=nfiles+1
                    lucmd=ludisk+nfiles
                    inquire(file=fileID(:lth),exist=fexist)
                    if (fexist) open(lucmd,file=fileID(:lth), &
                    status='old',access='sequential', &
                    form='formatted',iostat=ioerr)
                    if (( .NOT. fexist) .OR. ioerr /= 0) then
                    !                                               *** open error
                        write (luttyo,1030) fileID(:lth)
                        nfiles=nfiles-1
                        if (nfiles == 0) then
                            lucmd=luttyi
                        else
                            lucmd=lucmd-1
                        end if
                    
                    else
                        files(nfiles)=fileID
                    
                    !                     call catchc( ihltcmd )
                    
                    end if
                end if
            !                              *** File identification not specified
            else
                write (luttyo,1020)
            end if
        
        !----------------------------------------------------------------------
        !        RESET command
        !----------------------------------------------------------------------
        else if (token == 'RESET') then
            call nlsinit
        
        !----------------------------------------------------------------------
        !        SEARCH command
        !----------------------------------------------------------------------
        else if (token == 'SEARCH') then
            call srchc(line)
        
        !----------------------------------------------------------------------
        !        SERIES command
        !----------------------------------------------------------------------
        else if (token == 'SERIES') then
            call series(line)
        
        !----------------------------------------------------------------------
        !        SPECTRA command (alias SPEC)
        !        specify multiple experimental data sets
        !----------------------------------------------------------------------
        
        else if (token == 'SPECTRA' .OR. token == 'SPEC') then
            call spectra(line)
        
        !----------------------------------------------------------------------
        !        SPHERICAL command
        !----------------------------------------------------------------------
        else if (token == 'SPHERICAL' .OR. token == 'SPHER') then
            call convtc(line,2)
        
        !----------------------------------------------------------------------
        !        STATUS command
        !----------------------------------------------------------------------
        else if (token == 'STATUS') then
            call statc(line,luttyo)
        
        
        !----------------------------------------------------------------------
        !        TEMPERATURE command (alias TEMP)
        !----------------------------------------------------------------------
        !	 else if (token.eq.'TEMPERATURE' .or. token.eq.'TEMP')
        !     #		then
        !	    tempvar = .true.
        
        !----------------------------------------------------------------------
        !        TILT command	Not needed yet.
        !----------------------------------------------------------------------
        !	 else if (token.eq.'TILT') then
        !	   psivar = .true.
        
        !----------------------------------------------------------------------
        !        UNIFORM command
        !----------------------------------------------------------------------
        else if (token == 'UNIFORM') then
            uniflg=.true.
        
        !----------------------------------------------------------------------
        !        VARY command
        !----------------------------------------------------------------------
        else if (token == 'VARY') then
            call varyc(line)
        
        !----------------------------------------------------------------------
        !        Unknown command
        !----------------------------------------------------------------------
        else
            write(luttyo,1040) token(:lth)
        end if
    
    !----------------------------------------------------------------------
    !    No more lines (getlin returned .false. or it received ^d).
    !    Close current input unit; if there are no open macro files,
    !    stop program.  Interactive usage will not stop as getlin reads
    !    from keyboard by default, in bach mode, the end of the last macro
    !    file will terminate here.
    !----------------------------------------------------------------------
    
    else
        if (nfiles == 0) then
            write(luttyo,1000)
            stop 'end of program NLSPMC'
        else
            close(lucmd)
            nfiles=nfiles-1
            if (nfiles == 0) then
                lucmd=luttyi
            !               call uncatchc(ihltcmd)
            else
                lucmd=lucmd-1
            end if
        end if
    end if
    go to 25

!----------------------------------------------------------------------
!     Exit program
!----------------------------------------------------------------------

    9999 continue
    stop

!## format statements ###############################################

    1000 format(//,2x,70('#'),//)
    1004 format(/20x,'*** MACRO execution halted by user ***')
    1010 format(25x,'PROGRAM : NLSPMC',/,25x,'---------------',//)

    1011 format(/,'Program defaults to single spectra, single component.', &
    /,'You must first specify other choice with commands:',/, &
    'spectra # (number of data sets), components # (# in ', &
    'each spectra).',//,'Further, temperature variation', &
    ' is not allowed ',/,'unless first specified with temperature', &
    ' command.',/ )
    1020 format('*** File identification must be specified ***'/)
    1021 format('*** Log file is not open ***')
    1030 format('*** Error opening or reading file ''',a,''' ***'/)
    1040 format('*** Unknown command : ''',a,''' ***')
    1050 format('*** Cannot open ''',a,''': more than',i2, &
    ' read files ***')
    1060 format(a)
    END PROGRAM
