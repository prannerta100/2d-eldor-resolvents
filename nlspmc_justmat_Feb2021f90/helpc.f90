! NLSPMC VERSION (VERSION 1.0)  2/5/99
!----------------------------------------------------------------------
!                  =========================
!                      subroutine helpc
!                  =========================

!     Main program for a nonlinear least-squares fit using an
!     EPRP-family slow-motional calculation. The options in running
!     this program are too numerous to detail here. Read the manual...
!     or better yet, wait until the movie comes out.  DB

!----------------------------------------------------------------------
    subroutine helpc(line)
    implicit none

    include 'stdio.inc'

    character hlptxt*132,line*80, cat1*40, cat2*40, hlpcat*40
    integer :: ibar,iblk,lth1,lth2,nlines,LINES
    logical :: found1,found2,kywrd1,kywrd2,match1,noncmd
    parameter(LINES=23)

    call gettkn( line, cat1, lth1 )
    call gettkn( line, cat2, lth2 )
    if (lth1 /= 0) call touppr(cat1,lth1)
    if (lth2 /= 0) call touppr(cat2,lth2)
    open (ludisk,file='/afs/msc/home/freed/sanghyuk/bin/nlshlp.txt', &
    status='old',access='sequential',err=10)

    found1=.false.
    found2=.false.
    nlines=0

!  Search through lines in the help text file

    1 read (ludisk,'(a)',end=10,err=11) hlptxt
    ibar=1
    2 if (hlptxt(ibar:ibar) /= '|' .AND. ibar < 132) then
        ibar=ibar+1
        go to 2
    end if

    kywrd1=hlptxt(1:1).eq.'*'
    kywrd2=hlptxt(1:1).eq.'>'
    noncmd=hlptxt(1:1).eq.' '
    hlpcat=hlptxt(2:ibar-1)
    hlptxt=hlptxt(ibar+1:)
    ibar=ibar-2

!     Find last nonblank character in help text line

    iblk=132
    3 if (hlptxt(iblk:iblk) == ' ') then
        iblk=iblk-1
        go to 3
    end if

!     If help text entry represents a major category, check the
!     first keyword specified in the help command (if any) against it

    if (kywrd1) then
        match1=.false.
        if (lth1 == 0) then
            call linchk(nlines)
            write(luttyo,1000) hlpcat(:ibar),hlptxt(:iblk)
        else if (cat1(:lth1) == hlpcat(:lth1)) then
            call linchk(nlines)
            write(luttyo,1000) hlpcat(:ibar),hlptxt(:iblk)
            match1=.true.
            found1=.true.
        end if
    
    !     If help text entry represents a subcategory, check the
    !     second keyword specified in the help command (if any) against it

    else if (kywrd2) then
        if (match1 .AND. lth2 == 0) then
            call linchk(nlines)
            write(luttyo,1004) hlpcat(:ibar),hlptxt(:iblk)
        else if (match1 .AND. cat2(:lth2) == hlpcat(:lth2)) then
            call linchk(nlines)
            write(luttyo,1004) hlpcat(:ibar),hlptxt(:iblk)
            found2=.true.
        end if
    
    else if (noncmd .AND. lth1 /= 0) then
        if (cat1(:lth1) == hlpcat(:lth1)) then
            call linchk(nlines)
            write(luttyo,1000) hlpcat(:ibar),hlptxt(:iblk)
            found1=.true.
        end if
    end if
    go to 1

    10 if ((lth1 /= 0 .AND. .NOT. found1) .OR. &
    (lth2 /= 0 .AND. .NOT. found2)) &
    write (luttyo,1001) cat1(:lth1),cat2(:lth2)

    close(ludisk)
    write (luttyo,*)
    return

    11 write (luttyo,1003) hlpcat(:ibar),hlptxt(:iblk)
    close(ludisk)
    write (luttyo,*)
    return

    1000 format(a,t22,a)
    1001 format('*** No help available for ''',a,' ',a,''' ***')
    1002 format('*** File ''nlshlp.txt'' not available ***')
    1003 format('*** Error reading file ''nlshlp.txt'' ***')
    1004 format(2x,a,t24,a)
    end subroutine helpc


    subroutine linchk( nlines )
    implicit none
    include 'stdio.inc'

    integer :: nlines,MXLINES
    character dummy*1
    parameter(MXLINES=20)

    if (nlines == 0) write (luttyo,1001)
    nlines=nlines+1
    if (nlines > MXLINES) then
        write (luttyo,1000)
        read (luttyo,'(a)') dummy
        nlines=1
    end if
    return
    1000 format('...press <RETURN> to continue...')
    1001 format(/15x,' *** NLSPMC on-line help ***')
    end subroutine linchk
