! NLSPMC Version 1.0 2/5/99
!----------------------------------------------------------------------
!             I/O UTILITY ROUTINES FOR NLS COMMAND INTERPRETER

!  Contains the subroutines:
!       getlin  -- issue a prompt and retreive a line from command stream
!       gettkn  -- return a token from the given line
!       touppr  -- converts a string to uppercase
!       ftoken  -- returns a real number
!       itoken  -- returns an integer number
!       indtkn  -- returns an index (parenthesized number)
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!                    =========================
!                       function GETLIN
!                    =========================
!----------------------------------------------------------------------
    function getlin( line )
    implicit none
    logical :: getlin

    include 'stdio.inc'
    include 'limits.inc'

    integer :: ioerr
    character*(LINELG) line

!      call wpoll

    if(lucmd == luttyi) call lprmpt
    read (lucmd,1001,iostat=ioerr) line

!   -- Echo line if required

    if (ioerr == 0 .AND. lucmd /= luttyi .AND. luecho /= 0) &
    write(luecho,1001) line
    if (luecho /= 0 .AND. luout == lulog) write(lulog,1001) line
    getlin=ioerr.eq.0
    return

    1001 format(a)
    end function getlin

!------------------------------------------------------------------------
!                         ===================
!                          subroutine GETTKN
!                         ===================

!  Written for free-form input of parameters for slow-motional
!  calculations.  Returns a token consisting of nonseparator
!  characters (separators are space, tab, and ',') with all control
!  characters filtered out.
!  Special cases:  '(', ')', '=', '*', and end-of-line, which are
!  returned as single-character tokens.

! ------------------------------------------------------------------------

    subroutine gettkn(line,token,lth)
    implicit none
    include 'limits.inc'
    integer :: lth


    character line*(LINELG),token*(WORDLG),chr*1

    integer :: i,j,ichr,ichar

! *** Function definitions

    logical :: issepr,is1tok,isctrl,istab
    isctrl(chr) = ichar(chr).lt.32
    istab(chr) = ichar(chr).eq.9
    issepr(chr) = chr.eq.' '.or. ichar(chr).eq.9 .or. chr.eq.',' &
    .or.chr.eq.';'
    is1tok(chr) = chr.eq.'('.or.chr.eq.')'.or.chr.eq.'*' &
    .or.chr.eq.'='
! ***

!------------------------------------------
!  Find the next non-whitespace character
!------------------------------------------
    i=0
    2 i=i+1
    3 chr=line(i:i)

!     -------------------------
!      skip control characters
!     -------------------------
    if (isctrl(chr) .AND. .NOT. istab(chr)) then
        line(i:)=line(i+1:)
        go to 3
    end if

    if (issepr(chr) .AND. i < LINELG) goto 2

          
    if (i >= LINELG) then
        lth=0
        token=' '
        return
    end if

!     -----------------------------------
!     Check for single-character tokens
!     -----------------------------------
    if (is1tok(chr)) then
        token=chr
        lth=1
        line=line(i+1:)
        return
    end if

!----------------------------------------------------------------
! Place the next continuous string of characters in the token
! (stop at whitespace, punctuation, and single-character tokens)
!----------------------------------------------------------------
    j=i
    4 j=j+1
    5 chr=line(j:j)

!     -----------------------
!     Skip control characters
!     -----------------------
    if (isctrl(chr) .AND. .NOT. istab(chr)) then
        line(j:)=line(j+1:)
        go to 5
    end if

    if ( issepr(chr) .OR. is1tok(chr) ) then
        token=line(i:j-1)
        lth=j-i
        line=line(j:)
        return
    else
        go to 4
    end if
    end subroutine gettkn

!----------------------------------------------------------------------
!                     =========================
!                       subroutine UNGETT
!                     =========================
!  Replaces given token at the beginning of the given line
!  (Oh for the string functions of C..)
!----------------------------------------------------------------------
    subroutine ungett(token,lth,line)
    implicit none
    include 'limits.inc'
    character line*(LINELG),tmplin*(LINELG),token*(WORDLG)
    integer :: lth
    if (lth > 0 .AND. lth < LINELG) then
        tmplin=line
        line=token(:lth) // tmplin
    end if
    return
    end subroutine ungett
          
!----------------------------------------------------------------------
!                    =========================
!                       subroutine TOUPPR
!                    =========================
! Converts string argument to uppercase
!----------------------------------------------------------------------

    subroutine touppr(string,lth)
    implicit none
    include 'limits.inc'
    character string*(WORDLG),chr
    integer :: i,ich,ichar,lth

! *** Function definition
!      logical islc
!      islc(chr) = ichar(chr).ge.97 .and. ichar(chr).le.122
! ***
!----------------------------------------------------------------------

    do i=1,lth
        chr=string(i:i)
        ich=ichar(chr)
        if (ich >= 97 .AND. ich <= 122) string(i:i)=char( ich-32 )
    end do
    return
    end subroutine touppr

!----------------------------------------------------------------------
!                    =========================
!                       subroutine FTOKEN
!                    =========================
!  Decodes a token into a floating point number
!----------------------------------------------------------------------
    function ftoken( token,lth,val )
    implicit none
    include 'limits.inc'
    character token*(WORDLG),tkn*(WORDLG),tmptkn*(WORDLG),chr*1
    integer :: i,lth,idot,ibrk
    double precision :: val
    logical :: ftoken

! *** Function definitions --- these don't work with every FORTRAN
!     implementation.

    logical :: isdot,isexp,isdig
    isdot(chr)=chr.eq.'.'
    isexp(chr)=(chr .eq. 'd' .or. chr .eq. 'D') .or. &
    (chr .eq. 'e' .or. chr .eq. 'E')
    isdig(chr)=chr.ge.'0' .and. chr.le.'9'
! ***
!----------------------------------------------------------------------

    tkn=token
    idot=0
    ibrk=0

!----------------------------------------------------------------------
!     Find where a '.' is needed in the string
!     (this is to overcome the implied decimal used by some compilers)

!     Also, check for illegal characters (this is for FORTRAN compilers
!     that don't return an error when non-numeric characters
!     are encountered in the read.)
!----------------------------------------------------------------------
    do 10 i=1,lth
        chr=tkn(i:i)
        if (isdot(chr))  then
            idot=i
        else if (isexp(chr)) then
            ibrk=i
        else if ( .NOT. isdig(chr) .AND. chr /= '-') then
            go to 13
        end if
    10 END DO

    if (idot == 0) then
        if (ibrk == 0) then
            tkn=tkn(:lth)//'.'
        else
            tmptkn=tkn(ibrk:)
            tkn=tkn(:ibrk-1) // '.' // tmptkn
        end if
        lth=lth+1
    end if
     
    read(tkn,1000,err=13) val
    ftoken=.true.
    return

    13 ftoken=.false.
    return

    1000 format(bn,f20.10)
    end function ftoken


!----------------------------------------------------------------------
!                    =========================
!                       function ITOKEN
!                    =========================
!----------------------------------------------------------------------
    function itoken( token,lth,ival )
    implicit none
    include 'limits.inc'
    character*(WORDLG) token
    integer :: lth,ival
    logical :: itoken

!..................................................
!      read (token,1000,err=13) ival
! 12   itoken=.true.
!      return

! 13   itoken=.false.
!      return

! 1000 format(bn,i20)
!..................................................

    double precision :: fval
    logical :: ftoken
    external ftoken
    itoken=ftoken(token,lth,fval)
    ival=fval
    return
    end function itoken

!----------------------------------------------------------------------
!                       =========================
!                          function ITRIM
!                       =========================

!   Returns position of first blank in a string
!----------------------------------------------------------------------

    function itrim( string )
    implicit none
    integer :: itrim,j,lth
    character*(*) string
    lth=len(string)
    do 1 j=1,lth
        if (string(j:j) == ' ') go to 2
    1 END DO
    itrim=lth
    return

    2 itrim=j-1
    return
    end function itrim

!----------------------------------------------------------------------
!                    =========================
!                      subroutine INDTKN
!                    =========================

!     Looks for a secondary index specified by the series of tokens
!     '(' <n,m> { ')' } or '(' '*,*' { ')' }.  The information is
!     returned in the variables spectid and compid.  In the case of
!     no parameters or *, zeros are returned meaning all sites or
!     spectra.  Return -1 on error.

!----------------------------------------------------------------------
    subroutine indtkn(line,spectid,compid)
    implicit none
    include 'limits.inc'
    include 'stdio.inc'
    include 'miscel.inc'

    character line*(LINELG),token*(WORDLG)
    integer :: ival,lth,spectid,compid
    logical :: wldcrd
    logical :: itoken
    external itoken

!######################################################################

    call gettkn(line,token,lth)

!----------------------------------------------------------------------
!     Look for parenthesis indicating a second index will be specified
!----------------------------------------------------------------------
    if ( token /= '(' ) then		! got ()?
        spectid=0	! no
        compid=0
        call ungett(token,lth,line)	! restore last token to line
        return
    else	! yes, parse the info
        call gettkn(line,token,lth)
    
    !----------------------------------------------------------------------
    !     Check for a valid index: '*' or an integer in range
    !----------------------------------------------------------------------
        wldcrd=token.eq.'*'
    
    !                                   *** Empty '()': indices are 0
        if (token == ')') then
            spectid=0
            compid=0
            return
        else
        !                                   *** Wildcard: return 0
            if (wldcrd) then
                spectid=0
            else
            !                                   *** Check for legal number
            
                if (itoken(token,lth,ival)) then
                    spectid=ival
                else
                !                                            *** Illegal index
                
                    write(luttyo,1000) token(:lth)
                    spectid=-1
                    compid=-1
                    return
                end if
            end if
        ! now go for second index:
            call gettkn(line,token,lth)
            wldcrd=token.eq.'*'
            if (wldcrd) then
                compid=0		! second entry is *
                call gettkn(line,token,lth)
                if (token == ')') return
                call ungett(token,lth,line)
                return	! if not ), put it back and return
            else
                if (itoken(token,lth,ival)) then
                    compid=ival
                    call gettkn(line,token,lth)
                    if (token == ')') return
                    call ungett(token,lth,line)
                    return  ! if not ), put it back and return
                else
                    if (token == ')') then
                        compid=0
                        return
                    else
                        write(luttyo,1000) token(:lth)
                        spectid=-1
                        compid=-1
                        return
                    end if
                end if
            end if
        end if
    end if

    1000 format('*** Illegal index: ''',a,''' ***')
    end subroutine indtkn

