! NLSPMC Version 1.0 2/5/99
!----------------------------------------------------------------------
!                    =========================
!                      subroutine COMPS
!                    =========================

! comp <ncomps>

!      ncomps  : number of components for all spectra

!  NOTE: must apply to all spectra, ncomps defaults to 1


! Special rules:

!----------------------------------------------------------------------
    subroutine comps ( line )
    implicit none
    character(80) :: line

    include 'limits.inc'
    include 'simparm.inc'
    include 'parmequ.inc'
    include 'parms.inc'
    include 'lpnam.inc'
    include 'stdio.inc'
    include 'miscel.inc'

    integer :: ival,lth
    character token*30
    logical :: itoken

    call gettkn(line,token,lth)
!                                        *** No value given
    if (lth == 0) then
        write(luttyo,*) 'Number of components expected'
        return
    end if

    if (itoken(token,lth,ival)) then
    ! check max ival
        ncomps=ival
        return
    else
        write(luttyo,1001) token(:lth)
        return
    end if
    1001 format('*** number of components expected: ''',a,''' ***')
    end subroutine comps
