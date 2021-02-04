! NLSPMC Version 1.0 2/5/99
!----------------------------------------------------------------------
!                    =========================
!                      subroutine SPECTRA
!                    =========================

! spec <nspectra>

!      nspectra  : number of spectra

!  NOTE: defaults to 1


!----------------------------------------------------------------------
    subroutine spectra ( line )
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
        write(luttyo,*) 'Number of spectra expected'
        return
    end if

    if (itoken(token,lth,ival)) then
        if (ival > MXSPEC) then
            write(luttyo,1002) MXSPEC
            1002 format('*** Maximum of ''',i2,''' spectra allowed ***')
            nspectra = MXSPEC
            return
        else
            nspectra = ival
            return
        end if
    else
        write(luttyo,1001) token(:lth)
        return
    end if
    1001 format('*** number of spectra expected: ''',a,''' ***')
    end subroutine spectra
