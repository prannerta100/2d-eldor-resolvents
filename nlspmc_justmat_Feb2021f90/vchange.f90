!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                        function VCHANGE
!                    =========================

! check if the difference between two FP variables is significant
! .true. means they are different.  Must be double precision v1,v2.
!----------------------------------------------------------------------
    function vchange(v1,v2)

    implicit none
    logical :: vchange
    double precision :: v1,v2
    include 'rndoff.inc'

    vchange=.false.
    if (abs((v1+v2)/2.0D0) > rndoff) then
        vchange=(abs(v1-v2)/(abs(v1+v2)/2.0D0) .gt. 0.000001*abs(v1+v2))
    end if
    return
    end function vchange
