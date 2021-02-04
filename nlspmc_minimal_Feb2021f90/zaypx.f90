!  NLSPMC VERSION 1.0   2/5/99
!**********************************************************************

!             complex double precision vector scale and add
!             ---------------------------------------------

!        This subroutine will update a complex double precision
!        vector Y by the formula,

!                         Y=x+a*Y ,

!        where a is a complex double precision constant and X is a
!        complex double presion vector.

!       Includes:
!               nlsdim.inc
!               rndoff.inc

!       Uses:

!**********************************************************************

    subroutine zaypx(x,y,scale,ndim)

    include 'limits.inc'
    include 'rndoff.inc'

    integer :: ndim
    complex*16 scale

    complex*16 x,y
    dimension x(mxdim),y(mxdim)

    integer :: iel
    complex*16 czero
    parameter (czero=(0.0D0,0.0D0))

!######################################################################

    do 10 iel=1,ndim
        y(iel)=x(iel)+scale*y(iel)
    ! djs        if (abs(y(iel)).lt.rndoff) y(iel)=czero
    10 END DO

    return
    end subroutine zaypx
