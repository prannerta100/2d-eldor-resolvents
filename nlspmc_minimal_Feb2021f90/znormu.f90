!  VERSION 1.0  (NLSPMC version)   2/5/99
!**********************************************************************

!                                ZNORMU
!                                ======

!       ZNORMU is a complex double precision function for calculating
!       the Euclidean pseudo-norm of a complex vector, i.e. the square
!       root of the sum of the squares of the components, NOT the
!       square root of the sum of the moduli of the components.

!       Calling Parameters:
!       ------------------

!       vectors:

!               v       : a NDIM dimensional complex vector whose
!                         Euclidean pseudo-norm is desired.

!       scalars:

!               ndim    : number of elements in the vector v.


!       Local Variables:
!       ---------------

!       vectors:

!       scalars:

!               ir      : index for loop over components of v
!               amp     : amplitude of pseudo-norm of v.
!               phase   : phase of pseudo-norm of v.
!               tempr   : temporarily stores the real part of the
!                         pseudo-norm when checking for roundoff.
!               tempi   : temporarily stores the real part of the
!                         pseudo-norm when checking for roundoff.
!               accum   : holds the partial sums in loop over the
!                         in the loop over the components of v.

!       Notes:
!       ------
!               1)  The size of the vector v is determined by the
!                   parameter MXDIM in the include file stddim.inc.
!                   In the calling routine, the vector v should be
!                   dimensioned as complex*16 v(mxdim), or
!                   alternatively as real*8 v(2,mxdim).

!               2)  The square root of the final value of accum is
!                   not calculated using the usual Fortran square root
!                   function.  The methods used here allows one to
!                   take the square root of a number on the negative
!                   real axis.

!               3)  The standard include files are used to specify
!                   the dimensions of arrays, the unit roundoff error
!                   of the machine, etc.

!       written by DJS 5-12-87

!       Includes:
!               nlsdim.inc
!               rndoff.inc

!       Uses:

!**********************************************************************

    function znormu(v,ndim)

    complex*16 znormu

    include 'limits.inc'
    include 'rndoff.inc'

    integer :: ndim
    double precision :: v
    dimension v(2,mxdim)

    integer :: ir
    double precision :: amp,phase,accv,accr,acci

!######################################################################

    accv=0.0D0
    accr=0.0D0
    acci=0.0D0

    do 10 ir=1,ndim
        accv=accv+v(1,ir)*v(1,ir)+v(2,ir)*v(2,ir)
        accr=accr+v(1,ir)*v(1,ir)-v(2,ir)*v(2,ir)
        acci=acci+2.0D0*v(1,ir)*v(2,ir)
    10 END DO

    amp=sqrt(sqrt(accr*accr+acci*acci))
    phase=0.5D0*datan2(acci,accr)

    if (amp/sqrt(accv) < rndoff) amp=0.0D0

    accr=amp*cos(phase)
    acci=amp*sin(phase)

    znormu=dcmplx(accr,acci)

    return
    end function znormu
