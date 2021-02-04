!       NLSPMC VERSION 2/5/99
!*********************************************************************

!       This complex double precision function returns the dot
!       product of the two input vectors without complex
!       conjugation.

!       Calling Parameters:
!       ------------------

!       vectors:
!               x    : input vector
!               y    : input vector

!       scalars:
!               ndim : dimension of the vectors x and y


!       Local Variables:
!       ---------------
!               ir   : row counter
!               accr : real part of dot product
!               acci : imaginary part of dot product

!       Notes:
!       -----
!               It is assumed that the input vectors are stored in a
!               double precision complex array or as a real double
!               precision array dimensioned as 2 X mxdim in the
!               calling routine.


!       Includes:
!               nlsdim.inc
!               rndoff.inc

!       Uses:

!       written by DJS 26-AUG-87

!*********************************************************************

    function zdotu2(x,y,ndim)

    complex*16 zdotu2

    include 'limits.inc'
    include 'rndoff.inc'

    integer :: ndim
    double precision :: x,y
    dimension x(2,mxdim),y(2,mxdim)

    integer :: ir
    double precision :: accx,accy,accr,acci,scale

!######################################################################

    accx=0.0D0
    accy=0.0D0
    accr=0.0D0
    acci=0.0D0

    do 10 ir=1,ndim
        accx=accx+x(1,ir)*x(1,ir)+x(2,ir)*x(2,ir)
        accy=accy+y(1,ir)*y(1,ir)+y(2,ir)*y(2,ir)
        accr=accr+x(1,ir)*y(1,ir)-x(2,ir)*y(2,ir)
        acci=acci+x(1,ir)*y(2,ir)+x(2,ir)*y(1,ir)
    10 END DO

    scale=sqrt(accx*accy)

    if (dabs(accr)/scale < rndoff) accr=0.0D0
    if (dabs(acci)/scale < rndoff) acci=0.0D0

    zdotu2=dcmplx(accr,acci)

    return
    end function zdotu2
