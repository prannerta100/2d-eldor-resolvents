! NLSPMC Version 1.0 2/5/99
    function dpmpar(i)
    implicit none
    double precision :: dpmpar
    integer :: i
!     **********

!     function dpmpar

!     this function provides double precision machine parameters

!     the function statement is

!       double precision function dpmpar(i)

!     where

!       i is an integer input variable set to 1, 2, or 3 which
!         selects the desired machine parameter. if the machine has
!         t base b digits and its smallest and largest exponents are
!         emin and emax, respectively, then these parameters are

!         dpmpar(1) = b**(1 - t), the machine precision,

!         dpmpar(2) = b**(emin - 1), the smallest magnitude,

!         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.

!     Argonne National Laboratory. MINPACK project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

!     *******
!     Adapted for use on IBM 3090  David Budil May 1991

!        t = 56, emin=-128, emax=127
!     **********

    real*8 :: dmach(3)
    integer*2 :: imach(12)
    equivalence(imach,dmach)

    data imach / 13328, 0, 0, 0, &
    16, 0, 0, 0, &
    32767, -1, -1, -1 /

    dpmpar = dmach(i)
    return
    end function dpmpar
