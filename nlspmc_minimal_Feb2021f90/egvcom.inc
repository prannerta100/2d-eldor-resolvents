!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!     Declarations of common block to hold eigenvalues and eigenvectors

!     NOTE: the file 'nlsdim.inc' must be included before this file.

!*********************************************************************

    complex*16 egvalx,egvalz,egvecx,egvecz
    common	/egvcom/  egvalx(MXDIM,MXORT,MXSITE), &
    egvalz(MXDIM,MXORT,MXSITE)
    common /bigg/ egvecx(MXDIM,MXEGV,MXORT,MXSITE)
    common /bigg2/ egvecz(MXDIM,MXEGV,MXORT,MXSITE)
!         common /bigg/ egvecx(MXDIM,256,MXORT,MXSITE)
