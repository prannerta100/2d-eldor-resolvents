!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!     Declarations of common block to hold starting vector for off-
!     diagonal & diagonal space.

!     NOTE: the file 'nlsdim.inc' must be included before this file.

!*********************************************************************

    complex*16 stvo,stvd
    common	/stvcom/  stvo(mxdim),stvd(mxdim,5)
