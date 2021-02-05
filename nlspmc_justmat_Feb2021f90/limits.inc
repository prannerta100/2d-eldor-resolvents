!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************
!                    =========================
!                      include file LIMITS
!                    =========================

!  Defines array dimensioning parameters for NLS programs and other
!  numerical values

! LINELG  Length of line read
! WORDLG  Length of words
! NAMELG  Length of filenames
! NFPRM   Number of floating point parameters for simulation
! NIPRM   Number of integer point parameters for simulation
! MXSPEC  Maximum number of experimental spectra to fit
! MXSITE  Maximum number of components of each spectra
! MXVAR   Maximum number of simultanious variables for least squares
! MXORT   Maximum MOMD orientations
! MXDIM   Maximum dimension of matrix
! MMXDIM  Max combined size of basis sets = MXDIM*MXSITE
! MXSTEP  Maximum number CG steps ??
! MXEL    Maximum matrix elements
! MXEGV   Maximum number of eigenvalues
! MXFILE  Maximum number of data files
! MXNPT   Maximum number of points along one of the axis
!         NOTE: THIS SHOULD BE A POWER OF 2
! MXSPT   Maximum number of points in an individual spectrum
!              (can be mxnpt*mxnpt maximum)
! MXPT    Maximum total number of data points
!              (can be MXSPT*MXSPEC*MXSITE maximum) - store all spec.
! MXSPH   Maximum number of tensors that may be specified in spherical
!         form
! PG added MXDFT
! MXDFT      max array size of ip, we take a DFT of ip

!*********************************************************************

    integer :: LINELG,WORDLG,NAMELG,NFPRM,NIPRM,MXSPEC,MXSITE,MXVAR, &
    MXORT,MXDIM,MMXDIM,MXSTEP,NVPRM, &
    MXEL,MXEGV,MXFILE, &
    MXNPT,MXSPT,MXPT, &
    MXSPH,NALIAS,NSYMTR,NSYMBL,MXDFT,MXPOL

! later remove namelg and set all to 40.  see call to gettkn
! added 5 site weights as variable parameters.  Don't vary them
! all!
! mxdim was 800,
! mxsite was 3
! mxort was 20
! mxegv was 1024 ywc
! mxstep was 1024
! mxpt was 163840
! MAXDIM was 1700
!      parameter (LINELG=80,WORDLG=40,NAMELG=40,NFPRM=52,NIPRM=25,
!     #	MXSPEC=10,MXSITE=2,MXVAR=15,MXORT=40,MXDIM=5000,
!     #	MMXDIM=MXDIM*MXSITE,MXSTEP=1012,
!     #  NVPRM=41,MXEL=50000,MXEGV=128,MXFILE=20,
!     #           MXNPT=512,MXSPT=40000,MXPT=16384,
!     #           MXSPH=3,NALIAS=9,NSYMTR=3,NSYMBL=5,MXDFT=32,
!     #           MXPOL=10)
    parameter (LINELG=80,WORDLG=40,NAMELG=40,NFPRM=52,NIPRM=25, &
    MXSPEC=10,MXSITE=1,MXVAR=15,MXORT=16,MXDIM=50050, &
    MMXDIM=MXDIM*MXSITE,MXSTEP=1024, &
    NVPRM=41,MXEL=500000,MXEGV=256,MXFILE=20, &
    MXNPT=128,MXSPT=16384,MXPT=163840, &
    MXSPH=3,NALIAS=9,NSYMTR=3,NSYMBL=5,MXPOL=10,MXDFT=32)

!     mxpol=20

!*********************************************************************

