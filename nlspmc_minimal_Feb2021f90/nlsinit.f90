! NLSPMC Version 1.0 2/5/99
!----------------------------------------------------------------------
!                    =========================
!                       subroutine NLSINIT
!                    =========================

!  Extract initialization duties from nlspc.f

!     Initializes the following:
!       NLS parameter arrays
!       NLS convergence criteria
!       miscellaneous
!----------------------------------------------------------------------

    subroutine nlsinit
    implicit none

!      include 'nlsdim.inc'
!      include 'iterat.inc'
!      include 'nlsnam.inc'
!      include 'stdio.inc'
!      include 'expdat.inc'
!      include 'parcom.inc'
!      include 'eprprm.inc'
!      include 'prmeqv.inc'
!      include 'lmcom.inc'
    include 'stdio.inc'
    include 'limits.inc'
    include 'names.inc'
    include 'parms.inc'
!      include 'parmequ.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'lmcomm.inc'
    include 'basis.inc'
    include 'miscel.inc'

    integer :: ix,ix2,ix3

!      tempvar = .false.		! do we vary temperature in series?
!      psivar = .false.		! do we vary psi (why this variable?)
    nspectra = 1		! One simulation is default, unrelated
!  to presence or absence of data
    ncomps = 1		! one component is default

    ihltcmd=0

    nfiles=0			! # of experimental data files
    lucmd=luttyi
    luout=luttyo
    luecho=luttyo
    nspc=0
    nptot=0
    nptotcw=0
    nprm=0
    nser=0
    nsparm=0
    negegv=0
    ldelta=1
    kdelta=1

!----------------------------------------
!     Initialize parameter arrays
!----------------------------------------
    do 22 ix=1,NFPRM
        do 22 ix2=1,MXSPEC
            do 22 ix3=1,MXSITE
                fparm(ix,ix2,ix3)=0.0d0
                ixx(ix,ix2,ix3)=0
    22 END DO
! special values:
    do 222 ix3=1,MXSITE
        diaflg=.false.
        matrxok(ix3)=0	! need matrix!
        do 222 ix2=1,MXSPEC
            fparm(IDWGT,ix2,ix3)=1.0d0  ! spectrum weighting
            fparm(ISWGT,ix2,ix3)=1.0	! site weighting default
            fparm(ICGTOL,ix2,ix3)=1.0D-4
            fparm(ISHIFTR,ix2,ix3)=0.1D0
        ! parm(ISHIFTR,ix2,ix3)=2.0D0
    222 END DO

    do 23 ix=1,NIPRM
        do 23 ix2=1,MXSPEC
            do 23 ix3=1,MXSITE
                iparm(ix,ix2,ix3)=0
    23 END DO

    do 233 ix2=1,MXSPEC
        do 233 ix3=1,MXSITE
            iparm(INORT,ix2,ix3)=1
            iparm(ILEMX,ix2,ix3)=8
            iparm(ILEMX+1,ix2,ix3)=5
            iparm(ILEMX+2,ix2,ix3)=4
            iparm(ILEMX+3,ix2,ix3)=2
            iparm(ILEMX+4,ix2,ix3)=2
            iparm(IITYPE,ix2,ix3)=1
            iparm(ISPTST,ix2,ix3)=0
        ! basis set information
            iparm(INDIMO,ix2,ix3)=0
            iparm(INDIMD,ix2,ix3)=0
        !	  basindx(ix2+(ix3-1)*MXSPEC)=0
        !	  triindx(ix2+(ix3-1)*MXSPEC)=0
            do 233 ix=1,5
                basinfo(ix,ix2,ix3)=0
    233 END DO
    nbasis=0		! no basis sets yet.
    pidptr(1)=1	! starting location for basis storage
    dpidptr(1)=1	! starting location for basis storage
    basisptr=0	! no basis set mapped into jeq1 etc.

    aflag=.false.
    gflag=.false.
    rotflg=.false.
    potflg=.false.
    uniflg=.false.

!  -- Set initial values for convergence criteria

    xtol=1.0d-6
    ftol=1.0d-6
    gtol=1.0d-8
    factor=1.0d2
    maxev=100
    maxitr=10
    maxfit=1
    itrace=0
    idebug=0

!  -- Set initial values for line search parameters

    sstep=0.05d0
    stol=1.0d-3
    return
    end subroutine nlsinit
