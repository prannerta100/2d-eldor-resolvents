!----------------------------------------------------------------------
!                       =========================
!                          common block LMCOM
!                       =========================
!     Work arrays for Levenberg-Marquardt least squares fitting
!     procedure and parameters controlling termination criteria.

!     NOTE: File 'nlsdim.inc' should be included before this file
!           in the program source code.

!----------------------------------------------------------------------
    double precision :: corr,diag,fjac,fvec,gnvec,gradf,qtf, &
    work1,work2,work3,work4,xxx,ftol,xtol,gtol,factor
    integer :: ipvt,maxev,maxitr,maxfit,nprint,n,m,itrace,mode, &
    info,jacobi
    common /lmcom/ fjac(mxpt,MXVAR),fvec(mxpt),xxx(MXVAR), &
    diag(MXVAR),qtf(MXVAR),corr(MXVAR,MXVAR), &
    work1(MXVAR),work2(MXVAR),work3(MXVAR), &
    work4(mxpt),gnvec(MXVAR),gradf(MXVAR), &
    ftol,gtol,xtol,factor, &
    ipvt(MXVAR),n,m,nprint,maxev,maxitr,maxfit, &
    itrace,jacobi,mode,info

    integer :: nflmpr,nilmpr
    parameter (nflmpr=4,nilmpr=4)

    double precision :: flmprm(nilmpr)
    integer :: ilmprm(nilmpr)
    equivalence (flmprm(1),ftol),(ilmprm(1),maxev)

