!----------------------------------------------------------------------
!                    =======================
!                         common PARMS
!                    =======================

!     Common block containing a list of quantities associated with
!     each of the fitting parameters being varied (i.e. each parameter
!     in the x array for LMDER-family programs).

!     NOTE: the file 'limits.inc' must be included before this file.

! prmax  : Vector containing upper limits for each parameter in x array
! prmin  : Vector containing lower limits for each parameter in x array
! prscl  : Vector containing desired absolute accuracy for each parameter
!          in x array
! xfdstp : Vector containing size of forward-difference step for each
!          parameter in x array
! xerr   : Vector containing uncertainty estimates for each parameter in x
! ibnd   : Flag for boundaries imposed on each parameter:
!          0=none, 1=minimum, 2=maximum, 3=both
! ixpr   : Index of each variable parameter into the fepr array in /eprprm/
! nser   : Number of spectra in the series
! nsparm : Number of series parameters set
!          (should be 8 or 9 before data command)
! nprm   : Number of parameters being varied
! matrxok: Define what is needed to give spectra.  0=need matrix,
!		1=have eval,evec, 2=have xspec, 3=have spec.
!----------------------------------------------------------------------

    double precision :: fparm,prmax,prmin,prscl,xfdstp,xerr, &
    stol,sstep,smn,smx,stval,fnorm
    integer :: iparm,idebug,ievec,iter,ihltcmd,ihltfit, &
    ixx,ibnd,ixpr,nser,nsparm,nprm,negegv

    logical :: gflag,aflag,rotflg,potflg,uniflg,dataOK
    integer :: matrxok
    character(6) :: tag
    common /parms/ fparm(NFPRM,MXSPEC,MXSITE),prmax(MXVAR), &
    prmin(MXVAR), &
    prscl(MXVAR),xfdstp(MXVAR),xerr(MXVAR), &
    stol,sstep,smn,smx,stval,fnorm, &
    iparm(NIPRM,MXSPEC,MXSITE),idebug,ievec, &
    gflag,aflag,rotflg,potflg,uniflg, &
    iter,ihltcmd,ihltfit, &
    ixx(NFPRM,MXSPEC,MXSITE),ibnd(MXVAR), &
    ixpr(MXVAR),nser,nsparm,nprm,negegv, &
    matrxok(MXSITE), &
    dataOK,tag(MXVAR)

