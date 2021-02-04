!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       ================
!                       SUBROUTINE CSVAL
!                       ================

!     This subroutine calculates the distinct eigenvalues of a given
!     tridiagonal matrix.  The identification test for good versus
!     spurious eigenvalues is done according to the procedure in
!     Cullum & Willoughby's book.

!       Arguments:
!           alpha    diagonal T matrix elements
!           beta     extradiagonal T matrix elements
!           v1,v2    work spaces
!           gr,gc    work spaces
!           vs       work space, on exit = eigenvalues of T-matrix sorted
!                    in the order of decreasing weighting factor
!           w        weighting factor
!           ngood    # of good eigenvalues (passed via eprdat.inc common)
!           ndiag    dimension of T-matrix used in the calculation of the
!                    eigenvalues (passed via eprdat.inc common)
!           ndone    dimension of T-matrix generated in Lanczos routine
!           ierr     error flag for failure in cmtqli routine

!     modified from CSLEVAL.F by SHL 17-AUG-92

!     Includes :
!               nlsdim.inc
!               stdio.inc
!               eprprm.inc
!               rndoff.inc
!     Uses :
!               cmtqli.f
!               compev.f
!               lump.f
!               comgap.f
!               isoev.f
!               inverr.f

!*********************************************************************

    subroutine csval(alpha,beta,v1,v2,gr,gc,vs,w,ndone,ierr)

    implicit none

    include 'limits.inc'
    include 'stdio.inc'
    include 'simparm.inc'
    include 'rndoff.inc'

    integer ::  mxinit,mev,ndis,loop,ndone
    integer ::  i,j,niso,it,iwrite,l,m,mmb,ierr
    integer ::  mp(mxstep),mp2(mxstep)

    real*8 ::  g(mxstep),gg(mxstep),gtp,ggtp

    double precision ::  gr(mxstep),gc(mxstep)
    double precision ::  gaptol,ttol,epsm,reltol,evmax,t0,t1
    double precision ::  scale1,scale2,sputol,contol,multol
    double precision ::  one,zero
    parameter (one=1.0D0,zero=0.0D0)

    complex*16  alpha(mxstep),beta(mxstep),vs(mxstep),w(mxstep)
    complex*16  v1(mxstep),v2(mxstep)
    complex*16  betam,wt,egv
    complex*16  czero,ci
    parameter (czero=(0.0D0,0.0D0),ci=(0.0D0,1.0D0))

!#####################################################################

    if (sptst == 1) go to 99

!---------------------------------------------------------------------
!     Eigenvalue calculation without test for spurious eigenvalues
!---------------------------------------------------------------------

    if (ndone > mxstep) then
    end if
!      stop
    do 82 i=1,ndone
        vs(i)=alpha(i)
        v1(i)=beta(i)
    82 END DO

    call cmtqli(ndone,vs,v1,w,1,ierr)

    if (ierr /= 0) then
        write(lulog,1000)
        if (lulog /= luttyo) write(luttyo,1000)
        ierr=-abs(ierr)
        return
    end if

    ngood=ndone
    ndiag=ndone

    return

!---------------------------------------------------------------------
!     Eigenvalue calculation with test for spurious eigenvalues
!---------------------------------------------------------------------

    99 epsm=100.0D0*rndoff
    scale1=5.0D2
    scale2=5.0D0
    gaptol=1.0D-7
    reltol=1.0D-8
    mxinit=10
    mev=ndone

    betam=beta(mev)
    multol=500.D0*dfloat(mev+1000)*epsm
    sputol=multol

!---------------------------------------------------------------------
!     calculate the eigenvalues and test for the spurious values
!---------------------------------------------------------------------

    call compev(alpha,beta,v1,v2,vs,w,gr,multol,sputol,mp,mp2, &
    mev,ndis,ierr)

    if (ierr /= 0) then
        write(lulog,1000)
        if (lulog /= luttyo) write(luttyo,1000)
        ierr=-abs(ierr)
        return
    end if

    if (ndis == 0) then
        write (lulog,1010)
        if (lulog /= luttyo) write (luttyo,1010)
        ierr=-abs(ierr)
        return
    end if

    evmax=gr(ndis)
    loop=ndis

!---------------------------------------------------------------------

    call lump(vs,v1,gr,w,reltol,sputol,scale2,mp,mp2,loop)

    ndis=loop

!---------------------------------------------------------------------
!     calculate mingaps for distinct t(1,mev) eigenvalues.
!---------------------------------------------------------------------

    call comgap(vs,gr,gg,mp,mp2,ndis)

!     set convergence critierion

    ttol=epsm*evmax
    contol=cdabs(betam)*1.D-10

    beta(mev)=betam

!---------------------------------------------------------------------

    call isoev(vs,gr,gg,gaptol,sputol,scale1,mp,ndis,niso)
!      stop

    if (niso == 0) go to 370

!---------------------------------------------------------------------
!     error estimates for isolated good T-eigenvalues
!---------------------------------------------------------------------

    it=mxinit

    do 330 i=1,ndone
        gr(i)=zero
        gc(i)=zero
    330 END DO

    gr(1)=one
    gc(1)=one

    call inverr(alpha,beta,v1,v2,vs,epsm,gr,gc,g,gg, &
    mp,mp2,mev,ndis,niso,it)

!      write(luttyo,340) contol

    do 350 i=1,niso
        if (abs(g(i)) > contol) then
            write (lulog,1020)
            if (lulog /= luttyo) write (luttyo,1020)
            go to 370
        end if
    350 END DO

!      write(luttyo,360) contol

    370 ngood=0
    do 390 i=1,ndis
        if (mp(i) == 0) go to 390
        ngood=ngood+1
        mp(ngood)=mp(i)
        vs(ngood)=vs(i)
        w(ngood)=w(i)
        g(ngood)=g(i)
        gg(ngood)=gg(i)
    390 END DO

!---------------------------------------------------------------------
!     sort goodev in the order of decreasing weighting factor
!---------------------------------------------------------------------

!      write (luttyo,440) ngood,niso,ndis

    do 450 l=1,ngood-1
        t0=cdabs(w(l))
        do 450 m=l+1,ngood
            t1=cdabs(w(m))
            if (t1 > t0) then
                t0=t1
                wt=w(m)
                w(m)=w(l)
                w(l)=wt
                egv=vs(m)
                vs(m)=vs(l)
                vs(l)=egv
                gtp=g(m)
                g(m)=g(l)
                g(l)=gtp
                ggtp=gg(m)
                gg(m)=gg(l)
                gg(l)=ggtp
            end if
    450 END DO

    ndiag=mev

    beta(mev)=betam

!     end of loop on different size T-matrices allowed.

!      write(luttyo,750) ngood,ndiag
!      write(lulog,750) ngood,ndiag

    return

!=====================================================================
!     format statements
!=====================================================================

    340 format(' convergence is tested using the convergence tolerance', &
    e13.4)
    360 format(' all computed error estimates were less than',e15.4/ &
    ' therefore procedure terminates')
    440 format(/i6,' good T-eigenvalues have been computed'/ &
    i6,' of these are isolated'/ &
    i6,'=number of distinct T-eigenvalues computed'/)
    750 format(/,2x,'# of good eigenvalues : ',i4,/,2x,'Actual ', &
    'Lanczos steps used in the eigenvalue estimation : ',i4)
    1000 format(2x,'on return from CMTQLI, error flag was not zero'/)
    1010 format(2x,'on return from COMPEV no distinct eigenvalues ', &
    'were found'/)
    1020 format(2x,'inverse iteration did not converge'/)

    end subroutine csval
