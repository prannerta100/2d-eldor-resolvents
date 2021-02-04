! NLSPMC Version 1.0 2/5/99
!*********************************************************************

!                       ================
!                       SUBROUTINE ISOEV
!                       ================

!     Subroutine ISOEV uses tmingaps to label the isolated good
!     T-eigenvalues that are very close to spurious ones.
!     Error estimates will not be computed for these T-eigenvalues.

!     on entry and exit:
!        vs contains the computed distinct eigenvalues of t(1,mev)
!        gr(k) = |vs(k)|, k = 1,ndis, gr(k).le.gr(k+1)
!        gg(k) = min(j.ne.k,|vs(k)-vs(j)|)   mingap
!        mp contains the corresponding T-multiplicities
!        ndis = number of distinct T-eigenvalues
!        gaptol = relative gap tolerance set in main

!     on exit:
!        mp(j) is not changed except that  mp(j)=-1, if mp(j)=1,
!        and a spurious T-eigenvalue is too close.

!        if mp(i)=-1 that simple good T-eigenvalue will be skipped
!        in the subsequent error estimate computations in inverr
!        that is, we compute error estimates only for those good
!        T-eigenvalues with mp(j)=1.

!     includes :

!        nlsdim.inc
!        eprprm.inc

!*********************************************************************

    subroutine isoev(vs,gr,gg,gaptol,sputol,scale1,mp,ndis,niso)

    implicit none

    include 'limits.inc'
    include 'simparm.inc'

    complex*16 vs(mxstep),t0
    double precision ::  sputol,gaptol,scale1,temp,tol,tj,dgap,one
    double precision ::  gr(mxstep)
    real :: gg(mxstep)
    integer ::  i,j,ndis,niso,mp(mxstep)

!---------------------------------------------------------------------

    one=1.0D0
    dgap=scale1*sputol
    niso=0
    do 40 j=1,ndis
        if (mp(j) /= 1) go to 40
        tj=gr(j)
        t0=vs(j)
        tol=dmax1(dgap,gaptol*tj)
    !     tol=dmax1(one,tj)*gaptol
    
    !     vs(j) is next simple good T-eigenvalue
    
        niso=niso+1
        if (abs(gg(j)) > tol) go to 40
        i=j
        10 i=i-1
        if (i < 1) go to 20
        if (tj-gr(i) > tol) go to 20
        if (mp(i) /= 0) go to 10
        temp=cdabs(t0-vs(i))
        if (temp > tol) go to 10
        mp(j)=-mp(j)
        niso=niso-1
        go to 40
        20 i=j
        30 i=i+1
        if (i > ndis) go to 40
        if (gr(i)-tj > tol) go to 40
        if (mp(i) /= 0) go to 30
        temp=cdabs(t0-vs(i))
        if (temp > tol) go to 30
        mp(j)=-mp(j)
        niso=niso-1
    40 END DO

    return

    end subroutine isoev
