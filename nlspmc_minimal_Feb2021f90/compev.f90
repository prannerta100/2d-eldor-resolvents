!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       =================
!                       SUBROUTINE COMPEV
!                       =================

!     Subroutine COMPEV calls subroutine CMTQLI to compute the
!     T-eigenvalues.  COMPEV then applies the T-multiplicity and
!     spurious tests to the computed T-eigenvalues.

!     on return from COMPEV:
!        ndis=number of distinct eigenvalues of T(1,mev)
!        vs=distinct T-eigenvalues in increasing order of magnitude
!        gr(k)=|vs(k)|, k=1,ndis,  gr(k).le.gr(k+1)
!        mp=T-multiplicities of the T-eigenvalues in vs
!        mp(i)=(0,1,mi), mi>1, i=1,ndis  means:
!           (0)  vs(i) is spurious
!           (1)  vs(i) is simple and good
!           (mi) vs(i) is T-multiple and is therefore not only good but
!                also a converged good T-eigenvalue.

!     includes:
!           nlsdim.inc
!           eprprm.inc
!     uses:
!           cmtqli.f

!*********************************************************************

    subroutine compev(alpha,beta,v1,v2,vs,w,evmag,multol,sputol, &
    mp,t2flag,mev,ndis,ierr)

    implicit none

    include 'limits.inc'
    include 'simparm.inc'

    complex*16 alpha,beta,vs,v1,v2,w,w2,eval,ctemp,wtemp,czero
    dimension alpha(mxstep),beta(mxstep),vs(mxstep),v1(mxstep), &
    v2(mxstep),w(mxstep),w2(mxstep)
    parameter (czero=(0.0D0,0.0D0))
    double precision :: evmag(mxstep)
    double precision :: temp,dgap,tol,delmin
    double precision :: multol,sputol,evalr,evalc
    integer :: i,j,mev,mev1,imin,index,jp1,isort,ndis,ierr,mp,t2flag
    dimension mp(mxstep),t2flag(mxstep)

!---------------------------------------------------------------------

    mev1=mev - 1

    do 50 j=1,mev
        vs(j)=alpha(j)
        v1(j)=beta(j)
    50 END DO

    v1(mev)=czero

    isort=0
    call cmtqli(mev,vs,v1,w,isort,ierr)
    if (ierr /= 0) return

!     T-eigenvalues are in vs in increasing order of magnitude

    do 100 j=1,mev
        evmag(j)=cdabs(vs(j))
    100 END DO

    multol=multol*evmag(mev)
    sputol=sputol*evmag(mev)
    tol=1000.0d0*sputol

!----------------------------------------------------------------------
!     T-multiplicity determination
!----------------------------------------------------------------------

    j=0
    ndis=0
    do 150 i=1,mev
        t2flag(i)=0
    150 END DO

    160 j=j+1

    if (j <= mev) then
        if (t2flag(j) == 1) then
            go to 160
        else
            wtemp=w(j)
            ctemp=vs(j)
            eval=ctemp
            temp=evmag(j)
            ndis=ndis+1
            index=1
            t2flag(j)=1
            i=j
            170 i=i+1
            if (i > mev) go to 180
            if (t2flag(i) == 1) go to 170
            dgap=evmag(i)-temp
            if (dgap > multol) go to 180
            dgap=cdabs(eval-vs(i))
            if (dgap > multol) go to 170
        
        !     T-multiplicity increases
        
            index=index+1
            ctemp=ctemp+vs(i)
            wtemp=wtemp+w(i)
            t2flag(i)=1
            go to 170
        
        !     T-multiplicity for vs(ndis) has been determined
        
            180 vs(ndis)=ctemp/dble(index)
            w(ndis)=wtemp/dble(index)
            mp(ndis)=index
            go to 160
        end if
    end if

!----------------------------------------------------------------------
!     T(2,mev) eigenvalue calculation
!----------------------------------------------------------------------

    do 210 j=1,mev1
        jp1=j+1
        v2(j)=alpha(jp1)
        v1(j)=beta(jp1)
    210 END DO

    v1(mev1)=czero

    isort=0
    call cmtqli(mev1,v2,v1,w2,isort,ierr)
    if (ierr /= 0) return

    do 250 j=1,mev1
        evmag(j)=cdabs(v2(j))
    250 END DO

!----------------------------------------------------------------------
!     test for the spurious eigenvalues
!----------------------------------------------------------------------

    do 280 i=1,mev1
        t2flag(i)=0
    280 END DO

!     go through the eigenvalues of t2-hat.  find the closest eigenvalue
!     of T(1,mev).  if it is T-multiple go on.  if it is simple declare
!     spurious whenever delmin < sputol by setting mp(i)=0

    j=0
    290 j=j+1
    if (j > mev1) go to 390

    temp=evmag(j)
    eval=v2(j)
    evalr=temp+sputol
    evalc=temp-sputol
    delmin=2.D0*cdabs(vs(mev))
    imin=0

!     backward search

    i=j+1
    310 i=i-1
    if(i < 1) go to 320
    if(i > ndis) i=ndis

    temp=cdabs(vs(i))
    if (temp < evalc) go to 320
    if(mp(i) == 0) go to 310
    dgap=cdabs(vs(i)-eval)
    if (dgap >= delmin) go to 310
    delmin=dgap
    imin=i

    go to 310

!     forward search

    320 i=j
    330 i=i+1
    if(i > ndis) go to 340

    temp=cdabs(vs(i))
    if (temp > evalr) go to 340
    if(mp(i) == 0) go to 330
    dgap=cdabs(vs(i)-eval)
    if (dgap >= delmin) go to 330
    delmin=dgap
    imin=i

    go to 330

    340 continue

    if(imin == 0) go to 370

    if(delmin > sputol) go to 290
    if(mp(imin) > 1)  go to 290
    mp(imin)=0

    go to 290

    370 continue

    go to 290

    390 continue

    do 400 j=1,ndis
        evmag(j)=cdabs(vs(j))
    400 END DO

    return

    end subroutine compev
