!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       =================
!                       SUBROUTINE COMGAP
!                       =================

!     Subroutine COMGAP computes the minimal gaps between
!     distinct T-eigenvalues supplied.

!     on return :
!        |gg(k)| = min(j.ne.k,|vs(k)-vs(j)|)
!        mp2(k)=j index for minimum.

!        gg(k) < 0 means nearest neighbor is spurious.

!     includes :

!        nlsdim.inc

!*********************************************************************

    subroutine comgap(vc,va,gg,mp,ind,m)

    implicit none

    include 'limits.inc'

    complex*16 vc(mxstep),z
    double precision ::  va(mxstep),t0,t1,tu,tk
    real ::  gg(mxstep),gtemp
    integer ::  index,itemp,j,k,kk,km1,kp1,l,m
    integer ::  mp(mxstep),ind(mxstep)

!---------------------------------------------------------------------
!     va(k) = |vc(k)|  va(k) <= va(k+1)
!     gg(k) = min |vc(k)-vc(j)|  j .ne. k.
!---------------------------------------------------------------------

    tu=va(m)+va(m)
    k=0
    10 k=k+1
    if (k > m) return
    index=0
    t1=tu
    tk=va(k)
    z=vc(k)
    j=k
!     backwards
    20 j=j-1
    if (j < 1) go to 30
    t0=tk-va(j)
    if (t0 > t1) go to 30
    t0=cdabs(z-vc(j))
    if (t1 <= t0) go to 20
    t1=t0
    index=j
    go to 20
!     forwards
    30 j=k
    40 j=j+1
    if (j > m) go to 50
    t0=va(j)-tk
    if (t0 > t1) go to 50
    t0=cdabs(z-vc(j))
    if (t1 <= t0) go to 40
    t1=t0
    index=j
    go to 40
    50 ind(k)=index
    gg(k)=t1
    if(mp(index) == 0)  gg(k)=-gg(k)
    go to 10

    end subroutine comgap
