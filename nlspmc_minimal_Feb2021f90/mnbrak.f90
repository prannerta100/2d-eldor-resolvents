!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                         =======================
!                            subroutine MNBRAK
!                         =======================

!   Modified from the routine in NUMERICAL RECIPES to restrict the search
!   region according to the user-specified range.  The error flag (ierr=1)
!   is returned if it does not find the bracket for the minimum.  Note
!   that the region for the parabolic extrapolation to be valid is reduced
!   (glimit=10.0).

!   Error flag of negative value means the function "func" failed.

!   Given a function func, and given distinct initial points ax and bx,
!   routine searches in the downhill direction (defined by the function as
!   evaluated at the initial points) and returns new points ax, bx, cx
!   that bracket a minimum of the function.  Also returned are the function
!   values at the three points, fa, fb, and fc.
!----------------------------------------------------------------------

    subroutine mnbrak(ax,bx,cx,fa,fb,fc,func,smin,smax,ierr)

    implicit none

    include 'limits.inc'
    include 'parms.inc'
!      include 'iterat.inc'

    double precision :: ax,bx,cx,fa,fb,fc,smin,smax,func
    external func

    integer :: ierr,idir,i
    double precision :: dum,fu,r,q,u,ulim,bnd

    double precision :: glimit,gold,tiny
    parameter (gold=1.618034d0,glimit=10.,tiny=1.d-20)

!  -- Rearrange the two points so that a to b is always "downhill"


    bnd=smax
    idir=1

    fa=func(ax,ierr)
    if (ierr < 0 .OR. ihltcmd /= 0) return
! loop on sites
!      do 101 i=1,ncomps
!        seteval(i)=.true.	! got eigenvalues
! 101    continue
    fb=func(bx,ierr)
    if (ierr < 0 .OR. ihltcmd /= 0) return
    if(fb > fa)then
        bnd=smin
        idir=-1
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
    endif

    cx=bx+gold*(bx-ax)

    if ( idir*(cx-bnd) > 0.)            cx=bnd

    fc=func(cx,ierr)
    if (ierr < 0 .OR. ihltcmd /= 0) return

    1 if(fb >= fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
        if ( idir*(u-bnd) > 0.)          u=bnd
        ulim=bx+glimit*(cx-bx)
        if ( idir*(ulim-bnd) > 0.)       ulim=bnd
        if((bx-u)*(u-cx) > 0.)then
        !                                    ** parabolic u between b and c
            fu=func(u,ierr)
            if (ierr < 0 .OR. ihltcmd /= 0) return
            if(fu < fc)then
            !                                      * minimum between b and c
                ax=bx
                fa=fb
                bx=u
                fb=fu
                ierr=0
                return
            else if(fu > fb)then
            !                                      * minimum between a and u
                cx=u
                fc=fu
                ierr=0
                return
            endif
        !                                      * default magnification
            if (cx == bnd) go to 99
            u=cx+gold*(cx-bx)
            if ( idir*(u-bnd) > 0.)       u=bnd
            fu=func(u,ierr)
            if (ierr < 0 .OR. ihltcmd /= 0) return
        else if((cx-u)*(u-ulim) > 0.)then
        !                                    ** parabolic u between c and ulim
            if (cx == bnd) go to 99
            fu=func(u,ierr)
            if (ierr < 0 .OR. ihltcmd /= 0) return
            if(fu < fc)then
                bx=cx
                cx=u
                u=cx+gold*(cx-bx)
                if ( idir*(u-bnd) > 0.)    u=bnd
                fb=fc
                fc=fu
                fu=func(u,ierr)
                if (ierr < 0 .OR. ihltcmd /= 0) return
            endif
        else if((u-ulim)*(ulim-cx) >= 0.)then
            if (cx == bnd) go to 99
            u=ulim
            fu=func(u,ierr)
            if (ierr < 0 .OR. ihltcmd /= 0) return
        else
            if (cx == bnd) go to 99
            u=cx+gold*(cx-bx)
            if ( idir*(u-bnd) > 0.)       u=bnd
            fu=func(u,ierr)
            if (ierr < 0 .OR. ihltcmd /= 0) return
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        go to 1
    endif
    ierr=0
    return
    99 ierr=1
    return
    end subroutine mnbrak
