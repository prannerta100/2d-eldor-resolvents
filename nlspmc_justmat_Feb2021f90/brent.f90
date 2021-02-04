!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                         =======================
!                            function BRENT
!                         =======================

!   Modified from the routine in NUMERICAL RECIPES for the non-linear
!   least squares iteration.  Notice that it uses the absolute tolerance
!   instead of the fractional tolerence.

!   Error flag of negative value means the function "func" failed.

!   Given a function func, and given a bracketing triplet of abscissas ax,
!   bx, cx (such that bx is between ax and cx, and func(bx) is less than
!   both func(ax) and func(cx), this routine isolates the minimum to
!   a ABSOLUTE precision of about tol using Brent's method.

!   Includes :
!       stdio.inc, iterat.f
!----------------------------------------------------------------------

    function brent(ax,bx,cx,func,tol,xmin,prmID,ierr)

    implicit none

    include 'limits.inc'
    include 'stdio.inc'
    include 'parms.inc'

    character prmID*6
    double precision :: brent,ax,bx,cx,func,tol,xmin,cgold
    double precision :: a,b,d,e,etemp,fu,fv,fw,fx,r,q,p,tol1,tol2,u,v,w, &
    x,xm
    integer :: iterr,itmax,ierr
    parameter (itmax=100,cgold=.3819660d0)

    external func

    tol1=tol
    tol2=2.*tol1

    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.
    fx=func(x,ierr)
    write(*,*)'brent1 call func, x= ',x,fx
    if (ierr < 0 .OR. ihltcmd /= 0) return
!      seteval=.true.	! got eigenvalues - not used anymore
    fv=fx
    fw=fx
    do 11 iterr=1,itmax
        xm=0.5*(a+b)
    
        if (iterr == 1) then
            write (luout,1000) prmID
            if (luout /= luttyo) write (luttyo,1000) prmID
        end if
    
        write (luout,1010) iterr,fx,x
        if (luout /= luttyo) write (luttyo,1010) iterr,fx,x
    
        if(abs(x-xm) <= (tol2-.5d0*(b-a))) goto 3
        if(abs(e) > tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.*(q-r)
            if(q > 0.) p=-p
            q=abs(q)
            etemp=e
            e=d
            if(abs(p) >= abs(.5*q*etemp) .OR. p <= q*(a-x) .OR. &
            p >= q*(b-x)) goto 1
            d=p/q
            u=x+d
            if(u-a < tol2 .OR. b-u < tol2) d=sign(tol1,xm-x)
            goto 2
        endif
        1 if(x >= xm) then
            e=a-x
        else
            e=b-x
        endif
        d=cgold*e
        2 if(abs(d) >= tol1) then
            u=x+d
        else
            u=x+sign(tol1,d)
        endif
        fu=func(u,ierr)
        write(*,*)'brent2 call func, x= ',x,fu
        if (ierr < 0 .OR. ihltcmd /= 0) return
        if(fu <= fx) then
            if(u >= x) then
                a=x
            else
                b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
        else
            if(u < x) then
                a=u
            else
                b=u
            endif
            if(fu <= fw .OR. w == x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if(fu <= fv .OR. v == x .OR. v == w) then
                v=u
                fv=fu
            endif
        endif
    11 END DO
    pause 'Brent exceed maximum iterations.'
    3 xmin=x
    brent=fx
    write (luout,1020)
    if (luout /= luttyo) write (luttyo,1020)
    return

! ### format statements ##############################################

    1000 format(/,28('-'),/'Iter',5x,'RmsDv',6x,a,/28('-'))
    1010 format(i4,4x,g10.4,2x,g24.18)
    1020 format(28('-')/)

    end function brent
