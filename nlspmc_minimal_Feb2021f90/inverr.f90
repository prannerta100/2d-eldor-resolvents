!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       =================
!                       SUBROUTINE INVERR
!                       =================

!     Subroutine INVERR computes error estimates for computed isolated
!     good T-eigenvalues in vs and writes these eigenvalues and
!     estimates to file.  By definition a good T-eigenvalue is
!     isolated if its closest neighbor is also good, or if its closest
!     neighbor is spurious but that neighbor is far enough away.  So
!     in particular, we will compute estimates for any good
!     T-eigenvalue that is in a cluster of good T-eigenvalues.

!     uses inverse iteration on T(1,mev) solving the equation
!     (T - x*I)v2 = right-hand side for each such good T-eigenvalue x.

!     program refactors T-x*I on each iteration of inverse iteration.
!     typically only one iteration is needed per T-eigenvalue x.

!     on entry and exit:

!        mev = order of T
!        n = order of original matrix A
!        vs = computed distinct eigenvalues of T(1,mev)
!        mp = T-multiplicity of each T-eigenvalue in vs.
!             mp(i) = -1 means vs(i) is a good T-eigenvalue but that
!             it is sitting close to a spurious T-eigenvalue.
!             mp(i) = 0 means vs(i) is spurious.  Estimates are
!             computed only for those T-eigenvalues with mp(i) = 1.
!             Flagging was done in subroutine isoev prior to entering
!             INVERR.
!        niso = number of isolated good T-eigenvalues contained in vs
!        ndis =  number of distinct T-eigenvalues in vs

!     In program:

!        iter = maximum number of inverse iteration steps allowed for
!             each x.  iter = it on entry.
!        gr,gc = arrays of dimension at least mev+niso.  used to
!             store randomly-generated right-hand side.  This is not
!             regenerated for each x. g is also used to store error
!             estimates as they are computed for later printout.
!        v1,v2 = work spaces used in the factorization of T(1,mev).
!             At the end of the inverse iteration computation for x,
!             v2 contains the unit eigenvector of T(1,mev) corresponding
!             to x.  v1 and v2 must be of dimension at least mev.

!     on exit:

!        gg(j) = minimum gap in T(1,mev) for each vs(j), j=1,ndis
!        g(i) = |betam|*|v2(mev)| = error estimate for isolated good
!             T-eigenvalues, where i = 1,niso  and  betam = beta(mev)
!             T(1,mev) corresponding to ith isolated good T-eigenvalue.
!             If for some x it.gt.iter then the error estimate in g
!             is marked with a - sign.
!        v2 = isolated good T-eigenvalues
!        v1 = minimal T-gaps for the T-eigenvalues in v2.
!             These are constructed for write-out purposes only and not
!             needed elsewhere in the program.

!     includes :

!        nlsdim.inc
!        eprprm.inc

!*********************************************************************

    subroutine inverr(a,b,v1,v2,vs,eps,gr,gc,g,gg,mp, &
    intc,mev,ndis,niso,it)

    implicit none

    include 'limits.inc'
    include 'simparm.inc'

    integer ::  i,ii,iso,igood,ispur,it,iter
    integer ::  j,jev,mev,mp1,mm1,ndis,ng,niso
    integer ::  mp(mxstep), intc(mxstep)

    real :: g(mxstep),gg(mxstep)

    double precision :: est,estr,estc,sum,xu,norm,tsum,gsum,gap
    double precision ::  eps,eps3,eps4,zero,one,gr(mxstep),gc(mxstep)

    complex*16  a(mxstep),b(mxstep),v1(mxstep),v2(mxstep)
    complex*16  u,z,x,ratio,betam,temp,czero,csum,znormu,vs(mxstep)

    parameter (one=1.0D0,zero=0.0D0,czero=(0.0D0,0.0D0))

    external znormu

!#######################################################################

!     initialization and parameter specification

    ng=0
    niso=0
    iter=it
    mp1=mev+1
    mm1=mev-1
    betam=b(mev)
    b(mev)=czero

!     calculate scale and tolerances
    tsum=zero
    do 30 i=1,mev
        tsum=tsum+cdabs(a(i))+cdabs(b(i))
    30 END DO

    eps3=eps*tsum
    eps4=dble(mev)*eps3

!     generate scaled right-hand side
    gsum=zero
    do 60 i=1,mev
        gsum=gsum+dabs(gr(i))+dabs(gc(i))
    60 END DO
    gsum=eps4/gsum

    do 70 i=1,mev
        gr(i)=gsum*gr(i)
        gc(i)=gsum*gc(i)
    70 END DO

!     loop on isolated good T-eigenvalues in vs (mp(i)=1) to
!     calculate corresponding unit eigenvector of T(1,mev)

    do 200 jev=1,ndis
        if (mp(jev) == 0) go to 200
        ng=ng+1
        if (mp(jev) /= 1) go to 200
        it=1
        niso=niso+1
        x=vs(jev)
    
    !     initialize right hand side for inverse iteration
    !     and the flag on which rows are interchanged
        do 80 i=1,mev
            intc(i)=0
            v2(i)=dcmplx(gr(i),gc(i))
        80 END DO
    
    !-----------------------------------------------------------------------
    !      triangular factorization
    !-----------------------------------------------------------------------
    
        90 continue
        u=a(1)-x
        z=b(1)
    
        do 110 i=2,mev
            if (cdabs(b(i-1)) <= cdabs(u)) then
                v1(i-1)=z/u
                v2(i-1)=v2(i-1)/u
                v2(i)=v2(i)-b(i-1)*v2(i-1)
                ratio=b(i-1)/u
                u=a(i)-x-z*ratio
                z=b(i)
            else
                ratio=u/b(i-1)
                intc(i)=1
                v1(i-1)=a(i)-x
                u=z-ratio*v1(i-1)
                z=-ratio*b(i)
                temp=v2(i-1)
                v2(i-1)=v2(i)
                v2(i)=temp-ratio*v2(i)
            end if
        110 END DO
    
        if (cdabs(u) == zero) u=dcmplx(eps3,eps3)
    
    !     back substitution
        v2(mev)=v2(mev)/u
        do 130 ii=1,mm1
            i=mev-ii
            if (intc(i+1) /= 1) then
                v2(i)=v2(i)-v1(i)*v2(i+1)
            else
                v2(i)=(v2(i)-v1(i)*v2(i+1)-b(i+1)*v2(i+2))/b(i)
            end if
        130 END DO
    
    !-----------------------------------------------------------------------
    !      tests for convergence of inverse iteration
    !-----------------------------------------------------------------------
    
        norm=cdabs(v2(mev))
        do 140 ii=1,mm1
            i=mev-ii
            norm=norm+cdabs(v2(i))
        140 END DO
    
        if (norm >= one) go to 160
    
        it=it+1
        if (it > iter) go to 160
    
        xu=eps4/norm
    
        do 150 i=1,mev
            intc(i)=0
            v2(i)=v2(i)*xu
        150 END DO
    
        go to 90
    !     another inverse iteration step
    
    !     inverse iteration finished
    !     normalize computed T-eigenvector : v2=v2/||v2||
        160 continue
    
        csum=znormu(v2,mev)
    
        do 170 ii=1,mev
            v2(ii)=v2(ii)/csum
        170 END DO
    
    !     save error estimate for later output
        est=cdabs(betam)*cdabs(v2(mev))
        estr=dabs(dreal(v2(mev)))
        estc=dabs(dimag(v2(mev)))
        gsum=cdabs(betam)
        if (it > iter) est=-est
        g(niso)=est
    
    200 END DO

!     end error estimate loop on isolated good T-eigenvalues.
!     Generate distinct mingaps for T(1,mev).  This is useful as an
!     indicator of the goodness of the inverse iteration estimates.
!     Transfer isolated good T-eigenvalues and corresponding tmingaps
!     to v2 and v1 for output purposes only.

    iso=0
    do 210 j=1,ndis
        if (mp(j) /= 1) go to 210
        iso=iso+1
        gr(iso)=gg(j)
        v2(iso)=vs(j)
    210 END DO
    if(niso == 0) go to 270

    ispur=0
    i=0
    do 260 j=1,ndis
        if(mp(j) /= 0) go to 240
        ispur=ispur+1
        go to 260
        240 if(mp(j) /= 1) go to 260
        i=i+1
        igood=j-ispur
    260 END DO

!     restore b(mev+1)=betam
    270 b(mev)=betam

    return

    end subroutine inverr
