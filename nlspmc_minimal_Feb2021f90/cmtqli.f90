!  VERSION 1.0  (NLSPMC version)   2/5/99
!**********************************************************************

!       This subroutine will diagonalize a complex symmetric
!       tridiagonal matrix using the QL algorithm with implicit
!       shifts.  It also accumulates the first row of the eigenvector
!       matrix which is required for calculating ESR spectra from
!       the eigenvalues of the stochastic Liouville matrix.  The
!       algorithm used here is a slight modification of the algorithm
!       discussed CMTQL1 in "Lanczos Algorithms for Large Symmetric
!       Eigenvalue Computations", vol. 2, J. Cullum and R. Willoughby,
!       Birkhauser, 1985.

!       written by DJS 26-NOV-87

!       Includes:
!               nlsdim.inc
!               rndoff.inc

!       Uses:

!**********************************************************************

    subroutine cmtqli(n,d,e,z,isort,ierr)

    include 'limits.inc'
    include 'rndoff.inc'

    integer :: n,ierr
    complex*16 d,e,z
    dimension d(mxstep),e(mxstep),z(mxstep)

    integer :: i,j,l,m,mml,isort
    double precision :: temp,t0,t1,eps
    complex*16 b,c,f,g,p,r,s,w

    complex*16 czero,cone
    parameter (czero=(0.0D0,0.0D0),cone=(1.0D0,0.0D0))

!######################################################################

    eps=100.0D0*rndoff

!----------------------------------------------------------------------
!       intialize first row of eigenvector matrix
!----------------------------------------------------------------------

    z(1)=cone
    do 10 i=2,n
        z(i)=czero
    10 END DO

    ierr=0
    if (n == 1) go to 180

!======================================================================
!               loop over eigenvalues
!======================================================================

    e(n)=czero

    do 140 l=1,n
    
        j=0
    
    !----------------------------------------------------------------------
    !               find a small subdiagonal matrix element
    !----------------------------------------------------------------------
    
        20 continue
    
        do 30 m=l,n-1
            temp=cdabs(d(m))+cdabs(d(m+1))
            if (cdabs(e(m)) <= temp*rndoff) go to 40
        30 END DO
    
        m=n
    
        40 p=d(l)
        if (m /= l) then
            if (j == 100) go to 170
            j=j+1
        
        !----------------------------------------------------------------------
        !               form shift
        !----------------------------------------------------------------------
        
            g=(d(l+1)-p)*0.5D0
            t0=cdabs(g)
            t1=cdabs(e(l))
        
            if (t0 <= t1) then
                w=g/e(l)
                r=cdsqrt(cone+w*w)
                t0=cdabs(w+r)
                t1=cdabs(w-r)
                if (t1 <= t0) then
                    g=d(m)-p+e(l)/(w+r)
                else
                    g=d(m)-p+e(l)/(w-r)
                end if
            else
                w=e(l)/g
                r=cdsqrt(cone+w*w)
                t0=cdabs(cone+r)
                t1=cdabs(cone-r)
                if (t1 <= t0) then
                    g=d(m)-p+w*e(l)/(cone+r)

                else
                    g=d(m)-p+w*e(l)/(cone-r)
                end if
            end if
        
            s=cone
            c=-cone
            p=czero
            mml=m-l
        
            do 90 i=m-1,l,-1
                f=s*e(i)
                b=-c*e(i)
                t0=cdabs(g)
                t1=cdabs(f)
                if (t1 <= t0) then
                    w=f/g
                    r=cdsqrt(cone+w*w)
                    e(i+1)=g*r
                    c=cone/r
                    s=w*c
                else
                    w=g/f
                    r=cdsqrt(cone+w*w)
                    e(i+1)=f*r
                    s=cone/r
                    c=w*s
                end if
                temp=1.0D0+cdabs(w)**2
                t0=dsqrt(temp)
                t1=cdabs(r)
            
                if (t1 <= eps*t0) then
                    ierr=-l
                    go to 180
                else
                    ierr=0
                end if
            
            !----------------------------------------------------------------------
            !               finish up loop over i
            !----------------------------------------------------------------------
            
                g=d(i+1)-p
                r=(d(i)-g)*s+2.0D0*c*b
                p=s*r
                d(i+1)=g+p
                g=b-c*r
            
            !----------------------------------------------------------------------
            !               keep track of effect of rotations on first row
            !----------------------------------------------------------------------
            
                w=z(i+1)
                z(i+1)=s*z(i)+c*w
                if (cdabs(z(i+1)) < rndoff) z(i+1)=czero
                z(i)=c*z(i)-s*w
                if (cdabs(z(i)) < rndoff) z(i)=czero
            
            !----------------------------------------------------------------------
            !               end of loop to restore to tridiagonal form
            !----------------------------------------------------------------------
            
            90 END DO
        
        !----------------------------------------------------------------------
        !               go back to find another small off-diagonal element
        !----------------------------------------------------------------------
        
            d(l)=d(l)-p
            e(l)=g
            e(m)=czero
        
            go to 20
        end if
    
    !----------------------------------------------------------------------
    !               end of loop over eigenvalues
    !----------------------------------------------------------------------
    
    140 END DO

    do 145 l=1,n
        z(l)=z(l)*z(l)
    145 END DO

!----------------------------------------------------------------------
!      sort in the order of increasing magnitude of the eigenvalues
!      or decreasing order of weighting factor
!----------------------------------------------------------------------

!                           *** sort by eigenvalue ***
    if (isort == 0) then
        do 150 l=1,n-1
            t0=cdabs(d(l))
            do 150 m=l+1,n
                t1=cdabs(d(m))
                if (t1 < t0) then
                    t0=t1
                    w=d(m)
                    d(m)=d(l)
                    d(l)=w
                    w=z(m)
                    z(m)=z(l)
                    z(l)=w
                end if
        150 END DO
    else
    !                           *** sort by weighting factor ***
        do 160 l=1,n-1
            t0=cdabs(z(l))
            do 160 m=l+1,n
                t1=cdabs(z(m))
                if (t1 > t0) then
                    t0=t1
                    w=d(m)
                    d(m)=d(l)
                    d(l)=w
                    w=z(m)
                    z(m)=z(l)
                    z(l)=w
                end if
        160 END DO
    end if

    go to 180

    170 ierr=l

!----------------------------------------------------------------------
!       return to caller
!----------------------------------------------------------------------

    180 return
    end subroutine cmtqli
