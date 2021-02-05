!  VERSION 1.0  (NLSPMC version)   2/5/99
!**********************************************************************

!                         SUBROUTINE : CSCG
!                         -----------------

!       (C)omplex (S)ymmetric (C)onjugate (G)radient Algorithm

!       This subroutine attempts a solution of an Ax=b type problem
!       where A is a known sparse complex symmetric matrix, b is a
!       known vector and the approximate solution vector is obtained
!       via a simple complex symmetric conjugate gradient procedure.
!       The algorithm used here is based on algorithm 10.3-1 in
!       "Matrix Computations", first edition, by G. Golub and C.
!       Van Loan, Johns Hopkins Univ. Press, 1983.  The "Lanczos"
!       tridiagonal matrix in constructed using formula 10.2-14 in
!       the reference above.

!       It is assumed in this routine that the right hand vector is
!       normalized and the solution vector is initialized at the start
!       of the routine.  The matrix and index arrays are passed to
!       this routine through the inclusion of the file eprmat.inc.

!       Modified by SHL :  Lanczos vector is stored in array 'lnzvec'.
!       They are used in the conversion of the tridiagonal eigenvector
!	into the eigenvector in the original ritz basis.  The role of
!	the old cgltri.f routine (calculation of the Lanczos tridiagonal
!	from the quantities which arise naturally in the course of the
!	CG algorithm) is also incorporated in this routine to resolve
!	the sign ambiguity that comes in when one takes sqrt(complex #).
!	This happens twice here.
!	  i)  converting beta(CG) into beta(LA)
!	  ii) normalizing the residual vector
!	As long as the sign is consistent between these, the resulting
!	eigenvectors and eigenvalues are the same, even though the sign
!	of eigenvectors and beta's of the tridiagonal matrix may be
! 	different from that of LA.

!       arguments
!       ---------

!             b     : right hand vector
!             ndim  : number of rows and columns in the matrix
!             mxcgs : maximum number of conjugate gradient steps allowed
!             cgtol : maximum residual allowed for termination
!             shift : a complex shift value used to avoid extraneous
!                     divisions by zero
!             ter,r,p,w : work spaces
!      returns
!      -------
!             x     : approximate solution vector
!             lnzvec: Lanczos vectors that defines the transformation
!                     from the original to the tridiagonal matrix
!             al    : diagonal elements of tridiagonal matrix
!             bl    : off-diagonal elements of tridiagonal matrix
!             ter   : b-A*x
!             ndone : number of cg steps executed (positive if
!                     converged, negative otherwise)
!             terror: true error of the solution vector estimated
!                     Termination of the CG steps is based on this
!                     true error.

!       written by DJS 20-SEP-87
!	modified by SHL 22-JUL-92

!       Includes:
!               nlsdim.inc
!               parcom.inc
!               stdio.inc
!               rndoff.inc
!       Uses:
!               zdotu2.f
!               zaypx.f
!               scmvm.f
!               zaxpy2.f

!**********************************************************************

    subroutine cscg(b,ndim,mxcgs,cgtol,shift, &
    x,lnzvec,al,bl,ter,r,p,w,ndone,terror,trmin,ntrmin)

    include 'limits.inc'
    include 'parms.inc'
    include 'stdio.inc'
    include 'rndoff.inc'
    include 'eprmat.inc'

    integer :: ndim,mxcgs,ndone,ntrmin
    double precision :: cgtol,error,terror,trmin
!      real(kind=8) r8_normal_ab,mu,sigma
    integer(kind=4) seed
    complex*16 alpha,alph1,alph2,beta,x,b,lnzvec,al,bl
    dimension x(mxdim),b(mxdim),lnzvec(mxdim,mxstep), &
    al(mxstep),bl(mxstep)

    integer :: i,nstpcg,k,l
    double precision :: zero,one,amp,phase,tr,ti,trhalf,tihalf,norm,flg
    complex*16 rho1,rho2,rho1h,rho2h,czero,ci,shift,bubba
    parameter (zero=0.0D0,one=1.0D0,czero=(0.0D0,0.0D0), &
    ci=(0.0D0,1.0D0))

    complex*16 ter,p,r,w
    dimension ter(mxdim),p(mxdim),r(mxdim),w(mxdim)

    complex*16 zdotu2
    external zdotu2


!######################################################################

!      write(luttyo,*) ndim,mxcgs,cgtol,shift
    if (idebug /= 0) write (ludeb,1010)

!----------------------------------------------------------------------
!    Initialize r,x and p
!----------------------------------------------------------------------


!     Print zdiag(1,:),zdiag(2,:)
!      write(luttyo,*) "In cscg.f, about to start the CSCG calculation,
!     # printing zdiag(1:2,:) now"
!      do 69 i=1,ndim
!        write(luttyo,*) zdiag(1,i),zdiag(2,i)
! 9    continue

    if (idebug /= 0) then
        do 9890 i=1,ndim
            write(ludeb,*) b(i)
        9890 END DO
    end if

    do 10 i=1,ndim
        r(i)=b(i)
        p(i)=b(i)
    10 END DO

    rho1=zdotu2(r,r,ndim)

    alph1=one
    rho2=one
    rho2h=one


!======================================================================
!     begin loop over CG steps
!======================================================================

    nstpcg=0
    100 nstpcg=nstpcg+1

!----------------------------------------------------------------------
!     calculate beta
!----------------------------------------------------------------------

    if (nstpcg == 1) then
        beta=czero
    else
        rho2=rho1
        rho2h=rho1h
        rho1=zdotu2(r,r,ndim)
        beta=rho1/rho2
    
    end if

!----------------------------------------------------------------------
!   calculate the pseudo norm of r (sqrt(rho1))
!----------------------------------------------------------------------

    tr=dreal(rho1)
    ti=dimag(rho1)
    amp=sqrt(sqrt(tr*tr+ti*ti))
    if (amp > rndoff) then
        phase=0.5D0*datan2(ti,tr)
        trhalf=amp*cos(phase)
        tihalf=amp*sin(phase)
        if (abs(trhalf) < rndoff) trhalf=zero
        if (abs(tihalf) < rndoff) tihalf=zero
    else
        trhalf=zero
        tihalf=zero
    end if
    rho1h=trhalf+ci*tihalf

!----------------------------------------------------------------------
!     store the residual vectors to get the Lanczos vector
!----------------------------------------------------------------------

    do 20 i=1,ndim
        lnzvec(i,nstpcg)=r(i)/rho1h
    20 END DO

!----------------------------------------------------------------------
!     update p  ( p <- r+beta*p )
!----------------------------------------------------------------------

    if (nstpcg /= 1) call zaypx(r,p,beta,ndim)

!----------------------------------------------------------------------
!     calculate w ( w <- A*p )
!----------------------------------------------------------------------

    call scmvm(p,w,ndim)

!----------------------------------------------------------------------
!     add diagonal shift term to matrix-vector product
!----------------------------------------------------------------------

    call zaxpy2(p,w,shift,ndim)

!----------------------------------------------------------------------
!     calculate alpha and convert into LA tridiagonal form
!----------------------------------------------------------------------

    alph2=alph1
    alph1=rho1/zdotu2(p,w,ndim)

    al(nstpcg)=one/alph1+beta/alph2-shift

    bl(nstpcg)=-beta*rho2h/(alph2*rho1h)

!----------------------------------------------------------------------
!     update x ( x <- x+alpha*p )
!----------------------------------------------------------------------

    call zaxpy2(p,x,alph1,ndim)

!----------------------------------------------------------------------
!     update r ( r <- r-alpha*w )
!----------------------------------------------------------------------

    call zaxpy2(w,r,-alph1,ndim)

!----------------------------------------------------------------------
!     calculate error and terror & check for convergence
!----------------------------------------------------------------------

    call scmvm(x,ter,ndim)
    call zaxpy2(x,ter,shift,ndim)

    terror=zero
    do 30 i=1,ndim
        ter(i)=b(i)-ter(i)
        tr=dreal(ter(i))
        ti=dimag(ter(i))
        terror=terror+tr*tr+ti*ti
    30 END DO

    terror=dsqrt(terror)

    error=dsqrt(abs(rho1))

!      write(luttyo,*) 'terror=',terror,', error=',error
    if (error < 1.D-17) go to 200
    if (terror < cgtol) go to 200
! record minimum so far found
    if (terror < trmin)then
        trmin=terror
        ntrmin=nstpcg
    end if

    if (idebug /= 0) then
        if ( (nstpcg/10)*10 == nstpcg ) &
        write (ludeb,1020) nstpcg,error,terror
    end if

    if (nstpcg < mxcgs) go to 100

!======================================================================
!     end of loop over CG steps
!======================================================================

    200 continue	! exit loop

    write (ludeb,1020) nstpcg,error,terror

!   shift the index of beta by 1 to match convention

    do 40 i=1,nstpcg-1
        bl(i)=bl(i+1)
    40 END DO
    bl(nstpcg)=czero

!----------------------------------------------------------------------
!     return error code if not converged
!----------------------------------------------------------------------

    if (terror <= cgtol) then
    !        ndone=nstpcg-1
        ndone=nstpcg   ! test this version, re neg. real parts.
    ! this choice provides more stable operaition and maybe no neg reals!
    else
    !        ndone=-nstpcg+1
        ndone=-nstpcg
    end if



!            write(luttyo,*) "LNZVEC^T LNZVEC"
!            norm=0.0d0
!            do 3654 k=1,nstpcg
!                do 3655 l=1,nstpcg
!                    bubba=czero
!                    do 3656 i=1,ndim
!                        bubba=bubba+lnzvec(i,k)*lnzvec(i,l);
! 3656               continue
!                    write(luttyo,*) k,l,bubba
!                  flg=0.0q0
!                    if(k.eq.l) flg=1.0d0
!                    norm=norm+abs(bubba-flg)**2
! 3655           continue
! 3654       continue
!           write(luttyo,*) "L^t L norm square = ",norm

!      mu=0.0d0
!      sigma=1.0d0
!      seed=12345
!      write(luttyo,*) "randn(mu=0,sigma=1)",r8_normal_ab(mu,sigma,seed)
!      write(luttyo,*) "randn(mu=0,sigma=1)"
!      write(luttyo,*) "Lanczos Vector 1"
!      do 3660 k=1,ndim
! 660     write(luttyo,*) lnzvec(k,1)


!----------------------------------------------------------------------
!     return to calling routine
!----------------------------------------------------------------------


    return



!=====================================================================
!     format statements
!=====================================================================

    1010 format(' ** CG calculation **'/,3x,'step',5x,'error',3x, &
    'true error')
    1020 format(2x,i4,3x,g12.5,2x,g12.5)

    end subroutine cscg
