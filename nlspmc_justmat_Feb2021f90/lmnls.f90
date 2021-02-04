! NLSPMC Version 1.0 2/5/99
!----------------------------------------------------------------------
!                    =========================
!                      subroutine LMNLS
!                    =========================
!     (L)evenberg-(M)arquardt (N)onlinear (L)east (S)quares

! This is a modification of the original lmder subroutine from the
! MINPACK subroutine library. It uses a Levenberg-Marquardt nonlinear
! least-squares algorithm modified to carry out a local optimization
! constrained to lie within a "trust region" defined by a step bound
! delta using scaling of the variables.

! For a description of the trust region approach for least squares
! problems, see J.E. Dennis and R.B. Schnabel, Numerical Methods for
! Unconstrained Optimization and Nonlinear Equations, Prentice-Hall,
! Englewood Cliffs, NJ (1983), sections 6.4, 7.1, and 10.2.

! Modified to use limits on the variables as set in the calling
! program.  x when called is the physical real parameter, convert it
! to an infinite range parameter for real range defined by prmin,
! prmax.  Convert back to physical parameters on call to simulation
! program.

!----------------------------------------------------------------------
    subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, &
    maxfev,maxitr,diag,scale,factor,nprint,itrace, &
    jacobi,info,nfev,njev,ipvt,qtf,gnvec,gradf, &
    wa1,wa2,wa3,wa4)
    integer :: m,n,ldfjac,maxfev,maxitr,nprint,info,istep,nfev,njev, &
    itrace,jacobi
    integer :: ipvt(n)
    double precision :: ftol,xtol,gtol,factor
    double precision :: x(n),fvec(m),fjac(ldfjac,n),diag(n),scale(n), &
    qtf(n),gnvec(n),gradf(n),wa1(n),wa2(n),wa3(n), &
    wa4(m)
    external fcn

!----------------------------------------------------------------------

!     The purpose of LMDER is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the Levenberg-Marquardt algorithm. The user must provide a
!     subroutine which calculates the functions and the Jacobian.

!     The subroutine statement is

!       subroutine lmnls(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                        maxfev,maxitr,diag,scale,factor,nprint,info,
!                        nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)

!     where

!       FCN is the name of the user-supplied subroutine which
!         calculates the functions and the Jacobian. FCN must
!         be declared in an external statement in the user
!         calling program, and should be written as follows:

!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         double precision x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         If iflag=1 calculate the functions at x and
!         return this vector in fvec. Do not alter fjac.
!         If iflag=2 calculate the Jacobian at x and
!         return this matrix in fjac. Do not alter fvec.
!         ----------
!         return
!         end

!         The value of IFLAG should not be changed by FCN unless
!         the user wants to terminate execution of LMDER.
!         In this case set iflag to a negative integer.

!       M is a positive integer input variable set to the number
!         of functions.

!       N  is a positive integer input variable set to the number
!          of variables. N must not exceed M.

!       X is an array of length N. On input X must contain
!         an initial estimate of the solution vector. On output X
!         contains the final estimate of the solution vector.

!       FVEC is an output array of length M which contains
!         the functions evaluated at the output X.

!       FJAC is an output M by N array. the upper N by N submatrix
!         of FJAC contains an upper triangular matrix R with
!         diagonal elements of nonincreasing magnitude such that

!                T     T           T
!               P *(JAC *JAC)*P = R *R,

!         where P is a permutation matrix and JAC is the final
!         calculated Jacobian. column j of P is column IPVT(j)
!         (see below) of the identity matrix. The lower trapezoidal
!         part of FJAC contains information generated during
!         the computation of R.

!       LDFJAC is a positive integer input variable not less than M
!         which specifies the leading dimension of the array FJAC.

!       FTOL is a nonnegative input variable. Termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most FTOL.
!         Therefore, FTOL measures the relative error desired
!         in the sum of squares.

!       XTOL is a nonnegative input variable. Termination
!         occurs when the relative error between two consecutive
!         iterates is at most XTOL. Therefore, XTOL measures the
!         relative error desired in the approximate solution.

!       GTOL is a nonnegative input variable. Termination
!         occurs when the cosine of the angle between FVEC and
!         any column of the Jacobian is at most GTOL in absolute
!         value. therefore, GTOL measures the orthogonality
!         desired between the function vector and the columns
!         of the Jacobian.

!       MAXFEV is a positive integer input variable. Termination
!         occurs when the number of calls to FCN with IFLAG=1
!         has reached MAXFEV.

!       SCALE is an array of length N containing multiplicative scale
!         factors for each of the variables in X. If an element of SCALE
!         is non-positive, it will be reset internally to unity.
!         Positive entries in the SCALE array will be retained as
!         user-specified scaling factors for the trust-region search
!         of the algorithm. The step for the Ith parameter will be scaled
!         by SCALE(I) times the norm of the Ith column of the Jacobian.
!         The default value for all parameters is unity (i.e., the
!         column norms of the Jacobian will be used).

!         NB: This convention differs from the original
!         specifications of LMDER in MINPACK.

!       FACTOR is a positive input variable used in determining the
!         initial trust region bound. This bound is set to the product of
!         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else
!         to FACTOR itself. In most cases FACTOR should lie in the
!         interval (.1,100.). 100 is a generally recommended value.

!       NPRINT is an integer input variable that enables controlled
!         printing of iterates if it is positive. In this case,
!         fcn is called with IFLAG=0 at the beginning of the first
!         iteration and every NPRINT iterations thereafter and
!         immediately prior to return, with X, FVEC, and FJAC
!         available for printing. FVEC and FJAC should not be
!         altered. If NPRINT is not positive, no special calls
!         of FCN with IFLAG=0 are made.

!       INFO is an integer output variable. If the user has
!         terminated execution, INFFO is set to the (negative)
!         value of IFLAG. See description of FCN. Otherwise,
!         INFO is set as follows:

!         INFO=0  Improper input parameters.

!         INFO=1  Both actual and predicted relative reductions
!                   in the sum of squares are at most FTOL.

!         INFO=2  Relative error between two consecutive iterates
!                   is at most XTOL.

!         INFO=3  conditions for INFO=1 and INFO=2 both hold.

!         INFO=4  The cosine of the angle between FVEC and any
!                   column of the Jacobian is at most GTOL in
!                   absolute value.

!         INFO=5  number of calls to FCN with IFLAG=1 has
!                   reached MAXFEV.

!         INFO=6  FTOL is too small. No further reduction in
!                   the sum of squares is possible.

!         INFO=7  XTOL is too small. No further improvement in
!                   the approximate solution X is possible.

!         INFO=8  GTOL is too small. FVEC is orthogonal to the
!                   columns of the Jacobian to machine precision.

!       NFEV is an integer output variable set to the number of
!         calls to FCN with IFLAG=1.

!       NJEV is an integer output variable set to the number of
!         calls to NJEV with IFLAG=2.

!       IPVT is an integer output array of length N. IPVT
!         defines a permutation matrix p such that JAC*P=Q*R,
!         where JAC is the final calculated Jacobian, Q is
!         orthogonal (not stored), and R is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of P is column IPVT(j) of the identity matrix.

!       QTF is an output array of length N which contains
!         the first N elements of the vector (Q transpose)*FVEC.

!       WA1, WA2, and WA3 are work arrays of length N.

!       WA4 is a work array of length M.

!     Subprograms called

!       User-supplied ...... FCN

!       MINPACK-supplied ... DPMPAR,ENORM,LMPAR,QRFAC

!       FORTRAN-supplied ... DABS,DMAX1,DMIN1,DSQRT,MOD

!     Argonne National Laboratory. MINPACK Project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

!----------------------------------------------------------------------
    integer :: i,iflag,j,l,ld

    double precision :: actred,delta,dirder,epsmch,fnorm1,gnorm, &
    one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio, &
    sum,temp,temp1,temp2,xnorm,zero
    double precision :: dpmpar,enorm
    data one,p1,p5,p25,p75,p0001,zero &
    /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/

! ----added for EPR NLS
    include 'limits.inc'
    include 'parms.inc'
!      include 'iterat.inc'

!     epsmch is the machine precision.

    epsmch=dpmpar(1)

! on entry, convert physical xxx into infinite range lmnls variables.

    call mapxxx(x,n,-1)
    info=0
    iflag=0
    nfev=0
    njev=0

!----------------------------------------------------------------------
!     Check the input parameters for errors.
!----------------------------------------------------------------------
    if (n <= 0 .OR. m < n .OR. ldfjac < m &
     .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero &
     .OR. maxfev <= 0 .OR. maxitr <= 0 .OR. factor <= zero) &
    go to 300

!**********************************************************************

!----------------------------------------------------------------------
!     Evaluate the function at the starting point
!     and calculate its norm.
!----------------------------------------------------------------------
    iflag=1
    call mapxxx(x,n,1)	! map to physical variables
    call fcn(m,n,x,fvec,fjac,ldfjac,iflag)	! pfun
    if (ihltcmd /= 0) return
    call mapxxx(x,n,-1)	! map to lmnls variables
    nfev=1
!      seteval=.true.	! got eigenvalues. - not used anymore.
    if (iflag < 0) go to 300
    fnorm=enorm(m,fvec)

!----------------------------------------------------------------------
!     Initialize Levenberg-Marquardt parameter and iteration counter
!----------------------------------------------------------------------
    par=zero
    iter=1

!----------------------------------------------------------------------
! ********* Beginning of the outer loop *******************************
!----------------------------------------------------------------------

    30 continue
!----------------------------------------------------------------------
!        Calculate the Jacobian matrix (iflag=2)
!        and signal fcn to output it if necessary (iflag=3)
!----------------------------------------------------------------------
    iflag=2
    if (jacobi /= 0) iflag=3
    call mapxxx(x,n,1)	! map to physical variables
    call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    if (ihltcmd /= 0) return
    call mapxxx(x,n,-1)	! map to lmnls variables

!       call ftest1(fjac(1,1),16384,'in lmnls, fjac(1)',17)
!       call ftest1(fjac(1,2),16384,'in lmnls, fjac(2)',17)
    njev=njev+1
    nfev=nfev+n
    if (iflag < 0) go to 300
!----------------------------------------------------------------------
!        If requested, call fcn to enable printing of iterates
!----------------------------------------------------------------------
    if (nprint > 0) then
        iflag=0
        if (mod(iter-1,nprint) == 0)then
            call mapxxx(x,n,1)	! map to physical variables
            call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            if (ihltcmd /= 0) return
            call mapxxx(x,n,-1)	! map to lmnls variables
        end if
        if (iflag < 0) go to 300
    end if
!----------------------------------------------------------------------
!        Compute the QR factorization of the Jacobian
!----------------------------------------------------------------------
    call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
! ***** wa2 here is large?
!----------------------------------------------------------------------
!        On the first iteration, set each non-positive element of the
!        SCALE scaling array according to the norms of the columns of
!        the initial Jacobian
!----------------------------------------------------------------------
    if (iter == 1) then
        do 50 j=1, n
            if (scale(j) <= zero) then
                diag(j)=wa2(j)
            else
                diag(j)=wa2(j)/scale(j)
            end if
            if (diag(j) == zero) diag(j)=one
        50 END DO
    
        if (itrace /= 0) then
            write(itrace,1000) (tag(j),j=1,n)
            write(itrace,1001) (wa2(j),j=1,n)
            write(itrace,1002) (diag(j),j=1,n)
        end if
    
    !----------------------------------------------------------------------
    !        On the first iteration, calculate the norm of the scaled x
    !        and initialize the step bound delta
    !----------------------------------------------------------------------
        do 70 j=1, n
            wa3(j)=diag(j)*x(j)
        70 END DO
        xnorm=enorm(n,wa3)
        delta=factor*xnorm
        if (delta == zero) delta=factor
        if (itrace /= 0) write(itrace,1003) xnorm,delta,factor
    end if

!----------------------------------------------------------------------
!        Form (Q transpose)*fvec and store the first n components in
!        QtF.
!----------------------------------------------------------------------
    do 90 i=1, m
        wa4(i)=fvec(i)
    90 END DO
    do 130 j=1, n
    
        if (fjac(j,j) /= zero) then
            sum=zero
            do 100 i=j, m
                sum=sum + fjac(i,j)*wa4(i)
            100 END DO
            temp=-sum/fjac(j,j)
            do 110 i=j, m
                wa4(i)=wa4(i) + fjac(i,j)*temp
            110 END DO
        end if
    
        fjac(j,j)=wa1(j)
        qtf(j)=wa4(j)
    130 END DO

!----------------------------------------------------------------------
!        Compute the norm of the scaled gradient.
!----------------------------------------------------------------------
    gnorm=zero
    if (fnorm /= zero) then
        do 160 j=1, n
            l=ipvt(j)
            if (wa2(l) /= zero) then
                sum=zero
                do 140 i=1, j
                    sum=sum + fjac(i,j)*(qtf(i)/fnorm)
                140 END DO
                gnorm=dmax1(gnorm,dabs(sum/wa2(l)))
            end if
        160 END DO
    end if

!----------------------------------------------------------------------
!        Test for convergence of the gradient norm
!----------------------------------------------------------------------
    if (gnorm <= gtol) info=4
    if (info /= 0) go to 300

!----------------------------------------------------------------------
!        Rescale diag array
!----------------------------------------------------------------------
    do 180 j=1,n
        if (scale(j) <= zero) then
            temp=wa2(j)
        else
            temp=wa2(j)/scale(j)
        end if
        diag(j)=dmax1(diag(j),temp)
    180 END DO


!----------------------------------------------------------------------
!  ******** Beginning of the inner loop ******************************
!----------------------------------------------------------------------
    istep=1
    200 continue
!----------------------------------------------------------------------
!           Determine the Levenberg-Marquardt parameter.
!----------------------------------------------------------------------
    call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2, &
    wa3,wa4,gnvec,gradf)
!----------------------------------------------------------------------
!           Store the direction p and X + p. Calculate the norm of p.
!----------------------------------------------------------------------

    do 210 j=1, n
        wa1(j)=-wa1(j)
        wa2(j)=x(j) + wa1(j)
        wa3(j)=diag(j)*wa1(j)
    210 END DO
    pnorm=enorm(n,wa3)

    if (itrace /= 0 .AND. istep == 1) then
        write(itrace,1004) iter,fnorm,(tag(j),j=1,n)
        write(itrace,1005) (x(j),j=1,n)
        write(itrace,1006) (diag(j),j=1,n)
        if (par /= zero) write(itrace,1007) (gradf(j),j=1,n)
        write(itrace,1008) (-gnvec(j),j=1,n)
    end if

!----------------------------------------------------------------------
!        On the first iteration, adjust the initial trust region bound
!        to the size of the initial step.
!----------------------------------------------------------------------
    if (iter == 1) then
        if (itrace /= 0) then
            write(itrace,1010) pnorm
            if (delta > pnorm) write (itrace,1011)
        end if
        delta=dmin1(delta,pnorm)
    end if

!----------------------------------------------------------------------
!           Evaluate the function at x + p and calculate its norm.
!----------------------------------------------------------------------
    iflag=1
    call mapxxx(wa2,n,1)	! map to physical variables
    call fcn(m,n,wa2,wa4,fjac,ldfjac,iflag)
    if (ihltcmd /= 0) return
    call mapxxx(wa2,n,-1)	! map to lmnls variables
    nfev=nfev+1
    if (iflag < 0) go to 300
    fnorm1=enorm(m,wa4)

    if (itrace /= 0) write (itrace,1012) istep,par, &
    delta,fnorm1*fnorm1
    istep=istep+1

!----------------------------------------------------------------------
!           Compute the scaled actual reduction.
!----------------------------------------------------------------------
    actred=-one
    if (p1*fnorm1 < fnorm) actred=one - (fnorm1/fnorm)**2

!----------------------------------------------------------------------
!           Compute the scaled predicted reduction and
!           the scaled directional derivative.
!----------------------------------------------------------------------
    do 230 j=1,n
        wa3(j)=zero
        l=ipvt(j)
        temp=wa1(l)
        do 220 i=1, j
            wa3(i)=wa3(i) + fjac(i,j)*temp
        220 END DO
    230 END DO
    temp1=enorm(n,wa3)/fnorm
    temp2=(dsqrt(par)*pnorm)/fnorm
    prered=temp1**2 + temp2**2/p5
    dirder=-(temp1**2 + temp2**2)

!----------------------------------------------------------------------
!           Compute the ratio of the actual to the predicted
!           reduction.
!----------------------------------------------------------------------
    ratio=zero
    if (prered /= zero) ratio=actred/prered

!----------------------------------------------------------------------
!           Update the step bound.
!----------------------------------------------------------------------
    if (ratio <= p25) then
    
    !----------------------------------------------------------------------
    !             If actual reduction is too much smaller than the predicted
    !             reduction (i.e. actred/prered ratio is too small)
    !             the function is not well-approximated by a quadratic
    !             equation. Reduce the size of the trust region by a
    !             factor of 0.1 to 0.5 and increase the L-M parameter.
    !----------------------------------------------------------------------
    
        if (actred >= zero) temp=p5
        if (actred < zero) temp=p5*dirder/(dirder + p5*actred)
        if (p1*fnorm1 >= fnorm .OR. temp < p1) temp=p1
        if (p1*delta > pnorm) then
            if (itrace /= 0) write(itrace,1013) one/p1,pnorm/p1
            delta=pnorm/p1
        endif
        delta=delta*temp
        par=par/temp
    
    else
    
    !----------------------------------------------------------------------
    !             If ratio of actual to predicted reduction is close to 1,
    !             the quadratic model is a good approximation to the function,
    !             and we can try increasing the trust region by a factor of
    !             two to see if a better solution is available. Otherwise,
    !             the size of the step bound is left unchanged.
    !----------------------------------------------------------------------
    
        if (par == zero .OR. ratio >= p75) then
            delta=pnorm/p5
            par=p5*par
            temp=one/p5
        else
            temp=one
        end if
    end if

    if (itrace /= 0) write (itrace,1014) ratio,temp
    if (itrace /= 0) write (itrace,1015) (wa2(j)-x(j),j=1,n)

!----------------------------------------------------------------------
!           Test for successful iteration.
!----------------------------------------------------------------------
    if (ratio >= p0001) then
    !----------------------------------------------------------------------
    !           Successful iteration. Update X, FVEC, and their norms.
    !----------------------------------------------------------------------
        do 270 j=1, n
            x(j)=wa2(j)
            wa2(j)=diag(j)*x(j)
        270 END DO
        do 280 i=1, m
            fvec(i)=wa4(i)
        280 END DO
        xnorm=enorm(n,wa2)
        fnorm=fnorm1
        iter=iter+1
    end if
!----------------------------------------------------------------------
!           Tests for convergence.
!----------------------------------------------------------------------
    info=0
    if (dabs(actred) <= ftol .AND. prered <= ftol &
     .AND. p5*ratio <= one) info=1
    if (delta <= xtol*xnorm) info=info + 2
    if (info /= 0) go to 300
!----------------------------------------------------------------------
!           Tests for termination and stringent tolerances.
!----------------------------------------------------------------------
    if (nfev >= maxfev) info=5
    if (iter >= maxitr) info=6
    if (dabs(actred) <= epsmch .AND. prered <= epsmch &
     .AND. p5*ratio <= one) info=7
    if (delta <= epsmch*xnorm) info=8
    if (gnorm <= epsmch) info=9
    if (info /= 0) go to 300
!----------------------------------------------------------------------
!           End of the inner loop. Repeat if iteration unsuccessful.
!----------------------------------------------------------------------
    if (ratio < p0001) go to 200
!----------------------------------------------------------------------
!        End of the outer loop.
!----------------------------------------------------------------------
    go to 30
    300 continue
!----------------------------------------------------------------------
!     Termination, either normal or user-imposed.
!----------------------------------------------------------------------
    if (iflag < 0) then
        info=10
        nprint=0
    end if
    iflag=0
    call mapxxx(x,n,1)	! map to physical variables
    if (nprint > 0) call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    return

! ##### format statements for trace printout ########################

     
    1000 format(65('#')/10x,'TRACE OF LEVENBERG-MARQUARDT MINIMIZATION', &
    /65('#')// &
    'INITIAL VALUES:'/16x,10(3x,a6,3x))
    1001 format('Col norms of J:',10(2x,g10.4)/)
    1002 format(9x,'Scale:',10(2x,g10.4))
    1003 format(/10x,'Scaled X norm: ',g11.5/5x,'Trust region bound:', &
    g11.5,'  =(',g9.3,'*Xnorm)')
    1004 format(/65('#')/'Iteration',i3,': Chi-Sq=',g11.5//12x, &
    10(3x,a6,3x))
    1005 format(7x,'Xvec:',10(2x,g10.4))
    1006 format(6x,'Scale:',10(2x,g10.4))
    1007 format(3x,'Gradient:',10(2x,g10.4))
    1008 format(1x,'G-N vector:',10(2x,g10.4))
    1010 format(/3x,'Initial step size=',g10.4)
    1011 format(3x,'(TR bound has been reduced to the initial ', &
    'step size)'/)
    1012 format(5x,65('-')/5x,'Step',i3,'; LMpar=',g10.4,' TR bound=', &
    g10.4,'; Chi-Sq=',g12.5)
    1013 format(5x,'(TR bound reduced to ',f4.1,'*step length =',g11.5,')')
    1014 format(5x,'Actual/predicted Chi-Sq reduction=',g10.4, &
    '; TR scaled by ',f4.1/5x,65('-'))
    1015 format(12x,10(2x,g10.4))
    2000 format(/6x,'*** Jacobian calculation : ',i3,' ***'/6x, &
    'trans(J)*fvec',15x,'trans(J)*J'/)
    2001 format(6x,g12.5,5x,6(g12.5,1x))
    end subroutine lmnls
