! NLSPMC Version 1.0 2/5/99
    subroutine lmpar(n,r,ldr,ipvt,diag,qtf,delta,par,x,sdiag,wa1, &
    wa2,gnvec,gradf)
    integer :: n,ldr
    integer :: ipvt(n)
    real*8 :: delta,par
    real*8 :: r(ldr,n),diag(n),qtf(n),x(n),sdiag(n),wa1(n), &
    wa2(n),gnvec(n), gradf(n)
!     **********

!     Subroutine lmpar

!     Given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta,
!     the problem is to determine a value for the parameter
!     par such that if x solves the system

!         A*x = f ,     sqrt(par)*D*x = 0 ,

!     in the least squares sense, and dxnorm is the Euclidean
!     norm of D*x, then either par is zero and

!           (dxnorm-delta) .le. 0.1*delta ,

!     or par is positive and

!           abs(dxnorm-delta) .le. 0.1*delta .

!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization, with column pivoting, of A. That is, if
!     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!     columns, and R is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then lmpar expects
!     the full upper triangle of R, the permutation matrix P,
!     and the first n components of (Q transpose)*b. On output
!     lmpar also provides an upper triangular matrix S such that

!            t   t                   t
!           P *(A *A + par*D*D)*P = S *S .

!     S is employed within lmpar and may be of separate interest.

!     Only a few iterations are generally needed for convergence
!     of the algorithm. If, however, the limit of 10 iterations
!     is reached, then the output par will contain the best
!     value obtained so far.

!     The subroutine statement is

!       subroutine lmpar(n,r,ldr,ipvt,diag,qtf,delta,par,x,sdiag,
!                        wa1,wa2)

!     where

!       n is a positive integer input variable set to the order of R.

!       R is an n by n array. On input the full upper triangle
!         must contain the full upper triangle of the matrix R.
!         On output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix S.

!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.

!       ipvt is an integer input array of length n which defines the
!         permutation matrix P such that A*P = Q*R. Column j of P
!         is column ipvt(j) of the identity matrix.

!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix D.

!       qtf is an input array of length n which must contain the first
!         n elements of the vector (Q transpose)*f (f is the vector
!         containing the function values to be minimized).

!       delta is a positive input variable which specifies an upper
!         bound on the Euclidean norm of D*x.

!       par is a nonnegative variable. On input par contains an
!         initial estimate of the Levenberg-Marquardt parameter.
!         On output par contains the final estimate.

!       x is an output array of length n which contains the least
!         squares solution of the system A*x = b, sqrt(par)*D*x = 0,
!         for the output par.

!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix S.

!       gnvec is an output array of length n which contains the
!         Gauss-Newton vector corresponding to the input matrix R.

!       gradf is an output array of length n which contains the
!         gradient (steepest descent) vector

!       wa1 and wa2 are work arrays of length n.

!     Subprograms called

!       MINPACK-supplied ... dpmpar,enorm,qrsolv

!       Fortran-supplied ... dabs,dmax1,dmin1,dsqrt

!     Argonne National Laboratory. MINPACK project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

!     **********
    integer :: i,iter,j,jm1,jp1,k,l,nsing
    real*8 :: dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001, &
    sum,temp,zero
    real*8 :: dpmpar,enorm
    data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/

!     dwarf is the smallest positive magnitude.

    dwarf = dpmpar(2)

!----------------------------------------------------------------------
!     Compute and store in x the Gauss-Newton direction. If the
!     Jacobian is rank-deficient, obtain a least squares solution.
!     Store the scaled Gauss-Newton direction in gnvec.

!     The Gauss-Newton direction is the solution of the problem given
!     above for par=0. The solution to the above least-squares problem
!     in this case is equivalent to finding the vector

!           +        -1 T         +
!     x = -J f  = - R  Q f (i.e. J  is the pseudoinverse of the Jacobian)

!     See exercise 24 in chapter 6 of Dennis and Schnabel.
!----------------------------------------------------------------------
    nsing = n
    do 10 j = 1, n
        wa1(j) = qtf(j)
        if (r(j,j) == zero .AND. nsing == n) nsing = j - 1
        if (nsing < n) wa1(j) = zero
    10 END DO

    if (nsing >= 1) then
        do 40 k = 1, nsing
            j = nsing - k + 1
            wa1(j) = wa1(j)/r(j,j)
            temp = wa1(j)
            jm1 = j - 1
            if (jm1 >= 1) then
                do 20 i = 1, jm1
                    wa1(i) = wa1(i) - r(i,j)*temp
                20 END DO
            end if
        40 END DO
    end if

    do 60 j = 1, n
        l = ipvt(j)
        x(l) = wa1(j)
        gnvec(l)=x(l)
    60 END DO

!----------------------------------------------------------------------
!     Initialize the iteration counter. If the (scaled) Gauss-Newton
!     step size is within 10% of the trust-region bound, terminate the
!     search (par will be returned as zero).
!----------------------------------------------------------------------
    iter = 0
    do 70 j = 1, n
        wa2(j) = diag(j)*x(j)
    70 END DO
    dxnorm = enorm(n,wa2)
    fp = dxnorm - delta
    if (fp <= p1*delta) go to 220

!----------------------------------------------------------------------
!     If the Jacobian is not rank deficient, the Gauss-Newton
!     step provides a lower bound, parl, for the zero of
!     the function. Otherwise set this bound to zero.
!----------------------------------------------------------------------
    parl = zero
    if (nsing >= n) then
        do 80 j = 1, n
            l = ipvt(j)
            wa1(j) = diag(l)*(wa2(l)/dxnorm)
        80 END DO
        do 110 j = 1, n
            sum = zero
            jm1 = j - 1
            if (jm1 >= 1) then
                do 90 i = 1, jm1
                    sum = sum + r(i,j)*wa1(i)
                90 END DO
            end if
            wa1(j) = (wa1(j) - sum)/r(j,j)
        110 END DO
        temp = enorm(n,wa1)
        parl = ((fp/delta)/temp)/temp
    end if

!----------------------------------------------------------------------
!     Calculate the steepest descent direction of the residuals at
!     point where fvec has been evaluated.

!                T       T T
!     Grad(f) = J f = P R Q f

!     This provides an upper bound, paru, for the zero of the function.
!----------------------------------------------------------------------
    do 140 j = 1, n
        sum = zero
        do 130 i = 1, j
            sum = sum + r(i,j)*qtf(i)
        130 END DO
        l = ipvt(j)
        wa1(j) = sum/diag(l)
        gradf(j) = wa1(j)
    140 END DO
    gnorm = enorm(n,wa1)
    paru = gnorm/delta
    if (paru == zero) paru = dwarf/dmin1(delta,p1)

!----------------------------------------------------------------------
!     If the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.
!----------------------------------------------------------------------
    par = dmax1(par,parl)
    par = dmin1(par,paru)
    if (par == zero) par = gnorm/dxnorm

!----------------------------------------------------------------------
!     *** Beginning of an iteration.
!----------------------------------------------------------------------
    150 continue
    iter = iter + 1

!----------------------------------------------------------------------
!        Evaluate the function at the current value of par.
!----------------------------------------------------------------------
    if (par == zero) par = dmax1(dwarf,p001*paru)
    temp = dsqrt(par)
    do 160 j = 1, n
        wa1(j) = temp*diag(j)
    160 END DO
    call qrsolv(n,r,ldr,ipvt,wa1,qtf,x,sdiag,wa2)
    do 170 j = 1, n
        wa2(j) = diag(j)*x(j)
    170 END DO
    dxnorm = enorm(n,wa2)
    temp = fp
    fp = dxnorm - delta

!----------------------------------------------------------------------
!        If the function is small enough, accept the current value
!        of par. Also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.
!----------------------------------------------------------------------
    if (dabs(fp) <= p1*delta &
     .OR. parl == zero .AND. fp <= temp &
     .AND. temp < zero .OR. iter == 10) go to 220

!----------------------------------------------------------------------
!        Compute the Newton correction.
!----------------------------------------------------------------------
    do 180 j = 1, n
        l = ipvt(j)
        wa1(j) = diag(l)*(wa2(l)/dxnorm)
    180 END DO
    do 210 j = 1, n
        wa1(j) = wa1(j)/sdiag(j)
        temp = wa1(j)
        jp1 = j + 1
        if (n >= jp1) then
            do 190 i = jp1, n
                wa1(i) = wa1(i) - r(i,j)*temp
            190 END DO
        end if
    210 END DO
    temp = enorm(n,wa1)
    parc = ((fp/delta)/temp)/temp

!----------------------------------------------------------------------
!        Depending on the sign of the function, update parl or paru.
!----------------------------------------------------------------------
    if (fp > zero) parl = dmax1(parl,par)
    if (fp < zero) paru = dmin1(paru,par)

!----------------------------------------------------------------------
!        Compute an improved estimate for par.
!----------------------------------------------------------------------
    par = dmax1(parl,par+parc)

!----------------------------------------------------------------------
!        *** End of an iteration.
!----------------------------------------------------------------------
    go to 150
    220 continue

!----------------------------------------------------------------------
!     Termination.
!----------------------------------------------------------------------
    if (iter == 0) par = zero

    return
    end subroutine lmpar
