! NLSPMC Version 1.0 2/5/99
    subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
    integer :: m,n,lda,lipvt
    integer :: ipvt(lipvt)
    logical :: pivot
    real*8 :: a(lda,n),rdiag(n),acnorm(n),wa(n)
!     **********

!     subroutine qrfac

!     This subroutine uses Householder transformations with column
!     pivoting (optional) to compute a QR factorization of the
!     m by n matrix A. That is, qrfac determines an orthogonal
!     matrix Q, a permutation matrix P, and an upper trapezoidal
!     matrix R with diagonal elements of nonincreasing magnitude,
!     such that A*P = Q*R. The Householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form

!                           T
!           i - (1/u(k))*u*u

!     where u has zeros in the first k-1 positions. The form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding LINPACK subroutine.

!     The subroutine statement is

!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

!     where

!       m is a positive integer input variable set to the number
!         of rows of a.

!       n is a positive integer input variable set to the number
!         of columns of a.

!       a is an m by n array. On input a contains the matrix for
!         which the qr factorization is to be computed. On output
!         the strict upper trapezoidal part of A contains the strict
!         upper trapezoidal part of R, and the lower trapezoidal
!         part of a contains a factored form of Q (the non-trivial
!         elements of the u vectors described above).

!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.

!       pivot is a logical input variable. If pivot is set true,
!         then column pivoting is enforced. If pivot is set false,
!         then no column pivoting is done.

!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix P such that A*P = Q*R.
!         Column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.

!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.

!       rdiag is an output array of length n which contains the
!         diagonal elements of r.

!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.

!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.

!     Subprograms called

!       MINPACK-supplied ... dpmpar,enorm

!       Fortran-supplied ... dmax1,dsqrt,min0

!     Argonne National Laboratory. MINPACK project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

!     **********
    integer :: i,j,jp1,k,kmax,minmn
    real*8 :: ajnorm,epsmch,one,p05,sum,temp,zero
    real*8 :: dpmpar,enorm
    data one,p05,zero /1.0d0,5.0d-2,0.0d0/

!     epsmch is the machine precision.

    epsmch = dpmpar(1)

!----------------------------------------------------------------------
!     Compute the initial column norms and initialize several arrays.
!----------------------------------------------------------------------
    do 10 j = 1, n
        acnorm(j) = enorm(m,a(1,j))
        rdiag(j) = acnorm(j)
        wa(j) = rdiag(j)
        if (pivot) ipvt(j) = j
    10 END DO

!----------------------------------------------------------------------
!     Reduce A to R with Householder transformations.
!----------------------------------------------------------------------
    minmn = min0(m,n)
    do 110 j = 1, minmn
        if (pivot) then
        
        !----------------------------------------------------------------------
        !        Bring the column of largest norm into the pivot position.
        !----------------------------------------------------------------------
            kmax = j
            do 20 k = j, n
                if (rdiag(k) > rdiag(kmax)) kmax = k
            20 END DO
        
            if (j /= kmax) then
                do 30 i = 1, m
                    temp = a(i,j)
                    a(i,j) = a(i,kmax)
                    a(i,kmax) = temp
                30 END DO
            
                rdiag(kmax) = rdiag(j)
                wa(kmax) = wa(j)
                k = ipvt(j)
                ipvt(j) = ipvt(kmax)
                ipvt(kmax) = k
            end if
        end if
    
    !----------------------------------------------------------------------
    !        Compute the Householder transformation to reduce the
    !        j-th column of a to a multiple of the j-th unit vector.
    !----------------------------------------------------------------------
        ajnorm = enorm(m-j+1,a(j,j))
    
        if (ajnorm /= zero) then
            if (a(j,j) < zero) ajnorm = -ajnorm
            do 50 i = j, m
                a(i,j) = a(i,j)/ajnorm
            50 END DO
        
            a(j,j) = a(j,j) + one
        
        !----------------------------------------------------------------------
        !        Apply the transformation to the remaining columns
        !        and update the norms.
        !----------------------------------------------------------------------
            jp1 = j + 1
            if (jp1 <= n) then
                do 90 k = jp1, n
                    sum = zero
                    do 60 i = j, m
                        sum = sum + a(i,j)*a(i,k)
                    60 END DO
                
                    temp = sum/a(j,j)
                    do 70 i = j, m
                        a(i,k) = a(i,k) - temp*a(i,j)
                    70 END DO
                
                    if (pivot .AND. rdiag(k) /= zero) then
                        temp = a(j,k)/rdiag(k)
                        rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
                        if (p05*(rdiag(k)/wa(k))**2 <= epsmch) then
                            rdiag(k) = enorm(m-j,a(jp1,k))
                            wa(k) = rdiag(k)
                        end if
                    end if
                90 END DO
            
            end if
        end if
    
        rdiag(j) = -ajnorm
    110 END DO

    return

!     Last card of subroutine qrfac.

    end subroutine qrfac
