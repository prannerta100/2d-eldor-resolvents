    SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
    BETA, C, LDC )
!     .. Scalar Arguments ..
    CHARACTER(1) ::        TRANSA, TRANSB
    INTEGER ::            M, N, K, LDA, LDB, LDC
    COMPLEX*16         ALPHA, BETA
!     .. Array Arguments ..
    COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - COMPLEX*16      .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
    LOGICAL ::            LSAME
    EXTERNAL           LSAME
!     .. External Subroutines ..
    EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
    INTRINSIC          DCONJG, MAX
!     .. Local Scalars ..
    LOGICAL ::            CONJA, CONJB, NOTA, NOTB
    INTEGER ::            I, INFO, J, L, NCOLA, NROWA, NROWB
    COMPLEX*16         TEMP
!     .. Parameters ..
    COMPLEX*16         ONE
    PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
    COMPLEX*16         ZERO
    PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!     B  respectively are to be  transposed but  not conjugated  and set
!     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!     and the number of rows of  B  respectively.
!
    NOTA  = LSAME( TRANSA, 'N' )
    NOTB  = LSAME( TRANSB, 'N' )
    CONJA = LSAME( TRANSA, 'C' )
    CONJB = LSAME( TRANSB, 'C' )
    IF( NOTA )THEN
        NROWA = M
        NCOLA = K
    ELSE
        NROWA = K
        NCOLA = M
    END IF
    IF( NOTB )THEN
        NROWB = K
    ELSE
        NROWB = N
    END IF
!
!     Test the input parameters.
!
    INFO = 0
    IF(      ( .NOT. NOTA                 ) .AND.  &
    ( .NOT. CONJA                ) .AND.  &
    ( .NOT. LSAME( TRANSA, 'T' ) )      )THEN
        INFO = 1
    ELSE IF( ( .NOT. NOTB                 ) .AND.  &
        ( .NOT. CONJB                ) .AND.  &
        ( .NOT. LSAME( TRANSB, 'T' ) )      )THEN
        INFO = 2
    ELSE IF( M  < 0               )THEN
        INFO = 3
    ELSE IF( N  < 0               )THEN
        INFO = 4
    ELSE IF( K  < 0               )THEN
        INFO = 5
    ELSE IF( LDA < MAX( 1, NROWA ) )THEN
        INFO = 8
    ELSE IF( LDB < MAX( 1, NROWB ) )THEN
        INFO = 10
    ELSE IF( LDC < MAX( 1, M     ) )THEN
        INFO = 13
    END IF
    IF( INFO /= 0 )THEN
        CALL XERBLA( 'ZGEMM ', INFO )
        RETURN
    END IF
!
!     Quick return if possible.
!
    IF( ( M == 0 ) .OR. ( N == 0 ) .OR.  &
    ( ( ( ALPHA == ZERO ) .OR. ( K == 0 ) ) .AND. ( BETA == ONE ) ) ) &
    RETURN
!
!     And when  alpha.eq.zero.
!
    IF( ALPHA == ZERO )THEN
        IF( BETA == ZERO )THEN
            DO 20, J = 1, N
                DO 10, I = 1, M
                    C( I, J ) = ZERO
                    CONTINUE
                10 END DO
                CONTINUE
            20 END DO
        ELSE
            DO 40, J = 1, N
                DO 30, I = 1, M
                    C( I, J ) = BETA*C( I, J )
                    CONTINUE
                30 END DO
                CONTINUE
            40 END DO
        END IF
        RETURN
    END IF
!
!     Start the operations.
!
    IF( NOTB )THEN
        IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B + beta*C.
        !
            DO 90, J = 1, N
                IF( BETA == ZERO )THEN
                    DO 50, I = 1, M
                        C( I, J ) = ZERO
                        CONTINUE
                    50 END DO
                ELSE IF( BETA /= ONE )THEN
                    DO 60, I = 1, M
                        C( I, J ) = BETA*C( I, J )
                        CONTINUE
                    60 END DO
                END IF
                DO 80, L = 1, K
                    IF( B( L, J ) /= ZERO )THEN
                        TEMP = ALPHA*B( L, J )
                        DO 70, I = 1, M
                            C( I, J ) = C( I, J ) + TEMP*A( I, L )
                            CONTINUE
                        70 END DO
                    END IF
                    CONTINUE
                80 END DO
                CONTINUE
            90 END DO
        ELSE IF( CONJA )THEN
        !
        !           Form  C := alpha*conjg( A' )*B + beta*C.
        !
            DO 120, J = 1, N
                DO 110, I = 1, M
                    TEMP = ZERO
                    DO 100, L = 1, K
                        TEMP = TEMP + DCONJG( A( L, I ) )*B( L, J )
                        CONTINUE
                    100 END DO
                    IF( BETA == ZERO )THEN
                        C( I, J ) = ALPHA*TEMP
                    ELSE
                        C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                    END IF
                    CONTINUE
                110 END DO
                CONTINUE
            120 END DO
        ELSE
        !
        !           Form  C := alpha*A'*B + beta*C
        !
            DO 150, J = 1, N
                DO 140, I = 1, M
                    TEMP = ZERO
                    DO 130, L = 1, K
                        TEMP = TEMP + A( L, I )*B( L, J )
                        CONTINUE
                    130 END DO
                    IF( BETA == ZERO )THEN
                        C( I, J ) = ALPHA*TEMP
                    ELSE
                        C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                    END IF
                    CONTINUE
                140 END DO
                CONTINUE
            150 END DO
        END IF
    ELSE IF( NOTA )THEN
        IF( CONJB )THEN
        !
        !           Form  C := alpha*A*conjg( B' ) + beta*C.
        !
            DO 200, J = 1, N
                IF( BETA == ZERO )THEN
                    DO 160, I = 1, M
                        C( I, J ) = ZERO
                        CONTINUE
                    160 END DO
                ELSE IF( BETA /= ONE )THEN
                    DO 170, I = 1, M
                        C( I, J ) = BETA*C( I, J )
                        CONTINUE
                    170 END DO
                END IF
                DO 190, L = 1, K
                    IF( B( J, L ) /= ZERO )THEN
                        TEMP = ALPHA*DCONJG( B( J, L ) )
                        DO 180, I = 1, M
                            C( I, J ) = C( I, J ) + TEMP*A( I, L )
                            CONTINUE
                        180 END DO
                    END IF
                    CONTINUE
                190 END DO
                CONTINUE
            200 END DO
        ELSE
        !
        !           Form  C := alpha*A*B'          + beta*C
        !
            DO 250, J = 1, N
                IF( BETA == ZERO )THEN
                    DO 210, I = 1, M
                        C( I, J ) = ZERO
                        CONTINUE
                    210 END DO
                ELSE IF( BETA /= ONE )THEN
                    DO 220, I = 1, M
                        C( I, J ) = BETA*C( I, J )
                        CONTINUE
                    220 END DO
                END IF
                DO 240, L = 1, K
                    IF( B( J, L ) /= ZERO )THEN
                        TEMP = ALPHA*B( J, L )
                        DO 230, I = 1, M
                            C( I, J ) = C( I, J ) + TEMP*A( I, L )
                            CONTINUE
                        230 END DO
                    END IF
                    CONTINUE
                240 END DO
                CONTINUE
            250 END DO
        END IF
    ELSE IF( CONJA )THEN
        IF( CONJB )THEN
        !
        !           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
        !
            DO 280, J = 1, N
                DO 270, I = 1, M
                    TEMP = ZERO
                    DO 260, L = 1, K
                        TEMP = TEMP + &
                        DCONJG( A( L, I ) )*DCONJG( B( J, L ) )
                        CONTINUE
                    260 END DO
                    IF( BETA == ZERO )THEN
                        C( I, J ) = ALPHA*TEMP
                    ELSE
                        C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                    END IF
                    CONTINUE
                270 END DO
                CONTINUE
            280 END DO
        ELSE
        !
        !           Form  C := alpha*conjg( A' )*B' + beta*C
        !
            DO 310, J = 1, N
                DO 300, I = 1, M
                    TEMP = ZERO
                    DO 290, L = 1, K
                        TEMP = TEMP + DCONJG( A( L, I ) )*B( J, L )
                        CONTINUE
                    290 END DO
                    IF( BETA == ZERO )THEN
                        C( I, J ) = ALPHA*TEMP
                    ELSE
                        C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                    END IF
                    CONTINUE
                300 END DO
                CONTINUE
            310 END DO
        END IF
    ELSE
        IF( CONJB )THEN
        !
        !           Form  C := alpha*A'*conjg( B' ) + beta*C
        !
            DO 340, J = 1, N
                DO 330, I = 1, M
                    TEMP = ZERO
                    DO 320, L = 1, K
                        TEMP = TEMP + A( L, I )*DCONJG( B( J, L ) )
                        CONTINUE
                    320 END DO
                    IF( BETA == ZERO )THEN
                        C( I, J ) = ALPHA*TEMP
                    ELSE
                        C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                    END IF
                    CONTINUE
                330 END DO
                CONTINUE
            340 END DO
        ELSE
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
            DO 370, J = 1, N
                DO 360, I = 1, M
                    TEMP = ZERO
                    DO 350, L = 1, K
                        TEMP = TEMP + A( L, I )*B( J, L )
                        CONTINUE
                    350 END DO
                    IF( BETA == ZERO )THEN
                        C( I, J ) = ALPHA*TEMP
                    ELSE
                        C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                    END IF
                    CONTINUE
                360 END DO
                CONTINUE
            370 END DO
        END IF
    END IF
!
    RETURN
!
!     End of ZGEMM .
!
    END
