!----------------------------------------------------------------------
!                    =========================
!                       include file LMTXT
!                    =========================
!  Contains definitions of strings identifying return codes from
!  the subroutine LMDER in MINPACK
!----------------------------------------------------------------------
    character(32) :: minerr(0:10)

    data minerr / 'Bad input parameters', &
    'Chi-squared convergence', &
    'X vector convergence', &
    'Chi-squared/X vector convergence', &
    'Gradient convergence', &
    'Max function evaluations reached', &
    'Max iterations reached', &
    'FTOL too small', &
    'XTOL too small' , &
    'Zero Jacobian or GTOL too small', &
    'Terminated internally'/
