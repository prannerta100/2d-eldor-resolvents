!	VERSION	1.0	2/5/99
!***********************************************************************

!	define the largest value of the	L quantum number allowed
!	in the basis set

!	Notes:
!		This file is included in the program lbll.f for	the
!		purpose	of verifying that the parameters input by the
!		user are consistent with the arrays dimensioned	in
!		matrix element calculation subroutine, matrll.f, which
!		also includes this file.

!	written	by DJS 11-SEP-87

!***********************************************************************

    integer	mxlval
    parameter (mxlval=150)

