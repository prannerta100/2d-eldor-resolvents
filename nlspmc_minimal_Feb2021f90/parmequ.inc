!----------------------------------------------------------------------
!                    =========================
!                      include file PARMEQU
!                    =========================
!   Defines equivalence of arrays fepr and iepr with common block
!   /eprprm/ contained in sims.inc

!   fepr    is a vector corresponding to all the user-settable floating
!           point EPRLL parameters
!   iepr    is a vector corresponding to all the user-settable integer
!           parameters

!   NOTE: both nlsdim.inc and eprprm.inc must be included before
!   this file
!----------------------------------------------------------------------
    double precision :: fepr(NFPRM)
    integer :: iepr(NIPRM)
    equivalence (fepr(1),gib),(iepr(1),in2)
