! NLSPMC Version 1.0 2/5/99
!---------------------------------------------------------------------
!                      =========================
!                          subroutine ORDRPR
!                      =========================

! This routine calculates the molecular order parameters corresponding
! to a specified axial orienting potential. The potential is
! expressed as

!                                L   L
!    U(phi,theta,0) = Sum   Sum c   D   (phi,theta,0)
!                      L     K   K   0K

! where D^L_{0K} is a Wigner rotation matrix element (equivalent to
! the spherical harmonic Y^L_K), c^L_K are user-specified coefficients,
! and the sums are restricted to L=0 to 4 and K=0 to L (for even L,K only)

! The coefficients are specified in the cf array in the order
! c20, c22, c40, c42, c44. The routine calculates an order parameter for
! each nonzero coefficient specififed. The order parameter is the average
! of <D^L_{0K}(phi,theta)>  weighted by the value of exp(-U(phi,theta))
! and integrated over angles theta and phi.


! Includes:
!   rndoff.inc

! Uses:
!    ccrints.f  (Integration routines ccrint ccrin1)
!    ftht       (Function returning theta-dependent part of U)
!    fphi       (Function returning phi-dependent part of U)

!    The last two functions are included in this file.
!    Written by David E. Budil
!    Cornell Department of Chemistry, March 1992

!----------------------------------------------------------------------

    subroutine ordrpr(cf,d)
    implicit none
    double precision :: cf(5),d(5)

    double precision :: zero,one,small
    parameter(zero=0.0D0,one=1.0D0,small=1.0D-16)

    include 'rndoff.inc'

    integer :: i,npt,nup,id
    double precision :: fz,pfn

    integer :: kount
    double precision :: acc,c,d20,d40,bd22,bd42,bd44
    common/potprm/acc,c(5),kount
    common/func/d20,d40,bd22,bd42,bd44

    double precision :: ftht
    external ftht
!......................................................................
    acc=0.001
    npt=0

    do 10 kount=1,5
        c(kount)=cf(kount)
        if (dabs(c(kount)) > rndoff) npt=kount
    10 END DO

    if (npt > 0) then
    
    !---------------------------------------------
    ! Integrate unweighted potential function ftht
    !---------------------------------------------
        call ccrint(zero,one,acc,small,pfn,nup,ftht,id)
    
    !-----------------------------------------------
    ! Integrate potential weighted by D20, D22, etc.
    !-----------------------------------------------
        do 25 kount=1,npt
            call ccrint(zero,one,acc,small,fz,nup,ftht,id)
            d(kount)=fz/pfn
        25 END DO
    endif
    return
    end subroutine ordrpr

!----------------------------------------------------------------------
!                    =========================
!                          function FTHT
!                    =========================

!     Contains theta-dependence of orienting pseudopotential
!----------------------------------------------------------------------

    function ftht(ctht)
    implicit none
    double precision :: ftht,ctht

    integer :: nup,id
    double precision :: ctht2,stht2,result

    integer :: kount
    double precision :: acc,c,d20,d40,bd22,bd42,bd44
    common/potprm/acc,c(5),kount
    common/func/d20,d40,bd22,bd42,bd44

    double precision :: fphi
    external fphi

!---------------------------------------------------------
!  definition of some constants
!    a22 = sqrt(3/2), a42 = sqrt(5/8),  a44 = sqrt(35/8)
!---------------------------------------------------------
    double precision :: pi,a22,a42,a44,one,zero,small
    parameter (pi=3.14159265358979d0, &
    a22=1.22474487139159d0, &
    a42=0.790569415042095d0, &
    a44=1.04582503316759d0 )
    parameter(one=1.0D0,zero=0.0D0,small=1.0D-16 )

!......................................................................

    ctht2=ctht*ctht
    stht2=one-ctht2
    d20 =1.5d0*ctht2-0.5d0
    d40 =( (4.375d0*ctht2)-3.75d0)*ctht2+0.375d0
    bd22=a22*stht2
    bd42=a42*stht2*(7.0d0*ctht2-one)
    bd44=a44*stht2*stht2

    call ccrin1(zero,pi,acc,small,result,nup,fphi,id)
    ftht=result
    return
    end function ftht
                                   
!----------------------------------------------------------------------
!                    =========================
!                          function FPHI
!                    =========================

!     Contains phi-dependence of orienting pseudopotential
!----------------------------------------------------------------------
    function fphi(phi)
    implicit none
    double precision :: fphi,phi,c2phi,c4phi

    integer :: kount
    double precision :: acc,c,d20,d40,bd22,bd42,bd44
    common/potprm/acc,c(5),kount
    common/func/d20,d40,bd22,bd42,bd44

    double precision :: one,two
    parameter (one=1.0D0,two=2.0D0)

    c2phi=dcos(phi+phi)
    c4phi=two*c2phi*c2phi - one

    fphi=dexp(  c(1)*d20 &
    + c(2)*bd22*c2phi &
    + c(3)*d40 &
    + c(4)*bd42*c2phi &
    + c(5)*bd44*c4phi )

    if(kount == 0) return
    if(kount == 1) fphi=d20*fphi
    if(kount == 2) fphi=bd22*fphi*c2phi
    if(kount == 3) fphi=d40*fphi
    if(kount == 4) fphi=bd42*fphi*c2phi
    if(kount == 5) fphi=bd44*fphi*c4phi
    return
    end function fphi
