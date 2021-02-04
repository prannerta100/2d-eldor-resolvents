! NLSPMC VERSION (VERSION 1.0) 2/5/99
!----------------------------------------------------------------------
!            =====================================================
!             subroutines TOCART, TOSPHR, and TOAXIL, and TCHECK
!            =====================================================

! These subroutines interconvert among the three possible representations
! of the various tensors (in 3-space, including magnetic and diffusion
! tensors) used in EPR calculations, namely:
!   (1) Cartesian representation
!   (2) Spherical representation
!   (3) Axial (actually also in Cartesian coordinates)

! The Cartesian representation is the familiar set of (x,y,z) components.
! If the Cartesian components of a tensor M are Mx, My, and Mz, then
! the "spherical" components M1, M2, and M3 are given by

!   M1=(Mx+My+Mz)/3   (isotropic component)
!   M2=Mz-(Mx+My)/2   (axial component)
!   M3=Mx-My          (rhombic component)

! Note that these components differ from the "true" second-order
! spherical tensor components, M(0,0), M(2,0), and M(2,2), respectively,
! by constant factors. They are expressed in this form so that they may
! be more directly correlated with experimental ESR spectra. These
! transformations are carried out by the TOSPHR routine.

! The third representation, "axial" is really in Cartesian coordinates
! as well, but affords the user the option of maintaining axial tensor
! symmetry while varying the Cartesian components.

! The transformation from spherical to cartesian components, the inverse
! of that given above, is carried out by TOCART and given by

!   Mx=M1 - M2/3 + M3/2
!   My=M1 - M2/3 - M3/2
!   Mz=M1 + 2*M2/3
!----------------------------------------------------------------------
    subroutine tocart( t, iflg )
    implicit none
    double precision :: two,three,t(3),t1,t2,t3
    integer :: x,y,z,iflg,CARTESIAN,SPHERICAL,AXIAL

    parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
    parameter (two=2.0D0,three=3.0D0,x=1,y=2,z=3)

    if (iflg == CARTESIAN) return
    if (iflg == AXIAL) then
        t1=t(x)
        t2=t(y)
        t3=t(z)
        t(x)=t1
        t(y)=t1
        t(z)=t3
    else if (iflg == SPHERICAL) then
        t1=t(x)
        t2=t(y)
        t3=t(z)
        t(x)=t1 - t2/three + t3/two
        t(y)=t1 - t2/three - t3/two
        t(z)=t1 + two*t2/three
    end if
    return
    end subroutine tocart


    subroutine tosphr( t, iflg )
    implicit none
    double precision :: two,three,t(3),tx,ty,tz,t1,t2,t3
    integer :: iflg,CARTESIAN,SPHERICAL,AXIAL

    parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
    parameter (two=2.0D0,three=3.0D0)

    if (iflg == SPHERICAL) return
    if (iflg == AXIAL) then
        tx=t(1)
        ty=t(1)
        tz=t(3)
    else if (iflg == CARTESIAN) then
        tx=t(1)
        ty=t(2)
        tz=t(3)
    end if
    t(1)=(tx+ty+tz)/three
    t(2)=tz-(tx+ty)/two
    t(3)=(tx-ty)
    return
    end subroutine tosphr


    subroutine toaxil( t, iflg )
    implicit none
    include 'stdio.inc'
    include 'rndoff.inc'
    double precision :: two,t(3),tmp
    integer :: x,y,z,iflg,CARTESIAN,SPHERICAL,AXIAL

    parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
    parameter (two=2.0D0)

    if (iflg == AXIAL) return
    if (iflg == SPHERICAL) call tocart( t,iflg )
    if (abs(t(1)-t(2)) > rndoff) write(luout,*)'caution ', &
    'information lost, tensor was not axial before conversion'
    tmp=t(2)
    t(1)=(t(1)+t(2))/two
    t(2)=0.0D0
    return
    end subroutine toaxil

!----------------------------------------------------------------------
!                     =========================
!                        function TCHECK
!                     =========================

!  Checks whether the component of the g, A, or R (T) tensors specified
!  by the ix index is consistent with the previous setting of the tensor
!  mode flags. (igflg, iaflg, and irflg) TCHECK returns .true. if
!  it is consistent (or if ix does not specify one of these tensors)
!  and .false. if an inconsistency is detected.  We require all sites
!  and spectra to have the same mode.
!  If a nonzero logical unit number is specified, any error messages will
!  be directed to that unit.

!         ix specifies tensor symmetry in the following way:
!             ix < -100: axial mode
!      -100 < ix < 0   : spherical mode
!             ix > 0   : cartesian mode

!         Mode flag =  0 indicates that mode has not yet been set
!                      1 indicates that Cartesian mode has been set
!                      2 indicates that spherical mode has been set
!                      3 indicates that axial mode has been set

!----------------------------------------------------------------------

    function tcheck(ix,token,lumsg)
    implicit none
    logical :: tcheck

    include 'limits.inc'
    include 'simparm.inc'
    include 'parms.inc'
    include 'lpnam.inc'
    include 'miscel.inc'
    include 'stdio.inc'

    integer :: ispect,ix,ixa,ixf,lumsg,mflag,CARTESIAN,SPHERICAL,AXIAL
    integer :: isite,ispec
    parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)

    character token*(*)

!......................................................................

    if (ix > 0) ispect=CARTESIAN
    if (ix < -100) ispect=AXIAL
    if (ix > -100 .AND. ix < 0) ispect=SPHERICAL

    ixa=abs(mod(ix,100))
    ixf=(ixa-IGXX)/3		! 0,1 or 2

!     --- Return .true. for all parameters that have no aliases

    if (ixa < IGXX .OR. ixf > mxsph-1) then
        tcheck=.true.
        return
    end if

    mflag=iparm(IIGFLG+ixf,1,1)    ! get the flag for this parameter

! Check if mode is not yet set

    if (mflag == 0) then
        if (lumsg /= 0) write (lumsg,1003) symstr(ispect),token(1:1)
    
    ! Check if tensors are specified as Cartesian when another mode is set
    
    else if (ispect == CARTESIAN .AND. &
        (mflag == SPHERICAL .OR. mflag == AXIAL)) then
        if (lumsg /= 0) write(lumsg,1004) token,symstr(mflag)
        tcheck=.false.
        return
    
    ! Check if tensors are specified as spherical when another mode is set
    
    else if (ispect == SPHERICAL .AND. &
        (mflag == CARTESIAN .OR. mflag == AXIAL)) then
        if (lumsg /= 0) write(lumsg,1004) token,symstr(mflag)
        tcheck=.false.
        return
    
    ! Check if tensors are specified as axial when another mode is set
    
    else if (ispect == AXIAL .AND. &
        (mflag == CARTESIAN .OR. mflag == SPHERICAL)) then
        if (lumsg /= 0) write(lumsg,1004) token,symstr(mflag)
        tcheck=.false.
        return
    end if

!  Set tensor mode flags according to type of tensor that has been
!  specified in all sites and spectra.

    do 30 isite=1,ncomps
        do 30 ispec=1,nspectra
            iparm(IIGFLG+ixf,ispec,isite)=ispect
    30 END DO
    tcheck=.true.
    return

! ### format statements ########################################

    1003 format(' *** ',a,' form assumed for ',a,' tensor ***')
    1004 format(' *** Cannot modify ''',a,''': ',a,' form', &
    ' has been specified ***')
    end function tcheck
