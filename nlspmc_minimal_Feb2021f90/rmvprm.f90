!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                       subroutine RMVPRM
!                    =========================

! Remove a parameter from the list of parameters being varied for nonlinear
! least-squares. Also maintain the lists of
!    (1) ibnd    : boundary flag for variable parameter
!    (2) prmin   : minimum for variable parameter
!    (3) prmax   : maximum for variable parameter
!    (4) prscl   : desired accuracy for given parameter
!    (5) xfdstp  : Size of forward-differences step
!    (6) tag     : Name of the parameter (character*6)
!    (7) ixx     : index of each element in the fparm array into the x vector

! ixr identifies variable, with axis identification,
! ident is the the location in xxx of the variable to be removed.
! Note the xxx array itself is set elsewhere on calling LM routine

! Includes:

!      limits.inc
!      parms.inc
!      simparm.inc
!      stdio.inc
!      miscel.inc

! Spectrum and site information is checked against current variable
! parameter specification.

!----------------------------------------------------------------------
    subroutine rmvprm(ident)
    implicit none
    integer :: varid,ident,icomp,ispec

    include 'limits.inc'
    include 'parms.inc'
    include 'simparm.inc'
    include 'stdio.inc'
    include 'miscel.inc'

    integer :: j,j1
    character(6) :: idparm

!######################################################################

    if (nprm == 0) then
        write (luttyo,1000)
        1000 format('*** No parameters are being varied ***')
        return	! if no variable parameters, return
    end if

    idparm=tag(ident)		! get ascii parameter name
    write (luttyo,1001) idparm,nprm-1
    1001 format('*** Parameter ''',a,''' fixed: ', &
    i3,' variable parameters remaining ***')

!  Delete it and move up elements below it in the list.
!  (Note that move loop will not execute if we are removing last element;
!  i.e. if ident=nprm)
!------------------------------------------------------------------------

    varid=ixpr(ident)		! identity of the variable
    do 1 j=ident,nprm-1
        j1=j+1
        ixpr(j)=ixpr(j1)
        prmin(j)=prmin(j1)
        prmax(j)=prmax(j1)
        prscl(j)=prscl(j1)
        xfdstp(j)=xfdstp(j1)
        ibnd(j)=ibnd(j1)
        tag(j)=tag(j1)
    1 END DO
    do 2 icomp=1,ncomps
        do 2 ispec=1,nspectra
        ! zero ixx pointer for parameter being fixed
            if (ixx(varid,ispec,icomp) == ident) then
                ixx(varid,ispec,icomp)=0
            end if
    2 END DO
! decrement ixx pointer for moved parms
    do 3 j=1,nprm-1
        do 3 icomp=1,ncomps
            do 3 ispec=1,nspectra
                if (ixx(ixpr(j),ispec,icomp) > ident) &
                ixx(ixpr(j),ispec,icomp)=ixx(ixpr(j),ispec,icomp)-1
    3 END DO

    nprm=nprm-1
    return

! #### format statements ###########################################

    end subroutine rmvprm
