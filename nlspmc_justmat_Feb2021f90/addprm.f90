!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                       subroutine ADDPRM
!                    =========================

! Add a parameter to list of parameters being varied for nonlinear
! least-squares. Also maintain the lists of
!    (1) ibnd    : boundary flag for variable parameter
!    (2) prmin   : minimum for variable parameter
!    (3) prmax   : maximum for variable parameter
!    (4) prscl   : desired accuracy for given parameter
!    (5) xfdstp  : Size of forward-differences step
!    (6) tag     : Name of the parameter (character*6)
!    (7) ixx     : index of each element in the fparm array into the x vector

! The values of these parameters for the given variable parameter
! are passed to ADDPRM.

! ix - returned by ipfindc
! ibd - ibnd   : Flag for boundaries imposed on each parameter:
!          0=none, 1=minimum, 2=maximum, 3=both
! prmn,prmx,prsc,step - min, max, accuracy, fd-step
! prmID - parameter prefix
! ixrspec,ixrsite - specify spectra and sites to vary via indtkn.

!----------------------------------------------------------------------
    subroutine addprm(ixr,ibd,prmn,prmx,prsc,step,prmID, &
    ixrspec,ixrsite)
    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'lmcomm.inc'
    include 'stdio.inc'
    include 'rndoff.inc'
    include 'miscel.inc'

    integer :: ixr,ibd,ispec,isite,specptr,siteptr,iparmptr, &
    ixrspec,ixrsite
    double precision :: prmn,prmx,prsc,step
    character(6) :: prmID
    integer :: i,ixabs,j,k,k1,nmov

    logical :: tcheck
    external tcheck

!######################################################################

!----------------------------------------------------------------------
! See if there is enough room in the parameter list
!----------------------------------------------------------------------
    if (nprm+1 > MXVAR) then
        write(luout,1001) MXVAR
        if (luout /= luttyo) write(luttyo,1001) MXVAR
        return
    end if
! Check axis specification for consistency
    if ( .NOT. tcheck(ixr,prmID,luttyo)) then
        write(*,*)'axis information inconsistent in addprm'
        return
    end if

!----------------------------------------------------------------------
! Search parameter list to see if parameter is already there
! Must check specific site/spectrum information here.
!----------------------------------------------------------------------

!      ixabs=abs(mod(ixr,100))
    specptr=1                 ! get the old information
    siteptr=1
    if (ixrspec /= 0) specptr=ixrspec         ! get first spec
    if (ixrsite /= 0) siteptr=ixrsite         ! get first site
    iparmptr=iabs(mod(ixr,100))
! consider all sites, spec:
    do 30 isite=1,ncomps
        do 30 ispec=1,nspectra
        ! check matching sites with requested site
            if ((isite == ixrsite .OR. ixrsite == 0) .AND. &
            (ispec == ixrspec .OR. ixrspec == 0)) then
            ! if vary previously selected
                if (ixx(iparmptr,ispec,isite) /= 0) then
                    write(luttyo,1000) prmID
                    1000 format('*** Parameter ''',a,''' already varied ***')
                    return
                end if
            end if
    30 END DO
! so add parameter to list
! We use ixpr to keep track of the parameters in the list of
! variable paremeters.

!      do 10 i=1,nprm
!                                      *** parameter already in list
!         if (ixpr(i).eq.ixabs) then
!            write(luout,1000) prmID
!            if (luout.ne.luttyo) write(luttyo,1000) prmID
! 1000       format('*** Parameter ''',a,''' already being varied ***')
!            return
!         else if (ixpr(i).gt.ixabs) then
!            go to 14
!         end if
! 10   continue

!----------------------------------------------------------------------
! Reached end of list: parameter may be added to the end
!----------------------------------------------------------------------
    nprm=nprm+1

!----------------------------------------------------------------------
! Found proper location for new parameter: move top of list up
!----------------------------------------------------------------------
! 14   nmov=nprm-i+1
!      do 16 j=1,nmov
!         k=nprm-j+1
!         k1=k+1
!         x(k1)=x(k)
!         ixpr(k1)=ixpr(k)
!         prmin(k1)=prmin(k)
!         prmax(k1)=prmax(k)
!         prscl(k1)=prscl(k)
!         xfdstp(k1)=xfdstp(k)
!         ibnd(k1)=ibnd(k)
!         tag(k1)=tag(k)
!         ixx(ixpr(k1))=ixx(ixpr(k1))+1
! 16   continue
!      nprm=nprm+1

!----------------------------------------------------------------------
!  Add new parameter(s) to proper location list
!----------------------------------------------------------------------
    18 xxx(nprm)=fparm(iparmptr,specptr,siteptr)
    ixpr(nprm)=iparmptr
    do 20 isite=1,ncomps
        do 20 ispec=1,nspectra
        ! check matching sites with requested site
            if ((isite == ixrsite .OR. ixrsite == 0) .AND. &
            (ispec == ixrspec .OR. ixrspec == 0)) then
                ixx(iparmptr,ispec,isite)=nprm
            ! check that all fparm values are equal.
                if ((abs(fparm(iparmptr,specptr,siteptr)-fparm(iparmptr, &
                ispec,isite))) > 0.00001*abs(fparm( &
                iparmptr,specptr,siteptr)))then
                    write(luttyo,1110) prmID
                    write(*,*)'values ',fparm(iparmptr,specptr,siteptr), &
                    fparm(iparmptr,ispec,isite)
                    1110 format('*** Parameter ''',a,''' not consistent ***')
                end if
            end if
    20 END DO

    prmin(nprm)=prmn
    prmax(nprm)=prmx
    prscl(nprm)=prsc
    xfdstp(nprm)=step*xxx(nprm)
    if(dabs(xfdstp(nprm)) < 1000.0d0*rndoff) xfdstp(nprm)=step
    ibnd(nprm)=ibd
    tag(nprm)=prmID

    write(luout,1003) nprm
    if (luout /= luttyo) write(luttyo,1003) nprm
    return

! #### format statements ###########################################

    1001 format('*** Limit of',i2,' variable parameters exceeded ***')
    1003 format(i3,' parameters now being varied')
    end subroutine addprm
