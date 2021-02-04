! NLSPMC VERSION (VERSION 1.0) 2/5/99
!----------------------------------------------------------------------
!                         ====================
!                          subroutine CONVTC
!                         ====================

!     Subroutine to convert the tensor named on the passed command line
!     according to the symmetry specified by iflg. The following are
!     valid symmetries:
!         1:  Cartesian
!         2:  Spherical
!         3:  Axial
!     See source file tensym.f for an explanation of these symmetries.

!     Uses:
!         gettkn, touppr (in strutl.f)
!         ipfind.f
!         rmvprm.f
!         tocart, tosphr, toaxil (in tensym.f)

!     Includes
!         limits.inc
!         simparm.inc
!         parms.inc
!         stdio.inc
!         miscel.inc
!         names.inc
!         lpnam.inc

! 7/31/98 - modified to not use ipfind but rather to just for A, G or R.
! 	    Previous version didn't handle G correctly, (saw GIB rather
!	    than GXX)
! 	Symmetries apply to all spectra, all components.

!----------------------------------------------------------------------

    subroutine convtc( line,iflg )
    implicit none
    include 'limits.inc'
    include 'simparm.inc'
    include 'parms.inc'
    include 'stdio.inc'
    include 'miscel.inc'
    include 'names.inc'
    include 'lpnam.inc'

    character line*(LINELG), token*(WORDLG)
    integer :: i,iflg,ix,ixa,ixf,ixten,j,k,lth,lu,ispec,icomp,ident
    integer :: CARTESIAN,SPHERICAL,AXIAL
    parameter (CARTESIAN=1,SPHERICAL=2,AXIAL=3)

    integer :: ipfind
    external ipfind

!######################################################################

    if (iflg < CARTESIAN .OR. iflg > AXIAL) iflg=CARTESIAN

    call gettkn(line,token,lth)
    call touppr(token,lth)

!                                          *** No tensor specified
    if (lth == 0) then
        write (luttyo,1000)
        return
    end if

! Check for specific tensors:

    if (token(:1) == 'A') then
        ixa=IAXX
    elseif (token(:1) == 'G') then
        ixa=IGXX
    elseif (token(:1) == 'R') then
        ixa=IDX
    else
        write (luout,1001) token(1:1)
        if (luout /= luttyo) write (luttyo,1001) token(1:1)
        return
    end if

!----------------------------------------------------------------------
!     Tensor found: Check existing symmetry: a) all equal? b) set?
!----------------------------------------------------------------------
!      else
    ixf=IIGFLG+(ixa-IGXX)/3	! 0,1 or 2 added to IIGFLG
    ixten=IGXX+3*(ixf-IIGFLG)	! 0,3 or 6 added to IGXX
!----------------------------------------------------------------------
! check existing symmetry, insist equal for all sites, spectra:
!----------------------------------------------------------------------
    i=iparm(ixf,1,1)
    do 11 icomp=1,ncomps
        do 11 ispec=1,nspectra
            if (iparm(ixf,ispec,icomp) /= i) then
                write(luttyo,*)'in convtc, error, multiple symmetries '
                return
            end if
    11 END DO
!----------------------------------------------------------------------
!        Symmetry not set yet: set it and exit
!----------------------------------------------------------------------
    if (iparm(ixf,1,1) == 0) then
        do 20 icomp=1,ncomps
            do 20 ispec=1,nspectra
                iparm(ixf,ispec,icomp)=iflg
        20 END DO
        write(luttyo,1002) token(1:1),symstr(iflg)
        return
    
    !----------------------------------------------------------------------
    !        Symmetry same as the one specified: exit
    !----------------------------------------------------------------------
    else if (iparm(ixf,1,1) == iflg) then
        write(luttyo,*)'in convtc, symmetry unchanged, &
        nothing changed'
        return
    end if

!----------------------------------------------------------------------
!        Here when want new symmetry:
!        Remove any tensor components of any symmetry from
!        the list of variable parameters
!----------------------------------------------------------------------
    do 10 i=0,2	!all components of tensor to be changed
        j=ixten+i
        do 30 icomp=1,ncomps
            do 30 ispec=1,nspectra
                if (ixx(j,ispec,icomp) /= 0) then	! if we vary this one
                    ident=ixx(j,ispec,icomp)		! get xxx index
                    call rmvprm(ident)       ! remove it
                ! rmvprm when asked to remove a variable, removes all specs/comps
                ! where that variable applies, so this call should only happen for the
                ! first occasion.
                end if
        30 END DO
    10 END DO

!----------------------------------------------------------------------
!        Now...convert the tensor symmetry!
!----------------------------------------------------------------------
    do 40 icomp=1,ncomps
        do 40 ispec=1,nspectra
            if (iflg == 2) then
                call tosphr( fparm(ixten,ispec,icomp), &
                iparm(ixf,ispec,icomp) )
            else if (iflg == 3) then
                call toaxil( fparm(ixten,ispec,icomp), &
                iparm(ixf,ispec,icomp) )
            else
                call tocart( fparm(ixten,ispec,icomp), &
                iparm(ixf,ispec,icomp) )
            end if
            iparm(ixf,ispec,icomp)=iflg
        !           end if
    40 END DO
    write (luout,1003) token(1:1),symstr(iflg)
    if (luout /= luttyo) write (luttyo,1003) token(1:1),symstr(iflg)

!----------------------------------------------------------------------
!        Report new values for site 1 spectrum 1
!----------------------------------------------------------------------
    lu=luout
    300 continue
    write(lu,*)'in convtc ',ncomps,iflg
    do 400 icomp=1,ncomps
        write (lu,*)'For component # ',icomp
        if (ixf == IIGFLG) then
        !                                                  * g tensor
            if (iflg == 1) then
                write (lu,1011) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            else if (iflg == 2) then
                write (lu,1012) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            else
                write (lu,1013) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            end if
        else if (ixf == IIGFLG+1) then
        !                                                  * hf tensor
            if (iflg == 1) then
                write (lu,1021) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            else if (iflg == 2) then
                write (lu,1022) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            else
                write (lu,1023) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            end if
        else
        !                                                  * r tensor
            if (iflg == 1) then
                write (lu,1031) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            else if (iflg == 2) then
                write (lu,1032) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            else
                write (lu,1033) icomp,(fparm(ixten+i,1,icomp),i=0,2)
            end if
        end if
    400 END DO
    if (lu /= luttyo) then
        lu=luttyo
        go to 300
    end if

    return

    1000 format('*** No tensor specified ***')
    1001 format('*** Unknown tensor: ''',a,''' ***')
    1002 format('*** ',a1,' tensor set for all spec,comps to ',a,' ***')
    1003 format('*** ',a1,' tensor converted to ',a,' ***')
    1011 format('spec 1, site ',i2,': gxx,gyy,gzz = ',2(f9.6,','),f9.6)
    1012 format('spec 1, site ',i2,': g1,g2,g3 = ',2(f9.6,','),f9.6)
    1013 format('spec 1, site ',i2,': gprp,grhm,gpll = ',2(f9.6,','),f9.6)
    1021 format('spec 1, site ',i2,': Axx,Ayy,Azz = ',2(f9.6,','),f9.6)
    1022 format('spec 1, site ',i2,': A1,A2,A3 = ',2(f9.6,','),f9.6)
    1023 format('spec 1, site ',i2,': Aprp,Arhm,Apll = ',2(f9.6,','),f9.6)
    1031 format('spec 1, site ',i2,': log(Rx,Ry,Rz) = ',2(f9.6,','),f9.6)
    1032 format('spec 1, site ',i2,': log(Rbar,N,Nxy) = ',2(f9.6,','), &
    f9.6)
    1033 format('spec 1, site ',i2,': log(Rprp,Rrhm,Rpll) = ', &
    2(f9.6,','),f9.6)
    end subroutine convtc
