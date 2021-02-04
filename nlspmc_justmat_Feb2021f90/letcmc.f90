!  VERSION 1.0  NLSPMC   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                      subroutine LETCMC
!                    =========================
! Interprete the command:
! let <name>{(index)} {, <name>,...} = <value> {, <value> ... }

! let only allows setting of variables that apply to a specific site
! and all spectra.  To input variables applying to one or more spectra,
! in a multiple series, use the series command.  In the case of only
! one of a series, the let may be used or the series command.

!      name    : name of the parameter to be assigned a value
!      index   : site or spectral index of parameter
!      value   : value to be assigned to the given parameter

!   Up to 10 variables may be assigned values in a single let command

! Parse the input line into arrays (leftid, leftspect and leftsite)
! containing the name and site information in the form: 0=all
! sites/spectra, else=specific site/spectrum.

! Then parse into the final arrays fparm(variable,spectrum,site) and
! iparm(variable,spectrum,site).  Based on the CW version with mods.

! Later for each spectrum,site, we copy the current values into xxx
! for fitting, or into fepr and iepr in sim2d for simulation.

! Errors in various checks for consistency result in return without
! further processing of line.

!----------------------------------------------------------------------
    subroutine letcmc(line)
    implicit none

!      include 'limits.inc'
!	integer mxdim
!	parameter mxdim=1200

    include 'limits.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'lpnam.inc'
    include 'miscel.inc'
    include 'stdio.inc'

!  maximum number of variables in command line:
    character*(LINELG) line
    integer :: mxlft,specptr,siteptr,varyinfo,iparmptr,isite,ispec
    parameter (mxlft=10)
    logical :: vchange
    external vchange

    integer :: ixs,ival,instp,vspc
!      integer ileft,irt,ival,ixf,ixr,ixs,ixri,ixt,lth,left(mxlft)
!  ileft,irt = counters for parameter names and values
!  leftid,leftspec,leftsite = hold identification # for parameter,
!   spectrum and site info (0=all spectra or all sites).

    integer :: leftid(mxlft),leftspec(mxlft),leftsite(mxlft)
!      integer ileft,irt,lth,left(mxlft),leftx(mxlft),llth
    integer :: ileft,irt,lth,llth,ixr,ixrspec,ixrsite
    double precision :: fval

    character token*(WORDLG),prmID*6
!      logical rhs,logflg,serprm
    logical :: rhs,logflg

    integer :: isfind,varspec(mxlft)
    logical :: ftoken,itoken,tcheck
    external ftoken,itoken,tcheck,isfind

!######################################################################

    rhs=.false.
    logflg=.false.
! change starting point for irt and ileft counters.
    ileft=0
    irt=0

!----------------------------------------------------------------------
! Get the next parameter name (or value for right-hand-side of statement)
! line = line of text to be parsed, token = returned token, lth = length
!----------------------------------------------------------------------
    1 call gettkn(line,token,lth)
    lth=min(lth,6)
    call touppr(token,lth)

!----------------------------------------------------------------------
!  Check for error conditions
!----------------------------------------------------------------------
    if (lth == 0) then		! expect a value or parameter
        write (luttyo,1000)
        1000 format('*** not enough values ***')
        return
    end if

    if (token == ')' .OR. token == '(') then
    
    ! If (, not expected at this point.
    

        write(luttyo,1014)
        1014 format('*** Unexpected ( or ) ***')
        return
    end if

    if (token == 'LOG') then
        if(rhs) then
            logflg=.true.
            go to 1
        else
            write(luttyo,1015)
            1015 format('*** Log before =? ***')
        end if
        return
    end if

!----------------------------------------------------------------------
!  Check for '=' sign in assignment
!----------------------------------------------------------------------
    if (token == '=') then
        rhs=.true.
        go to 1
    end if

!----------------------------------------------------------------------
!  Left-hand side: assign parameter index and site information
!----------------------------------------------------------------------

    if ( .NOT. rhs) then	! if lefthand side:
    
    ! get parameter id # and save in leftid.
    
        if(ileft > mxlft-1) then
            write (luttyo,1016) token(:lth)
            1016 format('*** Too many parmeters in let''',a, ''' ***')
            return
        end if
    ! get the parameter id # including axis information:
    !         leftid(ileft+1)=ipfind(token,lth)	! # coded in eprprm.inc
        call ipfind(leftid(ileft+1),token,lth,varspec(ileft+1))
        if (leftid(ileft+1) == 0) then
            write (luttyo,1003) token(:lth)	! name not found
            1003 format('*** ''',a,''' is not a parameter ***')
            go to 1
        end if
    
        call indtkn(line,leftspec(ileft+1),leftsite(ileft+1))
        if(leftspec(ileft+1) == -1 .OR. leftsite(ileft+1) == -1) then
            write(luttyo,*)'error in indtkn '
            return
        end if
    
        ileft=ileft+1

        go to 1			! go get another
    end if


!----------------------------------------------------------------------
!  Right-hand side: assign value to parameter specified in leftid array
!----------------------------------------------------------------------

    if (rhs) then
        irt=irt+1
        if(irt > ileft) Then
            write(luttyo,1017)
            1017 format('*** too many right hand side values ***')
            return
        end if
        ixr=leftid(irt)		! pointer to parameter name
        ixrspec=leftspec(irt)		! spectra info, 0=all spectra
        ixrsite=leftsite(irt)		! site info, 0=all site
    
    !   --- EPR parameter assignment.  See ipfind.f for meaning of values.
    !  Note, axial and spherical tensor elements are pointed to by IGXX,
    !  while fp and integer are pointed to by ixr directly.
    
        if (ixr < -100 .AND. ixr > -200) then		! axial tensor
            prmID=alias2(-99-(IGXX+ixr))
        else if (ixr < 0) then				! spherical
            prmID=alias1(1-(IGXX+ixr))
        else if (ixr > 100) then			! integer
            prmID=parnam(ixr-100)
        else if (ixr > 0 .AND. ixr < 100) then	! floating point
            prmID=parnam(ixr)
        else
            write(luttyo,1018)ixr
            1018 format('*** Error, parm unknown, ixr ***''',i8,''' ***')
            return
        end if
    
        if (ixr < 100) then		! floating point.
        ! get fval if fp:
            if ( .NOT. ftoken(token,lth,fval)) then	! assigns fval
                write(luttyo,1004) token(:lth)
                1004 format('*** Real value expected: ''',a,''' ***')
                return
            end if
        ! got fval, store it:
            if (logflg) then		! take log if necessary
                if (fval > 0.0d0) then
                    fval=dlog10(fval)
                else
                    write(luttyo,1011) fval
                    1011 format(' *** Illegal log argument:',g12.5,' ***')
                    return
                end if
                logflg=.false.
            end if
        !           -------------------------------
        !           Assign the floating point value
        !           -------------------------------
        !  Check for consistent mode of tensors that depend of axis system.
        !  Set mode if previously undefined.  Compare ixr with mode flag.
        
            if (tcheck(ixr,prmID,luout)) then	! tensym.f, if valid axis
            
            !                 ---------------------------------------------
            !                 Issue a warning if:
            !                    (1) the parameter is being set for a specific
            !                        site/spectrum when it is being varied for
            !                        *all* sites/spectra
            !                    (2) the parameter is being set for all sites/spectra
            !                        when its current value is different for the
            !                        currently defined sites/spectra
            !                 ---------------------------------------------
            !  for now, we only allow variation of parameters that apply to all
            !  spectra.  Just check the vary information, don't change it.
            
                specptr=1			! get the old information
                siteptr=1
                if (ixrspec /= 0) specptr=ixrspec		! get first spec
                if (ixrsite /= 0) siteptr=ixrsite		! get first site
                iparmptr=iabs(mod(ixr,100))
                varyinfo=ixx(iparmptr,specptr,siteptr)  ! vary this?
            !	      mattrisp=varspec(iparmptr)
            ! consider all sites, spec:
                do 30 isite=1,ncomps
                    do 30 ispec=1,nspectra
                        if ((isite == ixrsite .OR. ixrsite == 0) .AND. &
                        (ispec == ixrspec .OR. ixrspec == 0)) then
                            if (ixx(iparmptr,specptr,siteptr) /= &
                            varyinfo) then
                                write(luttyo,1021) prmID
                                1021 format(' *** Parm variation inconsistent:', &
                                a,' ***')
                            else
                            ! if parameter is changed, update "have matrix" information:
                            ! Drop basinfo(2 information, replace with matrxok.
                            !                      if (basinfo(2,ispec,isite) .ge. mattrisp .and.
                            !     #			(abs(fparm(iparmptr,ispec,isite)-fval)) .gt.
                            !     #			1.d-8*abs(fval)) then
                            !                        basinfo(2,ispec,isite)= mattrisp-1
                            !                      end if
                                if (vchange(fval,fparm(iparmptr,ispec,isite))) &
                                matrxok(isite)= &
                                min(matrxok(isite),specid(iparmptr))
                                fparm(iparmptr,ispec,isite)=fval
                            end if
                        end if
                30 END DO
            
            else
                write(luout,*)'Tcheck failure in letcmc'
                return
            end if		! end if tcheck
        
        !  --- EPR parameter integer assignment
        
        else if (ixr > 100) then		! integer parameter
        ! do we have a symbol?
            ixs=isfind(token,lth)		! zero if not symbol
            if (ixs > 0) then
                ival=symval(ixs)
            else if ( .NOT. itoken(token,lth,ival)) then  ! integer->ival
            ! error result (itoken=.false.):
                write(luttyo,1005) token(:lth)
                1005 format('*** Integer value expected: ''',a,''' ***')
                return
            end if
        
        ! have ival, do some checks here:
        
        
        ! check # of CG steps against maximum
        
            call ipfind(instp,'NSTEP',5,vspc)  ! get pointer to 'NSTEP'
            if (ixr == instp) then	! if the keyword is NSTEP:
                if (ival > MXSTEP-2) then
                    write(luttyo,*)'nstep too big, set to ',MXSTEP-2
                    ival=MXSTEP-2
                end if
            end if
        
        
        
        
        ! set for specified sites, spectra:
        
            do 40 isite=1,ncomps
                do 40 ispec=1,nspectra
                ! if this spectra/site is to be changed:
                    if ((isite == ixrsite .OR. ixrsite == 0) .AND. &
                    (ispec == ixrspec .OR. ixrspec == 0)) then
                        iparm(ixr-100,ispec,isite)=ival
                    end if
            40 END DO
        end if
    end if	! end of if(rhs)

!    --- Return if all assignments have been made

    if(irt == ileft) then
        return
    end if
    go to 1
! ###### format statements ########################################

    1001 format('*** variable name expected ***')
    1002 format('*** no ''='' specified ***')
    1012 format(' *** Series parameters in let command : series ', &
    'assumed ***')
    end subroutine letcmc
