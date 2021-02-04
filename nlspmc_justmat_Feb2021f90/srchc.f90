!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                       subroutine SRCHC
!                    =========================

!   search  <parameter> {spectrum,site} { tol <tol> step <step>
!                         min <smn> max <smx> start <val> }

!    spectrum,site : specify which spectra or sites to apply parameter to
!    tol     : Convergence tolerance for fitting parameter
!    step    : Size of first step to take away from initial parameter value
!    smn,smx : minimum and maximum that defines the search range
!    val     : starting value of the Brent's search procedure
!              (Note that the bracketing procedure is not performed
!               if this option is used.  Therefore the used have to
!               make sure that there exists a minimum chi-squared within
!               the range.)

!----------------------------------------------------------------------
    subroutine srchc( line )
    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'parms.inc'
    include 'lpnam.inc'
    include 'stdio.inc'
    include 'miscel.inc'
!      include 'iterat.inc'

    logical :: vchange,gotit
    external vchange
    character line*(LINELG)
    integer :: ixp1p,i,itmp,ierr,lth,ibrckt,ixxsv(NFPRM,MXSPEC,MXSITE)
    double precision :: ax,bx,cx,fa,fb,fc,fval,xmin,dummy
    character token*(WORDLG),prmID*6

    integer :: nkeywd,ipprm
    parameter(nkeywd=5)

    character(8) :: keywrd(nkeywd)
    data keywrd / 'TOL','STEP','MIN','MAX','START'/

    logical :: ftoken,prmsOK,tcheck
    integer :: isite,ispec,specinfo,siteinfo
    integer :: nprmsv,varspec
    double precision :: brent,p1pfun
    external ftoken,prmsOK,tcheck,brent,p1pfun

!######################################################################

! test all basis sets here.
    setbas=.true.
    gotit=.true.
    do 1 isite=1,ncomps
    !        seteval(isite)=.false.	 ! make no assumption about having matrix
        do 1 ispec=1,nspectra
            setbas=setbas .and. (basinfo(1,ispec,isite) .ne. 0)
            if (basinfo(1,1,isite) /= basinfo(1,ispec,isite)) then
                write(luttyo,*)'Error, multiple basis sets not allowed'
                write(luttyo,*)'except for different sites.'
                return
            end if
    1 END DO
    if ( .NOT. setbas) then
        write (luout,1200)
        if (luout /= luttyo) write (luttyo,1200)
        1200 format(/,'** Basis sets not properly set **')
        return
    end if

!  -- Check whether any starting parameters have been input before
!     proceeding with the FIT command

    if ( .NOT. prmsOK(luttyo) ) then	! this checks all site/spect
        write(luttyo,1040)
        return
    end if

    ierr=0
    ibrckt=0

!----------------------------------------------------------------------
! Get the name of the parameter
!----------------------------------------------------------------------
    call gettkn(line,token,lth)
    lth=min(lth,6)
    call touppr(token,lth)

!      ixp1p=ipfind(token,lth)
    call ipfind(ixp1p,token,lth,varspec)

!----------------------------------------------------------------------
!     Check whether parameter may be varied
!----------------------------------------------------------------------
    if (ixp1p == 0 .OR. ixp1p > NFPRM) then
        write(luttyo,1002) token(:lth)
        return
    end if

    if (ixp1p < -100) then
        prmID=alias2( -99-(IGXX+ixp1p) )
    else if (ixp1p < 0) then
        prmID=alias1( 1-(IGXX+ixp1p) )
    else
        prmID=parnam(ixp1p)
    end if
! get spec/site information:
    call indtkn(line,specinfo,siteinfo)
    isite=siteinfo
    ispec=specinfo
    if (specinfo == 0) ispec=1
    if (siteinfo == 0) isite=1
    if (specinfo /= 0 .AND. nspectra > 1) then
        write(luttyo,*)'Error, search must apply to all spectra.'
        write(luttyo,*)'Search may apply to one site.'
        return
    end if


!  --- Check whether proper symmetry has been specified for
!      tensor components - all spec,sites have same symmetry.

    if ( .NOT. tcheck(ixp1p,prmID,luttyo)) return
    ixp1p=iabs( mod(ixp1p,100) )

! Set default search range

    smn=fparm(ixp1p,ispec,isite)-1.0d0	! ispec must = 1
    smx=fparm(ixp1p,ispec,isite)+1.0d0
    stval=fparm(ixp1p,ispec,isite)
!					! stol and sstep set in nlsinit.
    if (ixp1p == IC20) smn=1.0d-2

!----------------------------------------------------------------------
!  Look for a keyword
!----------------------------------------------------------------------

    14 call gettkn(line,token,lth)
    lth=min(lth,8)

!------------------------------------------------
!************************************************
! No more keywords: call the line search routine
!************************************************
!------------------------------------------------
    if (lth == 0) then
    
    !         seteval=.false.
    
    ! Set two initial parameter values for the parameter to be searched
    
        ax=fparm(ixp1p,ispec,isite)
        bx=fparm(ixp1p,ispec,isite)+sstep
    
        if (smn > ax) smn=ax-1.0d0
        if (smx < ax) smx=ax+1.0d0
    
    !  Swap indices for the search parameter into the first elements
    !  of the parameter for the NLS x-vector
    !  (This is so x can be used by the pfun routine as in a NLS procedure)
    
    ! proceed as follows: copy ixx, the index identifying which parameters
    ! are in xxx, into temporary array (don't modify xxx), set all ixx
    ! entries = 0, input search variable id into ixx(1) at all applicable
    ! sites and spectra.  Check all sites in fvec are equal, then continue.
    ! At end of search undo these steps.  Be sure to reset this on return.
    
        nprmsv=nprm
        itmp=ixpr(1)
        ixpr(1)=ixp1p
        nprm=1
        do 10 ispec=1,nspectra
            do 10 isite=1,ncomps
            ! loop over all parms
                do 10 ipprm=1,NFPRM
                    ixxsv(ipprm,ispec,isite)=ixx(ipprm,ispec,isite)
                    ixx(ipprm,ispec,isite)=0
                ! if this parmeter, site and spectra are to be searched, set it as xxx(1):
                    if (ixp1p == ipprm .AND. &
                    (isite == siteinfo .OR. siteinfo == 0) .AND. &
                    (ispec == specinfo .OR. specinfo == 0)) then
                        ixx(ipprm,ispec,isite) = 1	! set this variable
                    end if
        10 END DO
    
        if (ibrckt == 0) then
            write (luout,1100) smn,ax,smx
            if (luout /= luttyo) write (luttyo,1100) smn,ax,smx
            call mnbrak(ax,bx,cx,fa,fb,fc,p1pfun,smn,smx,ierr)
            if (ierr < 0 .OR. ihltcmd /= 0) go to 202	! bad pfun return
            write (luout,1105) cx,bx,ax,fc,fb,fa
            if (luout /= luttyo) write (luttyo,1105) cx,bx,ax,fc,fb,fa
            if(ierr /= 0) then
                write(luout,*)'error returned from mnbrak ',ierr
                if (luout /= luttyo) write (luttyo,*) &
                'error returned from mnbrak ',ierr
                go to 202
            end if
        else
            ax=smn
            bx=stval
            cx=smx
            write (luout,1100) smn,stval,smx
            if (luout /= luttyo) write (luttyo,1100) smn,stval,smx
        end if
    
        if (ierr == 0) then
            dummy=brent(ax,bx,cx,p1pfun,stol,xmin,prmID,ierr)
            if (ierr < 0 .OR. ihltcmd /= 0) go to 202
        else
            gotit=.false.
            write (luout,1110)
            if (luout /= luttyo) write (luttyo,1110)
            if (ax > bx) then
                write (luout,1120) cx,ax,fc,fa
                if (luout /= luttyo) write (luttyo,1120) cx,ax,fc,fa
            else
                write (luout,1120) ax,cx,fa,fc
                if (luout /= luttyo) write (luttyo,1120) ax,cx,fa,fc
            end if
        end if
    ! put the parameters back, put the new found value into fparm
        nprm=nprmsv
        ixpr(1)=itmp
        do 20 ispec=1,nspectra
            do 20 isite=1,ncomps
                do 20 ipprm=1,NFPRM
                    ixx(ipprm,ispec,isite)=ixxsv(ipprm,ispec,isite)
                    if (ixp1p == ipprm .AND. &
                    (isite == siteinfo .OR. siteinfo == 0) .AND. &
                    (ispec == specinfo .OR. specinfo == 0)) then
                    ! put back new found minimum, No need to change basinfo, was set when
                    ! spectrum was calculated
                    !                   if (basinfo(2,ispec,isite) .ge. varspec .and.
                    !     #               (abs(fparm(iparmptr,ispec,isite)-fval)) .gt.
                    !     #               1.d-8*abs(fval)) then
                    !                     basinfo(2,ispec,isite)= varspec-1
                    !                   end if
                    
                    ! if there is a parameter change, set matrixok for new calculation
                    
                        if (gotit .AND. vchange(xmin, &
                        fparm(ixp1p,ispec,isite))) then
                            matrxok(isite)=min(matrxok(isite),specid(ixp1p))
                            fparm(ixp1p,ispec,isite)=xmin
                        end if
                    end if
        20 END DO
        return
    ! error return, reset ixx first:
        202 nprm=nprmsv
        ixpr(1)=itmp
        do 201 ispec=1,nspectra
            do 201 isite=1,ncomps
                do 201 ipprm=1,NFPRM
                    ixx(ipprm,ispec,isite)=ixxsv(ipprm,ispec,isite)
        201 END DO
        return
    end if

!------------------------------
! Check keyword list (lth .ne. 0)
!------------------------------
    call touppr(token,lth)
    do 15 i=1,nkeywd
        if (token(:lth) == keywrd(i)(:lth)) goto 16
    15 END DO
!                                        *** Unrecognized keyword
    write (luttyo,1000) token(:lth)
    return

!----------------------------------------------------------------------
!  Keyword found: for keywords requiring an argument, convert
!  next token and assign appropriate value
!----------------------------------------------------------------------
    16 call gettkn(line,token,lth)
!                                        *** No value given
    if (lth == 0) then
        write(luttyo,1003) keywrd(i)
        return
    end if

    if (ftoken(token,lth,fval)) then
    !                                          *** TOL keyword
        if (i == 1) then
            stol=fval
        !                                          *** STEP keyword
        else if (i == 2) then
            sstep=fval
        !                                          *** MIN keyword
        else if (i == 3) then
            smn=fval
        !                                          *** MIN keyword
        else if (i == 4) then
            smx=fval
        !                                          *** START keyword
        else if (i == 5) then
            ibrckt=1
            stval=fval
        !                                      *** Illegal numeric value
        end if
    else
        write(luttyo,1001) token(:lth)
    end if
    go to 14

! ##### Format statements ########################################

    1000 format('*** Unrecognized SEARCH keyword: ''',a,''' ***')
    1001 format('*** Numeric value expected: ''',a,''' ***')
    1002 format('*** ''',a,''' is not a variable parameter ***')
    1003 format('*** No value given for ''',a,''' ***')
    1040 format('*** Initial parameter values are required before', &
    ' a fit ***')
    1100 format(/4x,'Initial Range  : ',g13.6,2x,g13.6,2x,g13.6)
    1105 format(4x,'Search Bracket : ',g13.6,2x,g13.6,2x,g13.6, &
    /4x,'Rms Deviation  : ',g13.6,2x,g13.6,2x,g13.6)
    1110 format('*** No minimum found within the valid range ***')
    1120 format(4x,'Last Bracket   : ',g13.6,2x,g13.6, &
    /4x,'Rms Deviation  : ',g13.6,2x,g13.6)
    end subroutine srchc
