!  VERSION 1.0  (NLSPMC version)   2/5/99
!**********************************************************************
!                    =========================
!                          file CMDS
!                    =========================

!  Contains the following routines:
!    varyc(line)    Interprets given line as an NLS "vary" command
!    fixc(line)     Interprets given line as an NLS "fix" command
!    basis(line)    Interprets given line as an NLS "basis" command
!    series(line)   Interprets given line as an NLS "series" command
!    fitc(line)     Interprets given line as an NLS "fit" command
!    parc(line,lu)  Interprets line as NLS "parms" command and prints
!                    parameter values to the logical unit specified
!    statc(line,lu) Interprets line as NLS "status" command

!    Includes
!      limits.inc
!      simparm.inc
!      datas.inc
!      parcom.inc
!      lpnam.inc
!      lmcomm.inc
!      stdio.inc

!    Uses
!      ipfind
!      gettkn
!      ftoken
!      itoken

!----------------------------------------------------------------------
!                    =========================
!                      subroutine VARYC
!                    =========================

! vary <name> {  {sprectrum,site} minimum <minval> maximum <maxval>
!                scale <scl> fdstep <step> }

!      name    : name of the parameter to be varied
!      sprectrum,site : parameters or * or absent, identifying same.
!      minval  : Minimum permissible value for parameter
!      maxval  : Maximum permissible value for paramete
!      scale   : Factor by which to scale search vector for this parameter
!      step    : Relative step size for forward-differences approximation

!----------------------------------------------------------------------
    subroutine varyc(line)
    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'parmequ.inc'
    include 'parms.inc'
    include 'lpnam.inc'
    include 'stdio.inc'
    include 'miscel.inc'


    character*(LINELG) line
    integer :: i,ibd,ix,lth,specinfo,siteinfo
    double precision :: fval,prmn,prmx,prsc,step
    character token*(WORDLG),prmID*6

    integer :: nkeywd
    parameter(nkeywd=4)
    character(8) :: keywrd(nkeywd)
    data keywrd /'MINIMUM','MAXIMUM','SCALE','FDSTEP'/

    double precision :: one,zero
    parameter (one=1.0d0,zero=0.0d0)

    integer :: varspec
    logical :: ftoken
    external ftoken

!----------------------------------------------------------------------
! Get the name of the parameter
!----------------------------------------------------------------------
    call gettkn(line,token,lth)
    lth=min(lth,6)
    call touppr(token,lth)

    call ipfind(ix,token,lth,varspec)

!----------------------------------------------------------------------
!     Check whether parameter may be varied
!----------------------------------------------------------------------
    if (ix == 0 .OR. ix > nvprm) then
        write(luttyo,1002) token(:lth)
        return
    end if

    if (ix < -100) then
        prmID=alias2( -99-(IGXX+ix) )
    else if (ix < 0) then
        prmID=alias1( 1-(IGXX+ix) )
    else
        prmID=parnam(ix)
    end if

! get parenthisized information

    call  indtkn(line,specinfo,siteinfo)

! set defaults

    prmn=zero
    prmx=zero
    prsc=one
    step=1.0D-5
    ibd=0

!----------------------------------------------------------------------
!  Look for a keyword
!----------------------------------------------------------------------
    14 call gettkn(line,token,lth)
    lth=min(lth,8)

!------------------------------------------------
! No more keywords: add the parameter and exit
!------------------------------------------------
    if (lth == 0) then
    !         call addprm(ix,ibd,prmn,prmx,prsc,step,prmID)
        if (specinfo /= 0 .AND. nspectra > 1) then
            write(luttyo,*)'Vary parameter must apply to all spectra.'
            write(luttyo,*)'but may apply to only one site if you wish.'
            return
        end if
    
        call addprm(ix,ibd,prmn,prmx,prsc,step,prmID,specinfo,siteinfo)
        return
    end if

!------------------------------
! Check keyword list
!------------------------------
    call touppr(token,lth)
    do 15 i=1,nkeywd
        if (token(:lth) == keywrd(i)(:lth)) goto 16
    15 END DO
!                                        *** Unrecognized keyword
    write (luttyo,1000) token(:lth)
    return

!----------------------------------------------------------------------
!  Keyword found: convert next token and assign appropriate value
!----------------------------------------------------------------------
    16 call gettkn(line,token,lth)
!                                        *** No value given
    if (lth == 0) then
        write(luttyo,1003) keywrd(i)
        return
    end if

    if (ftoken(token,lth,fval)) then
    !                                        *** MINIMUM keyword
        if (i == 1) then
            prmn=fval
            if (mod(ibd,2) == 0) ibd=ibd+1
        !                                        *** MAXIMUM keyword
        else if (i == 2) then
            prmx=fval
            ibd=ibd+2
        !                                        *** SCALE keyword
        else if (i == 3) then
            prsc=fval
        !                                        *** FDSTEP keyword
        else if (i == 4) then
            step=fval
        end if
    !                               *** Illegal real value
    else
        write(luttyo,1001) token(:lth)
    end if

    go to 14

! ###### format statements ########################################

    1000 format('*** Unrecognized VARY keyword: ''',a,''' ***')
    1001 format('*** Real value expected: ''',a,''' ***')
    1002 format('*** ''',a,''' is not a variable parameter ***')
    1003 format('*** No value given for ''',a,''' ***')

    end subroutine varyc


!----------------------------------------------------------------------
!                    =========================
!                      subroutine FIXC
!                    =========================

!  Removes the given parameter name from the list of variable parameters.

! fix <name>

!      name    : name of the parameter to be removed
!----------------------------------------------------------------------
    subroutine fixc(line)
    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'parmequ.inc'
    include 'parms.inc'
!      include 'lpnam.inc'
    include 'stdio.inc'
    include 'miscel.inc'

    character*(LINELG) line
    integer :: i,j,lth,specinfo,siteinfo,isite,ispec,specptr,siteptr
    character token*(LINELG),prmID*6

    integer :: ixr,ixabs,iprm,varspec,ident
    logical :: ftoken,okparm
    external ftoken

!----------------------------------------------------------------------
! Get the name of the parameter
!----------------------------------------------------------------------
    call gettkn(line,token,lth)
    lth=min(lth,6)
    if (lth <= 0) return

    call touppr(token,lth)

    if (token(:lth) == 'ALL') then
    !                                           ** Fix all parameters
        nprm=0
        do 1 iprm=1,nvprm
            do 1 isite=1,ncomps
                do 1 ispec=1,nspectra
                ! zero parameter being fixed, decrement location for moved parms
                    ixx(iprm,ispec,isite)=0
        1 END DO
        return
    end if

    call ipfind(ixr,token,lth,varspec)
    ixabs=abs(mod(ixr,100))	! token id without axis info

!----------------------------------------------------------------------
!     Check whether parameter may be varied
!----------------------------------------------------------------------
    if (ixabs == 0 .OR. ixabs > nvprm) then
        write(luttyo,1002) token(:lth)
        return
    end if

! get parenthisized information

    call  indtkn(line,specinfo,siteinfo)

! fix all relevent sites

    okparm=.false.
    specptr=1                 ! get the old information
    siteptr=1
    if (specinfo /= 0) specptr=specinfo         ! get first spec
    if (siteinfo /= 0) siteptr=siteinfo         ! get first site
! consider all sites, spec:
    do 30 isite=1,ncomps
        do 30 ispec=1,nspectra
        ! check matching sites with requested site
            if ((isite == siteinfo .OR. siteinfo == 0) .AND. &
            (ispec == specinfo .OR. specinfo == 0)) then
            ! if vary previously selected
                if (ixx(ixabs,ispec,isite) /= 0) then
                    ident=ixx(ixabs,ispec,isite)
                    okparm=.true.
                    call rmvprm(ident)
                ! all sites/spectra should be fixed with the first call to rmvprm.
                end if
            end if
    30 END DO
    if ( .NOT. okparm) then
        if (specinfo == 0 .AND. siteinfo == 0) then
            write(luttyo,1003) token(:lth)
        else if (specinfo == 0) then
            write(luttyo,1004) token(:lth),siteinfo
        else if (siteinfo == 0) then
            write(luttyo,1005) token(:lth),specinfo
        else
            write(luttyo,1006) token(:lth),specinfo,siteinfo
        end if
    end if
    return

! ###### format statements ########################################

    1002 format('*** ''',a,''' is not a variable parameter ***')
    1003 format('*** ''',a,'(*,*)',''' is not a variable parameter ***')
    1004 format('*** ''',a,' *,',i1,''' is not a variable parameter ***')
    1005 format('*** ''',a,i2,',*',''' is not a variable parameter ***')
    1006 format('*** ''',a,'(',i2,',',i1,')', &
    ''' is not a variable parm ***')

    end subroutine fixc


!----------------------------------------------------------------------
!                    =========================
!                      subroutine BASIS
!                    =========================

! Generates a basis set to be used in simulation.
! The program just reads the index file if it exists.
! If the index file does not exist or the filename is not specified,
! the full basis vectors within the maximum allowed indices (lemx,
! lomx,kmx,mmx,ipnmx) are generated.  The eigenvector calculation
! with the full basis vectors is very costly especially with the
! Rutishauser algorithm, and the pruned basis set should be used
! as much as possible.

! basis {<name> {site,spectrum} print }

!      name    : name (in full) of the basis index file
!      print   : prints out the indices on the screen

!----------------------------------------------------------------------
    subroutine basis(line,lu)
    implicit none

    include 'limits.inc'
    include 'names.inc'
    include 'simparm.inc'
    include 'parms.inc'
    include 'parmequ.inc'
    include 'basis.inc'
    include 'miscel.inc'
    include 'stdio.inc'
    character*(LINELG) line

    integer :: i,n,lth,lu,ierr
    character token*(WORDLG)

    integer :: ipar,iprt,specinfo,siteinfo,isite,ispec,specptr,siteptr
    external ipar

    iprt=0

! get parenthisized information, if any:

    call  indtkn(line,specinfo,siteinfo)

! store basis id # in basinfo

    if (nbasis < nspectra*ncomps) then
        nbasis=nbasis+1
    else		! Too many basis sets
        if (lu /= 0) write (lu,1000)
        if (lu /= luttyo) write (luttyo,1000)
        return
    end if

    specptr=1                 ! get the old information
    siteptr=1
    if (specinfo /= 0) specptr=specinfo         ! get first spec
    if (siteinfo /= 0) siteptr=siteinfo         ! get first site

! Set up lemx etc for calculation from this site:
! (Presumes this has been entered before the basis command)

    do 28 i=1,NFPRM
        fepr(i)=fparm(i,specptr,siteptr)
    28 END DO
    do 4 i=1,NIPRM
        iepr(i)=iparm(i,specptr,siteptr)
    4 END DO

!----------------------------------------------------------------------
! Get the name of the basis set if any:
!----------------------------------------------------------------------

    10 call gettkn(line,ixname,lth)

    if (lth /= 0) then
        if ((ixname(:lth) == 'print') .OR. (ixname(:lth) == 'PRINT')) then
            iprt=1
            go to 10
        end if
    end if
! set up mjqe1 etc:
    if (lth /= 0) then
        call fbasis(ixname,0,ierr)	! read file
    else
        call fbasis(ixname,1,ierr)	! full basis
        if (lu /= 0) write (lu,1001)
        if (lu /= luttyo) write (luttyo,1001)
    end if

    if (ierr /= 0) then
        if (ierr == 1) then
            if (lu /= 0) then
                write (lu,1002) ixname(:lth)	! file not exist
                write (lu,1001)
            end if
            if (lu /= luttyo) then
                write (luttyo,1002) ixname(:lth)
                write (luttyo,1001)
            end if
        else if (ierr == 2) then 	! dimension exceeded
            if (lu /= 0) then
                write (lu,1003) mxdim
                write (lu,1100) ndimoss(nbasis),ndimdss(nbasis)
            end if
            if (lu /= luttyo) then
                write (luttyo,1003) mxdim
                write (luttyo,1100) ndimoss(nbasis),ndimdss(nbasis)
            end if
            return
        else if (ierr == 3) then	! mag parms not set
            if (lu /= 0) write (lu,1004)
            if (lu /= luttyo) write (luttyo,1004)
            return
        end if
    end if

!----------------------------------------------------------------------
!     Obtain MTS index for this basis
!----------------------------------------------------------------------

    lemx=-1
    lomx=-1
    kmx=-1
    mmx=-1
    ipnmx=-1
! loop over this basis set:
    do 210 i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1
        n=ml1(i)
        if (ipar(n) == 1) then
            if (n > lemx) lemx=n
        else
            if (n > lomx) lomx=n
        end if
        n=mk1(i)
        if (n > kmx) kmx=n
        n=abs(mm1(i))
        if (n > mmx) mmx=n
        n=abs(mpi1(i))
        if (n > ipnmx) ipnmx=n
    210 END DO

    if (lomx == -1) lomx=0
    if (kmx == -1) kmx=0
    if (mmx == -1) mmx=0
    if (ipnmx == -1) ipnmx=0

    write (lu,1100) ndimoss(nbasis),ndimdss(nbasis)
    write (lu,1110) lemx,lomx,kmx,mmx,ipnmx
    if (lu /= luttyo) then
        write (luttyo,1100) ndimoss(nbasis),ndimdss(nbasis)
        write (luttyo,1110) lemx,lomx,kmx,mmx,ipnmx
    end if

! store the lemx etc in relevent sites/spectra

    do 40 isite=1,ncomps
        do 40 ispec=1,nspectra
        ! Fill matching sites with correct contents:
            if ((isite == siteinfo .OR. siteinfo == 0) .AND. &
            (ispec == specinfo .OR. specinfo == 0)) then
                iparm(ILEMX,ispec,isite)=lemx
                iparm(ILEMX+1,ispec,isite)=lomx
                iparm(ILEMX+2,ispec,isite)=kmx
                iparm(ILEMX+3,ispec,isite)=mmx
                iparm(ILEMX+4,ispec,isite)=ipnmx
                matrxok(isite)=0	! need new matrix!
            ! store basis info for this site/spec
                basinfo(1,ispec,isite)=nbasis
            ! store basis size information
                iparm(INDIMO,ispec,isite)=ndimoss(nbasis)
                iparm(INDIMD,ispec,isite)=ndimdss(nbasis)
            
            end if
    40 END DO

    call gettkn(line,token,lth)
    if ((iprt /= 1) .AND. (lth <= 0)) return

    write(lu,1200)
    do 100 i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1
        write(lu,1210) i-pidptr(nbasis)+1,ml1(i),mjk1(i),mk1(i), &
        mjm1(i),mm1(i),mpi1(i),mqi1(i)
    100 END DO

!      write(lu,1220)
!      do 110 i=dpidptr(nbasis),dpidptr(nbasis)+ndimdss(nbasis)-1
!         write(lu,1210) i-dpidptr(nbasis)+1,mdl1(i),mdjk1(i),mdk1(i),
!     #                  mdjqe1(i),mdm1(i),mdpi1(i),mdqi1(i)
! 110  continue

    return

! ###### format statements ############################################

    1000 format('** Too many basis sets requested **')
    1001 format('** Full basis set is generated **')
    1002 format(' Index file : ',a,' does not exist')
    1003 format(' Dimension of matrix too large : maximum dimension = ', &
    i4)
    1004 format(' magnetic parameters are not properly set')
    1100 format(' Dimension of matrix: ndimoss = ',i4,' ndimdss = ',i4)
    1110 format(' MTS lemx,lomx,kmx,mmx,ipnmx = ',i3,4(',',i2))
    1200 format(/,15x,'BASIS SET : off-diagonal space',//,2x, &
    'element',7x,'L  jK   K  jqM  M  pI  qI',/,2x,42('-'))
    1210 format(4x,i5,4x,'|',6(i3,','),i3,'  >')
    1220 format(/,15x,'BASIS SET : diagonal space',//,2x,'element',7x, &
    'L  jK   K  jqM  M  pI  qI',/,2x,42('-'))

    end subroutine basis


!----------------------------------------------------------------------
!                    =========================
!                       subroutine SERIES
!                    =========================

! Interpret a "series" command from the given line.
! The command is of the form:

!    series <name> {=} <value>, <value>, ...

! For now allowed series parameters apply to all sites.

!         name  : Name of a variable parameter
!         value : List of values assumed by the named parameter
!                 for each spectrum in the series
!    Note that nspectra (number of series calculations)
!    values are expected, if fewer, the last is repeated.

!----------------------------------------------------------------------
    subroutine series( line )
    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'miscel.inc'
    include 'stdio.inc'

    character line*(LINELG)
    logical :: serprm
    integer :: itmp,lth,lprm,ix,i
    double precision :: fval,serval(MXSPEC)
    character token*(WORDLG),prmID*(6),defnm*(6)

    logical :: ftoken
    external ftoken

!----------------------------------------------------------------------
!     Get name of series parameter
!----------------------------------------------------------------------
    call gettkn(line,prmID,lth)
    lth=min(lth,6)
!                               *** No series variable specified
    if (lth == 0) then
        write(luttyo,1002)
        return
    end if

!-----------------------------------------------------
! Check if the specified parameter is series variable.
!-----------------------------------------------------
    call touppr(prmID,lth)

    serprm = (prmID(:lth).eq.'IEXP').or.(prmID(:lth).eq.'ICOMB') &
    .or.(prmID(:lth).eq.'NPT1').or.(prmID(:lth).eq.'NPT2') &
    .or.(prmID(:lth).eq.'INIT1').or.(prmID(:lth).eq.'INIT2') &
    .or.(prmID(:lth).eq.'STEPT1').or.(prmID(:lth).eq.'STEPT2') &
    .or.(prmID(:lth).eq.'TFIX').or.(prmID(:lth).eq.'WEIGHT')

    lprm=lth
    if ( .NOT. serprm) then
        write (luttyo,1003) prmID(:lth)
        return
    end if

    itmp=0

!----------------------------------------------------------------------
!         Get a list of values for the series parameter
!----------------------------------------------------------------------
    10 call gettkn(line,token,lth)
    if (lth /= 0) then		! have a value or =
    
        if (token(:lth) == '=') goto 10
    
        if (ftoken(token,lth,fval)) then	! if is a fp #
            if (itmp >= nspectra) then
                write (luttyo,1005) nspectra
                1005 format('*** SERIES may not have more than',i2, &
                ' values: did you set nspectra? ***')

                return
            end if
            itmp=itmp+1
            serval(itmp)=fval
        !                                *** Illegal real number
        else
            write(luttyo,1006) token(:lth)
            1006 format('*** ',a,' : Not a real value ***')
        end if
    
        go to 10
    else
    
    !------------------------------
    ! No more values--map them.
    !------------------------------
    !                                   ** Set nser here **
    ! if this is the first call to series, set nser
        if (nsparm == 0) nser=nspectra
    
    ! if later calls, fill in if list is short
        if (itmp < nser) then
            write(luttyo,1004)
            do 11 ix=itmp+1,nspectra
                serval(ix)=serval(itmp)
            11 END DO
        end if
    
    !--------------------------------------------------
    ! Store current series variables into proper array
    !--------------------------------------------------
    
        if (prmID(:lprm) == 'IEXP') then
            do 20 ix=1,nser
            ! diaflg indicates ANY of the experiments are type greater than
            ! or equal to 3.
                if (serval(ix) >= 3) diaflg= .TRUE. 
                siexp(ix)=nint(serval(ix))
            20 END DO
        else if (prmID(:lprm) == 'ICOMB') then
            do 21 ix=1,nser
                sicomb(ix)=nint(serval(ix))
            21 END DO
        else if (prmID(:lprm) == 'NPT1') then
            do 22 ix=1,nser
                snpt1(ix)=nint(serval(ix))
            22 END DO
        else if (prmID(:lprm) == 'NPT2') then
            do 23 ix=1,nser
                snpt2(ix)=nint(serval(ix))
            23 END DO
        else if (prmID(:lprm) == 'INIT1') then
            do 24 ix=1,nser
                sinit1(ix)=serval(ix)
            24 END DO
        else if (prmID(:lprm) == 'INIT2') then
            do 25 ix=1,nser
                sinit2(ix)=serval(ix)
            25 END DO
        else if (prmID(:lprm) == 'STEPT1') then
            do 26 ix=1,nser
                sstept1(ix)=serval(ix)
            26 END DO
        else if (prmID(:lprm) == 'STEPT2') then
            do 27 ix=1,nser
                sstept2(ix)=serval(ix)
            27 END DO
        else if (prmID(:lprm) == 'TFIX') then
            do 28 ix=1,nser
                stfix(ix)=serval(ix)
            28 END DO
        else if (prmID(:lprm) == 'WEIGHT') then
            do 29 ix=1,nser
                datwgt(ix)=serval(ix)
            29 END DO
        end if
    !                * update # of series parameters set so far *
        nsparm=nsparm+1
    !                * initialize some arrays, after have snpt1, snpt2
    ! This assumes we will be reading experimental data sets.  This is not
    ! necessary in the case of no data, but we don't know that at this time.
    
        if (nsparm >= 9) then		! got all series parms
            ix=1
            defnm='simu'
            do 30 i=1,nspectra
                ixsp(i)=ix
                ndata(i)=snpt1(i)*snpt2(i)
                ix=ix+ndata(i)
                if (i < 10)write(defnm(5:5),'(i1)') i
                if (i > 10)write(defnm(5:6),'(i2)') i
                dataid(i)=defnm
            30 END DO
        end if
    
        return
    
    end if


! ###### format statements ######################################

    1001 format('*** nser NOT set: nser=1 assumed ***')
    1002 format('*** No SERIES variable specified ***')
    1003 format('*** ',a,' : Not a SERIES variable ***')
    1004 format('*** Not enough SERIES values : last repeated ***')

    end subroutine series


!----------------------------------------------------------------------
!                    =========================
!                       subroutine FITC
!                    =========================

!   fit      { trace xtol <xtol> ftol <ftol> gtol <ftol>
!              maxfun <mxf> maxitr <mxi> maxfit <mxfit> bound <factor> }

!          trace:  Specifies that a ".trc" file should be produced for
!                  the fit
!          xtol:   Convergence tolerance for scaled fitting parameters
!          ftol:   Convergence tolerance for chi-squared
!          gtol:   Convergence tolerance for gradient of chi-squared with
!                  respect to the fitting parameters
!          mxf:    Maximum number of function calls allowed
!          mxi:    Maximum number of iterations allowed
!          mxfit:  Maximum number of fit commands allowed until convergence
!          factor: Factor defining initial step bound used in parameter search

!----------------------------------------------------------------------
    subroutine fitc( line )
    implicit none

    include 'limits.inc'
    include 'lmcomm.inc'
    include 'parms.inc'
!      include 'iterat.inc'
    include 'miscel.inc'
    include 'stdio.inc'

    character line*(LINELG)

    logical :: ftoken
    external ftoken
    integer :: i,lth,ifit,ispec,isite
    double precision :: fval
    character*(WORDLG) token

    integer :: nkeywd
    parameter(nkeywd=9)

    character(8) :: keywrd(nkeywd)
    data keywrd / 'FTOL', 'GTOL', 'XTOL', 'BOUND', 'MAXFUN', &
    'MAXITR', 'MAXFIT', 'TRACE', 'JACOBI' /

    logical :: prmsOK
    external prmsOK

!######################################################################

!  -- Check whether any starting parameters have been input before
!     proceeding with the FIT command

    if ( .NOT. prmsOK(luttyo) ) then
        write(luttyo,1040)
        1040 format('*** Initial parameter values are required before', &
        ' a fit ***')
        return
    end if
!  check basis set settings:
    do 1 isite=1,ncomps
        do 1 ispec=1,nspectra
            if (basinfo(1,1,isite) /= basinfo(1,ispec,isite)) then
                write(luttyo,*)'Error, multiple basis sets not allowed'
                write(luttyo,*)'except for different sites.'
                return
            end if
    1 END DO

    ifit=0
!----------------------------------------------------------------------
!  Look for a keyword
!----------------------------------------------------------------------

    14 call gettkn(line,token,lth)
    lth=min(lth,8)

!   ***** THIS IS IT!!! ************
!------------------------------------------------
! No more keywords: call the NLS fitting routine
!------------------------------------------------
    if (lth == 0) then
    
    !        call catchc( ihltcmd )
    
        15 call fitp
    
        if (ihltcmd /= 0) then
        !          call uncatchc( ihltcmd )
            return
        else
            ifit=ifit+1
            if ((ifit < maxfit) .AND. ((info == 0) .OR. (info > 3))) &
            go to 15
        end if
    
        return
    end if

!------------------------------
! Check keyword list
!------------------------------
    call touppr(token,lth)
    do 16 i=1,nkeywd
        if (token(:lth) == keywrd(i)(:lth)) goto 17
    16 END DO
!                                        *** Unrecognized keyword
    write (luttyo,1000) token(:lth)
    return

!----------------------------------------------------------------------
!  Keyword found: for keywords requiring an argument, convert
!  next token and assign appropriate value
!----------------------------------------------------------------------
    17 if (i >= 1 .AND. i <= 7) then
        call gettkn(line,token,lth)
    !                                        *** No value given
        if (lth == 0) then
            write(luttyo,1003) keywrd(i)
            return
        end if
    
        if (ftoken(token,lth,fval)) then
        !                                          *** FTOL keyword
            if (i == 1) then
                ftol=fval
            !                                          *** GTOL keyword
            else if (i == 2) then
                gtol=fval
            !                                          *** XTOL keyword
            else if (i == 3) then
                xtol=fval
            !                                          *** BOUND keyword
            else if (i == 4) then
                factor=fval
            !                                          *** MAXFUN keyword
            else if (i == 5) then
                maxev=int(fval)
            !                                          *** MAXITR keyword
            else if (i == 6) then
                maxitr=int(fval)
            !                                          *** MAXFIT keyword
            else if (i == 7) then
                maxfit=int(fval)
            end if
        !                                          *** Illegal numeric value
        else
            write(luttyo,1001) token(:lth)
        end if
    !                                          *** TRACE keyword
    else if (i == 8) then
        if (luout == luttyo) then
            write (luttyo,1050)
            itrace=0
        else
            itrace=1
        end if
    !                                          *** JACOBI keyword
    else if (i == 9) then
        jacobi=1
    end if

    go to 14

    1000 format('*** Unrecognized FIT keyword: ''',a,''' ***')
    1001 format('*** Numeric value expected: ''',a,''' ***')
    1003 format('*** No value given for ''',a,''' ***')
    1050 format('*** A log file must be opened before using TRACE ***')
    end subroutine fitc


!----------------------------------------------------------------------
!                    =========================
!                       subroutine PARC
!                    =========================

!  Prints out a list of parameter values on the given logical unit
!  number.
!  Add lulog to do parms conversion and output parms into log file.
!----------------------------------------------------------------------
    subroutine parc( line )
    implicit none

    include 'limits.inc'
    include 'names.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'stdio.inc'
    include 'miscel.inc'
    include 'rndoff.inc'

    character line*(LINELG),hdr*7,shdr*4,fileID*(WORDLG), &
    flname*(WORDLG)
    integer :: i,ioerr,jx,lth,lu,npotn
    double precision :: totwgt

    character(2) :: nonbr(3)
    data nonbr/ 'l', 'xy', 'zz'/

! parms outputs varaiables if they are spec/site = 1/1 or if they vary
! from those valuse for other spectra and sites.  This test may still
! output copies that are equal but it should limit the output.

    integer :: idxx,isp,istx,ispec,isite
    logical :: test
    test(idxx,isp,istx)=((isp.eq.1 .and. istx.eq.1) .or. ((isp.eq.1) &
    .and. (abs(fparm(idxx,isp,istx)-fparm(idxx,1,1)) .gt. &
    0.00001*fparm(idxx,1,1))) .or. ((isp.ne.1) .and. &
    (abs(fparm(idxx,isp,istx)-fparm(idxx,1,istx)) .gt. &
    0.00001*fparm(idxx,1,istx))))

!......................................................................

    call gettkn(line,fileID,lth)

!----------------------------------------------------------------------
!     No name specified: output to terminal
!----------------------------------------------------------------------
    if (lth == 0) then
        lu=luttyo
        hdr=' '
        shdr=' '
    else
    
    !----------------------------------------------------------------------
    !     Set for output of "let" commands to specified file
    !----------------------------------------------------------------------
        flname=fileID(:lth)//'.run'
        open(ludisk,file=flname,status='unknown', &
        access='sequential',form='formatted',iostat=ioerr)
        if (ioerr /= 0) then
            write (luout,3000) flname
            return
        end if
        lu=ludisk
    ! the idea of writing a file that can be read does not work
    !  for multiple spectra/sites.
        hdr='let '
        shdr='c '
        write (lu,1022) fileID(:lth)
        write (lu,1024) ixname
    end if

! Loop over spectra and sites, outputting only when information
! differes from spectra 1

    write(lulog,4007)
    write(lulog,4008)
    do 10 ispec=1,nspectra
        do 10 isite=1,ncomps
            write(lu,101)ispec,isite
            101 format(/,' For Spectrum ',i3,' Site ',i3,/)
            write(lulog,101)ispec,isite
            if(test(IB0,ispec,isite)) &
            write (lu,1020) hdr,fparm(IB0,ispec,isite)
        !----------------------------------------------------------------------
        !     g-tensor
        !----------------------------------------------------------------------
            if(test(IGXX,ispec,isite)) then
                if (iparm(IIGFLG,ispec,isite) == 3) then
                    write (lu,2001) hdr,(fparm(IGXX+i,ispec,isite),i=0,2)
                else if (iparm(IIGFLG,ispec,isite) == 2) then
                    write (lu,1001) hdr,(fparm(IGXX+i,ispec,isite),i=0,2)
                else
                    write (lu,1000) hdr,(fparm(IGXX+i,ispec,isite),i=0,2)
                end if
            end if
        
        !----------------------------------------------------------------------
        !     Nuclear spin and hyperfine tensor
        !----------------------------------------------------------------------
            if(test(IAXX,ispec,isite)) then
                write (lu,1002) hdr,iparm(IIN2,ispec,isite)
            end if
        
            if (iparm(IIN2,ispec,isite) /= 0) then
                if(test(IAXX,ispec,isite)) then
                    if (iparm(IIGFLG+1,ispec,isite) == 3) then
                        write(lu,2004) hdr,(fparm(IAXX+i,ispec,isite),i=0,2)
                    else if (iparm(IIGFLG+1,ispec,isite) == 2) then
                        write(lu,1004) hdr,(fparm(IAXX+i,ispec,isite),i=0,2)
                    else
                        write(lu,1003) hdr,(fparm(IAXX+i,ispec,isite),i=0,2)
                    end if
                end if
            end if
        
        !----------------------------------------------------------------------
        !     Diffusion tensor
        !----------------------------------------------------------------------
            if(test(IDX,ispec,isite)) then
                if (iparm(IIGFLG+2,ispec,isite) == 3) then
                    write(lu,2006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
                    write(lulog,2006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
                    write(lulog,4001) hdr,(10.d0**fparm(IDX+i,ispec,isite) &
                    ,i=0,2)
                else if (iparm(IIGFLG+2,ispec,isite) == 2) then
                    write(lu,1006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
                    write(lulog,1006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
                    write(lulog,4003) hdr,(10.d0**fparm(IDX+i,ispec,isite) &
                    ,i=0,2)
                else
                    write(lu,1005) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
                    write(lulog,1005) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
                    write(lulog,4002) hdr,(10.d0**fparm(IDX+i,ispec,isite) &
                    ,i=0,2)
                end if
            end if
        
        !----------------------------------------------------------------------
        !     Non-Brownian rotational model parameters
        !----------------------------------------------------------------------
            if(test(IDJF,ispec,isite) .OR. test(IPML,ispec,isite)) then
                write (lu,1007) hdr,iparm(IIPDF,ispec,isite)
                if (iparm(IIPDF,ispec,isite) == 1) then
                    do 2 i=0,2
                        write(lu,1008) hdr,nonbr(i+1),nonbr(i+1), &
                        iparm(IML+i,ispec,isite),fparm(IPML+i,ispec,isite)
                    2 END DO
                
                !   -- Anisotropic viscosity
                
                else if (iparm(IIPDF,ispec,isite) == 2) then
                    write (lu,1009) hdr,fparm(IDJF,ispec,isite), &
                    fparm(IDJFPRP,ispec,isite)
                end if
            
            !   --- Discrete jump motion
                if (iparm(IIPDF,ispec,isite) /= 2 .AND. &
                iparm(IIST,ispec,isite) /= 0 .AND. &
                npotn > 0) write (lu,1010) hdr,iparm(IIST,ispec,isite), &
                fparm(IDJF,ispec,isite)
            end if
        
        !----------------------------------------------------------------------
        !   Orienting potential
        !----------------------------------------------------------------------
            if(test(IC20,ispec,isite) .OR. test(IC20+1,ispec,isite) &
             .OR. test(IPSI,ispec,isite)) then
                npotn=0
                do 3 i=0,4
                    if (dabs(fparm(IC20+i,ispec,isite)) > rndoff) npotn=i+1
                3 END DO
                if (npotn > 0) then
                    write (lu,1011) hdr,iparm(INORT,ispec,isite), &
                    fparm(IPSI,ispec,isite)
                    write (lu,1012) hdr,(fparm(IC20+i,ispec,isite),i=0,4)
                    write (lulog,1011) hdr,iparm(INORT,ispec,isite), &
                    fparm(IPSI,ispec,isite)
                    write (lulog,1012) hdr,(fparm(IC20+i,ispec,isite),i=0,4)
                end if
            end if
        
        !   --- Heisenberg exchange
        
            if(test(IOSS,ispec,isite)) then
                if (dabs(fparm(IOSS,ispec,isite)) > rndoff) then
                    write (lu,1013) hdr,fparm(IOSS,ispec,isite)
                    write (lulog,1013) hdr,fparm(IOSS,ispec,isite)
                    write (lulog,4004) hdr,10.d0**fparm(IOSS,ispec,isite)
                end if
            end if
        
        !----------------------------------------------------------------------
        !     Diffusion tilt & molecular tilt
        !----------------------------------------------------------------------
            if(test(IALD,ispec,isite) .OR. test(IALM,ispec,isite)) then
                if ((dabs(fparm(IALD,ispec,isite)) > rndoff) .OR. &
                (dabs(fparm(IALD+1,ispec,isite)) > rndoff) .OR. &
                (dabs(fparm(IALD+2,ispec,isite)) > rndoff)) &
                write (lu,1014) hdr,(fparm(IALD+i,ispec,isite),i=0,2)
            
                if ((dabs(fparm(IALM,ispec,isite)) > rndoff) .OR. &
                (dabs(fparm(IALM+1,ispec,isite)) > rndoff) .OR. &
                (dabs(fparm(IALM+2,ispec,isite)) > rndoff)) &
                write (lu,1015) hdr,(fparm(IALM+i,ispec,isite),i=0,2)
            end if
        
        !----------------------------------------------------------------------
        !  Relaxation constants
        !----------------------------------------------------------------------
            if(test(IT2EDI,ispec,isite) .OR. test(IT2EDI+1,ispec,isite) .OR. &
            test(IT2EDI+2,ispec,isite) .OR. test(IT2EDI+3,ispec,isite) &
             .OR. test(IT2EDI+4,ispec,isite)) then
                write (lu,1016) hdr,(fparm(IT2EDI+i,ispec,isite),i=0,4)
                write (lulog,1016) hdr,(fparm(IT2EDI+i,ispec,isite),i=0,4)
                write (lulog,4005) hdr,(10.d0**fparm(IT2EDI+i,ispec, &
                isite),i=0,3,3)
                write (lulog,4006) hdr,(1.d0/(10.d0**fparm(IT2EDI+i, &
                ispec,isite))*1.d9,i=0,3,3)
            end if
        
        !----------------------------------------------------------------------
        !  Siteweight
        !----------------------------------------------------------------------
            if(test(ISWGT,ispec,isite))then
                totwgt=0.0d0
                do 31 i=1,ncomps
                    totwgt=totwgt+fparm(ISWGT,ispec,i)
                31 END DO
            !c
                write(lu,1116) hdr,fparm(ISWGT,ispec,isite)/totwgt
                write(lulog,1116) hdr,fparm(ISWGT,ispec,isite)/totwgt
            end if
        
        !----------------------------------------------------------------------
        !  Broadening width
        !----------------------------------------------------------------------
            if(test(IMWID,ispec,isite) .OR. test(IGIB,ispec,isite) .OR. &
            test(ILIB,ispec,isite) .OR. test(IHWID,ispec,isite)) then
                write (lu,1017) hdr,fparm(IMWID,ispec,isite), &
                fparm(IGIB,ispec,isite),fparm(ILIB,ispec,isite), &
                fparm(IHWID,ispec,isite)
                write (lulog,1017) hdr,fparm(IMWID,ispec,isite), &
                fparm(IGIB,ispec,isite),fparm(ILIB,ispec,isite), &
                fparm(IHWID,ispec,isite)
            end if
        
        !----------------------------------------------------------------------
        !  Basis set and pruning tolerance
        !----------------------------------------------------------------------
        
            if(test(ICGTOL,ispec,isite) .AND. nbasis /= 0) then
                write (lu,1018) hdr,(iparm(ILEMX+i,ispec,isite),i=0,4), &
                hdr,ndimoss(basinfo(1,ispec,isite)), &
                ndimdss(basinfo(1,ispec,isite))
                write (lu,1019) hdr,fparm(ICGTOL,ispec,isite), &
                iparm(INSTEP,ispec,isite)
            end if
        
        !----------------------------------------------------------------------
        !  Series parameters
        !----------------------------------------------------------------------
        
            if(ispec == 1 .AND. isite == 1) then  ! parms apply to all sites
                if (lu == ludisk) hdr='series '
                write (lu,1030) shdr,nser
                write (lu,1031) hdr,(siexp(i),i=1,nser)
                write (lu,1032) hdr,(sicomb(i),i=1,nser)
                write (lu,1033) hdr,(snpt1(i),i=1,nser)
                write (lu,1034) hdr,(snpt2(i),i=1,nser)
                write (lu,1035) hdr,(sinit1(i),i=1,nser)
                write (lu,1036) hdr,(sstept1(i),i=1,nser)
                write (lu,1037) hdr,(sinit2(i),i=1,nser)
                write (lu,1038) hdr,(sstept2(i),i=1,nser)
                write (lu,1039) hdr,(stfix(i),i=1,nser)
                write (lu,1040) hdr,(sfac(i),i=1,nser)
                write (lu,1140) hdr,(sratio(i),i=1,nser)
                if (lu == luttyo) write (lu,1041) hdr,(idepnd(i),i=1,nser)
            end if
    10 END DO
    write(lulog,4009)
    write(lulog,4007)

    return

!######################################################################

    1000 format(a,'gxx,gyy,gzz = ',2(f9.6,','),f9.6)
    1001 format(a,'g1,g2,g3 = ',2(f9.6,','),f9.6)
    2001 format(a,'gprp,grhm,gpll = ',2(f9.6,','),f9.6)
    1002 format(a,'in2 =',i2)
    1003 format(a,'Axx,Ayy,Azz = ',2(f9.4,','),f9.4)
    1004 format(a,'A1,A2,A3 = ',2(f9.4,','),f9.4)
    2004 format(a,'Aprp,Arhm,Apll = ',2(f9.4,','),f9.4)
    1005 format(a,'log(Rx,Ry,Rz) = ',2(f9.4,','),f9.4)
    1006 format(a,'log(Rbar,N,Nxy) = ',2(f9.4,','),f9.4)
    2006 format(a,'log(Rprp,Rrhm,Rpll) = ',2(f9.4,','),f9.4)
    1007 format(a,'ipdf = ',i1)
    1008 format(a,'m',a2,', pm',a2,' = ',i2,',',g10.3)
    1009 format(a,'djf,djfprp = ',g10.3,',',g10.3)
    1010 format(a,'ist,djf =',i3,',',g10.3)
    1011 format(a,'nort,psi = ',i3,',',f8.3)
    1012 format(a,'c20,c22,c40,c42,c44 = ',4(f7.4,','),f7.4)
    1013 format(a,'oss = ',g11.4)
    1014 format(a,'ald,bed,gad =',f7.2,',',f7.2,',',f7.2)
    1015 format(a,'alm,bem,gam =',f7.2,',',f7.2,',',f7.2)
    1016 format(a,'t2edi,t2ndi,t2efi,t1edi,t1ndi = ',4(f6.3,','),f6.3)
    1116 format(a,'Siteweight = ',f10.6)
    1017 format(a,'mwid,gib,lib,hwid = ',3(f8.4,','),f8.4)
    1018 format(a,'lemx,lomx,kmx,mmx,ipnmx = ',4(i3,','),i3, &
    /a,'ndimoss,ndimdss = ',i4,i5)
    1019 format(a,'cgtol,nstep = ',g10.3,i5)
    1020 format(a,'B0 =',f10.3)
    1022 format('log    ',a)
    1024 format('basis  ',a)
    1030 format(a,'*** Parameters for ',i2,' spectra in series ***')
    1031 format(a,'iexp = ',9(i3,','),i3)
    1032 format(a,'icomb = ',9(i3,','),i3)
    1033 format(a,'npt1 = ',9(i4,','),i4)
    1034 format(a,'npt2 = ',9(i4,','),i4)
    1035 format(a,'init1  = ',9(f6.1,','),f6.1)
    1036 format(a,'stept1 = ',9(f6.1,','),f6.1)
    1037 format(a,'init2  = ',9(f6.1,','),f6.1)
    1038 format(a,'stept2 = ',9(f6.1,','),f6.1)
    1039 format(a,'tfix = ',9(f7.1,','),f7.1)
    1040 format(a,'spect weight = ',9(f10.6,','),f10.6)
    1140 format(a,'data scale factor = ',9(f10.6,','),f10.6)
    1041 format(a,'idepnd = ',9(i3,','),i3)
    3000 format('*** Unable to open file ',a,' for parameter output ***')
    4001 format(a,'Rprp,Rrhm,Rpll [1/sec]= ',2(e9.4,','),e9.4)
    4002 format(a,'Rx,Ry,Rz [1/sec] = ',2(e9.4,','),e9.4)
    4003 format(a,'Rbar[1/sec],N,Nxy = ',2(e9.4,','),e9.4)
    4004 format(a,'oss [1/sec]= ',e9.4)
    4005 format(a,'t2edi,t1edi [1/sec]= ',1(e9.4,','),e9.4)
    4006 format(a,'T2,T1 [nsec]= ',1(f9.2,','),f9.2)
    4007 format('=======================================================')
    4008 format('************** PARAMETER conversion *************')
    4009 format('************** PARAMETER conversion END *************')

    end subroutine parc


!----------------------------------------------------------------------
!                    =========================
!                       subroutine STATC
!                    =========================

!  Prints out the information on variables on the given logical unit
!  number.

!----------------------------------------------------------------------
    subroutine statc( line,lu )
    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'lmcomm.inc'
    include 'miscel.inc'
!      include 'iterat.inc'

    character line*(LINELG)
    integer :: lu,i,isite,ispec

!----------------------------------------------------------------------
!     parameters for fit command
!----------------------------------------------------------------------

! only output spec 1 site 1:
    ispec=1
    isite=1
    write(lu,1000)
    write(lu,1001) xtol,ftol,gtol
    write(lu,1002) maxitr
    write(lu,1003) maxev
    write(lu,1004) factor

!----------------------------------------------------------------------
!     parameters for search command
!----------------------------------------------------------------------
    write(lu,1006)
    write(lu,1007) stol,sstep,smn,smx

!----------------------------------------------------------------------
!     variable informations
!----------------------------------------------------------------------
    if (nprm >= 1) then
        write(lu,1010) nprm
        write(lu,1011)
        do 10 i=1,nprm
            write(lu,1012) tag(i),fparm(ixpr(i),ispec,isite), &
            prscl(i),xfdstp(i)
        10 END DO
    else
        write(lu,1020)
    end if

    write(lu,1030) lemx,lomx,kmx,mmx,ipnmx
    write(lu,1032) ndimoss(basinfo(1,ispec,isite)), &
    ndimdss(basinfo(1,ispec,isite))
    write(lu,1040) fnorm/dsqrt(dfloat(nptot))

    return

!######################################################################

    1000 format(/,2x,'FIT parameters : ')
    1001 format(5x,'xtol = ',1p,e8.1,3x,'ftol = ',1p,e8.1,3x,'gtol = ', &
    1p,e8.1)
    1002 format(5x,'maximum iterations = ',i2)
    1003 format(5x,'maximum function evaluations = ',i3)
    1004 format(5x,'initial step bound in parameter search = ',f8.1)

    1006 format(/,2x,'SEARCH parameters : ')
    1007 format(5x,'tol = ',1p,e8.1,3x,'step = ',1p,e8.1,3x,'min = ', &
    1p,e8.1,3x,'max = ',1p,e8.1)

    1010 format(/,2x,'VARY parameters : ',i1,' variables')
    1011 format(2x,56('-'),/,6x,'VARIABLE',8x,'VALUE',10x,'SCALE',5x, &
    'JC STEP',/,2x,56('-'))
    1012 format(7x,a6,5x,g14.7,3x,f6.1,5x,1p,e8.1)

    1020 format(/,2x,'No variables are being varied.')
    1030 format(/,2x,'MTS lemx,lomx,kmx,mmx,ipnmx = ',i3,4(',',i2))
    1032 format(2x,'Dimension of matrix: ndimoss = ',i4,' ndimdss = ',i4)
    1040 format(/,2x,'Rms Deviation = ',g14.7,/)

    end subroutine statc
