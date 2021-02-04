!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                       subroutine DATAC
!                    =========================

! Interprets a "data" command from the given line. The format of
! the command is as follows:

! data <dataid> { ascii|binary scale|noscale shift|noshift store}

! ** Special Caution necessary for NLSPMC **
!   The order of <dataid> in the fit of series spectra SHOULD
!   match the order of parameters specified in series command.

!   dataid : base name for datafile and associated output files:
!            <dataid>.SPC -- datafile
!            <dataid>.FIT -- fit spectrum
!   ascii|binary : specifies the format of the datafile
!            default (ascii)
!   scale|noscale : specifies that NLSPMC program should scale the
!            corresponding spectrum to fit.  The flag idepnd are
!            set according to the choice.  Refer comments in sscale.f.
!            default (scale)
!   shift|noshift : specifies that NLSPMC program should shift the
!            corresponding spectrum to fit.  The flag sishft is
!            set according to the choice.  Refer comments in pfunnew.f
!            before calling xshft.  default (shift)
!   store|nostore :  specifies that the actual data (may have been
!            spline-smoothened) be stored in the file <datid>.nsp
!            in 2D plotting format.

!----------------------------------------------------------------------
    subroutine datac( line )
    implicit none
    character line*80,token*30

    include 'limits.inc'
    include 'datas.inc'
    include 'lmcomm.inc'
    include 'parms.inc'
    include 'names.inc'
    include 'wkspcm.inc'
    include 'stdio.inc'

    integer :: nkeywd
    parameter (nkeywd=8)

    integer :: i,iret,ispc,ival,ix,ixcw,j,lth,iform,istore,lu, &
    isp,ispi,mpts,mptscw
    double precision ::  inif2,res2,inif1,res1,amin,amax,scfac,tog

    character(8) :: keywrd(nkeywd)
    data keywrd /'ASCII','BINARY','SCALE','NOSCALE', &
    'SHIFT','NOSHIFT','STORE','NOSTORE'/

    character(20) :: wndoid(MXSPEC)
    logical :: itoken
    integer :: getdat
    external itoken,getdat

!----------------------------------------------------------------------
!  Check if the series parameters are properly set
!----------------------------------------------------------------------
          
    if (nsparm < 9) then
        write(luttyo,1100)nsparm
        return
    end if

!----------------------------------------------------------------------
!  Get the name of the datafile
!----------------------------------------------------------------------
    call gettkn(line,token,lth)

!----------------------------------------------------------------------
!  If (1) no data name was specified, or (2) the number of datafiles
!  already read in equals the number of spectra to be calculated (nser)
!  then reset the data buffer
!----------------------------------------------------------------------
    if (lth == 0 .OR. nspc >= nser) then
        nspc=0
        nptot=0
        write(luttyo,1020)
        if (lth == 0) return
    end if

    nspc=nspc+1
    dataid(nspc)=token

! ....initialization of default values

    iform=0
    ndata(nspc)=snpt1(nspc)*snpt2(nspc)
    idepnd(nspc)=0
    idepnd(nspc+1)=0
    sishft(nspc)=1
    istore=0

!----------------------------------------------------------------------
!     Look for a keyword
!----------------------------------------------------------------------
    5 call gettkn(line,token,lth)

    if (lth /= 0) then
        lth=min(lth,8)
        call touppr(token,lth)
        do 6 i=1,nkeywd
            if (token(:lth) == keywrd(i)(:lth)) go to 7
        6 END DO
    
        write (luttyo,1000) token(:lth)
        go to 5
    
    !----------------------------------------------------------------------
    !     Keyword found: assign appropriate value using next token
    !----------------------------------------------------------------------
    
    ! --- Non-argument keywords:
    !     ASCII, BINARY, SCALE, NOSCALE, SHIFT, NOSHIFT
    
    !                                         *** ASCII keyword
        7 if (i == 1) then
            iform=0
        !                                         *** BINARY keyword
        else if (i == 2) then
            iform=1
            write (luttyo,1030)
            iform=0
        !                                         *** SCALE keyword
        else if (i == 3) then
            idepnd(nspc)=0
        !                                         *** NOSCALE keyword
        else if (i == 4) then
            idepnd(nspc)=1
        !                                         *** SHIFT keyword
        else if (i == 5) then
            sishft(nspc)=1
        !                                         *** NOSHIFT keyword
        else if (i == 6) then
            sishft(nspc)=0
        !                                         *** STORE keyword
        else if (i == 7) then
            istore=1
        !                                         *** NOSTORE keyword
        else if (i == 8) then
            istore=0
        end if
    
        go to 5
    
    !----------------------------------------------------------------------
    !     Read in datafile (end of data command)
    !----------------------------------------------------------------------
    else
    
        call setdat( dataid(nspc) )
    
    !    --- Check whether there is enough storage for the new data
    
        ix=nptot+1
        ixcw=nptotcw+1
        if ( (nptot+ndata(nspc)) > MXPT ) then
            write(luttyo,1050) MXPT
            nspc=nspc-1
            return
        end if
    
        lu=0
    
        iret=getdat (dtname,iform,snpt1(nspc),snpt2(nspc),lu, &
        siexp(nspc),sicomb(nspc),sstept1(nspc),sstept2(nspc), &
        sratio(nspc),data(ix),cwdata(ixcw),rwsp,cwsp1,cwsp2, &
        cwsp3,c2wsp1,c2wsp2)
    
        if (iret /= 0) then
        !                               *** Error opening/reading datafile
            if (iret == -1) write(luttyo,1060) dtname(:lthdnm)
            if (iret == -2) write(luttyo,1070) dtname(:lthdnm)
            if (iret == -3) write(luttyo,1080) dtname(:lthdnm)
            if (iret == -4) write(luttyo,1090) dtname(:lthdnm)
            nspc=nspc-1
            return
        end if
    
        if (istore == 1) then
        !                               *** store the splined datafile
        !	    tog=dble(snpt1(nspc)-1)/dble(snpt1(nspc))
        !            inif1=-5.0d2*tog/sstept1(nspc)
        !            res1=-2.0d0*inif1/dble(snpt1(nspc)-1)
        !	    inif1=-5.0d2/sstept1(nspc)
        !	    tog=dble(snpt2(nspc)-1)/dble(snpt2(nspc))
        !            inif2=-5.0d2*tog/sstept2(nspc)
        !            res2=-2.0d0*inif2/dble(snpt2(nspc)-1)
        !	    inif2=-5.0d2/sstept2(nspc)
            amin=0.0d0
            amax=0.0d0
            do 45 j=ix,ix+ndata(nspc)-1
                if ( data(j) > amax ) amax=data(j)
            !               if ( data(j).lt.amin ) amin=data(j)
            45 END DO
        !            call wrfit( data(ix),snpt1(nspc),snpt2(nspc),
        !     #            inif2,res2,inif1,res1,amin,amax,nsname,lthdnm )
            call wrfit( data(ix),nspc,amin,amax,nsname )
        
        end if
    
        ixsp(nspc)=ix
        ixspcw(nspc)=ixcw
        sshft(nspc)=0.0d0
        wndoid(nspc)=dataid(nspc)
        wndoid(nspc)(10:10)=char(0)
        nptot=nptot+ndata(nspc)
        nptotcw=nptotcw+snpt2(nspc)
    
    !----------------------------------------------------------------------
    !        Scale the series spectra within reasonable range (smax=1000)
    !----------------------------------------------------------------------
    ! if we have all the data sets:
        if (nspc == nser) then
        
            isp=1
            100 if (uniflg) then

                ispi=isp
                ix=ixsp(isp)
                ixcw=ixspcw(isp)
                mpts=ndata(isp)
                mptscw=snpt1(isp)
            !                              * search next independent spectra
                50 isp=isp+1
                if (idepnd(isp) == 1) then
                    mpts=mpts+ndata(isp)
                    mptscw=mptscw+snpt1(isp)
                    go to 50
                end if
            else
                ispi=1
                isp=nspc+1
                ix=1
                mpts=nptot
                ixcw=1
                mptscw=nptotcw
            end if
        !                              * obtain scale factor
            amax=0.0d0
            do 52 j=ispi,isp-1
                if (sratio(j) > amax) amax=sratio(j)
            52 END DO
            scfac=1.0d3/amax
            do 54 j=ispi,isp-1
                sratio(j)=scfac
            54 END DO
        !                              * scale the spectra
            do 56 j=ix,ix+mpts-1
                data(j)=data(j)*scfac
            56 END DO
            do 58 j=ixcw,ixcw+mptscw-1
                cwdata(j)=cwdata(j)*scfac
            58 END DO
        
            if (isp > nspc) then
                dataOK=.true.
                return
            end if
        
            go to 100
        
        end if
    
    !######################################################################
    
    end if

    return

! #### format statements ########################################

    1000 format('*** Unrecognized DATA keyword: ''',a,''' ***')
    1020 format('*** Data buffer has been reset *** ')
    1030 format('*** Binary format for input data is not supported ***')
    1050 format('*** Maximum number of data points (',i4,') exceeded ***')
    1060 format(13x,'*** Error opening or reading datafile ''',a, &
    ''' ***')
    1070 format(13x,'*** Too many data points in datafile ''',a,''' ***')
    1080 format(13x,'*** Data pts outside input range in ''',a,''' ***')
    1090 format(13x,'*** Can not extract cw-equivalent spectrum ***')
    1100 format('*** Series parameters not properly set ***',i3)
    end subroutine datac
