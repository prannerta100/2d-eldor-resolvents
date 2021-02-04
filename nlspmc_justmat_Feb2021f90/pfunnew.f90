
!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                         =======================
!                            subroutine PFUN
!                         =======================

! Subroutine for interfacing EPRCGF spectral calculations with the
! MINPACK version of the Levenberg-Marquardt nonlinear least squares
! algorithm for fitting experimental spectra. This routine loads the
! necessary parameters into common /eprprm/ and then calls the SIM2D
! subroutine to calculate a spectrum or evaluate a Jacobian as required
! by the calling program.

! *** NLSPMC only ***
! For evaluation of a series of spectra (which share the eigenvalues and
! eigenvectors), the routine re-evaluates the spectrum using the saved
! eigenvalues and vectors.  Therefore any attempt to fit a series of
! parameters that affects the eigenvalues is forbidden (e.g. series in
! tilt angle or temperature).

! For evaluation of the Jacobian, the routine uses the forward-differences
! approximation to calculate the partial derivative of the spectrum with
! respect to each of the parameters being varied. The routine assumes
! that the function values at the point where the Jacobian is being
! evaluated are contained in fvec and the scaling factor for each spectrum
! in series option is contained in sfac ( common /expdat/ )

! In the Jacobian evaluation, the routine makes a special check for
! variation of the parameters relevent in series option (now gib,lib,
! hwid,mwid).  If line-broadening parameter (gib) is varied, it simply
! re-evaluates the spectrum using the saved time domain data (xspec).

! NOTE that in order to follow this strategy, the Jacobian should
! be calculated with respect to any line-broadening parameters before
! any other parameter is varied (which would change the time domain
! signal in the absence of Gaussian inhomogeneous broadening).

! In the future, one might hope to recognize when only the real or
! imaginary part of the matrix changes, and ultimately to utilize a
! perturbation approach to calculating the Jacobian.

! Modified calculation of fvec and fjac to include weighting option for
! total data sets.

! iflag=0 just print variables
! iflag=1 function evaluation
! iflag=2 or 3 Jacobian
! iflag=3 output Jacobian matrix

! Includes:
!    nlsdim.inc
!    stdio.inc
!    eprprm.inc
!    expdat.inc
!    parcom.inc
!    tdspec.inc
!    wkspcm.inc
!    iterat.inc

! Uses
!    sim2d
!    xshft
!    sscale
!----------------------------------------------------------------------

    subroutine pfun( totpts,n,xxx,fvec,fjac,ldfjac,iflag )
    implicit none
    logical :: vchange
    double precision :: enorm
    external vchange,enorm
    integer :: totpts,n,iflag,ldfjac,ispec,isite
    double precision :: xxx(n),fvec(totpts),fjac(ldfjac,n),xxxo(n)

    include 'limits.inc'
    include 'stdio.inc'
    include 'lpnam.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'tdspec.inc'
    include 'wkspcm.inc'
    include 'miscel.inc'
!      include 'iterat.inc'
! for storing the individual simulations:
    character(20) :: cchar
    data cchar/'1234567890abcdefghij'/
    character*(WORDLG) fnmtst/'simss??.tst'/
    real*8 :: amin,amax
    logical :: sitewgc

    integer :: i,icalc,ierr,npt1,npt2,isp,ixi,ixs,ixscw,j,k,ld
    double precision :: rmsdv,shift,resol,ovlap,xtemp,delx,ccc
    character dashes*132,jacnam*10
    integer :: jspec,jsite,ixprmx(MXVAR),iparmptr,ivar

    double precision :: zero,totwgt,siteweight,totwgt2(MXSPEC)
    parameter (zero=0.0D0)

!######################################################################

!----------------------------------------------------------------------
!    Print parameter values for the current iteration
!----------------------------------------------------------------------

    if (iflag == 0) then
        ld=min(132,12*n+16)
        do 1 i=1,ld
            dashes(i:i)='-'
        1 END DO
        if (iter == 1) then
            if (luout /= luttyo) then
                write (luttyo,1006) dashes(:ld)
                write (luttyo,1007) (tag(i),i=1,n)
                write (luttyo,1008) dashes(:ld)
            end if
            write (luout,1006)  dashes(:ld)
            write (luout,1007)  (tag(i),i=1,n)
            write (luout,1008)  dashes(:ld)
        end if
        rmsdv=fnorm/dsqrt(dfloat(nptot))
        sitewgc=.false.
        do 2 i=1,n
            iparmptr=ixpr(i)
            if (iparmptr == ISWGT) then
                sitewgc=.true.
                totwgt=0.0d0
                do 3 j=1,ncomps
                    totwgt=totwgt+fparm(ISWGT,1,j)
                3 END DO
                xxxo(i)=xxx(i)/totwgt
            else
                xxxo(i)=xxx(i)
            end if
        2 END DO
    !c
        if (sitewgc) then
            if (luout /= luttyo) &
            write (luttyo,1009) iter,rmsdv,(xxxo(i),i=1,n)
            write (luout,1009) iter,rmsdv,(xxxo(i),i=1,n)
        else
            if (luout /= luttyo) &
            write (luttyo,1009) iter,rmsdv,(xxx(i),i=1,n)
            write (luout,1009) iter,rmsdv,(xxx(i),i=1,n)
        end if
    !c
    !c         if (luout.ne.luttyo)
    !c     #                write (luttyo,1009) iter,rmsdv,(xxx(i),i=1,n)
    !c         write (luout,1009) iter,rmsdv,(xxx(i),i=1,n)
        if (ihltcmd /= 0) return
    
    !**********************************************************************
    !**********************************************************************
    ! IFLAG=1: Function evaluation
    !**********************************************************************
    !**********************************************************************
    
    else if (iflag == 1) then
    
    ! Here, the variables contained in xxx determine the new spectra to
    ! be calculated.  If there is nothing changed, we don't repeat the
    ! spectral calculation, if something is changed, we calculate only
    ! the necessary parts.  We require the same parameters to apply
    ! to all spectra.
    
        do 41 j=1,n	! for eack xxx variable,
        !           ixprmx(j)=0	! keep track of variable severity
            iparmptr=ixpr(j)	! get variable id set in addprm, srchc
        ! for all sites,spectra affected
            do 42 isite=1,ncomps
                do 43 ispec=1,nspectra
                ! if this variable applies to this site/spectrum:
                    if (ixx(iparmptr,ispec,isite) == j) then
                    !                 if (ispec.eq.1)write(*,*)'affecting site ',isite
                    ! check if this parameter is changed
                        if (vchange(xxx(j),fparm(iparmptr,ispec,isite))) then
                        ! if so, update fparm and basinfo to define simulation
                            matrxok(isite)=min(matrxok(isite),specid(iparmptr))
                        !                   if (ispec.eq.1)write(*,*)'changed to ',xxx(j)
                            fparm(iparmptr,ispec,isite)=xxx(j)
                        end if
                    end if
                    if (ixx(iparmptr,ispec,isite) /= ixx(iparmptr,1,isite)) &
                    then
                        stop
                    end if
                43 END DO
            42 END DO
        41 END DO
    
    ! Now matrxok tell if need new matrix etc.  After first matrix
    ! calculation, update this information if that matrix applies to more
    ! than one site.  Must apply to all sites in this version
    !     ...Loop over all spectra in a series and all sites....
    
    !       call ftest2(fparm(1,1,1),NFPRM,1,'pfun,fparm 1c   ',14)
    !       call ftest2(fparm(1,1,2),NFPRM,1,'pfun,fparm 2c   ',14)

        do 5 jsite=1,ncomps
            cursite=jsite	   ! used to keep track of basis, matrix etc.
            do 5 jspec=1,nspectra
            
            ! this is checked in datac:
            !          if(ixsp(nser).gt.mxpt)then
            !            write(*,*)'error in pfun, inc. mxpt',ixsp(nser),mxpt
            !            stop
            !          end if
                ixs=ixsp(jspec)	! pointer to this data set
                ixscw=ixspcw(jspec)	! this cw data set
            
            ! replace by new arrays recording what to do:
                icalc=2
                if(matrxok(jsite) >= 0 .AND. matrxok(jsite) <= 2) then
                    icalc = 2-matrxok(jsite)  ! 2=need matrix,1=use eval
                elseif (matrxok(jsite) == 4) then
                ! here if matrxok = 4
                    icalc=0	! here if vary site weighting
                end if
            
            !----------------------------------------------------------------------
            ! If desired, find the shift that will give optimal overlap between
            ! a spectrum and the data.  This is achieved by comparing cw magnitude
            ! spectrum with the cw-equivalent spectrum extracted from 2D data,
            ! which was obtain in data command and saved in cwdata array in common
            ! block /expdat/.
            ! if got dat and want shift, must do cw calculation.  Here for multiple
            ! components, do for each time and sum to give composite cw spec.
            ! In the case of varying site weight, redo shift, but don't do it for
            ! case of vary gib or lib.
            !----------------------------------------------------------------------
            !  Don't redo shift if varying gib.
            !  For site weighting, don't call sim2d, to preserve xspec, but check
            !  shift.
            
                if (dataOK .AND. (icalc /= 0) .AND. (matrxok(jsite) /= 4) .AND. &
                (sishft(jspec) /= 0)) then
                ! set this spectrum parameters                  ** cw spectrum
                    call setspc( jspec,jsite,npt1,npt2 )
                ! specify cw calculation
                    fparm(ISHFT,jspec,jsite)=0.0d0
                    iparm(IIEXP,jspec,jsite)=0
                    iparm(IICOMB,jspec,jsite)=1
                    resol=1.0d3/(sstept2(jspec)*float(npt2-1))
                    npt1=npt2
                    npt2=1	! specify CW calculation, result to
                
                    if(matrxok(jsite) /= 4) then
                    ! need to mod sim2d to use site basis
                        call sim2d( icalc,xspectr(ixs,jsite),cspectr(ixs,jsite), &
                        npt1,npt2,n,iflag,ixi,delx,jspec,jsite,ierr )
                    
                        if(icalc < 1)stop 'icalc=0 error in pfunnew'
                    
                    ! return here means have eval, evec, and matrix for this spec,
                    ! reset matrxok to indicate  this.
                    ! note 2D spec requires calculation of xspec if not present.
                    
                        matrxok(jsite)=max(1,matrxok(jsite))  ! have matrix for site.
                        if ((ihltcmd == 0) .AND. .FALSE. ) then
                            if (luout /= luttyo) write(luttyo,1004)
                            write(luout,1004)
                            iflag=-1
                            return
                        else if (ierr /= 0) then
                            if (luout /= luttyo) write(luttyo,1003)
                            write(luout,1003)
                            iflag=-1
                        end if
                    end if
                
                ! save the cw spectra for all sites,spectra. Later, sum with weights
                !  do this so that if changing only weights, we don't recalculate cwsim
                
                    do 6 j=1,npt1
                        ccwsim(ixscw-1+j,jsite)=cspectr(ixs-1+j,jsite)
                    6 END DO
                end if	! end CW calculation.
        5 END DO
    
    !----------------------------------------------------------------------
    !        Calculate shift term (shft).
    !        First, sum over sites with site weighting
    !----------------------------------------------------------------------
    
    ! we use absolute values of the weights normalized:
    
        do 71 jspec=1,nspectra
            if (dataOK .AND. (icalc /= 0) .AND. (sishft(jspec) /= 0)) then
                ixscw=ixspcw(jspec)	! this cw data set
                totwgt=0.0D0
                do 72 jsite=1,ncomps
                    siteweight=abs(fparm(ISWGT,jspec,jsite))
                    totwgt=totwgt+siteweight
                    do 8 j=1,npt1
                        if(jsite == 1) ccwsimtot(ixscw-1+j)=(0.0D0,0.0D0)
                        ccwsimtot(ixscw-1+j)=ccwsimtot(ixscw-1+j)+ &
                        ccwsim(ixscw-1+j,jsite)*siteweight
                        if(jsite == ncomps) cwsimtot(ixscw-1+j)= &
                        cdabs(ccwsimtot(ixscw-1+j))/totwgt
                    8 END DO
                
                !  Then compare this spectrum with CW data to determine best shift
                
                72 END DO
                call xshft(cwdata(ixscw),cwsimtot(ixscw),npt1, &
                cwsp1,cwsp2,shift,ovlap)
                sshft(jspec)=resol*shift	! save shift for setspc
            end if
        71 END DO
    
    !----------------------------------------------------------------------
    !        Calculate 2D spectrum with desired shift terms (sshft).
    !----------------------------------------------------------------------
    
    !       call ftest2(fparm(1,1,1),NFPRM,1,'pfun,fparm 1b   ',14)
    !       call ftest2(fparm(1,1,2),NFPRM,1,'pfun,fparm 2b   ',14)
    
        do 11 jsite=1,ncomps
        !          cursite=jsite	; used to keep track of basis, matrix etc.
            do 10 jspec=1,nspectra
                ixs=ixsp(jspec)	! pointer to this data set
            !          ixscw=ixspcw(jspec)	! this cw data set
                icalc=2
                if(matrxok(jsite) >= 0 .AND. matrxok(jsite) <= 2) then
                    icalc = 2-matrxok(jsite)  ! 2=need matrix,1=use eval, 0 for gib
                elseif (matrxok(jsite) == 4) then
                ! here if matrxok = 4
                    icalc=0	! here if vary site weighting
                end if
                call setspc( jspec,jsite,npt1,npt2 )
                if(matrxok(jsite) /= 4) then
                
                !            time1=mclock()
                    call sim2d( icalc,xspectr(ixs,jsite),cspectr(ixs,jsite), &
                    npt1,npt2,n,iflag,ixi,delx,jspec,jsite,ierr )
                !       call ftest3(xspectr(1,1),16384,'pfun,xsptr 1b   ',14)
                !       call ftest3(xspectr(1,2),16384,'pfun,xsptr 2b   ',14)
                    matrxok(jsite)=max(1,matrxok(jsite))  ! have matrix for site.
                
                ! note don't set matrxok = 0 since it won't be correct for next spectrum
                
                    if ((ihltcmd == 0 .AND. .FALSE. )) then
                        if (luout /= luttyo) write(luttyo,1004)
                        write(luout,1004)
                        iflag=-1
                        return
                    else if (ierr /= 0) then
                        if (luout /= luttyo) write(luttyo,1003)
                        write(luout,1003)
                        iflag=-1
                        return
                    end if
                end if
            
            !               ... End of loop over spectra in a series and sites
            
            10 END DO
            matrxok(jsite)=4	! got everything for this site.  If
        !                            nothing is changed, next calc will be fast
        11 END DO
    
    ! Accumulate 2D spectra over sites, weighting appropriately
    ! We assume the site weights are identical for all spectra.
    ! This may not be true for temperature variation in the future.
    
        do 115 jspec=1,nspectra
            totwgt2(jspec)=0.0D0
            do 115 jsite=1,ncomps
                totwgt2(jspec)=totwgt2(jspec)+ &
                abs(fparm(ISWGT,jspec,jsite))
        115 END DO
        do 12 jsite=1,ncomps
            i=0
            do 12 jspec=1,nspectra
            !  for testing, write the site,spectrum simulations to files
                fnmtst(6:6)=cchar(jspec:jspec)
                fnmtst(7:7)=cchar(jsite:jsite)
                amin=0.0d0
                amax=0.0d0
                do 13 j=1,ndata(jspec)
                    spectr(j,jsite)=cdabs(cspectr(j,jsite))
                13 END DO
                call wrfit(spectr(1+(jspec-1)*ndata(jspec),jsite), &
                jspec,amin,amax,fnmtst)
            
                siteweight=abs(fparm(ISWGT,jspec,jsite))
                do 12 j=1,ndata(jspec)
                    i=i+1
                    if(jsite == 1)ctotspec(i)=(0.0D0,0.0D0)
                    ctotspec(i)=ctotspec(i)+cspectr(i,jsite)*siteweight
                    if(jsite == ncomps)totspec(i)=cdabs(ctotspec(i))/totwgt2( &
                    jspec)
                ! be careful, from this point on we work with totspec instead of spectr
        12 END DO
    
    !----------------------------------------------------------------------
    ! Calculate the scale factor(s) for each spectrum
    ! according to the dependencies (idepnd) that user specified in data
    ! command.  The subroutine returns the scaled spectra and the scale
    ! factors (sfac in common block /expdat/).
    !----------------------------------------------------------------------
    
        if (dataOK) call sscale(totspec)
    
    !----------------------------------------------------------------------
    ! Calculate the difference between the experimental and scaled
    ! calculated spectra and store the result in fvec.
    !----------------------------------------------------------------------
    
        i=0
        do 20 jspec=1,nspectra
            do 20 j=1,ndata(jspec)
                i=i+1
                if(dataOK) then
                    fvec(i)=(totspec(i)-data(i))*datwgt(jspec)
                !               else
                !                 fvec(i)=totspec(i)	! not used
                end if
        20 END DO
    
        if (idebug /= 0) then
            write (ludeb,2010)
            do 21 j=1,nser
                write (ludeb,2020) j,sfac(j),sshft(j)
            21 END DO
        end if
    
    !        ** Update X vector for the case where the real part of
    !           the eigenvalues are negative.
    ! this is no longer done.  RC 11/11/98.
    !         do 6 j=1,n
    !          if(xxx(j)-fparm(ixpr(j)) .gt. 1.0D-10) then
    !            if (luout.ne.luttyo) then
    !              write(luttyo,*)'remapping x: ',j,xxx(j),fparm(ixpr(j))
    !              write(luttyo,*)'This should not happen!'
    !            end if
    !              write(luout,*)'remapping x: ',j,xxx(j),fparm(ixpr(j))
    !              write(luout,*)'This should not happen!'
    !          end if
    !            xxx(j)=fparm(ixpr(j))
    ! 6       continue
    
    !**********************************************************************
    ! IFLAG >= 2: Jacobian evaluation by forward-difference approximation:
    !             fjac = (f(x0+dx)-f(x0))/dx
    !          where f(x0) = fvec+data, and dx is assumed to be positive
    ! We assume we have a previously calculated spectrum, for the current
    ! set of parameters.  Then each is varied by dx and a new spectrum
    ! calculated.
    
    !**********************************************************************
    
    else if (iflag >= 2) then
    
    !----------------------------------------------------------------------
    ! Loop over all parameters, introducing the forward-difference step
    ! into each parameter.  On entry, have spectra calculated
    ! from the current values of the parms and matrxok=4, save this info.
    ! As vary each parm, do only calculation necessary for that parameter.
    !----------------------------------------------------------------------
    
        do 30 ivar=1,n	        ! for eack xxx variable, increment it
            xtemp=xxx(ivar)		! save current value of xxx(i)
            xxx(ivar)=xxx(ivar)+xfdstp(ivar)	! introduce step into this parm
            ixi=ixpr(ivar) 	! get variable id set in addprm, srchc
            delx=xfdstp(ivar)	! save step size
        
        ! for all sites,spectra affected, map only this xxx into fparm.
        !  The others should not be changed.
        
            do 82 isite=1,ncomps
                do 83 ispec=1,nspectra
                ! if this variable applies to this site/spectrum:
                    if (ixx(ixi,ispec,isite) == ivar) then
                    ! update fparm to define simulation.  We know variable is changed.
                        matrxok(isite)=min(matrxok(isite),specid(ixi))
                        fparm(ixi,ispec,isite)=xxx(ivar)
                    end if
                    if (ixx(ixi,ispec,isite) /= ixx(ixi,1,isite)) &
                    then
                        write(*,*)'different parameters for different '
                        write(*,*)'spectra not allowed yet (in pfun).'
                        stop
                    end if
                83 END DO
            82 END DO
        
        !  for multiple variables, restore the previous one since we only
        !  vary one parameter at a time in computing the Jacobian
        
            if (ivar > 1) then
                do 821 isite=1,ncomps
                    do 831 ispec=1,nspectra
                    ! if this variable applies to this site/spectrum:
                        if (ixx(ixi,ispec,isite) == ivar-1) then
                        ! update fparm to define simulation.  We know variable is changed.
                            matrxok(isite)=min(matrxok(isite),specid(ixi))
                            fparm(ixi,ispec,isite)=xxx(ivar-1)
                        end if
                    831 END DO
                821 END DO
            end if
        
        !       call ftest2(fparm(1,1,1),NFPRM,1,'pfun,fparm 1a   ',14)
        !       call ftest2(fparm(1,1,2),NFPRM,1,'pfun,fparm 2a   ',14)
        
        !        ...Loop over all spectra in a series, same as do 10 loop
        
            do 28 jsite=1,ncomps
            !            cursite=jsite        ; used to keep track of basis, matrix etc.
                do 27 jspec=1,nspectra
                    ixs=ixsp(jspec)       ! pointer to this data set
                !             ixscw=ixspcw(jspec)   ! this cw data set
                    icalc=2
                    if(matrxok(jsite) >= 0 .AND. matrxok(jsite) <= 2) then
                        icalc = 2-matrxok(jsite)  ! 2=need matrix,1=use eval
                    elseif (matrxok(jsite) == 4) then
                    ! here if matrxok = 4
                        icalc=0	! here if vary site weighting
                    end if
                    call setspc( jspec,jsite,npt1,npt2 )
                    if(matrxok(jsite) /= 4) then
                    !               time1=mclock()
                        call sim2d( icalc,xspectr(ixs,jsite),cspectr(ixs,jsite), &
                        npt1,npt2,n,iflag,ixi,delx,jspec,jsite,ierr )
                    !       call ftest3(xspectr(1,1),16384,'pfun,xsptr 1a   ',14)
                    !       call ftest3(xspectr(1,2),16384,'pfun,xsptr 2a   ',14)
                        matrxok(jsite)=max(1,matrxok(jsite))  ! have matrix for site.
                    
                        if ((ihltcmd /= 0 .AND. .FALSE. )) then
                            if (luout /= luttyo) write(luttyo,1004)
                            write(luout,1004)
                            iflag=-1
                            return
                        else if (ierr /= 0) then
                            if (luout /= luttyo) write(luttyo,1003)
                            write(luout,1003)
                            iflag=-1
                            return
                        end if
                    end if
                
                !               ... End of loop over spectra in a series and sites
                
                27 END DO
                matrxok(jsite)=4      ! got everything for this site.  If
            !                            nothing is changed, next calc will be fast
            28 END DO
        
        ! sum over the sites with weighting.  Shift already in spectra.  Use
        ! previous scale factor.
        
            do 215 jspec=1,nspectra
                totwgt2(jspec)=0.0D0
                do 215 jsite=1,ncomps
                    totwgt2(jspec)=totwgt2(jspec)+ &
                    abs(fparm(ISWGT,jspec,jsite))
            215 END DO

            do 92 jsite=1,ncomps
                i=0
                do 92 jspec=1,nspectra
                    do 92 j=1,ndata(jspec)
                        i=i+1
                        if(jsite == 1)ctotspec(i)=(0.0D0,0.0D0)
                        ctotspec(i)=ctotspec(i)+cspectr(i,jsite)* &
                        abs(fparm(ISWGT,jspec,jsite))
                        if(jsite == ncomps)totspec(i)= &
                        cdabs(ctotspec(i))/totwgt2(jspec)
                    ! be careful, from this point on we work with totspec instead of spectr
            92 END DO
        !           call ftest1(totspec,i,'92 pfunnew   ',11 )
        
        !          ... Store function value before forward-differences step
        
            do 39 jspec=1,nspectra
                ixs=ixsp(jspec)
                do 25 j=1,ndata(jspec)
                    k=ixs+j-1
                ! since data and fvec are not changed here, use to calc. prev simulation
                ! at the original values of xxx.
                    fjac(k,ivar)=-(fvec(k)+data(k)*datwgt(jspec))
                    fjac(k,ivar)=(fjac(k,ivar)+sfac(jspec)*totspec(k)* &
                    datwgt(jspec))/delx
                25 END DO
            39 END DO
        !           call ftest2(fjac(1,ivar),mxpt,1,'39 pfunnew fjac  ',16 )
        
            xxx(ivar)=xtemp		! restore x
        
        ! Restore fparms array with this value, readying for next variable:
        
            do 820 isite=1,ncomps
                do 830 ispec=1,nspectra
                ! if this variable applies to this site/spectrum:
                    if (ixx(ixi,ispec,isite) == ivar) then
                    ! update fparm to define simulation.  We know variable is changed.
                        matrxok(isite)=min(matrxok(isite),specid(ixi))
                        fparm(ixi,ispec,isite)=xxx(ivar)
                    end if
                    if (ixx(ixi,ispec,isite) /= ixx(ixi,1,isite)) &
                    then
                        write(*,*)'different parameters for different '
                        write(*,*)'spectra not allowed yet (in pfun).'
                        stop
                    end if
                830 END DO
            820 END DO
        !                     *** end of loop over parameters
        30 END DO
    
    ! note here, the x values are as they were on entry to pfun, while
    !  the spectra are calculated from the last parameter incremented.
    !  This should not be a problem, since the fparm values do
    !  correspond to the matrix.
    
    !----------------------------------------------------------------------
    !     Optional output of Jacobian matrix for this iteration
    !----------------------------------------------------------------------
    
        if (iflag == 3) then
            write(jacnam,2000) iter
            open(unit=ludisk,file=jacnam,status='unknown', &
            access='sequential',form='formatted')
            do 40 i=1,totpts
                write (ludisk,2003) data(i),fvec(i),(fjac(i,j),j=1,n)
            40 END DO
            close(ludisk)
        end if
    
    end if

    return

! ### format statements ##############################################

    1003 format(/20x,'*** Error: fit or search procedure terminating', &
    ' ***'/)
    1004 format(/20x,'*** Fit halted by user ***')
    1006 format(/a)
    1007 format('Iter',5x,'RmsDv',10(6x,a))
    1008 format(a)
    1009 format(i4,4x,g10.4,2x,10(g12.6))
    2000 format('nlsjac.',i3.3)
    2003 format(10(g12.5,1x))
    2010 format(/' ** ',6x,'Spectrum',5x,'Scale Factor',5x, &
    'Shift Factor **')
    2020 format(13x,i1,7x,g14.7,3x,g14.7)

    end subroutine pfun

!----------------------------------------------------------------------
!                    =========================
!                       subroutine SETSPC
!                    =========================

!  Sets the series-dependent parameters in the fparm and iparm arrays.
!  These parameters are not allowed to change the eigenvalues.

!  Specifically, the parameters allowed in series option are

!    Integer :  iexp, icomb, npt1, npt2
!    Floating:  init1, stept1, init2, stept2, tfix, shft

!  on entry, ispc,isit specify the current spectrum,site.

!----------------------------------------------------------------------

    subroutine setspc( ispc,isit,npt1,npt2 )
    implicit none
    integer :: ispc,isit,npt1,npt2
    logical :: vchange
    external vchange
    logical :: chnge

    include 'limits.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'

    if (ispc <= 0 .OR. ispc > nser) return

    chnge=.false.	! are we changing anything?
    chnge=chnge.or.vchange(fparm(IINIT1,ispc,isit),sinit1(ispc))
    fparm(IINIT1,ispc,isit) = sinit1(ispc)
    chnge=chnge.or.vchange(fparm(ISTEPT1,ispc,isit), sstept1(ispc))
    fparm(ISTEPT1,ispc,isit) = sstept1(ispc)
    chnge=chnge.or.vchange(fparm(IINIT2,ispc,isit), sinit2(ispc))
    fparm(IINIT2,ispc,isit) = sinit2(ispc)
    chnge=chnge.or.vchange(fparm(ISTEPT2,ispc,isit),sstept2(ispc))
    fparm(ISTEPT2,ispc,isit) = sstept2(ispc)
    chnge=chnge.or.vchange(fparm(ITFIX,ispc,isit),stfix(ispc))
    fparm(ITFIX,ispc,isit) = stfix(ispc)
    chnge=chnge.or.vchange(fparm(ISHFT,ispc,isit),sshft(ispc))
    fparm(ISHFT,ispc,isit) = sshft(ispc)


    chnge=chnge.or.(iparm(IIEXP,ispc,isit).ne.siexp(ispc))
    iparm(IIEXP,ispc,isit) = siexp(ispc)
    chnge=chnge.or.(iparm(IICOMB,ispc,isit).ne.sicomb(ispc))
    iparm(IICOMB,ispc,isit) = sicomb(ispc)
! any changes here can use Eval,Evec already calculated:
    if(chnge) matrxok(isit)=min(matrxok(isit),1)

    npt1 = snpt1(ispc)
    npt2 = snpt2(ispc)

    return
    end subroutine setspc
