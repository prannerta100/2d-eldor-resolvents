!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                         SUBROUTINE FITP
!                    =========================
!     This routine carries out nonlinear least-squares analysis of
!  slow-motional 2D spectra with EPRESA calculation.   It uses a
!  modification of the Levenberg-Marquardt trust-region algorithm
!  available in MINPACK.

!  Includes:
!     nlsdim.inc
!     stdio.inc
!     nlsnam.inc
!     eprprm.inc
!     expdat.inc
!     parcom.inc
!     tdspec.inc
!     lmcom.inc
!     lmtxt.inc
!     iterat.inc

!  Uses
!     setdat
!----------------------------------------------------------------------

    subroutine fitp
    implicit none

    include 'limits.inc'
    include 'stdio.inc'
    include 'names.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'tdspec.inc'
    include 'lmcomm.inc'
    include 'miscel.inc'
    include 'lmtxt.inc'
!      include 'iterat.inc'

    logical :: fexist
    integer :: i,ifile,iflag,iret,ix,j,k,l,mrgn,nfev,njev, &
    npt1,npt2,nstot,iparmptr,isite,ispec
    double precision :: field,rmsdv,sigma,denom,inif1,res1, &
    inif2,res2,amax,amin,scfac
    character line*132,chr*2

    double precision :: enorm,ten
    external pfun,enorm
    parameter(ten=10.0D0)

!######################################################################

! set setbas from basinfo

    setbas=.true.
    do 1 isite=1,ncomps
    !        seteval(isite)=.false.	! make no assumption about matrix
        do 1 ispec=1,nspectra
            setbas=setbas .and. (basinfo(1,ispec,isite) .ne. 0)
    1 END DO
    if ( .NOT. setbas) then
        write (luout,1001)
        if (luout /= luttyo) write (luttyo,1001)
        1001 format(/,'** Basis Index set is not properly set **')
        return
    end if

!----------------------------------------------------------------------
!     Open trace file (if desired)
!----------------------------------------------------------------------
    if (itrace /= 0) then
        inquire(file=trname(:lthfnm),exist=fexist)
        open(lutrc,file=trname(:lthfnm),status='unknown', &
        access='sequential',form='formatted')
        if (fexist) then
        !                              * append the trace file
            5 read (lutrc,'(a2)',end=7) chr
            go to 5
        end if
        7 itrace=lutrc
    end if

    if (nspc /= nser) write(luttyo,1040)

!  Put correct values for the unknowns into the xxx vector.

    do 10 j=1,nprm	! for each parameter to be varied
        iparmptr=ixpr(j)	! identify variable
    ! for all sites, spectra, find one which is being varied:
        do 20 isite=1,ncomps
            do 20 ispec=1,nspectra
            ! check matching sites with requested site
                if (ixx(iparmptr,ispec,isite) == j) then
                    xxx(j)=fparm(iparmptr,ispec,isite)
                    cycle
                end if
        20 END DO
        write(*,*)'error in fitp, never should get here',j
        stop

    10 END DO

!----------------------------------------------------------------------
!     Call pfun once if only a single calculation is specified
!     or if there is no data file, implying a single calculation
!----------------------------------------------------------------------
    if (maxitr <= 1 .OR. maxev <= 1 .OR. nprm <= 0 .OR. &
     .NOT. dataOK) then
    !         seteval=.false.	! make no assumption about matrix
        if ( .NOT. dataOK) then
            nstot=nptot
            nptot=0
            do 11 j=1,nser
                nptot=nptot+ndata(j)
            11 END DO
        end if
        iflag=1
        call pfun(nptot,nprm,xxx,fvec,fjac,mxpt,iflag)
        if (iflag < 0 .OR. ihltcmd /= 0) return
    
        if (dataOK) then
            rmsdv=enorm(nptot,fvec)/dsqrt(dfloat(nptot))
            write(luout,1046) rmsdv
            if (luout /= luttyo) write(luttyo,1046) rmsdv
        else
            nptot=nstot
        end if
    
    else
    
    !======================================================================
    !     Call Marquardt-Levenberg nonlinear least squares routine
    !======================================================================
    !----------------------------------------------------------------------
    !     Search parameter list and set flags for different types of parameters
    !     : g-tensor, a-tensor, motional parameter, orienting potential
    !----------------------------------------------------------------------
        gflag=.false.
        aflag=.false.
        rotflg=.false.
        potflg=.false.
    !         seteval=.false.	! make no assumption about matrix
        do 12 i=1,nprm
            ix=ixpr(i)
            gflag=gflag .or. (ix.ge.IGXX .and. ix.le.IGZZ)
            aflag=aflag .or. (ix.ge.IAXX .and. ix.le.IAZZ)
            rotflg=rotflg .or. (ix.ge.IDX .and. ix.le.IPMZZ) &
            .or. (ix.eq.IDJF .or. ix.eq.IDJFPRP) &
            .or. (ix.eq.IOSS)
            potflg=potflg .or. (ix.ge.IC20 .and. ix.le.IC44)
        12 END DO
    
        nprint=1
    
        call lmnls( pfun,nptot,nprm,xxx,fvec,fjac,mxpt,ftol, &
        xtol,gtol,maxev,maxitr,diag,prscl,factor,nprint, &
        itrace,jacobi,info,nfev,njev,ipvt,qtf,gnvec,gradf, &
        work1,work2,work3,work4 )
    
    !----------------------------------------------------------------------
    !     calculate covariance matrix and estimate errors assuming
    !     normal distribution
    !----------------------------------------------------------------------
    
        if (ihltcmd /= 0) return
    !                                         *** MINPACK error return
        if (info == 0 .OR. info > 3) then
            write (luout,2034) minerr(info)
            if (luout /= luttyo) write (luttyo,2034) minerr(info)
            if (itrace /= 0)  write (lutrc,2034) minerr(info)
            return
        !                                         *** Normal return
        else
            fnorm=enorm(nptot,fvec)
            rmsdv=fnorm/dsqrt( dfloat(nptot) )
            sigma=fnorm/dsqrt( dfloat(nptot-nprm) )
            write(luout,1000)
            write(luout,2035) minerr(info),nfev,nprm,nptot
            write(luout,2036) rmsdv,fnorm*fnorm
            write(luout,1000)
        
            if (luout /= luttyo) then
                write(luttyo,2035) minerr(info),nfev,nprm,nptot
                write(luttyo,2036) rmsdv,fnorm*fnorm
            end if
            if (itrace /= 0) then
                write(lutrc,2035) minerr(info),nfev,nprm,nptot
                write(lutrc,2036) rmsdv,fnorm*fnorm
            end if
        end if
    
    !-----------------------------------------------------------------------
    !     Calculate residuals and covariance matrix
    !     COVAR calculates covariance matrix from the R matrix of QR
    !     factorization, which lmder stored in upper triangular form in fjac
    !     ipvt is permutation array returned by lmder
    !------------------------------------------------------------------------
    
        call covar( nprm,fjac,mxpt,ipvt,xtol*xtol,work1 )
    
    !----------------------------------------------------------------------
    !     Calculate correlation matrix from covariance and output it
    !----------------------------------------------------------------------
    
        line=' '
        mrgn=max(1,36-4*nprm)
        if (mrgn < 1) mrgn=1
        write (luout,2000)
        write(line(mrgn:),2002) (tag(i),i=1,nprm)
        write(luout,2001) line(:8*(nprm+1)+mrgn)
    
        do 23 i=1,nprm
            line=' '
            do 15 j=i,nprm
                if (i == j) xerr(i)=sigma*dsqrt( fjac(i,i) )
                denom = dsqrt( fjac(i,i) * fjac(j,j) )
                if (denom /= 0.0D0) corr(i,j) = fjac(i,j)/denom
                k=8*(j-1)+mrgn
                write(line(k:),2003) corr(i,j)
            15 END DO
            write (luout,2001) line(:8*(nprm+1)+mrgn)
        23 END DO
    
    !----------------------------------------------------------------------
    !     output final fit of parameters
    !----------------------------------------------------------------------
        if ( .NOT. (info == 0 .OR. info > 3)) then
            write (luout,1000)
            write (luout,2004)
            write (luout,2007) (tag(i),xxx(i),xerr(i),i=1,nprm)
            write (luout,*)
        end if
    
    !----------------------------------------------------------------------
    !     For fit of spherical tensor components, report Cartesian tensor
    !     and propagate uncertainty estimates
    !----------------------------------------------------------------------
    
    !       g-tensor
    
    !        if (gflag.and.igflg.ne.1) then
    !           write(luout,1060) gxx,gyy,gzz
    !        end if
    
    !       A-tensor
    
    !        if (aflag.and.iaflg.ne.1) then
    !           write(luout,1061) axx,ayy,azz
    !        end if
    
    !       R-tensor
    
    !       Calculate uncertainty in each rate from uncertainty in log(rate)
    
    !        if (rotflg) then
    
    !           if (irflg.ne.1) write(luout,1062) dx,dy,dz
    
    !      end if
    
    !         if (potflg)  call ordrpr(c20,ordpar)
    
    end if

!----------------------------------------------------------------------
!     Output fit (or calculated spectrum) for each datafile
!----------------------------------------------------------------------
    do 40 i=1,nser
    
    ! *** 2D-output format **
    
    !---------------------------------------------------------------------
    !     Preparation of graph (Frequencies in MHz)
    !---------------------------------------------------------------------
    
        amin=0.0d0
        amax=0.0d0
    !                                   * unscale the spectrum
        scfac=1.0d0
        if (dataOK) then
            scfac=1.0d0/sratio(i)
            do 45 j=1,ndata(i)
                k=ixsp(i)+j-1
                totspec(k)=(fvec(k)/datwgt(i)+data(k))*scfac
                if ( totspec(k) > amax ) amax=totspec(k)
            45 END DO
        else
        ! totspec already calculated, use it
            do 46 j=1,ndata(i)
                k=ixsp(i)+j-1
                if ( totspec(k) > amax ) amax=totspec(k)
            46 END DO
        end if
    
    !         npt1=snpt1(i)
    !         npt2=snpt2(i)
    !         inif2=-5.0d2*dble(snpt2(i)-1)/sstept2(i)/dble(snpt2(i))
    !         res2=-2.0d0*inif2/(snpt2(i)-1)
    !	 inif2=-5.0d2/sstept2(i)
    !         inif1=-5.0d2*dble(snpt1(i)-1)/sstept1(i)/dble(snpt1(i))
    !         res1=-2.0d0*inif1/(snpt1(i)-1)
    !	 inif1=-5.0d2/sstept1(i)
    !                              * obtain output file name
        call setdat (dataid(i))
    
    !         call wrfit( totspec(ixsp(i)),npt1,npt2,inif2,res2,
    !     #               inif1,res1,amin,amax,ftname,lthdnm )
        call wrfit( totspec(ixsp(i)),i,amin,amax,ftname )

    
    40 END DO

!----------------------------------------------------------------------
!     Close files and return to calling program
!----------------------------------------------------------------------
    if (itrace /= 0) close (lutrc)
    return

!# format statements ##################################################

    1000 format(/,2x,70('-'),/)
    1040 format(/10x,'*** Some datafiles are not read yet. ***', &
    /10x,'*** Fit is independent of the data. ***')
    1046 format(10x,'Single calculation : RMS deviation =',g14.7)
! 1060 format('Cartesian g-tensor = (',2(f10.7,','),f10.7,')')
! 1061 format('Cartesian A-tensor = (',2(f9.3,','),f9.3,')')
! 1062 format(10x,'Rx =',e13.6,'   Ry =',e13.6,'   Rz =',e13.6)
    2000 format(/24x,'*** Correlation Matrix ***'/)
    2001 format(a)
    2002 format(8(1x,a6,1x))
    2003 format(f8.4)
    2004 format(/9x,'*** Final Parameters ***',/, &
    7x,' Parameter      Value'/7x,26('-')/)
    2007 format(10x,a6,' = ',g13.7,' +/- ',g13.7)
    2034 format(/2x,'MINPACK could not find a solution: ', a )
    2035 format(/9x,'MINPACK completed: ', a/ &
    12x,'Function evaluations: ',i4/ &
    12x,'There were ',i2,' parameters and ',i6,' data points')
    2036 format(12x,'Rms Deviation = ',g13.6/12x,'chi-squared = ',g13.6/)

    end subroutine fitp
