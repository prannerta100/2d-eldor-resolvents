!  VERSION 1.0  (NLSPMC version)   2/5/99
!**********************************************************************

!                        =================
!                        subroutine PCHECK
!                        =================

!       This procedure replaces the function of program LBLF & TDIP for
!       performing initial computations based on the input parameters
!       for the EPRP 2D calculation programs.  It assumes the floating
!       point and integer variables in the fparm and iparm arrays have
!       been copied into common block /eprprm/ and checks to make sure
!       they are valid.

!       Returns ierr=0 if parameters are valid;  a negative number
!       identifying a fatal error, and a positive number issuing warning
!       messages. (it is reset to zero before returning to the calling
!       routine.)

!       lumsg specifies a logical unit for output of any error messages.
!       None are produced if lumsg=0.

!       Fatal parameter errors
!          1) zero gxx, gyy, or gzz
!          2) zero B0

!       Includes:
!               nlsdim.inc
!               nlsnam.inc
!               maxl.inc
!               eprprm.inc
!               rndoff.inc
!       Uses:
!               ipar.f
!               cd2km.f
!               tocart.f (coded in tensym.f)
!               fbasis.f

!**********************************************************************

    subroutine pcheck(lumsg,ierr)

    implicit none

    include 'limits.inc'
    include 'names.inc'
    include 'maxl.inc'
    include 'simparm.inc'
    include 'rndoff.inc'

    integer :: lumsg,ierr

    integer :: i,j,itmp,inlemx,inlomx,inkmx,inmmx,inpnmx
    double precision :: d2km(2,5,5),tmp,gmax,gmin,hmax

    double precision :: zero,half,dsq23,one,ten
    parameter (zero=0.0D0,half=0.5D0,one=1.0D0,ten=1.0D1)
    parameter (dsq23=0.816496580927726d0)

    logical :: axiala,axialg

    integer :: ipar
    external ipar

!######################################################################

    ierr=0

!----------------------------------------------------------------------
!     initialize spherical tensors
!----------------------------------------------------------------------
    do 10 j=1,5
        faa(j)=zero
        fgm(j)=zero
        do 9 i=1,2
            fam(i,j)=zero
            fgd(i,j)=zero
            fad(i,j)=zero
        9 END DO
    10 END DO

!----------------------------------------------------------------------
!     Check for zero field
!----------------------------------------------------------------------
    if (b0 < zero) then
        if (lumsg /= 0) write (lumsg,1002)
        ierr=2
        b0=abs(b0)
    elseif (b0 <= rndoff) then
        if (lumsg /= 0) write (lumsg,2001)
        ierr=-1
        return
    endif

!----------------------------------------------------------------------
!     Convert g, A, and R tensors to Cartesian form
!----------------------------------------------------------------------
    call tocart( gxx, igflg )
    call tocart( axx, iaflg )
    call tocart( dx,  irflg )

!  *** NOTE: in NLS version, diffusion rates are specified as
!            a power of 10

!     set criterion for weighting factor to store eigenvector

    tmp=7.99
    do 30 i=1,6
        if ((dx > tmp) .AND. (dy > tmp) .AND. (dz > tmp)) then
            go to 40
        else
            tmp=tmp-1.0d0
        end if
    30 END DO

!  Wtol eliminated.  RC 9/5/96

    40 wtol=1.0D-9

    dx=ten**dx
    dy=ten**dy
    dz=ten**dz

!----------------------------------------------------------------------
!     Check for zero values in g-tensor and zero A tensor
!----------------------------------------------------------------------
    if ( abs(gxx) < rndoff  .OR. &
    abs(gyy) < rndoff  .OR. &
    abs(gzz) < rndoff)   then
        if (lumsg /= 0) write (lumsg,2002)
        ierr=-2
        return
    endif

    if ( (in2 <= 0) .OR. &
    (dabs(axx) < rndoff .AND. &
    dabs(ayy) < rndoff .AND. &
    dabs(azz) < rndoff)     ) then
        in2=0
        axx=zero
        ayy=zero
        azz=zero
        if (lumsg /= 0) write (lumsg,1001)
        ierr=1
    endif

    gamman=zero

!----------------------------------------------------------------------
!     Calculate spherical components of tensors in magnetic frames
!----------------------------------------------------------------------
!                                        *** Electronic Zeeman ***
    g0=(gxx+gyy+gzz)/3.0D0
    fgm(1)=half*(gxx-gyy)*(b0/g0)
    fgm(3)=dsq23*(gzz-half*(gxx+gyy))*(b0/g0)
    fgm(5)=fgm(1)
    axialg = dabs(fgm(1)).lt.rndoff
!                                        *** Hyperfine ***
    a0=(axx+ayy+azz)/3.0D0
    faa(1)=half*(axx-ayy)

    faa(3)=dsq23*(azz-half*(axx+ayy))
    faa(5)=faa(1)
    axiala = dabs(faa(1)).lt.rndoff
!                                        *** Nuclear Zeeman ***
    zeen=zero

    gmin=dmin1( gxx, gyy, gzz )
    gmax=dmax1( gxx, gyy, gzz )
    hmax=dmax1( dabs(axx), dabs(ayy), dabs(azz) )

!----------------------------------------------------------------------
! Issue warning if high-field approximation has been violated
! (this is a very rough criterion for the high-field approx.!)
!----------------------------------------------------------------------
    if (b0 < 10.0D0*dmax1( hmax, (gmax-gmin)*b0/g0) ) then
        if (lumsg /= 0) write (lumsg,1003)
        ierr=3
    endif
!----------------------------------------------------------------------
!  Set ipt and cpot array according to the potential coefficients given
!----------------------------------------------------------------------
    ipt=0
    do 240 j=1,5
        do 245 i=1,5
            cpot(i,j)=zero
        245 END DO
    240 END DO

    if (dabs(c20) > rndoff) then
        ipt=ipt+1
        cpot(2,1)=c20
        lptmx=2
    endif
    if (dabs(c22) > rndoff) then
        ipt=ipt+1
        cpot(2,2)=c22
        lptmx=2
        kptmx=2
    endif
    if (dabs(c40) > rndoff) then
        ipt=ipt+1
        cpot(3,1)=c40
        lptmx=4
    endif
    if (dabs(c42) > rndoff) then
        ipt=ipt+1
        cpot(3,2)=c42
        lptmx=4
        kptmx=2
    endif
    if (dabs(c44) > rndoff) then
        ipt=ipt+1
        cpot(3,3)=c44
        lptmx=4
        kptmx=4
    endif

    if (lptmx >= 2) then
        lband=lptmx*2
    else
        lband=2
    end if

    kband=kptmx*2
    if (dabs(dx-dy) > rndoff) kband=kband+2

!----------------------------------------------------------------------
! Check for consistency of specified motion parameters with given model
!----------------------------------------------------------------------
    if (ipdf > 2) ipdf=2
    if (ipdf < 0) ipdf=0
    if (ist < 0) ist=0

!     Set exponents according to non-Brownian model

    pl=one
    pkxy=one
    pkzz=one
    if (ml == 1)  pl=half
    if (mxy == 1) pkxy=half
    if (mzz == 1) pkzz=half

!     *** Non-Brownian motion: must have at least one nonzero residence
!         time, and no potential specified

! NOTE: in NLS version, R*residence time products are specified as
!       powers of ten

    if (ipdf /= 0) then
        pml = ten**pml
        pmxy = ten**pmxy
        pmzz = ten**pmzz
    
        if ((ml == 0) .AND. (mxy == 0) .AND. (mzz == 0) ) then
            if (lumsg /= 0) write (lumsg,1004)
            ierr=4
            ipdf=0
        end if
    
    !  *** NOTE: in NLS version, djf, djfprp specified as powers of ten
    
        djf=ten**djf
        djfprp=ten**djfprp
    
        if (ipt > 0) then
            if (lumsg /= 0) write (lumsg,1005)
            ierr=5
            ipdf=0
        end if
    end if

!     *** Anisotropic viscosity model: must have potential and
!         no discrete jump motion specified

    if (ipdf == 2) then
        if (ipt == 0) then
            if (lumsg /= 0) write (lumsg,1006)
            ierr=6
            ipdf=0
        end if
    
        if (ist > 0) then
            if (lumsg /= 0) write (lumsg,1007)
            ierr=7
            ipdf=0
        end if
    end if

!     *** Brownian motion model: set residence times to zero

    if (ipdf == 1) then
        if (ml == 0)  pml =zero
        if (mxy == 0) pmxy=zero
        if (mzz == 0) pmzz=zero
    end if

!     *** MOMD calculation: must have potential

    if (nort > 1) then
        if (nort > MXORT) then
            if (lumsg /= 0) write (lumsg,*)'nort too big, reset to ', &
            MXORT
            nort=MXORT
        end if
        if (ipt < 1) then
            if (lumsg /= 0) write (lumsg,1015)
            ierr=15
        end if
    end if

!----------------------------------------------------------------------
!     Set director tilt flag
!----------------------------------------------------------------------
    ipsi0=0
    if (nort > 1) ipsi0=1
    if((abs(psi) > rndoff) .AND. (abs(psi-180.0d0) > rndoff)) ipsi0=1

!----------------------------------------------------------------------
!     Set diffusion tilt flag
!----------------------------------------------------------------------
    bed=abs( bed )
    itd=0
    if ( dabs(ald) > rndoff .OR. dabs(gad) > rndoff &
     .OR. bed > rndoff) itd=1

!----------------------------------------------------------------------
!     Get magnetic tilt angles
!----------------------------------------------------------------------
    bem = dabs( bem )
    if (axialg .AND. dabs(alm) > rndoff) then
        if (lumsg /= 0) write (lumsg,1013)
        ierr=13
        alm=zero
    end if

!     *** Only beta tilt angles allowed if all tensors are axial

    if (axialg .AND. axiala .AND. (dabs(alm) > rndoff .OR. &
    dabs(gam) > rndoff .OR. dabs(gad) > rndoff) ) then
        if (lumsg /= 0) write (lumsg,1014)
        ierr=14
        alm=zero
        gam=zero
        gad=zero
    end if

    itm=0
    if (dabs(alm) > rndoff .OR. bem > rndoff &
     .OR. dabs(gam) > rndoff)     itm=1

!----------------------------------------------------------------------
!     Set A tensor in magnetic (g-tensor) frame
!----------------------------------------------------------------------
    if (itm == 0) then
    !                           *** no tilt: copy tensor directly
        do 252 j=1,5
            fam(1,j)=faa(j)
        252 END DO
    else
    !                          *** transform A tensor into g axis system
    
        call cd2km(d2km,alm,bem,gam)
        do 255 i=1,5
            do 254 j=1,5
                fam(1,i)=fam(1,i)+d2km(1,i,j)*faa(j)
                fam(2,i)=fam(2,i)+d2km(2,i,j)*faa(j)
            254 END DO
        255 END DO
    end if

!----------------------------------------------------------------------
!     Set F_g and F_A tensors in the diffusion frame
!----------------------------------------------------------------------
    if (itd == 0) then
    !                           *** no tilt: copy tensors directly
        do 260 j=1,5
            fgd(1,j)=fgm(j)
            do 259 i=1,2
                fad(i,j)=fam(i,j)
            259 END DO
        260 END DO
    
    else
    !                    *** Transform A and g tensors into the diffusion frame
    
        call cd2km(d2km,ald,bed,gad)
        do 265 i=1,5
            do 264 j=1,5
                fgd(1,i)=fgd(1,i)+d2km(1,i,j)*fgm(j)
                fgd(2,i)=fgd(2,i)+d2km(2,i,j)*fgm(j)
                fad(1,i)=fad(1,i)+d2km(1,i,j)*fam(1,j) &
                -d2km(2,i,j)*fam(2,j)
                fad(2,i)=fad(2,i)+d2km(1,i,j)*fam(2,j) &
                +d2km(2,i,j)*fam(1,j)
            264 END DO
        265 END DO
    end if

!----------------------------------------------------------------------
!     Relaxation parameters & Heisenberg exchange parameter
!----------------------------------------------------------------------

! NOTE: in NLS version, these rates are given in powers of ten

    if (t2edi > rndoff) t2edi=ten**t2edi
    if (t2ndi > rndoff) t2ndi=ten**t2ndi
    if (t2efi > rndoff) t2efi=ten**t2efi
    if (t1edi > rndoff) t1edi=ten**t1edi
    if (t1ndi > rndoff) t1ndi=ten**t1ndi
    if (oss > rndoff) oss=ten**oss

!----------------------------------------------------------------------
!     Check basis set parameters
!----------------------------------------------------------------------
!                             ** basis already set **
    if (setbas) go to 500

    inlemx = lemx
    inlomx = lomx
    inkmx  = kmx
    inmmx  = mmx
    inpnmx = ipnmx

    if((lemx > mxlval)) lemx=mxlval
    if(ipar(lemx) /= 1) lemx=lemx-1
    if(ipar(lomx) /= -1) lomx=lomx-1
    if(lomx > lemx) lomx=lemx-1
    if(kmx > lemx) kmx=lemx
    if(mmx > lemx) mmx=lemx
    if(ipar(kmx) /= 1) kmx=kmx-1
    if(ipnmx > in2) ipnmx=in2
    if((ipsi0 == 0) .AND. (mmx > ipnmx+1)) mmx=ipnmx

    if(lemx < 0)  lemx=0
    if(lomx < 0)  lomx=0
    if(kmx < 0)   kmx=0
    if(mmx < 0)   mmx=0
    if(ipnmx < 0) ipnmx=0

    if (inlemx /= lemx .OR. &
    inlomx /= lomx .OR. &
    inkmx  /= kmx  .OR. &
    inmmx  /= mmx  .OR. &
    inpnmx /= ipnmx ) then
    
        if (lumsg /= 0) write (lumsg,1009) lemx,lomx,kmx,mmx,ipnmx
        ierr=9
    endif

!----------------------------------------------------------------------
!  Determine basis set using rules in M,I,I,M & Freed
!  (1) jm=1 (use only symmetric M combinations) for no nuclear Zeeman
!  (2) jkmn=1 (use only symmetric K combinations) if the magnetic
!          tensors in the diffusion frame are real-valued. This is the
!          case if alpha_m, gamma_m (magnetic tilt), and gamma_d
!          (diffusion tilt) are all zero.
!  (3) only even K values if there is no magnetic and diffusion tilt.
!  (4) only even L values no K values (kmx=0) in case of axial magnetic
!          tensors, axial potential, and no magnetic/diffusion tilt.
!----------------------------------------------------------------------
    jmmn=1
    if (abs(zeen) > rndoff) jmmn=-1

    jkmn=1
    do 270 i=1,5
        if (    dabs(fgd(2,i)) >= rndoff &
         .OR. dabs(fad(2,i)) >= rndoff ) jkmn=-1
    270 END DO

    if (itm == 0 .AND. itd == 0) then
        kdelta=2
    else
        kdelta=1
    end if

    if (axiala .AND. axialg .AND. kdelta == 2 .AND. kptmx == 0) &
    then
        ldelta=2
        lomx=0
        kmx=0
    else
        ldelta=1
    end if

!----------------------------------------------------------------------
!     check the calculation type and related information
!----------------------------------------------------------------------
    500 continue

    if (nstep > MXSTEP-2) then
        if (lumsg /= 0) write (lumsg,*) &
        'adjusting nstep to ',MXSTEP-2
        write(*,*)'adjusting nstep to ',MXSTEP-2
        nstep=MXSTEP-2
    end if

    if ((itype /= 1) .AND. (itype /= 4)) then
        if (lumsg /= 0) write (lumsg,1016)
        ierr=16
        itype=4
    end if
!                                   *** eprbf calculation ***
    if (itype == 2) then
        if (nstep == 0) then
            if(ndimo > 0) then
                if (lumsg /= 0) write (lumsg,1017)ndimo
                ierr=17
                nstep=ndimo
            else
                if (lumsg /= 0) write (lumsg,1170)
                ierr=-17
                return
            end if
        end if
    
        if (ptol < rndoff) then
            if (lumsg /= 0) write (lumsg,1018)
            ierr=18
            ptol=1.0d-2
        end if
    
        if ((abs(fieldi) < 1.0d0) .AND. (abs(fieldi) < 1.0d0)) then
            if (lumsg /= 0) write (lumsg,1019)
            ierr=19
            fieldi=-50.0d0
            fieldf=-50.0d0
        end if
    
        if (nfield < 2) then
            if (lumsg /= 0) write (lumsg,1020)
            ierr=20
            nfield=10
        end if
    
    end if
!                                   *** reset ierr flag ***
    if (ierr > 0) ierr=0
    return

!======================================================================

! Formats for warning and error messages
!----------------------------------------------------------------------
! Warning messages

    1001 format(' Nuclear spin in2=0 assumed')
    1002 format(' Positive B0 value assumed')
    1003 format(' Warning: high-field approximation may not apply for', &
    ' given parameters')
    1004 format(' No non-Brownian motion specified for l, xy, zz:', &
    ' ipdf=0 assumed')
    1005 format(' Nonzero potential with jump/free model:', &
    ' ipdf=0 assumed')
    1006 format(' Zero potential with anisotropic viscosity:', &
    ' ipdf=0 assumed')
    1007 format(' Discrete jumps with anisotropic viscosity:', &
    ' ipdf=0 assumed')
    1008 format(' lemx must be 48 or less with potential;', &
    ' lemx=48 assumed')
    1009 format(' Basis set parameters adjusted to ',4(i3,','), i3)
    1010 format(' Too many Lanczos/CG steps specified: nstep=',i3, &
    ' assumed ')
    1011 format(' Specified cgtol too small: cgtol=',g9.3,' assumed')
    1012 format(' Specified shiftr too small: shiftr=',g9.3,' assumed')
    1013 format(' Axial g-tensor: alm=0 assumed')
    1014 format(' Axial g, hf-tensors: alm=gam=gad=0 assumed')
    1015 format(' Zero potential with nort > 1 attempted')
    1016 format(' Calculation type: Rutishauser assumed')
    1017 format(' Zero CG steps: nstep=ndimo assumed',I8)
    1170 format(' Zero CG steps: fatel error in pcheck')
    1018 format(' Zero pruning tolerance for EPRBF: ptol=0.01 assumed')
    1019 format(' Zero field range for EPRBF: -50 -- +50 assumed')
    1020 format(' Zero field points for EPRBF: nfield=10 assumed')

!----------------------------------------------------------------------
! Fatal error messages

    2001 format(' Fatal parameter error: zero B0')
    2002 format(' Fatal parameter error: zero values in g-tensor')

!======================================================================

    end subroutine pcheck


!----------------------------------------------------------------------
!                      =========================
!                        subroutine PRMSOK
!                      =========================

!  Checks all the parameter blocks up to the currently defined
!  number of spectra and sites to see that:

!   (1) B0 is nonzero
!   (2) gxx,gyy,gzz are all nonzero (Cartesian tensor) or
!       at least one of g1,g2,g3 is nonzero (spherical tensor)
!   (3) Dimension of the matrix is valid
!   (4) nort is valid

!   Returns .false. if any of these conditions is violated for any
!   of the given sites.
!----------------------------------------------------------------------
    function prmsOK( lu )

    use basis        
    implicit none
    logical :: prmsOK
    integer :: lu

    include 'limits.inc'
    include 'simparm.inc'
    include 'datas.inc'
    include 'parms.inc'
    include 'miscel.inc'
!    include 'basis.inc'
    include 'rndoff.inc'

    integer :: j
    logical :: allOK,gOK,dimOK,dimdOK,ortOK

    integer :: CARTESIAN,SPHERICAL,AXIAL,isite,ispec
    parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)

    gOK=.true.
    allOK=.true.
    dimOK=.true.
    dimdOK=.true.
    ortOK=.true.
! total basis set storage check:
    dimOK = dimOK .and. ((pidptr(nbasis)+ndimoss(nbasis)-1.gt.0) &
    .and. (pidptr(nbasis)+ndimoss(nbasis)-1.le.MMXDIM))
    dimdOK = dimdOK .and. ((dpidptr(nbasis)+ndimdss(nbasis)-1.gt.0) &
    .and. (dpidptr(nbasis)+ndimdss(nbasis)-1.le.MMXDIM))

!  for all sites/spectra:

    do 10 ispec=1,nspectra
        do 10 isite=1,ncomps
            allOK = allOK .and. (dabs(fparm(IB0,ispec,isite)).gt. &
            rndoff)
        
        !   Check g-tensors. For Cartesian tensor, all elements must be nonzero.
        !   For spherical tensor, at least one element must be nonzero.
        
            if (iparm(IIGFLG,ispec,isite) == CARTESIAN) then
                gOK = gOK .and. (dabs(fparm(IGXX,ispec,isite)).gt. &
                rndoff .and. dabs(fparm(IGXX+1,ispec,isite)).gt. &
                rndoff .and. dabs(fparm(IGXX+2,ispec,isite)).gt.rndoff)
            else if (iparm(IIGFLG,ispec,isite) == SPHERICAL) then
                gOK = gOK .and. (dabs(fparm(IGXX,ispec,isite)).gt. &
                rndoff .or. dabs(fparm(IGXX+1,ispec,isite)).gt. &
                rndoff .or. dabs(fparm(IGXX+2,ispec,isite)).gt.rndoff)
            else if (iparm(IIGFLG,ispec,isite) == AXIAL) then
                gOK = gOK .and. (dabs(fparm(IGXX,ispec,isite)).gt. &
                rndoff .and. dabs(fparm(IGXX+1,ispec,isite)).gt.rndoff)
            end if
        
            allOK = allOK .and. gOK
        !  Individual basis size check:
            dimOK = dimOK .and. ((ndimoss(basinfo(1,ispec,isite)).gt.0) &
            .and. (ndimoss(basinfo(1,ispec,isite)).le.mxdim))
            dimdOK = dimdOK.and.((ndimdss(basinfo(1,ispec,isite)).gt.0) &
            .and. (ndimoss(basinfo(1,ispec,isite)).le.mxdim))
            allOK = allOK .and. dimOK .and. dimdOK
        
            ortOK = ortOK .and. ((iparm(INORT,ispec,isite).gt.0).and. &
            (iparm(INORT,ispec,isite).le.mxort))
            allOK = allOK .and. ortOK
    10 END DO

    if ( .NOT. allOK) write (lu,1000)
    if ( .NOT. gOK) write (lu,1001)
    if (diaflg) dimOK = dimOK .AND. dimdOK
    if ( .NOT. dimOK) write (lu,1002)
    if ( .NOT. ortOK) write (lu,1003)

    prmsOK=allOK
    return

    1000 format('*** B0 is zero for spectrum ***')
    1001 format('*** Zero g-tensor value set ***')
    1002 format('*** Dimension of the matrix is not valid ***')
    1003 format('*** Number of orientations for MOMD is not valid ***')
    end function prmsOK


!----------------------------------------------------------------------
!                      =========================
!                        subroutine EVALOK
!                      =========================

!  Checks if real part of all the eigenvalues is positive.

!----------------------------------------------------------------------
    function evalOK( ips,eval,nvar,lu )

    implicit none
    logical :: evalOK
    integer :: ips,nvar,lu,luold,i,j,nev,ii

    include 'limits.inc'
    include 'stdio.inc'
    include 'simparm.inc'
    include 'parms.inc'
    include 'rndoff.inc'

    double precision :: tol
    complex*16 eval(mxdim)

    evalOK=.true.
    tol=rndoff*1.0d9
    luold=lu
    nev=nevo
    if (ips == 0) nev=nevd

! RC - Modify code to set neg eigenvalues to zero and continue:

!      do 100 i=1,nev
!         if ( dreal(eval(i)).lt.-tol ) then
!            if (ips.eq.0) then
!               write (lu,1001) psi
!               if (lu.ne.luttyo) write (luttyo,1001) psi
!            else
!               write (lu,1002) psi
!               if (lu.ne.luttyo) write (luttyo,1002) psi
!            end if
!            write (lu,1003) (j,eval(j),j=1,nev)
!            if (nvar.gt.0) then
!               write (lu,1004)
!               write (lu,1005) (tag(j),fparm(ixpr(j)),j=1,nvar)
!            end if
!            go to 20
!         end if
! 100  continue
!      return

! 20   negegv=negegv+1
!      evalOK=.false.
!      lu=luold
!      return

!             write(*,*)(ii,eval(ii),ii=1,nev)
    do 100 i=1,nev
        if ( dreal(eval(i)) < -tol ) then
            if (ips == 0) then
                write (lu,2001) i,psi,dreal(eval(i)), &
                dimag(eval(i))
                if (lu /= luttyo) write (luttyo,2001) i,psi,dreal(eval(i)), &
                dimag(eval(i))
            !             write(*,*)(ii,eval(ii),ii=1,nev)
                eval(i)=(0.0d0,0.0d0)
            !             eval(i)=(0.0d0,dimag(eval(i)))
            else
                write (lu,2002) i,psi,dreal(eval(i)), &
                dimag(eval(i))
                if (lu /= luttyo) write (luttyo,2002) i,psi,dreal(eval(i)), &
                dimag(eval(i))
            !             write(*,*)(ii,eval(ii),ii=1,nev)
                eval(i)=(0.0d0,0.0d0)
            !             eval(i)=(0.0d0,dimag(eval(i)))
            end if
        end if
    100 END DO
    return

!## format statements ###############################################

    1001 format(/5x,'*** Negative real part in psi = ',f5.2, &
    ' of DIAGONAL eigenvalues ***')
    1002 format(/5x,'*** Negative real part in psi = ',f5.2, &
    ' of OFFDIAGONAL eigenvalues ***')
    2001 format(/5x,'*** Neg. real at psi = ',i7,f6.2, &
    ' of DIAG. egvs. *** ',2g15.5,/,'Forced to zero')
    2002 format(/5x,'*** Neg. real at psi = ',i7,f6.2, &
    ' of OFFDIAG. egvs. *** ',2g15.5,/,'Forced to zero')
    1003 format(8x,i4,4x,g14.7,2x,g14.7)
    1004 format(/8x,'EPR parameters')
    1005 format(11x,a,4x,g16.9)

    end function evalOK
