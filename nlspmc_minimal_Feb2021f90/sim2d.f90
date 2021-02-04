!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       ==================
!                       SUBROUTINE : SIM2D
!                       ==================

!     Subroutine version of MOMD program

!     This routine is intended for use with nonlinear least-squares
!     applications.  The routine calculates 2D-spectrum for the
!     parameters given in /eprprm/ using the eigenvalues and eigenvectors
!     already obtained by evcgf routine.

!     On Entry :

!        icalc : = 0  recalculate the spectrum using previous xspec.
!                     (valid when varying gib or lib parameters)
!                = 1  recalculate the spectrum using the previous
!                        eigenvalues and eigenvectors.
!                     (valid for series of 2D spectra with different
!                        experimental type, combination, mixing time
!                        and so on
!                      also valid when varying hwid parameter)
!                = 2  recalculate the eigenvalues and the spectrum
!        xspec : complex time-domain spectrum (used for icalc=0 case)

!     On Exit  :

!        xspec : complex time-domain spectrum (to be used with icalc=0)
!        cspec  : complex frequency spectrum
!        ierr  : error flag
!                      = 0   normal return
!                      < 0   critical error in one of
!                            pcheck ; wrong parameter specification
!                            matrx? ; wrong matrix generation
!                            cmtqli ; failure in QL-decomposition
!     Includes :
!               nlsdim.inc
!               stdio.inc
!               eprprm.inc
!               prmeqv.inc
!               parcom.inc
!               stvcom.inc
!               egvcom.inc

!     Uses :
!               evcgf.f
!               spcalc.f
!               convft.f

!*********************************************************************

    subroutine sim2d(icalc,xspec,cspec,npt1,npt2,nvar, &
    iflag,ixi,delx,ispec,isite,ierr)

    implicit none

    include 'limits.inc'
    include 'stdio.inc'
    include 'simparm.inc'
    include 'parmequ.inc'
    include 'parms.inc'
    include 'stvcom.inc'
    include 'basis.inc'
    include 'egvcom.inc'
    include 'miscel.inc'

    double precision :: zero,scal,pi
    complex*16 czero,ci
    parameter (zero=0.0D0,czero=(0.0D0,0.0D0),ci=(0.0d0,1.0d0))

    integer :: npt1,npt2,icalc,i,j,k,l,iort,bdw, &
    ntmp,nvar,iflag,ixi,ierr,ispec,isite
    integer :: tt1,tt2,tt3,nbindx
    double precision :: t1,t2,cspsi,amin,amax,delx,flg,norm

    complex*16 xspec(npt1,npt2),cspec(npt1,npt2),bubba
!      real*8 spec(npt1,npt2)

    logical :: evalOK
    external evalOK

!#####################################################################


    ierr=0
    pi=4.0d0*datan(1.0d0)

!----------------------------------------------------------------------
!     Load values from parameter array into /eprprm/
!----------------------------------------------------------------------

    10 continue

! Move the current ispec, isite parameters into the fepr and iepr
! arrays.

    do 2 i=1,nfprm
        fepr(i)=fparm(i,ispec,isite)
    !         write(*,*)'i,fepr ',i,fepr(i)
    2 END DO

! only copy part of eprprm, NIPRM=23 ???  Rest must not matter?
! same above too...

    do 4 i=1,niprm
        iepr(i)=iparm(i,ispec,isite)
    !         write(*,*)'i,iepr ',i,iepr(i)
    4 END DO
! record site and specturm being calculated:
!      specnum=ispec
!      sitenum=isite

! Move the basis set into position for this simulation if necessary:

    if (basinfo(1,ispec,isite) /= basisptr) then
        nbindx=basinfo(1,ispec,isite)
        basisptr=nbindx
        ndimo=ndimoss(nbindx)
        ndimd=ndimdss(nbindx)
        j=0
        tt1=pidptr(nbindx)
        tt2=pidptr(nbindx)+ndimoss(nbindx)-1
        do 41 i=tt1,tt2
        !        do 41 i=pidptr(nbindx),pidptr(nbindx)+ndimoss(nbindx)-1
            j=j+1
            jqe1(j)=mjqe1(i)
            pi1(j)=mpi1(i)
            qi1(j)=mqi1(i)
            l1(j)=ml1(i)
            jk1(j)=mjk1(i)
            k1(j)=mk1(i)
            jm1(j)=mjm1(i)
            m1(j)=mm1(i)
        41 END DO
        j=0
        tt1=dpidptr(nbindx)
        tt2=dpidptr(nbindx)+ndimdss(nbindx)-1
        do 42 i=tt1,tt2
        !        do 42 i=1+dpidptr(nbindx),dpidptr(nbindx)+ndimdss(nbindx)-1
            j=j+1
            djqe1(j)=mdjqe1(i)
            dpi1(j)=mdpi1(i)
            dqi1(j)=mdqi1(i)
            dl1(j)=mdl1(i)
            djk1(j)=mdjk1(i)
            dk1(j)=mdk1(i)
            djm1(j)=mdjm1(i)
            dm1(j)=mdm1(i)
            pp(j)=mpp(i)
        42 END DO
    ! rite the pulse propagator to pprop.txt
        open(unit=9,file='pprop.txt',status='unknown', &
        access='sequential')
        j=0
        do 669 i=tt1,tt2
            j=j+1
            write(9,*) pp(j)
        669 END DO
        close(unit=9)
                
        do 420 i=1,ndimoss(nbindx)
            pid(i)=mpid(i+pidptr(nbindx)-1)
        420 END DO

        open(unit=9,file='pid.txt',status='unknown', &
        access='sequential')
        do 421 i=1,ndimoss(nbindx)
            write(9,*) pid(i)
        421 END DO
        close(unit=9)

    end if
    if (icalc >= 1) then
    !                                    *** test validity of parameters
        if (basinfo(1,ispec,isite) == 0) then
            write(luout,1063)
            if (luout /= luttyo) write(luttyo,1063)
            1063 format(/5x,'** ERROR, basis not set in SIM2D **')
            return
        else
            setbas=.true.
        end if
        call pcheck(luttyo,ierr)
        if (ierr < 0) then
            write(luout,1000)
            if (luout /= luttyo) write(luttyo,1000)
            write(luout,1100)
            if (luout /= luttyo) write(luttyo,1100)
            return
        end if
    !                                    *** off-diagonal starting vector
    ! get stvo
    
        call stveco
    
    end if

!---------------------------------------------------------------------
!     calculate eigenvalues & eigenvectors, egval and egvec
!     Only do this once for each site.
!---------------------------------------------------------------------

    if (icalc >= 2) then
    !                        ===============================
    !                     === loop over orientations (MOMD) ===
    !                        ===============================
        do 100 iort=1,nort
        
            if (nort > 1) then
                cspsi=dfloat(iort-1)/dfloat(nort-1)
                psi=dacos(cspsi)*1.8D2/pi
            end if
        
            if (idebug /= 0) then
                write (ludeb,2000) iort,psi
                write (ludeb,2001)
            end if
        !                     *** solve for eigenvalues & vectors for pS=1
            call evcgf(1,stvo,nevo,egvalx(1,iort,isite), &
            egvecx(1,1,iort,isite),ierr) ! nevo returned is # of imp't evals.
        


            if (ierr < 0) then
                write(luout,1100)
                if (luout /= luttyo) write(luttyo,1100)
                return
            end if
        
            nev1(iort,isite)=nevo
            do 105 i=1,nevo
                egvalx(i,iort,isite)=egvalx(i,iort,isite)+ci*b0
            105 END DO
        !     #		egvalx(2,1,isite),egvecx(1,1,1,isite)
        !	stop
        
            if ( .NOT. evalOK(1,egvalx(1,iort,isite),nvar,luout) ) then
            
            !   ###                                                            ###
            !   ### This is a temporary fix for negative eigenvalues occuring  ###
            !   ### in function evaluation.  If this problem happens in OFF-   ###
            !   ### diagonal space, it USUALLY means the basis set being used  ###
            !   ### is not properly pruned yet and the program issues an error ###
            !   ### message and exit.  If the eigenvalues of the diagonal sub- ###
            !   ### space have negative real part, the first component of the  ###
            !   ### R tensor is increased by 1.0d-4 and the spectra are re-    ###
            !   ### calculated.  This seems to be related with the pruning     ###
            !   ### scheme for the diagonal basis and further development for  ###
            !   ### a new scheme is necessary to fix this problem.             ###
            !   ###                     July 19, 1993       Sanghyuk Lee       ###
            !   ###                                                            ###
            
                write (luout,1012)
                if (luout /= luttyo) write (luttyo,1012)
                ierr=-1
                return
            
            end if
        !                     *** solve for eigenvalues & vectors for pS=0
            if (diaflg .OR. (iexp > 2)) then
            
                if (idebug /= 0) write (ludeb,2011)
            
                call stvecd(egvalx(1,iort,isite),egvecx(1,1,iort,isite))
                call evcgf(0,stvd,nevd,egvalz(1,iort,isite), &
                egvecz(1,1,iort,isite),ierr)


            !            write(luttyo,*) "EVECZ^T*EVECZ at iort=1, isite=1, new"
            !            norm=0.0d0
            !            do 3657 i=1,nevd
            !                do 3658 j=1,nevd
            !                   bubba=czero!tmpry variable, stores the dot product
            !                    do 3659 k=1,ndimd
            !                        bubba=bubba+egvecz(k,i,1,1)*egvecz(k,j,1,1);
            ! 3659               continue
            !                    write(luttyo,*) i,j,bubba
            !                    flg=0.0d0
            !                    if(i.eq.j) flg=2.0d0
            !                    norm=norm+abs(bubba-flg)**2
            ! 3658           continue
            ! 3657       continue
            !            write(luttyo,*) "norm square = ",norm
                         
                stop












                if (ierr < 0) then
                    write(luout,1100)
                    if (luout /= luttyo) write(luttyo,1100)
                    return
                end if
            
                nev0(iort,isite)=nevd
            ! these conditions also eliminated by forcing negative eigenvalues
            ! to zero:
                if( .NOT. evalOK(0,egvalz(1,iort,isite),nvar,luout))then
                    write(luout,*)'evalOK error in egvalz, stopping'
                    if ( .TRUE. ) stop
                    if (iflag == 1) then
                        fparm(IDX,ispec,isite)=fparm(IDX,ispec,isite) &
                        +1.0d-4
                        write (luout,1014)
                        if (luout /= luttyo) write (luttyo,1014)
                    else
                        write (luout,1016)
                        fparm(ixi,ispec,isite)=fparm(ixi,ispec,isite) &
                        +9.0d0*delx
                        delx=delx*10.0d0
                        if (luout /= luttyo) write (luttyo,1016)
                    end if
                    go to 10
                end if
            
            end if
        
        100 END DO
    
    end if

!---------------------------------------------------------------------
!     calculate complex 2D spectrum in time domain
!---------------------------------------------------------------------

    if (icalc >= 1) then
    
        do 200 j=1,npt2
            do 200 i=1,npt1
                xspec(i,j)=czero
        200 END DO
    !                        ===============================
    !                     === loop over orientations (MOMD) ===
    !                        ===============================
        do 210 iort=1,nort
            scal=1.0d0
            if ((nort > 1) .AND. ((iort == 1) .OR. (iort == nort))) &
            scal=0.5d0
            nevo=nev1(iort,isite)
            nevd=nev0(iort,isite)
            call spcalc(scal,egvalx(1,iort,isite),egvalz(1,iort,isite), &
            egvecx(1,1,iort,isite),egvecz(1,1,iort,isite), &
            xspec,npt1,npt2)
        210 END DO
    !	stop

    
    end if

!---------------------------------------------------------------------
!     Calculate magnitude 2D spectrum with Gaussian and/or Lorentzian
!     inhomogeneous Broadening.  Don't modify xspec.
!---------------------------------------------------------------------

    call convft(xspec,cspec,npt1,npt2)
    return

!=====================================================================
!     format statements
!=====================================================================

    1000 format(/5x,'** Calculation attempted with illegal ', &
    'parameters **')
    1012 format(/5x,'** Basis set with looser pruning criterion is ', &
    'recommended **')
    1014 format(/8x,'First diffusion component is increased by 1e-4',/, &
    'This should not happen')
    1016 format(/8x,'Forward step size for current Jacobian evaluation ', &
    'is multiplied by 10'/)
    1100 format(/5x,'** Critical ERROR in SIM2D routine **')

    2000 format(/70('#'),/5x,'Orientation ',i2,'  : psi = ',f5.2, &
    /70('#'))
    2001 format(/'### OFF-DIAGONAL ###')
    2011 format(/'### DIAGONAL ###')

!=====================================================================

    end subroutine sim2d
