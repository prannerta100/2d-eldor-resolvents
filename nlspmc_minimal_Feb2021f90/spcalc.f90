!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       =================
!                       SUBROUTINE SPCALC
!                       =================

!     Subroutine version of EPRPS spectral calculation routine.

!     This routine is intended for use with nonlinear least-squares
!     applications.  The routine calculates 2D-spectrum for the
!     parameters given in /eprprm/ using the eigenvalues and eigenvectors
!     already obtained by evesa routine.

!     !!!  No nuclear modulation allowed here.  !!!

!     On Entry :

!        scal  :  scale factor for MOMD
!                 (1/2 for end points in Simpson's integration)
!        evalx,evalz,evecx,evecz   :  eigenvalues and eigenvectors
!                 for off-diagonal & diagonal subspace

!     On Exit  :
!        xspec :  complex 2D spectrum in time domain

!     Includes :
!               nlsdim.inc
!               eprprm.inc
!               stvcom.inc
!               wkspcm.inc
!               physcn.inc

!     Uses :
!               zgemul, zgemv   (ESSL library routines)
!               ZGEMV and ZGEMM from IMSL replace the ESSL libs.

!*********************************************************************


    subroutine spcalc(scal,evalx,evalz,evecx,evecz,xspec, &
    npt1,npt2)

    implicit none

    include 'omp_lib.h'
    include 'limits.inc'
    include 'simparm.inc'
    include 'stvcom.inc'
    include 'wkspcm.inc'
    include 'physcn.inc'
    include 'stdio.inc'
    include 'basis.inc'
!      include 'eprmat.inc'

    double precision :: zero,gtf,cfact,scal,twopi,omshft,om1
    parameter (zero=0.0d0)
    integer :: it3,mdft,ndft,mpol,nn,mxgrd
    complex*16 czero,ci,cone,tmpry
    parameter (czero=(0.0D0,0.0D0),ci=(0.0D0,1.0D0), &
    cone=(1.0d0,0.0d0),mxgrd=128)

    integer :: i,j,k,it1,it2,npt1,npt2,ntmp,ll
    double precision :: wid,t1,t2,f1,f2,ltfix, &
    linit1,lstept1,linit2,lstept2

    complex*16 sum1,sum2
    complex*16 evalx(mxegv),evalz(mxegv), &
    levalx(mxegv),levaly(mxegv),levalz(mxegv), &
    evecx(mxdim,mxegv),evecy(mxdim,mxegv), &
    evecz(mxdim,mxegv),xspec(npt1,npt2), &
    plsmix(mxegv,mxegv),w(mxdim),yr(mxdim)
    complex*16 spec_pS1(mxdim,mxgrd),ip(mxdft),xl(mxdim),xr(mxdim)
    complex*16 bubba,spec_pS1_cnorm(mxgrd)
    double precision :: flg,norm,ermin,error
    integer :: ierr,ndone,nmin
!        ermin,error,ierr,ndone,nmin,w
    double precision :: omT,lstepT,delta
    complex*16 cdft,sdft,corre,corim,sinint,cosint,endpts(8),fac1,fac2
!   cubic order correction, so 8 elements in the endpts array (NR 13.9)
    double precision :: cpol(mxpol),spol(mxpol),xpol(mxpol),a,b,c,s, &
    cerr,serr,corfac,en,cdftre,sdftre, &
    cpolt(mxpol),xpolt(mxpol)
    integer :: id,num,nthreads
!              fft, dftcor, polint to be called
! spec_pS1 now normalized spec_pS1, cnorms in spec_pS1(it1)


! Halal gyro ratio g0*betae/hbar stored in cfact; not sure if g0 ahould be multiplied or not
          

!#####################################################################

    twopi=8.0d0*datan(1.0d0)
    gtf=g0*betae/(hbar*twopi)*1.0d-6
    cfact=twopi*gtf
    omshft=twopi*shft

    ltfix=tfix*1.0d-3
    linit1=init1*1.0d-3
    lstept1=stept1*1.0d-3
    if (iexp >= 1) then
        linit2=init2*1.0d-3
        lstept2=stept2*1.0d-3
    end if

! RC modification 11/18/96 to allow T2 to vary with mixing time
! with exponential dependence.  To do this, the old variables
! hwid and mwid have been redefined.   The old equation was:
!      wid=hwid+(mwid*twopi/gtf)*ltfix
! In the revised equation below,  hwid becomes the increase in line
! broadening (Gauss) between Tmix=0 and Tmix=infinity in usec and
! mwid describes how fast that change takes place with Tmix.
! mwid is in micro seconds.  To provide for a constant
! wid, use t2edi which has a similar effect.
!       given in run file. [Gauss]
    if (abs(mwid) > 1.0D-10) then
        wid=abs(hwid)*(1.0d0-dexp(-ltfix/abs(mwid)))
    else
        wid=0.0d0
    endif

    do 10 i=1,nevo
        levalx(i)=(evalx(i)+wid)*cfact-ci*omshft
        levaly(i)=dconjg(levalx(i))
    10 END DO

    do 20 j=1,nevo
        do 20 i=1,ndimo
            evecy(i,j)=dconjg(evecx(i,j))
    20 END DO

    if (iexp >= 3) then
        do 30 i=1,nevd
            levalz(i)=evalz(i)*cfact
        30 END DO
    end if

!---------------------------------------------------------------------
!     calculate the projection of each eigenvector on starting vector
!---------------------------------------------------------------------

    if (icomb == 1) then
        do 50 i=1,nevo
            levaly(i)=levalx(i)
            do 55 j=1,ndimo
                evecy(j,i)=evecx(j,i)
            55 END DO
        50 END DO
    end if

    do 60 i=1,nevo
        sum1=czero
        sum2=czero
        do 65 j=1,ndimo
            sum1=sum1+evecx(j,i)*stvo(j)
            sum2=sum2+evecy(j,i)*stvo(j)
        65 END DO
        cwsp1(i)=sum1
        cwsp2(i)=sum2
    60 END DO
            

!---------------------------------------------------------------------
!     calculate constant part in 2-D calculation for the given
!     experimental type
!---------------------------------------------------------------------

!                 *** 2 Pulse exp. : COSY or SECSY ***

    if ((iexp == 1) .OR. (iexp == 2)) then
    
    !         call zgemul(evecy,mxdim,'t',evecx,mxdim,'n',
    !     #        plsmix,mxegv,nevo,ndimo,nevo)
    
        call zgemm('t','n',nevo,nevo,ndimo,cone,evecy,mxdim,evecx, &
        mxdim,czero,plsmix,mxegv)
    
    else if (iexp >= 3) then
    
        if ((iexp == 3) .OR. (iexp == 4)) then
        
        !                 *** 3 Pulse exp. : ELDOR, Echo-ELDOR ***
        
        !            do 3654 i=1,ndimd+1
        ! 3654           write(luttyo,*)(evecz(i,j),j=1,nevd+1)
            do 90 i=1,nevd
                levalz(i)=cdexp(-levalz(i)*ltfix)
            90 END DO


        !                    t
        !          spit out O  O for evecx, form not in matrix style: split out
        !          row 1, then row 2 and so on

        !            write(luttyo,*) "EVECZ^T*EVECZ old from spcalc"
        !            norm=0.0d0
        !            do 3657 i=1,nevo
        !                do 3658 j=1,nevo
        !                    bubba=czero
        !                    do 3659 k=1,ndimo
        !                        bubba=bubba+evecz(k,i)*evecz(k,j);
        ! 3659               continue
        !                    write(luttyo,*) i,j,bubba
        !                    flg=0.0d0
        !                    if(i.eq.j) flg=2.0d0
        !                    norm=norm+abs(bubba-flg)**2
        ! 3658           continue
        ! 3657       continue
        !            write(luttyo,*) "norm square = ",norm


        !            write(luttyo,*) "nevo,nevd=",nevo,nevd


        
        !            call zgemul(evecx,mxdim,'t',evecz,mxdim,'n',
        !     #                  c2wsp1,mxdim,nevo,ndimo,nevd)
        
            call zgemm('t','n',nevo,nevd,ndimo,cone,evecx,mxdim,evecz, &
            mxdim,czero,c2wsp1,mxdim)

        
        !            call zgemul(evecy,mxdim,'t',evecz,mxdim,'n',
        !     #                  c2wsp2,mxegv,nevo,ndimo,nevd)
        
            call zgemm('t','n',nevo,nevd,ndimo,cone,evecy,mxdim,evecz, &
            mxdim,czero,c2wsp2,mxegv)
        

            do 230 i=1,nevo
                do 230 j=1,nevo
                    sum1=czero
                    do 235 k=1,nevd
                        sum1=sum1+c2wsp2(i,k)*levalz(k)*c2wsp1(j,k)
                    235 END DO
                    plsmix(i,j)=sum1*icomb
            230 END DO
        
        else
        
        !                 *** 3 Pulse exp. : Stim.-SECSY ***
        
            do 240 i=1,nevo
                levaly(i)=cdexp(-levaly(i)*tfix)
            240 END DO
        
            do 245 i=1,nevo
                cwsp3(i)=cwsp2(i)*levaly(i)
            245 END DO
        
        !            call zgemul(evecy,mxdim,'t',evecz,mxdim,'n',
        !     #                  c2wsp1,mxdim,nevo,ndimo,nevd)
        
            call zgemm('t','n',nevo,nevd,ndimo,cone,evecy,mxdim,evecz, &
            mxdim,czero,c2wsp1,mxdim)
        !          ****print th evecz matrix to see whether it has ndimd rows;
        !              we will go till ndimd+1 to be really sure

            do 260 i=1,nevd
                sum1=czero
                do 265 j=1,nevo
                    sum1=sum1+c2wsp1(j,i)*cwsp3(j)
                265 END DO
                cwsp2(i)=sum1*icomb
            260 END DO
        
        !            call zgemul(evecz,mxdim,'t',evecx,mxdim,'n',
        !     #                  plsmix,mxegv,nevd,ndimo,nevo)
        
            call zgemm('t','n',nevd,nevo,ndimo,cone,evecz,mxdim,evecx, &
            mxdim,czero,plsmix,mxegv)

        
        end if
    
    end if

!=====================================================================
!     calculate 2-D spectrum in time-domain
!=====================================================================

!---------------------------------------------------------------------
!     step out t1 values
!---------------------------------------------------------------------

    do 300 it1=1,npt1
        t1=linit1+(it1-1)*lstept1
    
        if (iexp /= 5) then
            do 305 i=1,nevo
                c2wsp1(i,it1)=cdexp(-levaly(i)*t1)*cwsp2(i)
            305 END DO
        else
            do 307 i=1,nevd
                c2wsp1(i,it1)=cdexp(-levalz(i)*t1)*cwsp2(i)
            307 END DO
        end if
    300 END DO

    ntmp=nevo
    if (iexp == 5) ntmp=nevd

    if (iexp >= 1) then
    !         call zgemul(plsmix,mxegv,'t',c2wsp1,mxdim,'n',
    !     #               c2wsp2,mxegv,nevo,ntmp,npt1)
        call zgemm('t','n',nevo,npt1,ntmp,cone,plsmix,mxegv,c2wsp1, &
        mxdim,czero,c2wsp2,mxegv)
    else
    !                                      *** iexp = 0: CW 1D-FID ***
    !         call zgemv('t',nevo,npt1,cone,c2wsp1,mxdim,
    !     #              cwsp2,1,czero,plsmix,1)	! returns plsmix
    !						! = cone*c2wsp1'*cwsp2
    
        call zgemv('t',nevo,npt1,cone,c2wsp1,mxdim, &
        cwsp2,1,czero,plsmix,1)	! returns plsmix
    !						! = cone*c2wsp1'*cwsp2
    
        npt2=1
        go to 400
    end if

!---------------------------------------------------------------------
!     step out t2 values
!---------------------------------------------------------------------

!               *** t2 independent of t1 : no SECSY or echo-ELDOR ***

    if ((iexp == 1) .OR. (iexp == 3) .OR. (iexp == 5)) then
    
        do 330 it2=1,npt2
            t2=linit2+(it2-1)*lstept2
            if (iexp == 5) t2=ltfix+t2
        
            do 335 i=1,nevo
                c2wsp1(i,it2)=cwsp1(i)*cdexp(-levalx(i)*t2)
            335 END DO
        330 END DO
    
    !         call zgemul(c2wsp2,mxegv,'t',c2wsp1,mxdim,'n',
    !     #               plsmix,mxegv,npt1,nevo,npt2)
    
        call zgemm('t','n',npt1,npt2,nevo,cone,c2wsp2,mxegv,c2wsp1, &
        mxdim,czero,plsmix,mxegv)
    
    
    !               *** t2 dependent of t1 : SECSY or echo-ELDOR ***
    
    else if ((iexp == 2) .OR. (iexp == 4)) then
    
        do 350 it1=1,npt1
            t1=0.0d0
            do 360 it2=1,npt2
                t2=linit2+(it2-1)*lstept2
                t2=t2+t1
                do 365 i=1,nevo
                    c2wsp1(i,it2)=cwsp1(i)*cdexp(-levalx(i)*t2)
                365 END DO
            360 END DO
        
        !            call zgemv('t',nevo,npt2,cone,c2wsp1,mxdim,
        !     #                 c2wsp2(1,it1),1,czero,plsmix(it1,1),mxegv)
            call zgemv('t',nevo,npt2,cone,c2wsp1,mxdim, &
            c2wsp2(1,it1),1,czero,plsmix(it1,1),mxegv)
        350 END DO
    
    end if

!---------------------------------------------------------------------
!     Add spectrum
!---------------------------------------------------------------------

    if(iexp == 3) then

        call matrxo(ierr)
        do 4210 it1=1,npt1
        !            do 4220 it2=1,npt2 assuming that npt1=npt2, stept1=stept2; need to change
        !             WHAT ABOUT DEAD TIME????????
            om1=twopi*(-0.5d0+real(it1-1)/npt1)/(lstept1*cfact)
            do 3228 i=1,mxdim
                cwsp1(i)=czero
                cwsp2(i)=czero
                cwsp3(i)=czero
                cwsp4(i)=czero
                xl(i)=czero
                w(i)=czero
                do 3230 j=1,mxstep
                    c2wsp1(i,j)=czero
                3230 END DO
                cwsp5(i)=czero
            3228 END DO
        


            ermin=1.0d20
            call cscg(stvo,ndimo,nstep,cgtol,dcmplx(shiftr,shifti+om1+b0), &
        ! should be shifti+om1+b0 for pS=1
            xl,c2wsp1,cwsp1,cwsp2,cwsp3,cwsp4,cwsp5,w,ndone,error,ermin,nmin)
            if (ndone <= 0) then      ! not converged
                ndone=abs(ndone)
                write (luttyo,*) 'Off-diag CG error, om1=',om1
            else
            !                 write(lulog,*) "Yup!"
            end if

        !               write(luttyo,*) 'om1=',om1
        !               call scmvm(xl,y1,ndimo)
        !               write(luttyo,*) "checking whether cscg worked,
        !     #                           ie. Ax=b, in spcalc"
        !               tmpry=czero

            spec_pS1_cnorm(it1)=czero
            do 9989 j=1,ndimo
                spec_pS1_cnorm(it1)=spec_pS1_cnorm(it1)+ &
                spec_pS1(j,it1)*spec_pS1(j,it1)
            9989 END DO
            spec_pS1_cnorm(it1)=sqrt(spec_pS1_cnorm(it1))
            do 9990 j=1,ndimo
                spec_pS1(j,it1)=spec_pS1(j,it1)/spec_pS1_cnorm(it1)
            9990 END DO

            write(luttyo,*) "it1=",it1
            do 7995 j=1,ndimo
                write(luttyo,*) spec_pS1(j,it1)
            7995 END DO


        !               do 2990 ll=1,ndimo
        ! 990               spec_pS1(ll,it1)=xl(ll)

        ! 990              tmpry=tmpry+abs(stvo(ll)-y1(ll)-
        !     #                  dcmplx(shiftr,shifti+om1)*xl(ll))**2


        !               write(luttyo,*) 'check_norm_sq=',tmpry
        ! should be shifti+om1+b0 for off-diag space
        !              om2=twopi*(-0.5d0+real(it2-1)/npt2)/lstept2
        ! 220        continue
        4210 END DO

    ! print the spec_pS1 array, one vector at a time
    !      do 7991 it1=1,npt1
    !         do 7992 j=1,ndimd
    ! 992        write(luttyo,*) spec_pS1(j,it1)
    ! 991  continue
        open (unit=11,file='stvznop.txt',status='unknown', &
        access='sequential')
        do 2424 it1=1,npt1
            write(11,*) (spec_pS1(k,it1),k=1,ndimo)
        2424 END DO
        close(unit=11)
        open (unit=11,file='stvznop_norms.txt',status='unknown', &
        access='sequential')
        write(11,*) (spec_pS1_cnorm(it1),it1=1,npt1)
        close(unit=11)


    !                end if

         
    ! on to the 2d ELDOR calculation
    ! layout: xl, xr are the 2 arrays that have the off-diag om1, om2 cg sols
    ! y2 is the output of CG acting on x2, but in the diag space
    ! Note that we blew xl,xr (coming from spec_pS1) up to
    ! diagonal space by using the appropriate propagator
    !         lstepT=2.0d-3
    ! deciding a and b, see pg 125 of 2017-18 lab notebook;
    ! make sure that theta=Tm*(b-a)/M is less than pi


    ! defining parameters for RESOLVANTS
        NDFT=32768
        MDFT=83!99
        MPOL=6
        b=(twopi*1.0d0/ltfix)/cfact
        a=-b
        delta=(b-a)/MDFT
        write(luttyo,*) "a,b=",a,b

        call matrxd(ierr)
        do 5210 it1=1,1!npt1!should be npt1
        ! DEFINED THE STARTING VECTOR
            do 5220 it2=1,1!npt2!should be npt2
            ! T2 IS THE DOT PRODUCT GUY



            ! activate this chunk if many cores
            !c$omp parallel
            !c$omp& private (id,w,yr,xr,ermin)
            !c$omp& private (i,j,omT,nstep,ndone,error,nmin,it3)
            ! omp make threadprivate wkspcm:(cwsp1,cwsp2,cwsp3,cwsp4,cwsp5,c2wsp1)
            !!$omp threadprivate (/wkspcm/)
            
            !  Have each thread say hello.
            
            !c$OMP PARALLEL DO
            !      id = omp_get_thread_num ( )
            !      if (id.eq.0 .AND. it1.eq.1 .AND. it2.eq.1) then
            !         NTHREADS = OMP_GET_NUM_THREADS()
            !         write(luttyo,*) NTHREADS, " threads"
            !         num=(MDFT+1)/NTHREADS
            !      end if



            ! activate this chunk if only 1 core
                id=0
                NTHREADS=1
                num=(MDFT+1)/NTHREADS



            ! number of tasks for each core
            !  distribute work among cores; but everyone is using the same
            ! variables, so each array has to have an extra dimension telling
                      
                do 5230 it3=1+id*num,(id+1)*num

                    do 9999 i=1,ndimd
                        cwsp1(i)=czero
                        cwsp2(i)=czero
                        cwsp3(i)=czero
                        cwsp4(i)=czero
                        w(i)=czero
                        do 3231 j=1,mxstep
                            c2wsp1(i,j)=czero
                        3231 END DO
                        cwsp5(i)=czero
                        yr(i)=czero
                        xr(i)=spec_pS1(i,it2)
                    9999 END DO


                                   
                    ermin=1.0d20
                    omT=a+real(it3-1)*delta
                    call cscg(xr,ndimd,nstep,cgtol, &
                    dcmplx(shiftr,shifti+omT),yr,c2wsp1,cwsp1,cwsp2,cwsp3, &
                    cwsp4,cwsp5,w,ndone,error,ermin,nmin)
                ! should be shifti+omT for pS=0
                ! SHOULD I SUBTRACT THE SHIFT????????????????????????????????
                    if (ndone <= 0) then      ! not converged
                        ndone=abs(ndone)
                        if (it1 == 1 .AND. it2 == 1) then
                            write (luttyo,*) 'Diag CG error, omT=',omT
                        end if
                    !                   write (luttyo,*) 'ndone=',abs(ndone),', error=',error
                    !                   write (luttyo,*) 'CG min was, at step ',ermin,nmin
                    else
                    !                  write(lulog,*) "Yup!"
                    end if


                !                  do 7993 j=1,ndimd
                ! 993                 write(luttyo,*) y1(j)
                                      
                !                  take inner product between cscg ans and it2 soln
                !                  ip(nomT) is element #nomT of the inner product array
                    ip(it3)=czero
                    do 7993 j=1,ndimd
                        ip(it3)=ip(it3)+yr(j)*spec_pS1(j,it1)
                    7993 END DO


                ! need to define the following
                ! complex*16
                ! double precision omT
                ! integer it3,nomT
                                      
                5230 END DO
            !c$OMP END PARALLEL DO
            !c$omp end parallel
            !              take inverse Fourier transform
            !              endpts, MPOL, en, int nn, cpol(MPOL), spol(MPOL), cerr, serr,cdft,sdft,xpol(MPOL),
            !              fft, dftcor, polint,
            !              corfac,corre,corim
                do 12 j=MDFT+2,NDFT
                    ip(j)=czero
                12 END DO
                do 13 j=1,4
                    endpts(j)=ip(j)
                    endpts(j+4)=ip(MDFT-3+j)
                13 END DO
            !               call realft(data,NDFT,1)
                call fft(ip,NDFT)
            !              data(2)=0.
            !       what is this???? I am omitting. Don't see the point.
            ! it is ifft, not fft; so I redefine en accordingly
                en=real(NDFT)-ltfix*cfact*delta*real(NDFT)/twopi+1.
                nn=min(max(int(en-0.5*MPOL+1.),1),NDFT-MPOL+1)
                do 14 j=1,MPOL
                    cpol(j)=real(ip(nn))
                    spol(j)=aimag(ip(nn))
                    xpol(j)=nn
                    nn=nn+1
                !                 write(luttyo,*) 'xpol(j),cpol(j),spol(j)=',xpol(j),
                !     #                                              cpol(j),spol(j)
                14 END DO
            !               xpolt(1)=0.0d0
            !               xpolt(2)=1.0d0
            !               xpolt(3)=2.0d0
            !               xpolt(4)=3.0d0
            !               cpolt(1)=0.0d0
            !               cpolt(2)=1.0d0
            !               cpolt(3)=4.0d0
            !               cpolt(4)=9.0d0
            !               call polint(xpolt,cpolt,3,1.0d0,cdftre,cerr)
            !               write(luttyo,*) 'check polint: observed,expected=',
            !     #                         cdftre,',2.25',';err=',cerr
                call polint(xpol,cpol,MPOL,en,cdftre,cerr)
                call polint(xpol,spol,MPOL,en,sdftre,serr)
            ! cpol real part of the FFT, spol im part
            !c              write(*,*) "w,delta,a,b"
            !c              write(*,*) w,delta,a,b
            !c              write(*,*) "endpts"
            !              write(*,*) (endpts(j), j=1,8)
                cdft=cmplx(cdftre)
                sdft=cmplx(sdftre)
            ! later cdft and sdft could become complex, hence I used cdftre and sdftre
                call dftcor(ltfix*cfact,delta,a,b,endpts,corre,corim,corfac)
                cdft=cdft*corfac+corre
                sdft=sdft*corfac+corim
                c=delta*cos(ltfix*cfact*a)
                s=delta*sin(ltfix*cfact*a)
                cosint=c*cdft-s*sdft
                sinint=s*cdft+c*sdft
                xspec(it1,it2)=(cosint+cmplx(-aimag(sinint),real(sinint))) &
                *spec_pS1_cnorm(it1)*spec_pS1_cnorm(it2)


            ! end term correction factors
                fac1=czero
                do 79 i=1,ndimd
                    cwsp1(i)=czero
                    cwsp2(i)=czero
                    cwsp3(i)=czero
                    cwsp4(i)=czero
                    w(i)=czero
                    do 80 j=1,mxstep
                        c2wsp1(i,j)=czero
                    80 END DO
                    cwsp5(i)=czero
                    yr(i)=czero
                    xr(i)=spec_pS1(i,it2)
                79 END DO

                ermin=1.0d20
                omT=a
                call cscg(xr,ndimd,nstep,cgtol, &
                dcmplx(shiftr,shifti+omT),yr,c2wsp1,cwsp1,cwsp2,cwsp3, &
                cwsp4,cwsp5,w,ndone,error,ermin,nmin)
            ! should be shifti+omT for pS=0
            ! SHOULD I SUBTRACT THE
            ! SHIFT????????????????????????????????
                if (ndone <= 0) then      ! not converged
                    ndone=abs(ndone)
                    write (luttyo,*) 'Diag CG error fac1, omT=',omT
                    write (luttyo,*) 'ndone=',abs(ndone),', error=',error
                    write (luttyo,*) 'CG min was, at step ',ermin,nmin
                else
                !               write(lulog,*) "Yup!"
                end if
                do 81 i=1,ndimd
                    fac1=fac1+spec_pS1(i,it1)*yr(i)
                81 END DO
                fac1=-fac1*cmplx(cos(a*cfact*ltfix),sin(a*cfact*ltfix))* &
                spec_pS1_cnorm(it1)*spec_pS1_cnorm(it2)* &
                ci/(cfact*ltfix)


                fac2=czero
                do 89 i=1,ndimd
                    cwsp1(i)=czero
                    cwsp2(i)=czero
                    cwsp3(i)=czero
                    cwsp4(i)=czero
                    w(i)=czero
                    do 190 j=1,mxstep
                        c2wsp1(i,j)=czero
                    190 END DO
                    cwsp5(i)=czero
                    yr(i)=czero
                    xr(i)=spec_pS1(i,it2)
                89 END DO

                ermin=1.0d20
                omT=b
                call cscg(xr,ndimd,nstep,cgtol, &
                dcmplx(shiftr,shifti+omT),yr,c2wsp1,cwsp1,cwsp2,cwsp3, &
                cwsp4,cwsp5,w,ndone,error,ermin,nmin)
            ! should be shifti+omT for pS=0
            ! SHOULD I SUBTRACT THE
            ! SHIFT????????????????????????????????
                if (ndone <= 0) then      ! not converged
                    ndone=abs(ndone)
                    write (luttyo,*) 'Diag CG error fac2, omT=',omT
                    write (luttyo,*) 'ndone=',abs(ndone),', error=',error
                    write (luttyo,*) 'CG min was, at step ',ermin,nmin
                else
                !               write(lulog,*) "Yup!"
                end if
                do 91 i=1,ndimd
                    fac2=fac2+spec_pS1(i,it1)*yr(i)
                91 END DO
                fac2=fac2*cmplx(cos(b*cfact*ltfix),sin(b*cfact*ltfix))* &
                spec_pS1_cnorm(it1)*spec_pS1_cnorm(it2)* &
                ci/(cfact*ltfix)

                xspec(it1,it2)=xspec(it1,it2)+fac1+fac2
            !              spectrum=<it1|resolvant|it2> is my formulation
            ! remember to not take FFT again, or better to do an IFFT and then let
            ! the code do its thing to add broadening

            5220 END DO

        5210 END DO


    end if


    400 continue

!      do 410 j=1,npt2
!      do 410 i=1,npt1
! 410     xspec(i,j)=xspec(i,j)+scal*plsmix(i,j)

    if (mod(npt1,2) /= 0 .OR. mod(npt2,2) /= 0 ) then
        write(luttyo,*) "NO NO: don't flip ud/lr before fft"
        stop
    end if

! flip ud and flip lr
! I am assuming that npt1 and npt2 are divisible by 2
! any way, for fft, npt1, npt2 are supposed to be powers of 2
! (64,128,etc.)
    do 444 j=1,npt2
        do 444 i=2,npt1/2
            tmpry=xspec(i,j)
            xspec(i,j)=xspec(npt1+2-i,j)
            xspec(npt1+2-i,j)=tmpry
    444 END DO

    do 445 j=2,npt2/2
        do 445 i=1,npt1
            tmpry=xspec(i,j)
            xspec(i,j)=xspec(i,npt2+2-j)
            xspec(i,npt2+2-j)=tmpry
    445 END DO
        
    do 446 j=1,npt2
        do 447 i=1,npt1
            xl(i)=xspec(i,j)
        447 END DO
        call fft(xl,npt1)
        do 448 i=1,npt1
            xspec(i,j)=xl(i)
        448 END DO
    446 END DO

    do 449 i=1,npt1
        do 450 j=1,npt2
            xl(j)=xspec(i,j)
        450 END DO
        call fft(xl,npt2)
        do 451 j=1,npt2
            xspec(i,j)=xl(j)
        451 END DO
    449 END DO


    return

    1010 format(2x,' CG did not converge in ',i4,' steps : error = ', &
    g14.7)

    end subroutine spcalc
