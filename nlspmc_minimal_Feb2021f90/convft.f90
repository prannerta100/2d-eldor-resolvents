!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       =================
!                       SUBROUTINE CONVFT
!                       =================

!     This routine performs convolution with Gaussian inhomogeneous
!     broadening and subsequent FT on comlex 2D time-domain data to
!     generate frequency domain spectrum.  It is intended for use
!     with nonlinear least-squares applications.

!     On Entry :
!        xspec :  complex 2D spectrum in time domain

!     On Exit  :
!        cspec  :  complex spectrum in freq. domain

!     Includes :
!               nlsdim.inc
!               eprprm.inc
!               wkspcm.inc
!               physcn.inc

!     Uses :
!               fft.f
!               switch.f

!*********************************************************************

    subroutine convft(xspec,cspec,npt1,npt2)

    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'wkspcm.inc'
    include 'physcn.inc'
    include 'stdio.inc'

    double precision :: zero,gtf,cfact,twopi
    parameter (zero=0.0d0)

    integer :: i,j,npt1,npt2,k
    double precision :: gib0,lib0,tdabs,expnt,gconv,t1,t2, &
    linit1,lstept1,linit2,lstept2,frst,tdt

    complex*16 xspec(npt1,npt2),cspec(npt1,npt2),czero
    parameter (czero=(0.0d0,0.0d0))
    complex*16 fftchk(16)
!      real*8 spec(npt1,npt2)

!#####################################################################

    twopi=8.0d0*datan(1.0d0)
    gtf=g0*betae/(hbar*twopi)*1.0d-6
    cfact=twopi*gtf

    linit1=init1*1.0d-3
    lstept1=stept1*1.0d-3
    if (iexp >= 1) then
        linit2=init2*1.0d-3
        lstept2=stept2*1.0d-3
    end if

    gib0=dabs(gib*cfact/twopi)
    lib0=dabs(lib*cfact)

! check FFT
    do 1111 i=1,16
        fftchk(i)=cmplx(1.0d0,0.0d0)
    1111 END DO
          
    call fft(fftchk,16)

    write(luttyo,*) "FFT"
    do 1112 i=1,16
        write(luttyo,*) fftchk(i)
    1112 END DO



!---------------------------------------------------------------------
!     Convolution with Gaussian and/or Lorentzian inhomogeneous broad.
!---------------------------------------------------------------------

    400 if ((gib0 > 1.0d-13) .OR. (lib0 > 1.0d-13)) then
    
        gib0=twopi*twopi*gib0*gib0/2.0d0
        do 410 i=1,npt1
            t1=linit1+(i-1)*lstept1
        !                                     ** FID **
            if (iexp == 0) then
                c2wsp1(i,1)=xspec(i,1)*dexp(-gib0*t1*t1)*dexp(-lib0*t1)
                go to 410
            end if
        
            if ((iexp == 2) .OR. (iexp == 4) .OR. (iexp == 5)) t1=zero
            do 420 j=1,npt2
                t2=linit2+(j-1)*lstept2
                tdabs=dabs(t2+icomb*t1)
                expnt=gib0*tdabs*tdabs
                gconv=zero
                if (expnt <= 5.0d1)  gconv=dexp(-expnt)
                c2wsp1(i,j)=xspec(i,j)*gconv*dexp(-lib0*tdabs)
            420 END DO
        410 END DO
    
    else
        do 430 j=1,npt2
            do 430 i=1,npt1
                c2wsp1(i,j)=xspec(i,j)
        430 END DO
    end if
    if ((iexp == 2) .OR. (iexp == 4)) then
        ipt=0
        do 470 i=1,npt1
            tdt=init1+dble(i-1)*stept1
        !               * skip the region with asymmetric dead time
            if (ipt == 0) then
                if ((tdt >= init2) .OR. ((init2-tdt) < 1.0d-15)) &
                then
                    ipt=i
                else
                    go to 470
                end if
            end if
        !         * get the threshold index when t2 > t1 starts
            do 480 j=1,npt2
                frst=init2+(j-1)*stept2
                if ((frst >= tdt) .OR. ((tdt-frst) < 1.d-4)) &
                go to 440
            480 END DO
        !                    * copy FID and zero-fill
            440 do 450 k=j,npt2
                c2wsp1(i-ipt+1,k-j+1)=c2wsp1(i,k)
            450 END DO
            do 460 k=npt2-j+2,npt2
                c2wsp1(i-ipt+1,k)=czero
            460 END DO
        
        470 END DO
    end if
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!---------------------------------------------------------------------
!     Frequency spectrum by Fourier transform
!---------------------------------------------------------------------

    if (iexp >= 1) then
    !                                 **  along t2
        do 500 i=1,npt1
            do 510 j=1,npt2
                cwsp3(j)=c2wsp1(i,j)
            510 END DO
        ! divide first pt by 2 for baseline
            cwsp3(1)=cwsp3(1)/2.0d0
            call fft(cwsp3,npt2)
            call switch(cwsp3,npt2)
            do 520 j=1,npt2
                c2wsp1(i,j)=cwsp3(j)
            520 END DO
        500 END DO
    !                                 **  along t1
        do 550 j=1,npt2
            do 560 i=1,npt1
                cwsp3(i)=c2wsp1(i,j)
            560 END DO
        ! divide first pt by 2 for baseline
            cwsp3(1)=cwsp3(1)/2.0d0
            call fft(cwsp3,npt1)
            call switch(cwsp3,npt1)
            do 570 i=1,npt1
                cspec(i,j)=cwsp3(i)
            570 END DO
        550 END DO
    
    else
    
        do 580 i=1,npt1
            cwsp3(i)=c2wsp1(i,1)
        580 END DO
    ! divide first pt by 2 for baseline
        cwsp3(1)=cwsp3(1)/2.0d0
        call fft(cwsp3,npt1)
        call switch(cwsp3,npt1)
        do 590 i=1,npt1
        
        ! not a real CW spectrum, but something to compare with the diagonal
        ! slice from the 2D magnitude plot.  (iexp=0 case)
        
            cspec(i,1)=cwsp3(i)
        590 END DO
    
    end if

    return
    end subroutine convft
