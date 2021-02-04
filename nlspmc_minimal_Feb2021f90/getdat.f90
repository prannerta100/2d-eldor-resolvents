!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                     =========================
!                         function GETDAT
!                     =========================

! Inputs 2D experimental EPR spectra for nonlinear least squares
! fitting programs.

!      function getdat( filen,infmt,npt1,npt2,lumsg,iexp,icomb,
!                       stept1,stept2,spndat,cwdat,
!                       x1a,x2a,tmp1,tmp2,ydat,y2drv )

! Returns 0 for a successful read.

! Subroutine arguments:

!   filen       Name of data file
!   infmt       Specifies the input format of the datafile
!               0: ASCII format for stap or conp (DISSPLA 3D PLOT)
!               1: binary format (IBM in eagle, IEEE elsewhere)
!   npt1,npt2   Number of data points in 2D spectrum
!   lumsg       Logical unit for output of informational messages.
!               Set to 0 if no output is desired
!   iexp        Experimental type
!   icomb       Combination type
!   stept1      Step size in t1 specified in series option
!   stept2      Step size in t2 specified in series option
!   x1a,x2a,tmp1,tmp2,ydat,y2drv  Work spaces for spline

! Output:
!   smax        Maximum of the input data
!   spndat      2D splined spectrum
!   cwdat       cw-equivalent spectrum extracted from 2D spectrum

!----------------------------------------------------------------------

    function getdat( filen,infmt,npt1,npt2,lumsg,iexp,icomb, &
    spt1,spt2,smax,spndat,cwdat,x1a,x2a,tmp1, &
    tmp2,ydat,y2drv )
    implicit none

    include 'limits.inc'
    include 'stdio.inc'

    character(30) :: filen
    integer :: getdat,infmt,npt1,npt2,lumsg,iexp,icomb,i,j,k, &
    npt1i,npt2i,nend,ibuf
    double precision :: spndat(npt1,npt2),cwdat(npt2),spt1,spt2, &
    inif1,res1,inif2,res2,inif1i,res1i,inif2i,res2i,f1,f2, &
    smin,smax

    double precision :: x1a(mxegv),x2a(mxegv),tmp1(mxegv),tmp2(mxegv), &
    tmp3(mxegv),ydat(mxegv,mxegv),y2drv(mxegv,mxegv)
    integer :: imx,jmx
    double precision :: mxdat

!------ for binary input format

!      real wksp1(npt1,npt2)

!######################################################################

!----------------------------------------------------------------------
!     Input datafile according to specified format:
!----------------------------------------------------------------------
    if (lumsg /= 0) write(lumsg,1005) filen

!      if (infmt .eq. 1) then

!   .....format #1: binary single precision format : NOT ACTIVE OPTION

!         open (unit=ludisk,file=filen,status='old',
!     #         access='sequential',form='unformatted',err=98)
!         read (ludisk,err=98) npt2i,npt1i
!         read (ludisk,err=98) inif2i,res2i,inif1i,res1i
!         read (ludisk,err=98) wksp1
!         close (unit=ludisk)

!     ** convert single precision into double precision

!         do 10 i=1,npt1i
!         do 10 j=1,npt2i
!            ydat(i,j)=dble(wksp1(i,j))
! 10      continue
!      else

!  ......format #0 (default): standard ASCII format

    open (unit=ludisk,file=filen,status='old', &
    access='sequential',form='formatted',err=98)
    read (ludisk,*,err=98)
    read (ludisk,1110,err=98) npt2i
    read (ludisk,1110,err=98) npt1i
    read (ludisk,1120,err=98) inif2i,res2i,inif1i,res1i,smin,smax
    read (ludisk,1130,err=98) ((ydat(i,j),j=1,npt2i),i=1,npt1i)
    close (unit=ludisk)

!      end if

    if ((npt1i > mxegv) .OR. (npt2i > mxegv)) then
        getdat=-2
        return
    else
        go to 100
    end if

!     ** Abort processing if there was an input error
    98 close (unit=ludisk)
    getdat=-1
    return

!----------------------------------------------------------------------
!     Perform 2D natural bicubic spline interpolation to get desired
!     resolution & number of points
!     (Ref. Numerical Recipes : spline,splint,splie2,splin2 routines)
!----------------------------------------------------------------------

    100 inif1=-5.0d2*dble(npt1-1)/spt1/dble(npt1)
    res1=-2.0d0*inif1/dble(npt1-1)	! step size in f1
    inif1=-5.0d2/spt1
    inif2=-5.0d2*dble(npt2-1)/spt2/dble(npt2)
    res2=-2.0d0*inif2/dble(npt2-1)
    inif2=-5.0d2/spt2

    if ( ((inif1-inif1i) < -1.0d0) .OR. ((inif2-inif2i) < -1.0d0) ) &
    then
        getdat=-3
        return
    end if

    if ((npt1i /= npt1) .OR. (npt2i /= npt2) .OR. &
    (abs(inif1-inif1i) > 1.0d-3) .OR. &
    (abs(inif2-inif2i) > 1.0d-3) .OR. &
    (abs(res1-res1i) > 1.0d-3) .OR. &
    (abs(res2-res2i) > 1.0d-3) ) then	! do spline...
    
        do 110 j=1,npt1i
            x1a(j)=inif1i+dble(j-1)*res1i
        110 END DO
        do 112 j=1,npt2i
            x2a(j)=inif2i+dble(j-1)*res2i
        112 END DO
    
    !                         ** construct 2nd derivative table
        do 120 j=1,npt1i
            do 122 k=1,npt2i
                tmp1(k)=ydat(j,k)
            122 END DO
            call ncspln(x2a,tmp1,npt2i,tmp2)
            do 124 k=1,npt2i
                y2drv(j,k)=tmp2(k)
            124 END DO
        120 END DO
    
    !                         ** loop over f2 values
        do 130 i=1,npt2
            f2=inif2+dble(i-1)*res2
        
            do 132 j=1,npt1i
                do 134 k=1,npt2i
                    tmp1(k)=ydat(j,k)
                    tmp2(k)=y2drv(j,k)
                134 END DO
                call splnay(x2a,tmp1,tmp2,npt2i,f2,1.0d0,tmp3(j),1)
            132 END DO
        
            call ncspln(x1a,tmp3,npt1i,tmp2)
            call splnay(x1a,tmp3,tmp2,npt1i,inif1,res1, &
            spndat(1,i),npt1)
        130 END DO
    
    else
    !                 ** copy the spectrum
        do 200 j=1,npt1
            do 200 k=1,npt2
                spndat(j,k)=ydat(j,k)
        200 END DO
    end if

!----------------------------------------------------------------------
!     Extract cw-equivalent spectrum from 2D spetrum
!----------------------------------------------------------------------

    if ((iexp == 1) .OR. (iexp == 3)) then
    
    ! *** FID experiment ***
    
    !     The following restriction is only for obtaining quick working
    !     version of NLSPMC program.  One should modify those codes keeping
    !     in mind that correct information should be passed for the
    !     calculation of cw-spectrum.  (init2, spt2, npt2 are used in
    !     current version.  See calling sequence for xshft in pfunnew.f for
    !     detail.)
    !         1) inif1=inif2
    !         2) res1=res2
    !         3) npt1=npt2
    
        if ((npt1 /= npt2) .OR. (abs(res1-res2) > 1.0d-3) .OR. &
        (abs(inif1-inif2) > 1.0d-3)) then
            getdat=-4
            return
        end if
    
        if (icomb == 1) then
        !                              * positive diagonal
            do 210 i=1,npt2
                cwdat(i)=spndat(i,i)
            210 END DO
        else
        !                              * negative diagonal
        !            mxdat=0.0D0
        !            do 211 i=1,npt1
        !              do 211 j=1,npt2
        !                if(spndat(j,i).gt.mxdat) then
        !                  mxdat=spndat(j,i)
        !                  jmx=j
        !                  imx=i
        !                end if
        ! 11         continue
        
        ! modified to pick up expected max at (65,65) for 128x128 data file.
        
            cwdat(1)=0.0d0
            do 212 i=2,npt2
                cwdat(i)=spndat(npt2-i+2,i)
            212 END DO
        !        stop
        end if
    else if ((iexp == 2) .OR. (iexp == 4) .OR. (iexp == 5)) then
    
    ! *** ECHO experiment ***
    
        j=0
        220 j=j+1
        f1=inif1+(j-1)*res1
        if (f1 < 0.0d0) go to 220
        if (-(f1-res1) < f1) j=j-1
    
        do 222 i=1,npt2
            cwdat(i)=spndat(j,i)
        222 END DO
    end if

    getdat=0
    return

!### format statements ############################################

    1005 format(/5x,'Opening file ',a,/)
    1110 format(i10)
    1120 format(6e10.4)
    1130 format(8e10.4)

    end function getdat


!----------------------------------------------------------------------
!                       =======================
!                         subroutine NCSPLN
!                       =======================
! Perform a natural cubic spline fit. This is a modification of
! subroutine SPLINE from Numerical Recipes. Given the tabulated function
! of N points in arrays X and Y, return an array Y2 of length N which
! contains the second derivative of the interpolating function at the
! corresponding X values. The major difference from the original SPLINE
! routine is that the second derivatives at the both boundaries of
! the function are assumed to be zero (hence "natural cubic spline").
!----------------------------------------------------------------------
    subroutine ncspln(x,y,n,y2)
    implicit none
    integer :: nmax
    parameter(nmax=512)
    double precision :: zero,one,two,six
    parameter(zero=0.0d0,one=1.0d0,two=2.0D0,six=6.0d0)

    integer :: i,k,n
    double precision :: x(n),y(n),y2(n),u(nmax),sig,p

    y2(1)=zero
    u(1) =zero
    do 18 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1) + two
        y2(i)=(sig-one)/p
        u(i)=(six*( (y(i+1)-y(i))/(x(i+1)-x(i)) &
        -(y(i)-y(i-1))/(x(i)-x(i-1))  ) &
        /(x(i+1)-x(i-1)) &
        - sig*u(i-1))/p
    18 END DO
    y2(n)=zero
    do 19 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    19 END DO
    return
    end subroutine ncspln


!----------------------------------------------------------------------
!                         ====================
!                          subroutine SPLNAY
!                         ====================

! Modification of Numerical Recipes spline interpolation routine (SPLINT)
! which interpolates an entire array at once. Given the original N
! data points in XA,YA and the second derivative YA2 (calculated by SPLINE
! or NCSPLN) and the set of NS X values starting at X0 and spaced by DX,
! this routine returns the interpolated function values in Y.
! The original binary search in subroutine SPLINT has been replaced,
! since the interpolations start near the lower bound and are repeated
! for successively higher, but closely spaced x values.
!----------------------------------------------------------------------
    subroutine splnay(xa,ya,y2a,n,x0,dx,y,ns)
    implicit none
    integer :: i,klo,khi,n,ns
    double precision :: xa(n),ya(n),y2a(n),y(ns),x,x0,dx,h,a,b

    if (x0 < xa(1) .OR. x0+dx*(ns-1) > xa(n)) &
    pause 'SPLNAY: Bad interpolated X range.'
    x=x0-dx
    klo=1
    khi=2
    do 3 i=1, ns
        x=x+dx
        2 if (x > xa(khi) .AND. khi < n) then
            klo=khi
            khi=khi+1
            goto 2
        endif
        h=xa(khi)-xa(klo)
        if (h == 0.0D0) pause 'SPLNAY: Bad XA input.'
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y(i)=a*ya(klo) + b*ya(khi) + &
        ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.0D0
    3 END DO
    return
    end subroutine splnay
