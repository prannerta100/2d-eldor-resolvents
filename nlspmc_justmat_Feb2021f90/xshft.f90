!----------------------------------------------------------------------
!                      ==================
!                    VERSION 1.0 NLSPMC 2/5/99
!                        function XSHFT
!                      ==================

!   Given two arrays containing spectral data (presumably experimental)
!   and a calculated spectrum, both with evenly spaced x-values with
!   equal delta-x's, this function returns the amount by which the spectum
!   array must be shifted for maximum overlap with the data array, in
!   units of delta-x. This is accomplished by finding the maximum of
!   the correlation function of the two arrays using fast Fourier
!   transform methods as described by Press et al. (Numerical Recipes).

!   Inputs:
!     data     Data array
!     spctr    Array that is to be shifted to match data array
!     npt      Number of points in data and spctr
!              NOTE: For the FFT, this must be a power of two.
!              This routine does not check for this condition.
!     ans      Contains correlation function (should be dimensioned
!              as complex*16 or twice the size of data array)
!     fft      Work array dimensioned as complex*16(npt)

!     Shift is calculated by finding the maximum in the array containing
!     the correlation function and then doing a parabolic interpolation
!     with the two neighbor points if possible.

!     Function returns shift in units of delta-x (which need never
!     be specified). Note that a POSITIVE shift means that the data
!     lag the calculated spectrum; i.e. the calculation appears shifted
!     to the LEFT of the data.

!     Uses:
!       twofft   (from Numerical Recipes)
!       four1    (from Numerical Recipes)
!       realft   (from Numerical Recipes)

!----------------------------------------------------------------------
    subroutine xshft(data,spctr,npt,ans,fft,shift,ovlap)
    implicit none
    integer :: npt
    double precision :: data(npt),spctr(npt),ans(2*npt),fft(2*npt)

    integer :: i,ixmx,ii,ir,no2
    double precision :: a,b,c,ansr,amax,halfn,ovlap,tmp,xfrac,shift


!  Following code is from subroutine CORREL in Numerical Recipes, Press et al.

!  ----
    call twofft(data,spctr,fft,ans,npt)
    no2=npt/2
    halfn=float(no2)
    do 11 i=1,no2+1
        ii=i+i
        ir=ii-1
        ansr=ans(ir)
        ans(ir)=(fft(ir)*ansr+fft(ii)*ans(ii))/halfn
        ans(ii)=(fft(ii)*ansr-fft(ir)*ans(ii))/halfn
    11 END DO
    ans(2)=ans(npt+1)
    call realft(ans,no2,-1)
! ---

!  Note that correlation function is stored from ans(1) to ans(no2)
!  for zero and positive correlations, and from ans(npt) down to
!  ans(no2+1) for negative correlations. Rearrange the ans array
!  and find the maximum value. (Do not allow negative overlaps!)

    amax=0.0d0
! set a default value for ixmx in case routine below fails....
    ixmx=no2

    do 12 i=1,no2
        tmp=ans(i)
        if (tmp > amax) then
            ixmx=i+no2
            amax=tmp
        end if
        ans(i)=ans(no2+i)
        if (ans(i) > amax) then
            ixmx=i
            amax=ans(i)
        end if
        ans(no2+i)=tmp
    12 END DO
    shift=float(ixmx-1)-halfn

!--------------------------------------------------------------------------
!  Find ordinate and abcissa for extremum of a parabola through ixmx and
!  2 neighbor pts.
!--------------------------------------------------------------------------
    if (ixmx > 1 .AND. ixmx < npt) then
        a = ans(ixmx)
        b = 0.5d0*(ans(ixmx+1)-ans(ixmx-1))
        c = 0.5d0*(ans(ixmx+1)+ans(ixmx-1))-a
        xfrac= -b/(2.0d0*c)
        shift=shift+xfrac
        ovlap = a + b*xfrac + c*xfrac*xfrac
    else
        ovlap=ans(ixmx)
    end if

    return
    end subroutine xshft


!----------------------------------------------------------------------
!                    ======================
!                      subroutine TWOFFT
!                    ======================
!     From Numerical Recipes by Press et al.
!----------------------------------------------------------------------
    subroutine twofft(data1,data2,fft1,fft2,n)
    implicit none
    integer :: j,n,n2
    double precision :: data1(n),data2(n)
    complex*16 fft1(n),fft2(n),h1,h2,c1,c2

    double precision :: zero
    parameter(zero=0.0d0)

!......................................................................

    c1=dcmplx(0.5d0,zero)
    c2=dcmplx(zero,-0.5d0)
    do 11 j=1,n
        fft1(j)=dcmplx(data1(j),data2(j))
    11 END DO
    call four1(fft1,n,1)
    fft2(1)=dcmplx(dimag(fft1(1)),zero)
    fft1(1)=dcmplx(dreal(fft1(1)),zero)
    n2=n+2
    do 12 j=2,n/2+1
        h1=c1*(fft1(j)+dconjg(fft1(n2-j)))
        h2=c2*(fft1(j)-dconjg(fft1(n2-j)))
        fft1(j)=h1
        fft1(n2-j)=dconjg(h1)
        fft2(j)=h2
        fft2(n2-j)=dconjg(h2)
    12 END DO
    return
    end subroutine twofft

!----------------------------------------------------------------------
!                    ======================
!                      subroutine REALFT
!                    ======================
!     From Numerical Recipes by Press et al.
!----------------------------------------------------------------------
    subroutine realft(data,n,isign)
    implicit none
    integer :: i,i1,i2,i3,i4,isign,n,n2p3
    double precision :: c1,c2,h1r,h1i,h2r,h2i,wrs,wis,wr,wi,wpr,wpi, &
    wtemp,theta
    double precision :: data(2*n)

    theta=6.28318530717959d0/2.0d0/dble(n)
    c1=0.5d0
    if (isign == 1) then
        c2=-0.5d0
        call four1(data,n,+1)
    else
        c2=0.5d0
        theta=-theta
    endif
    wpr=-2.0d0*dsin(0.5d0*theta)**2
    wpi=dsin(theta)
    wr=1.d0+wpr
    wi=wpi
    n2p3=2*n+3
    do 11 i=2,n/2+1
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
    11 END DO
    if (isign == 1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
    else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n,-1)
    endif
    return
    end subroutine realft

!----------------------------------------------------------------------
!                    ======================
!                      subroutine FOUR1
!                    ======================
!     From Numerical Recipes by Press et al.
!----------------------------------------------------------------------
    subroutine four1(data,nn,isign)
    implicit none
    integer :: i,istep,j,m,mmax,n,nn,isign
    double precision :: wr,wi,wpr,wpi,wtemp,tempi,tempr,theta
    double precision :: data(2*nn)

    n=2*nn
    j=1
    do 11 i=1,n,2
        if(j > i)then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
        endif
        m=n/2
        1 if ((m >= 2) .AND. (j > m)) then
            j=j-m
            m=m/2
            go to 1
        endif
        j=j+m
    11 END DO
    mmax=2
    2 if (n > mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
            do 12 i=m,n,istep
                j=i+mmax
                tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
                tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
                data(j)=data(i)-tempr
                data(j+1)=data(i+1)-tempi
                data(i)=data(i)+tempr
                data(i+1)=data(i+1)+tempi
            12 END DO
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
        13 END DO
        mmax=istep
        go to 2
    endif
    return
    end subroutine four1
