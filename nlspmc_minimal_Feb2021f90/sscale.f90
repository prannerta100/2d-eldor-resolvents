!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                      subroutine SSCALE
!                    =========================

!     For spectra in series, this subroutine scales the simulated
!  spectrum so that it matches the experimental spectrum.  The scale
!  factor may have vary for spectra with idepnd=0.

!     In case more than two spectra should share the same scale
!  factor, the dependency flag for the first spectrum in that group
!  should be 0 and those for the rest of the spectra in that group
!  should be 1.  This situation happens when one want to fit Sc+ &
!  Sc- spectra from single experiment simutaneously or when one tries
!  to estimate We (T1 relaxation constant) from the intensities of
!  several ELDOR spectra obtained with the identical experimental
!  condition (e.g. amplifier gain, # of averages, same dead time,
!  etc.).  For this strategy to work, one needs to put those dependent
!  spectra in the contiguous space and only the first one should have
!  idepnd=0 in data command.

!  Includes:
!     nlsdim.inc
!     expdat.inc

!----------------------------------------------------------------------

    subroutine sscale(spct)

    implicit none

    include 'limits.inc'
    include 'datas.inc'

    double precision :: spct(mxpt)

    integer :: i,isp,ispi,ixs,mpts
    double precision :: asum,bsum

    double precision :: zero
    parameter (zero=0.0d0)

!######################################################################

    isp=1
    ispi=1
    ixs=ixsp(1)
    mpts=ndata(1)
!                   * search next independent spectra (idepnd=0)
    10 isp=isp+1
    if (idepnd(isp) == 1) then	! if 1, share this scale factor
        mpts=mpts+ndata(isp)
        go to 10
    end if
!                   * obtain scale factor from group of data sets
    asum=zero
    bsum=zero
    do 20 i=ixs,ixs+mpts-1
        asum=asum+spct(i)*spct(i)
        bsum=bsum+data(i)*spct(i)
    20 END DO

    do 30 i=ispi,isp-1
        sfac(i)=bsum/asum
    30 END DO
!                              * scale the spectra
    do 40 i=ixs,ixs+mpts-1
        spct(i)=spct(i)*sfac(ispi)
    40 END DO

    if (isp > nspc) return

    ispi=isp
    ixs=ixsp(isp)
    mpts=ndata(isp)
    go to 10

    end subroutine sscale
