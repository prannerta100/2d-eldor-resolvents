!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                       subroutine WRFIT
!                    =========================

! Writes the NLS fit 2D spectrum into the file.

!----------------------------------------------------------------------
!      subroutine wrfit( spec,npt1,npt2,inif2,res2,inif1,res1,
!     #                  amin,amax,name,lth )
    subroutine wrfit( spec,icnt,amin,amax,name )
    implicit none

    include 'limits.inc'
    include 'stdio.inc'
    include 'datas.inc'

    integer :: i,j,npt1,npt2,icnt
    double precision :: inif2,res2,inif1,res1,amin,amax
!c      double precision spec(snpt1(icnt),snpt2(icnt))
    double precision :: spec(128,128)

    character*(WORDLG) name
    character(30) :: xaxisf,yaxisf
    data xaxisf/'f2 (MHz)'/
    data yaxisf/'f1 (MHz)'/

!######################################################################

! the following is moved here from fitp.f

    npt1=snpt1(icnt)
    npt2=snpt2(icnt)
! calc amax if not done already
    if(amax < 1.0d-10) then
        do 10 j=1,npt2
            do 10 i=1,npt1
                if(spec(i,j) > amax) amax=spec(i,j)
        10 END DO
        amin=0.0d0
    end if

    inif1=-5.0D2/sstept1(icnt)
    inif2=-5.0D2/sstept2(icnt)
    res1=10.0D2/snpt1(icnt)/sstept1(icnt)
    res2=10.0D2/snpt2(icnt)/sstept2(icnt)

    open(unit=ludisk,file=name,status='unknown', &
    access='sequential',form='formatted')
    write (ludisk,3100) 0
    write (ludisk,3110) npt2,1,npt2,xaxisf
    write (ludisk,3110) npt1,1,npt1,yaxisf
    write (ludisk,3120) inif2,res2,inif1,res1,amin,amax
    write (ludisk,3120) ((spec(i,j),j=1,npt2),i=1,npt1)
    close (unit=ludisk)

    if(npt2 == 1) then
        open(unit=ludisk,file='cwspec',status='unknown', &
        access='sequential',form='formatted')
        write (ludisk,3130) (i,spec(i,1), &
        (spec(i+1,1)-spec(i-1,1))/2.0d0,i=2,npt1-1)
        close (unit=ludisk)
    end if
    return

! #### format statements ########################################

    3100 format(i10)
    3110 format(3i10,5x,a30)
    3120 format(8e10.4)
    3130 format(1i10,2x,e10.4,2x,e10.4,/)

    end subroutine wrfit
