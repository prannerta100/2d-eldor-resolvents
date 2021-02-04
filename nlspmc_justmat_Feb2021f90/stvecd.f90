!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       ===================
!                       SUBROUTINE : STVECD
!                       ===================

!     This subroutine calculates the starting vector for the diagonal
!     space, which is used to obtain the weighting factors of diagonal
!     eigenvectors.  It is intended for use with nonlinear lease-squares
!     applications.

!     On Entry :
!        evalx : eigenvalues of the off-diagonal space
!        evevx : eigenvectors of the off-diagonal space

!        Other parameters are passed through common block /eprprm/.

!     On Exit  :
!        stvd  : starting vector for diagonal space returned through
!                the common block /stvcom/

!     10-APR-93 Sanghyuk Lee

!     Includes :
!        nlsdim.inc
!        eprprm.inc
!        indexf.inc
!        stvcom.inc
!        wkspcm.inc

!*********************************************************************

    subroutine stvecd(evalx,evecx)

    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'basis.inc'
    include 'stvcom.inc'
    include 'wkspcm.inc'

    integer :: i,j,k,nt1
    double precision :: t1,cfact

    complex*16 czero
    parameter (czero=(0.0D0,0.0D0))

    complex*16 evalx,evecx,cnorm,znormu
    dimension evalx(mxegv),evecx(mxdim,mxegv)

    external znormu

!#####################################################################

    cfact=1.760844783D1

    do 10 i=1,nevo
        cwsp1(i)=czero
        do 20 j=1,ndimo
            cwsp1(i)=cwsp1(i)+evecx(j,i)*stvo(j)
        20 END DO
    10 END DO

    do 40 nt1=1,5
        t1=(init1+(nt1-1)*stept1)*cfact*1.0D-3
    
        do 50 i=1,nevo
            cwsp2(i)=cdexp(-evalx(i)*t1)
        50 END DO
    
        do 60 i=1,ndimo
            stvd(i,nt1)=czero
            do 70 j=1,nevo
                stvd(i,nt1)=stvd(i,nt1)+evecx(i,j)*cwsp2(j)*cwsp1(j)
            70 END DO
        60 END DO
    
    40 END DO
    do 100 nt1=1,5
        j=ndimd+1
        do 110 i=ndimo,1,-1
            if (pid(i) /= 1 .AND. pid(i) /= 2) then
                write(*,*)'pid error in stvecd ',i,pid(i)
                stop
            end if
            do 110 k=1,pid(i)
                j=j-1
                if(j > MXDIM)then
                    write(*,*)'in stvecd, MXDIM dimension exceeded ',MXDIM,j
                    stop
                end if
                stvd(j,nt1)=pp(j)*stvd(i,nt1)
        110 END DO
    
        if (nt1 == 5) then
            do 120 i=1,ndimd
                stvd(i,1)=0.3d0*stvd(i,1)+0.2d0*stvd(i,2)+ &
                0.2d0*stvd(i,3)+0.2d0*stvd(i,4)+0.1d0*stvd(i,5)
            120 END DO
        end if
    
    100 END DO

    cnorm=znormu(stvd(1,1),ndimd)
    do 140 i=1,ndimd
        stvd(i,1)=stvd(i,1)/cnorm
    140 END DO

    return
    end subroutine stvecd

