!  VERSION 1.1  (NLSPMC version)   2/5/99
!*********************************************************************

!                       ==================
!                       SUBROUTINE : EVCGF
!                       ==================

!     Subroutine version of EPRCGF program. (intended for NLSPMC)

!     It calculates the eigenvalues and eigenvectors for the general
!     complex matrix using Conjugate gradient algorithm.

!     On Entry :

!        Parameters are passed through the common block /simparm/.
!        Pulse propagator is passed through the common block /basis/.
!        Off-diagonal basis vectors are passed through the common
!               block /basis/.

!     On Exit  :

!        nev   : number of important eigenvalues
!        eval  : important eigenvalues
!        evec  : important eigenvectors
!        ierr  : error flag

!     Includes :
!               nlsdim.inc
!               stdio.inc
!               eprprm.inc
!               parcom.inc
!               eprmat.inc
!               indexf.inc
!               wkspcm.inc

!     Uses :
!               matrxo.f
!               matrxd.f
!               cscg.f
!               csval.f
!               csvec.f

!*********************************************************************

    subroutine evcgf(psval,stvx,nev,eval,evec,ierr)

    implicit none

    include 'limits.inc'
    include 'stdio.inc'
    include 'simparm.inc'
    include 'parms.inc'
    include 'eprmat.inc'
    include 'basis.inc'
    include 'wkspcm.inc'

    integer :: i,j,k,psval,ndim,nev,ierr,mxit,ival,ndone,nmin,ll
    double precision :: zr,zi,zero,bshft,wtolcg,error,errorv,ermin, &
    flg,norm
    complex*16 czero,cone,eval,evec,evtmp,stvx,y,w,cnorm,bubba,y1
    parameter (zero=0.D0,czero=(0.D0,0.D0),cone=(1.D0,0.D0))
    dimension eval(mxstep),evec(mxdim,mxegv),stvx(mxdim),y(mxdim), &
    w(mxdim),y1(mxdim)

!#####################################################################

!---------------------------------------------------------------------
!     calculate the matrix elements
!---------------------------------------------------------------------


    ierr=0
    if (psval /= 0) then
        ndim=ndimo
    !         write(luttyo,*) 'pS=1,nelim, nelre=',nelim,nelre
        call matrxo(ierr)
        write(*,*) '!nelreo ',nelre
        write(*,*) '!nelimo ',nelim
        write(*,*) '!ndimo ',ndim
    !         if (idebug.ne.0) then
    !                                ** store starting vector **
        write(*,*) "Hallo"
        open (unit=ludiskb,file='stvec.stvx',status='unknown', &
    !     #           access='sequential',form='unformatted')
        access='sequential')
        write (ludiskb,*) (stvx(i),i=1,ndim)
        close (unit=ludiskb)
    !                                ** store matrix elements **
        open (unit=ludiskb,file='matrx.mtxx',status='unknown', &
    !     #           access='sequential',form='unformatted')
        access='sequential')
        write (ludiskb,*) (zdiag(1,i),zdiag(2,i),i=1,ndim)
        write (ludiskb,*) (jzmat(i),i=1,ndim+1)
        write (ludiskb,*) (izmat(i),zmat(i),i=1,nelim)
        write (ludiskb,*) (kzmat(i),i=1,ndim+1)
        write (ludiskb,*) (izmat(mxel-i+1),zmat(mxel-i+1),i=1,nelre)
        close (unit=ludiskb)
    !         end if
    else
        ndim=ndimd
    !         write(luttyo,*) 'pS=0,nelim, nelre=',nelim,nelre
        call matrxd(ierr)
        write(*,*) '!nelred ',nelre
        write(*,*) '!nelimd ',nelim
        write(*,*) '!ndimd ',ndim

        write(*,*) "Hallo"
    !         if (idebug.ne.0) then
    !                                ** store starting vector **
    !            open (unit=ludiskb,file='stvec.stvz',status='unknown',
    !     #           access='sequential',form='unformatted')
    !     #           access='sequential')
    !            write (ludiskb,*) (stvx(i),i=1,ndim)
    !            close (unit=ludiskb)
    !                                ** store matrix elements **
        open (unit=ludiskb,file='matrx.mtxz',status='unknown', &
    !     #           access='sequential',form='unformatted')
        access='sequential')
        write (ludiskb,*) (zdiag(1,i),zdiag(2,i),i=1,ndim)
        write (ludiskb,*) (jzmat(i),i=1,ndim+1)
        write (ludiskb,*) (izmat(i),zmat(i),i=1,nelim)
        write (ludiskb,*) (kzmat(i),i=1,ndim+1)
        write (ludiskb,*) (izmat(mxel-i+1),zmat(mxel-i+1),i=1,nelre)
        close (unit=ludiskb)
    !         end if
    end if

    if (ierr /= 0) then
    !                  ** make ierr negative for critical error **
        ierr=-abs(ierr)
        write(luttyo,*)'Error in generating matrix elements, ierr',ierr
        return
    end if
    end subroutine evcgf
