!  VERSION 1.0  (NLSPMC version)   2/5/99
!**********************************************************************

!                       ===================
!                       SUBROUTINE : FBASIS
!                       ===================
!     *** COMBINATION OF FBASO & FBASD ***

!          This routine builds a list of the basis set indices for
!     both diagonal and off-diagonal spaces in common block /indexf/
!     for the truncation parameters set in /eprprm/.
!          If the file, <fileid>.ind exists, it simply reads the basis
!     set.  This basis would be the pruned basis obtained by running
!     program eprbf.  If the file does not exist, it generates full
!     basis set within the given MTS parameters.

!     Includes:
!        nlsdim.inc
!        eprprm.inc
!        basis.inc
!        stdio.inc

!     Uses:
!        ipar.f

!**********************************************************************

    subroutine fbasis(ixname,igen,ierr)

! igen=0 to read from file, else generate full.

    implicit none

    include 'limits.inc'
    include 'simparm.inc'
    include 'basis.inc'
    include 'stdio.inc'
    include 'miscel.inc'

    logical :: fexist
    character*(NAMELG) ixname
    integer :: lu,igen,ierrtmp

    integer :: lr,jkr,kr,krmx,jmr,mr,mrmx,iper,iqer,iqermn,iqermx, &
    ipnr,ipnrmx,ipnrmn,iqnr,iqnrmx,iqnrmn,nrow,i,j,ierr,ilm

    double precision :: one,sqrt2

    integer :: ipar
    external ipar

!######################################################################

    open (unit=9,file='ind_offdiag.indx',status='unknown', &
    access='sequential')
          
    ierr=0
    ierrtmp=0

    if (igen /= 0) go to 100	! generate full
! read from file:
    inquire(file=ixname,exist=fexist)
    if (fexist) then
    
    !----------------------------------------------------------------------
    !     Read basis set information from index file
    !----------------------------------------------------------------------
    
        write(*,*)'opening ixname ',ludiskb, ixname
        open (unit=ludiskb,file=ixname,status='old', &
        access='sequential',form='unformatted')
        write(*,*)'opened it'
        read (ludiskb) ndimoss(nbasis)
    ! fill next area:
        read (ludiskb) (mpi1(i),mqi1(i),ml1(i),mjk1(i),mk1(i),mjm1(i), &
        mm1(i),i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1)
        if (nbasis+1 <= nspectra*ncomps) then   ! ptr to next basis
            pidptr(nbasis+1)=pidptr(nbasis)+ndimoss(nbasis)
        end if
        close (unit=ludiskb)
        go to 200
    end if

    ierrtmp=1		! error, gen full basis (assumes have lr, etc.)

!----------------------------------------------------------------------
!     Check the magnetic parameters set in cmds.f
!----------------------------------------------------------------------

!  setbas should be false here since are asking for a basis.
    setbas=.false.
    ndimo=ndimoss(nbasis)
    100 call pcheck(luttyo,ierr)

    if (ierr < 0) then
        ierr=3
        return
    end if

    nrow=pidptr(nbasis)-1

!----------------------------------------------------------------------
!     *** loop over lr ***
!----------------------------------------------------------------------

    do 110 lr=0,lemx,ldelta
        if((ipar(lr) /= 1) .AND. (lr > lomx)) go to 110
    
    !----------------------------------------------------------------------
    !     *** loop over jkr ***
    !----------------------------------------------------------------------
    
        do 120 jkr=jkmn,1,2
        
        !----------------------------------------------------------------------
        !     *** loop over kr ***
        !----------------------------------------------------------------------
        
            krmx=min(kmx,lr)
            do 130 kr=0,krmx,kdelta
                if((kr == 0) .AND. (ipar(lr) /= jkr)) go to 130
            
            !----------------------------------------------------------------------
            !     *** loop over jmr ***
            !----------------------------------------------------------------------
            
                do 140 jmr=jmmn,1,2
                
                !----------------------------------------------------------------------
                !     *** loop over mr ***
                !----------------------------------------------------------------------
                
                    mrmx=min(mmx,lr)
                    do 150 mr=0,mrmx
                    
                    !----------------------------------------------------------------------
                    !     *** loop over ipnr ***
                    !----------------------------------------------------------------------
                    
                        ipnrmx=min(in2,ipnmx)
                        if (mr == 0) then
                            ipnrmn=0
                        else
                            ipnrmn=-ipnrmx
                        end if
                    
                        do 160 ipnr=ipnrmn,ipnrmx
                            if((mr == 0) .AND. (ipnr == 0) .AND. (ipar(lr) /= jmr)) &
                            go to 160
                            if((ipsi0 == 0) .AND. (ipnr /= mr)) go to 160
                        
                        !----------------------------------------------------------------------
                        !     *** loop over iqnr ***
                        !----------------------------------------------------------------------
                        
                            iqnrmx=in2-iabs(ipnr)
                            iqnrmn=-iqnrmx
                            do 170 iqnr=iqnrmn,iqnrmx,2
                            
                                nrow=nrow+1
                                mjqe1(nrow)=0
                                ml1(nrow)=lr
                                mjk1(nrow)=jkr
                                mk1(nrow)=kr
                                mjm1(nrow)=jmr
                                mm1(nrow)=mr
                                mpi1(nrow)=ipnr
                                mqi1(nrow)=iqnr
                                write(9,*) nrow,ml1(nrow),mjk1(nrow),mk1(nrow), &
                                mjm1(nrow),mm1(nrow),mpi1(nrow),mqi1(nrow)
                            
                            !----------------------------------------------------------------------
                            !     end loop over rows
                            !----------------------------------------------------------------------
                            
                            170 END DO
                        160 END DO
                    150 END DO
                140 END DO
            130 END DO
        120 END DO
    110 END DO

    close(unit=9)
    if (nbasis+1 <= nspectra*ncomps) then   ! ptr to next basis
        pidptr(nbasis+1)=nrow+1		! location of next basis
    end if
    ndimoss(nbasis)=nrow-pidptr(nbasis)+1	! size of this one


    200 continue

!      ndimo=nrow

!**********************************************************************
!     Basis index for diagonal space (FBASD)
!**********************************************************************

    one=1.0D0
    sqrt2=dsqrt(2.0D0)

!----------------------------------------------------------------------
!       convert pruned off-diagonal basis set into diagonal one.
!----------------------------------------------------------------------

    j=dpidptr(nbasis)-1
    do 210 i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1
        j=j+1
        mdl1(j)=ml1(i)
        mdjk1(j)=mjk1(i)
        mdk1(j)=mk1(i)
        mdjm1(j)=0
        mdm1(j)=mm1(i)
        mdjqe1(j)=mjm1(i)
        mdpi1(j)=mpi1(i)
        mdqi1(j)=mqi1(i)
        if ((mpi1(i) == 0) .AND. (mm1(i) == 0)) then
            mpid(i)=1
            mpp(j)=sqrt2
        else
            mpp(j)=one
        
            j=j+1
            mdl1(j)=ml1(i)
            mdjk1(j)=mjk1(i)
            mdk1(j)=mk1(i)
            mdjm1(j)=0
            mdm1(j)=-mm1(i)
            mdjqe1(j)=mjm1(i)
            mdpi1(j)=-mpi1(i)
            mdqi1(j)=mqi1(i)
            mpid(i)=2
            ilm=ipar(ml1(i)+mm1(i))
            if (ilm == mjm1(i)) then
                mpp(j)=one
            else
                mpp(j)=-one
            end if
        end if
    210 END DO
    do 211 i=pidptr(nbasis)+ndimoss(nbasis),j
        mpid(i)=-1	! test for problems
    211 END DO
    if (nbasis+1 <= nspectra*ncomps) then   ! ptr to next basis
        dpidptr(nbasis+1)=j+1		! location of next basis
    end if
    ndimdss(nbasis)=j-dpidptr(nbasis)+1	! length of dp

    if (ierrtmp == 1) ierr=1
    if ((pidptr(nbasis) > mmxdim) .OR. (dpidptr(nbasis) &
     > 2*mmxdim)) ierr=2

    return
    end subroutine fbasis
