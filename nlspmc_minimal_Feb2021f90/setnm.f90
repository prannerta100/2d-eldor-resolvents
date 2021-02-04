! NLSPMC Version 1.0 2/5/99
!----------------------------------------------------------------------
!                    =====================
!                     subroutine SETFIL
!                    =====================

! Given a generic name in <fileid>, this subroutine removes any
! file extension from <fileid> and determines the names for the
! various possible output files associated with the NLSL slow-
! motional fitting calculations by appending a file extension
! in the form '.xxx'. The extensions are

!      prname : <fileid>.PAR   Parameter file (formatted)
!      lgname : <fileid>.LOG   Log file for fitting results
!      trname : <fileid>.TRC   Trace of NLS steps
!      ixname : <fileid>.INDX  Basis index file
!           ** This default ixname is used only when the user
!              forgets to specify basis index filename using
!              "basis <name-in-full>" command.
!  Includes
!     nlsdim.inc
!     nlsnam.inc   Definition of names in common /nlsnam/

!----------------------------------------------------------------------
    subroutine setfil( fileid )
    implicit none

    include 'limits.inc'
    include 'names.inc'

    character(30) :: fileid

    integer :: iend,iroot
    character(1) :: chr
    external iroot

!######################################################################

    iend = iroot( fileid )
    if (iend < 1) then
        fileid = 'noname'
        iend = 6
    endif

    prname=fileid(:iend)//'.par'
    lgname=fileid(:iend)//'.log'
    dbname=fileid(:iend)//'.deb'
    trname=fileid(:iend)//'.trc'
    ixname=fileid(:iend)//'.indx'
    lthfnm=iend+4
    return
    end subroutine setfil


!----------------------------------------------------------------------
!                    =====================
!                     subroutine SETDAT
!                    =====================

! Given a generic name in <dataid>, this subroutine removes any
! file extension from <dataid> and determines the names for the
! various input and output files associated with the NLSL slow-
! motional fitting calculations by appending a file extension
! in the form '.xxx'. The extensions are

!      dtname : <dataid>.SPC   input 2D spectrum
!      nsname : <dataid>.NSP   splined 2D spectrum (used in NLS fit)
!      ftname : <dataid>.FIT   NLS fit to the 2D input spectrum

!  Includes
!     nlsdim.inc
!     nlsnam.inc   Definition of names in common /nlsnam/

!----------------------------------------------------------------------
    subroutine setdat( dataid )
    implicit none

    include 'limits.inc'
    include 'names.inc'

    character(30) :: dataid

    integer :: iend,iroot
    character(1) :: chr
    external iroot

!######################################################################

    iend = iroot( dataid )

    if (iend < 1) then
        dataid = 'noname'
        iend = 6
    endif

    dtname=dataid(:iend)//'.spc'
    nsname=dataid(:iend)//'.nsp'
    ftname=dataid(:iend)//'.fit'
    lthdnm=iend+4
    return
    end subroutine setdat


!----------------------------------------------------------------------
!                    =========================
!                       function IROOT
!                    =========================
!  Strips the dot and trailing file extension from a file name
!  and returns the length of the root name
!----------------------------------------------------------------------
    function iroot( fileid )
    implicit none

    include 'limits.inc'
    integer :: i,iroot
    character fileid*30,chr*1

    i=0
    1 i=i+1
    chr=fileid(i:i)
    if(chr /= '.' .AND. chr /= ' ' .AND. i < 30) goto 1

    fileid=fileid(:i-1)
    iroot=i-1
    return
    end function iroot
