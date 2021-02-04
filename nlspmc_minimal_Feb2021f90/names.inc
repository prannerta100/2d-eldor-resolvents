!----------------------------------------------------------------------
!                   =========================
!                      common block NLSNAM
!                   =========================

! Common block for all the filenames associated with a given
! fitting parameter file ID, or datafile ID.

!      prname : <fileid>.PAR   Parameter file (formatted)
!      lgname : <fileid>.LOG   Log file for fitting results
!      trname : <fileid>.TRC   Trace of NLS steps
!      dbname : <fileid>.DEB   Debugging aids
!      dtname : <dataid>.SPC   Input Data file
!      nsname : <dataid>.NSP   2D splined Data file (used in NLSP fit)
!      ftname : <dataid>.FIT   NLSP FIT
!      ixname : <index fileid> Specified in basid command

! NOTE: Files nlsdim.inc must be included BEFORE this file
! in the source code
!----------------------------------------------------------------------

    character*(NAMELG) prname,lgname,dtname,ftname,nsname,trname, &
    dbname,ixname
    integer :: lthfnm,lthdnm

    common /nlsnam/ lthfnm,lthdnm,prname,lgname,trname,dtname, &
    ftname,nsname,ixname,dbname

    integer :: nfiles
    character*(NAMELG) files
    common /filcom/ nfiles,files(MXSPEC)

