!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!     Declarations of common block to hold 2D spectra

!     NOTE: the file 'nlsdim.inc' must be included before this file.
!  Add MXSITE for multi components.  Note, multiple spectra are
!  indicated by MXPT=10*size of one spectrum.

!*********************************************************************

    double precision :: spectr,totspec
    complex*16 xspectr,cspectr,ctotspec
    common	/tdspec/  spectr(MXPT,MXSITE), &
    xspectr(MXPT,MXSITE),totspec(MXPT), &
    cspectr(MXPT,MXSITE),ctotspec(MXPT)
