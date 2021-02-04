!  NLSPMC VERSION 1.0 2/5/99
!*********************************************************************

!	declarations of parameters for new series and component fits.

!	Notes:
!	  1) The parameters define here are used as

!	    nspectra - specify number of spectra expected (SPECTRA)
!	    ncomps - specify number of components in each spectra (COMP)
!           basinfo - 1=basis set index, 2=have matrix and tridiag,
!			3=have spectra, 4=undefined yet
!           nbasis - current number of basis sets (max=MXSPEC*MXSITE)
!           ndimoss - off diagonal dimension
!           ndimdss - diagonal dimension
!           cursite - current site under consideration in simulation
!                     get from basinfo(1)
! c           setbass - got this basis set

! To use the basis set location information, use basinfo(index,spec,site)
!  to identify which basis (index=1) the given spec and site require.
! With this information, ndimoss and ndimdss can be used in conjunction with
! pidptr and dpidptr, in basis.inc, to locate the specific basis information.

!      logical setbass

          
    integer :: nspectra,ncomps,basinfo,nbasis, &
    ndimoss,ndimdss,cursite

! initflg is set true after program starts and user enters type of fit
!   information

    common /mulcomp/ nspectra,ncomps,basinfo(5,MXSPEC,MXSITE), &
    nbasis,ndimoss(MXSPEC*MXSITE),ndimdss(MXSPEC*MXSITE)

!     #		setbass(MXSPEC*MXSITE)

          
