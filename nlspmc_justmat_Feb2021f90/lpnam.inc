!----------------------------------------------------------------------
!                    =========================
!                       common block LPNAM
!                    =========================

! NOTE: the file nlsdim.inc must be included before this file
! specid tells what needs to be recalculated when a parameter is changed.
! 1=matrix, 3=Evalues,Evectors, 4=spectra only, 5=nothing
!----------------------------------------------------------------------
    character(6) :: parnam,alias1,alias2,iprnam
    character(10) :: symstr,symbol
    integer :: symval,specid
    common /lpnaa/ symval(NSYMBL),specid(NFPRM),parnam(NFPRM), &
    iprnam(NIPRM),alias1(NALIAS),alias2(NALIAS),symstr(NSYMTR), &
    symbol(NSYMBL)
         
