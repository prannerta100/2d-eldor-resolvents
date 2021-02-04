!  NLSPMC VERSION 1.0 2/5/99
!*********************************************************************

!	declarations of	standard I/O parameters

!	Notes:
!		1) The parameters define here are used as

!			luttyi : logical unit for terminal input
!			luttyo : logical unit for terminal output
!                       lulog  : logical unit for log file
!                       lutrc  : logical unit for trace file
!			ludisk : logical unit for other disk files
!			ludiskb : logical unit for unformatted disk files

!                       NOTE: The units that will be used for additional
!                       disk input depend on the MXSPEC parameter defined
!                       in "nlsdim.inc". If MXSPEC is greater than 1,
!                       units <ludisk+1>, <ludisk+2>, ... <ludisk+MXSPEC>
!                       are assumed to be available.

!		2) The appropriate values for these parameters are
!		   highly operating system and compiler	dependent.

!*********************************************************************

    integer :: luttyi,luttyo,ludisk,ludiskb,lulog,lutrc,ludeb
    parameter (ludeb=3,lulog=4,luttyi=5,luttyo=6,lutrc=7,ludisk=9, &
    ludiskb=8)

    integer :: lucmd,luecho,luout
    common /stdio/ lucmd,luecho,luout

    external mclock
    integer :: time1, time2, mclock
    real :: seconds
    common /atimes/seconds,time1,time2
