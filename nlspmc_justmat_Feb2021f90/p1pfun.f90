!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------
!                    =========================
!                        function P1PFUN
!                    =========================
!----------------------------------------------------------------------
    function p1pfun( parm,iflag )

    implicit none

    include 'limits.inc'
    include 'parms.inc'
    include 'datas.inc'
    include 'lmcomm.inc'

    integer :: iflag
    double precision :: p1pfun,parm,enorm
    external enorm

    iflag=1
    call pfun(nptot,1,parm,fvec,fjac,mxpt,iflag)
    if (ihltcmd /= 0) return
    if (iflag >= 0) p1pfun = enorm(nptot,fvec)/dsqrt(dfloat(nptot))
    return
    end function p1pfun
