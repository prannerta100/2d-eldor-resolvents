!  VERSION 1.0  (NLSPMC version)   2/5/99
!*********************************************************************

!                       ===============
!                       SUBROUTINE LUMP
!                       ===============

!     Subroutine LUMP combines the eigenvalues of T-matrix using the
!     relative tolerance reltol.  Lumping is necessary because it is
!     impossible to accurately precict the accuracy of the CMTQLI
!     routine.  LUMP combines T-eigenvalues that have slipped by the
!     tolerance that was used in the T-multiplicity tests.  In parti-
!     cular if for some j,

!          |eval(j)-eval(j-1)| < max(reltol*|eval(j)|,scale2*multol)

!     then these T-eigenvalues are combined.  multol is the tolerance
!     that was used in the T-multiplicity test in COMPEV.

!     If in a set of T-eigenvalues to be combined there is an eigenvalue
!     with lindex=1, then the value of the combined T-eigenvalues is set
!     equal to the value of that eigenvalue.  Note that if a spurious
!     T-eigenvalue is to be 'combined' with a good eigenvalue, then this
!     is done only by increasing the index, lindex, for that eigenvalue
!     numerical values of spurious T-eigenvalues are never combined with
!     those of good T-eigenvalues.

!     arguments :

!         vc(j) = jth distinct T-eigenvalue
!         va(j) = |vc(j)|, in order of increasing magnitude
!         lindex(j) = T-multiplicity of jth distinct T-eigenvalue
!         loop = number of distinct T-eigenvalues
!         value of reltol is 1.d-8.

!     Includes :
!               nlsdim.inc
!               eprprm.inc

!*********************************************************************

    subroutine lump(vc,v1,va,w,reltol,sputol,scale2,lindex, &
    tflag,loop)

    implicit none

    include 'limits.inc'
    include 'simparm.inc'

    complex*16 vc(mxstep),v1(mxstep),w(mxstep),csum,wsum,czero
    complex*16 wt(100)

    double precision ::  va(mxstep),reltol,sputol,scale2
    double precision ::  thold,th1,th2,dgap,one
    parameter (one=1.0D0,czero=(0.0D0,0.0D0))

    integer ::  i,icount,idif,in,indsum,ispur,j,jn,k,loop,nloop
    integer ::  lindex(mxstep),tflag(mxstep)

!-----------------------------------------------------------------------

    th2=scale2*sputol

    do 10 k=1,loop
        tflag(k)=0
    10 END DO

    nloop=0
    j=0

    20 j=j+1
    if (j > loop) go to 130
    if (tflag(j) == 1) go to 20
    nloop=nloop+1
    tflag(j)=1
    v1(1)=vc(j)
    wt(1)=w(j)
    icount=1
    jn=lindex(j)
    th1=reltol*va(j)
    thold=dmax1(th1,th2)
!     thold=reltol*dmax1(one,va(j))

    if (jn /= 0) then
        indsum=jn
        ispur=0
        csum=dfloat(jn)*vc(j)
        wsum=dfloat(jn)*w(j)
    else
        indsum=1
        ispur=1
        csum=czero
        wsum=czero
    end if

    if (j == loop) go to 70
    i=j
    50 i=i+1
    if (i > loop) go to 70
    if (tflag(i) == 1) go to 50
    dgap=va(i)-va(j)
    if (dgap >= thold) go to 70
    dgap=cdabs(vc(i)-vc(j))
    if (dgap >= thold) go to 50

!     lump vc(i) with vc(j)

    icount=icount+1
    tflag(i)=1
    v1(icount)=vc(i)
    wt(icount)=w(i)
    in=lindex(i)

    if (in == 0) then
        ispur=ispur+1
        indsum=indsum+1
    else
        indsum=indsum+in
        csum=csum+dfloat(in)*vc(i)
        wsum=wsum+dfloat(in)*w(i)
    end if
    go to 50

!     compute the 'combined' T-eigenvalue and the resulting
!     T-multiplicity

    70 continue

    if (icount == 1) indsum=jn

    idif=indsum-ispur

    if (icount == 1) then
        vc(nloop)=vc(j)
        va(nloop)=va(j)
        w(nloop)=w(j)
        lindex(nloop)=indsum
    else if (idif /= 0) then
        csum=csum/dfloat(idif)
        wsum=wsum/dfloat(idif)
        vc(nloop)=csum
        va(nloop)=cdabs(csum)
        w(nloop)=wsum
        lindex(nloop)=indsum
    else
        do 90 k=1,icount
            vc(nloop+k-1)=v1(k)
            va(nloop+k-1)=cdabs(v1(k))
            w(nloop+k-1)=wt(k)
        90 END DO
        nloop=nloop+icount-1
    end if

    go to 20

!     index j is finished

!     on return vc contains the distinct T-eigenvalues  va=|vc|
!     lindex contains the corresponding T-multiplicities

    130 continue
    loop=nloop

    return

    end subroutine lump
