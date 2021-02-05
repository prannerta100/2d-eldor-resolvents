!----------------------------------------------------------------------
!                    =========================
!                       common block INDEXF
!                    =========================

!  Common block containing a list of the basis set indices
!  for off-diagonal and diagonal subspace.

! pidptr - points to location within pi1 etc of basis set information
! dpidptr - points into dpi1 etc.
! basisptr - indicates which basis is currently mapped into jqe1 etc

!----------------------------------------------------------------------
  module basis  
    integer,dimension(:),allocatable,save:: jqe1,pi1,qi1,l1,jk1,k1,jm1,m1
    integer,dimension(:),allocatable,save:: mjqe1,mpi1,mqi1,ml1,mjk1,mk1,mjm1,mm1,&
                                            pidptr,dpidptr
    integer,save::                          basisptr

!    integer :: jqe1,pi1,qi1,l1,jk1,k1,jm1,m1
!    integer :: mjqe1,mpi1,mqi1,ml1,mjk1,mk1,mjm1,mm1, &
!    pidptr,dpidptr,basisptr
! allow for multiple sites, but not multiple spec basis sets.
!    common /indexo/ jqe1(MXDIM),pi1(MXDIM),qi1(MXDIM), &
!    l1(MXDIM),jk1(MXDIM),k1(MXDIM),jm1(MXDIM), &
!    m1(MXDIM)
!    common /mindexo/ mjqe1(MMXDIM),mpi1(MMXDIM),mqi1(MMXDIM), &
!    ml1(MMXDIM),mjk1(MMXDIM),mk1(MMXDIM),mjm1(MMXDIM), &
!    mm1(MMXDIM),pidptr(MXSPEC*MXSITE),basisptr

    integer,dimension(:),allocatable,save :: djqe1,dpi1,dqi1,dl1,djk1,dk1,djm1,dm1
    integer,dimension(:),allocatable,save :: mdjqe1,mdpi1,mdqi1,mdl1,mdjk1,mdk1,mdjm1,mdm1
!    common /indexd/ djqe1(2*MXDIM),dpi1(2*MXDIM),dqi1(2*MXDIM), &
!    dl1(2*MXDIM),djk1(2*MXDIM),dk1(2*MXDIM),djm1(2*MXDIM), &
!    dm1(2*MXDIM)
!    common/mindexd/mdjqe1(2*MMXDIM),mdpi1(2*MMXDIM),mdqi1(2*MMXDIM), &
!    mdl1(2*MMXDIM),mdjk1(2*MMXDIM),mdk1(2*MMXDIM), &
!    mdjm1(2*MMXDIM),mdm1(2*MMXDIM),dpidptr(MXSPEC*MXSITE)

!----------------------------------------------------------------------
!                    =========================
!                       common block PULSEP
!                    =========================

!  Common block containining the pulse propagator

!  pp(i)  : the actual pulse propagator element for the i-th basis
!           element in the diagonal space.
!  pid(j) : the number of the corresponding basis elements in the
!           diagonal space for the j-th basis element in off-diagonal
!           basis.  It can be either 1 or 2.

!----------------------------------------------------------------------

    integer,dimension(:),allocatable,save :: pid,mpid
    double precision,dimension(:),allocatable,save :: pp,mpp
!    common /pulsep/ pid(MXDIM),pp(2*MXDIM),mpid(MMXDIM),mpp(2*MMXDIM)
    save
  end module basis  
