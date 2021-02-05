!	VERSION	1.0	2/5/99
!*********************************************************************

!	Declarations of common block to hold matrix for EPRLL

!       These declarations include the additional jzmat and kzmat
!       arrays used with the updated versions of MATRLL and SCMVM
!       subroutines.

!	Notes:
!		1) The dimensions of arrays used declared here are
!		   determined by parameters declared in	the include
!		   file	stddim.inc, therefore the inclusion of the
!		   file	must follow the	inclusion of stddim.inc.
!		2) For more information	on the usage of	these
!		   arrays, see the header for the Lanczos algorithm
!		   subroutine (cslnzs.f) and the sparse	matrix-vector
!		   multiplication routine (scmvm.f).

!	written	by DJS 11-SEP-87

!*********************************************************************

    module eprmat
      integer,dimension(:),allocatable,save:: izmat,jzmat,kzmat
      double precision, dimension(:),allocatable,save:: zmat
!         common /scmat/ zmat(MXEL),zdiag(2,MXDIM),
!     #          izmat(MXEL),jzmat(MXDIM+1),kzmat(MXDIM+1)
      double precision, dimension(:,:),allocatable,save:: zdiag
!         allocate(zdiag(2,MXDIM),izmat(MXEL),jzmat(MXDIM+1),&
!          kzmat(MXDIM+1),zmat(MXEL)) !can't allocate in a module it seems
      save
    end module eprmat
