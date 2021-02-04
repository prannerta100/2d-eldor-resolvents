!  VERSION 1.0  (NLSPMC version)  2/5/99
!----------------------------------------------------------------------
!                    =========================
!                       common block EPRPRM
!                    =========================

!	Variable and common block definition for NLSPMC fitting of
!  experimental spectra


! IMPORTANT NOTES:

!   The order in which the parameters appear in /eprprm/ is critical to
!   the proper functioning of the NLSPMC program. If this common block is
!   to be changed, the following rules must be observed:

!   (1) The gib,lib,hwid,mwid parameters must appear before all other
!       parameters in the common block.  When these parameters are varied,
!       the spectrum may be recalculated from the existing complex time
!       domain data (xspectr) matrix without repeating the eigenvalue
!       and spectral calculation.  Thus, the Jacobian should be calculated
!       with respect to this parameter first.  This is most readily
!       accomplished by ordering the fitting parameters according to
!       their indices in the fepr array.

!   (2) To permit alias names to be used for the g,A, and diffusion
!       tensors, the parameters gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz should
!       appear contiguously in that order in the common block.

!   (3) The parameters igflg,iaflg,irflg should appear contiguously
!       in that order in the common block. (See function tcheck in file
!       tensym)

!       *** MORE CONTIGUOUS PARAMETERS FOR NLSPMC ***
!          2) c20, c22, c40, c42, c44
!          3) shiftr, shifti
!          4) fieldi, fieldf
!          5) ald, bed, gad
!          6) alm, bem, gam

!   (4) The order in which the parameters appear in /eprprm/ must
!       be consistent with the names defined in the block data section
!       of file "ipfind.f"

!   (5) The residence times tl,tkxy,tkzz have been replaced with the
!   	parameters pml, pmxy, pmzz. This is to facilitate fitting using
!	non-Brownian models where it is desirable to vary the residence
!	time and the diffusion rate constant together (i.e. keeping the
!       product, which is related to the rms jump angle, constant).

!   (6) The non-Brownian model flags have been changed to ml, mzz, and
!	mxy for perpendicular, parallel, and internal motions, respectively.
!       The flags may be set to 0, 1, and 2 for Brownian, free, and
!	jump diffusion, respectively.

! NOTE: Files limits.inc must be included BEFORE this file
! in the source code

! Access to the parameteris is given in prmeqv.inc by the line:
! equivalence (fepr(1),gib),(iepr(1),in2)

! For multiple sites, add weight term:  sitewgt.

!***********************************************************************

    double precision :: gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz, &
    pml,pmxy,pmzz,djf,djfprp,oss,psi,ald,bed,gad, &
    alm,bem,gam,c20,c22,c40,c42,c44,shft,gib,lib, &
    hwid,mwid,b0,gamman,fieldi,fieldf,cgtol, &
    shiftr,shifti,ptol,init1,stept1,init2,stept2, &
    tfix,datwg,t2edi,t2ndi,t2efi,t1edi,t1ndi, &
    a0,g0,zeen,pl,pkxy,pkzz,wtol, &
    faa,fgm,fam,fgd,fad,cxp,cpot,xlk,ordpar, &
    sitewgt

    logical :: diaflg,setbas

    integer	in2,ipdf,ist,itype,nort,igflg,iaflg,irflg, &
    ml,mxy,mzz,lemx,lomx,kmx,mmx,ipnmx,jkmn,jmmn, &
    nstep,nfield,iexp,icomb,sptst,ipt,itd,itm,ipsi0, &
    lband,kband,ldelta,kdelta,lptmx,kptmx,ndimo,ndimd, &
    nevo,nevd,neltot,nelv,nelre,nelim,nev1,nev0,ngood, &
    ndiag
!     #          ndiag,specnum,sitenum

    common /eprprm/ gib,lib,hwid,mwid,gxx,gyy,gzz,axx,ayy,azz, &
    dx,dy,dz,pml,pmxy,pmzz,djf,djfprp,oss,psi, &
    ald,bed,gad,alm,bem,gam,c20,c22,c40,c42,c44, &
    t2edi,t2ndi,t2efi,t1edi,t1ndi, &
    sitewgt, &
    b0,gamman,shft, &
    fieldi,fieldf,cgtol,shiftr,shifti,ptol, &
    init1,stept1,init2,stept2,tfix,datwg, &
    a0,g0,zeen,pl,pkxy,pkzz,wtol,faa(5),fgm(5),fam(2,5), &
    fgd(2,5),fad(2,5),cxp(6),cpot(5,5),xlk(5,5),ordpar(5), &
    in2,ipdf,ist,itype,nort,igflg,iaflg,irflg,ml,mxy,mzz, &
    lemx,lomx,kmx,mmx,ipnmx,jkmn,jmmn,nstep,nfield, &
    iexp,icomb,sptst,ipt,itd,itm,ipsi0,lband,kband,ldelta, &
    kdelta,lptmx,kptmx,ndimo,ndimd,nevo,nevd,neltot,nelv, &
    nelre,nelim,nev1(mxort,MXSITE),nev0(mxort,MXSITE),ngood, &
    ndiag,diaflg,setbas
!     #         ndiag,specnum,sitenum,diaflg,setbas

! *** The following constants identify the position of several
!     important parameters within the fepr (and fparm) arrays.
!     Also in the iepr array starting with in2 in the above common.
!     THESE CONSTANTS *MUST* BE REDEFINED IF THE PARAMETER ORDER
!     IN SIMPARM IS CHANGED!!!

    integer :: IGIB,ILIB,IHWID,IMWID,IGXX,IGZZ,IAXX,IAZZ,IDX,IDZ, &
    IPML,IPMXY,IPMZZ,IDJF,IDJFPRP,IOSS,IPSI,IALD,IBED, &
    IGAD,IALM,IGAM,IC20,IC44,IT2EDI,IT2NDI,IT2EFI, &
    IT1EDI,IT1NDI, &
    ISWGT, &
    IB0,IGAMMAN,ISHFT,IFIELDI,ICGTOL, &
    ISHIFTR,IPTOL,IINIT1,ISTEPT1,IINIT2,ISTEPT2,ITFIX,IDWGT, &
    IIN2,IIPDF,IIST,IITYPE,INORT,IIGFLG,IML,ILEMX,INSTEP, &
    INFIELD,IIEXP,IICOMB,ISPTST,INDIMO,INDIMD

    parameter(IGIB=1,ILIB=2,IHWID=3,IMWID=4,IGXX=5,IGZZ=7,IAXX=8, &
    IAZZ=10,IDX=11,IDZ=13,IPML=14,IPMXY=15,IPMZZ=16, &
    IDJF=17,IDJFPRP=18,IOSS=19,IPSI=20,IALD=21,IBED=22, &
    IGAD=23,IALM=24,IGAM=26,IC20=27,IC44=31,IT2EDI=32, &
    IT2NDI=33,IT2EFI=34,IT1EDI=35,IT1NDI=36, &
    ISWGT=37, &
    IB0=38, &
    IGAMMAN=39,ISHFT=40,IFIELDI=41,ICGTOL=43,ISHIFTR=44, &
    IPTOL=46,IINIT1=47,ISTEPT1=48,IINIT2=49,ISTEPT2=50, &
    ITFIX=51,IDWGT=52, &
    IIN2=1,IIPDF=2,IIST=3,IITYPE=4, &
    INORT=5,IIGFLG=6,IML=9,ILEMX=12,INSTEP=19,INFIELD=20, &
    IIEXP=21,IICOMB=22,ISPTST=23,INDIMO=24,INDIMD=25)
