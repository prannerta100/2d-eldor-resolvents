!  VERSION 1.0  (NLSPMC version)   2/5/99
!----------------------------------------------------------------------

!                         ===================
!                          subroutine MATRXD
!                         ===================

! *** MATRIX STORAGE SCHEME FOR CG ALGORITHM ***

!    This routine calculates the matrix elements diagonal subspace
! (pS=0, qS=+/- 1).  Nonsecular terms in the Hamiltonian is ignored
! due to the high field approximation.  This is done by not allowing
! the coupling between the spaces with different pS indices.  The only
! source of the coupling between diagonal and off-diagonal subspaces is
! pulse propagator.

!    qS-symmetrization, equivalent to M-symmetrization in off-diagonal
! space is performed for the basis set.  [Change the sign of qS also
! when the signs of pI and M change].  And M-symmetry is not used since
! it would not make the matrix block-diagonal.  In the absence of
! nuclear Zeeman, the block with jqS=1 is completely decoupled from the
! block with jqS=-1.  The pulse propagator does not couple different
! jqS and and jM.  In the presence of nuclear Zeeman, the block with
! jqS=-1 is included in the calculation.

!    This routine is based upon EPRLL routines ver.1.3., whose features are
!      1. improved matrix storage algorithm (store only half of the matrix)
!      2. support nonaxial diffusion tensor R (in case of Brownian motion)
!      3. uses pruning scheme by eprbl and tnll.  eprbl and tnll generates
!         an <fileid>.ind file which keeps only the basis elements with
!         significant projection on the starting vector.

!    The phenomenological relaxation processes are also included as
!    T1e, T2n, T1n.

!      T1e  : electronic T1 relaxation for the allowed ESR transitions
!      T2n  : orientationally invariant nuclear T2 relaxation
!             This defines the homogeneous line-width of NMR transition.
!      T1n  : nuclear T1 relaxation for the allowed NMR transitions
!             This parameter is implicitly calculated since the pseudo-
!             secular term is included in the calculation of the evolution
!             of the diagonal subspace (pS=0)
!   *Note that the sign of these relaxation constants are reversed
!    with respect to the time domain formalism of the relaxation prcesses
!    since we are calculating the frequency spectrum directly.

!      Original program by G. Moro with modifications by D. Budil
!      Further modification by Sanghyuk Lee to calculate the slow-motional
!              2D-FT spectra.

!       flags
!       -----

!             fld   : True if |lr-lc| <=2.
!                     Used as flag for bandwidth of Liouville operator
!                     matrix elements. The limit on the L-bandwidth is due
!                     to the restriction of spin operator tensors to
!                     ranks 0 and 2.

!             fldmd : True if |lr-lc| <=2 or mr=mc, false otherwise.
!                     Used as flag for possible non-zero matrix
!                     elements. The L-bandwidth restriction is
!                     explained above and the diffusion operator
!                     matrix elements can be found only when mr=mc.

!             flk0  : True if lr=lc and kr=kc.
!                     Used as a flag for determing the existence of
!                     non-zero Liouville (A-tensor) and diffusion
!                     (Heisenberg spin exchange) matrix elements.

!	      fkm0  : True if kr=kc and mr=mc

!       Includes :
!               nlsdim.inc
!               eprprm.inc
!               indexf.inc
!               eprmat.inc
!               maxl.inc
!               rndoff.inc
!               physcn.inc

!       Uses:
!               cd2km.f
!               anxlk.f
!               ipar.f
!               w3j.f

!----------------------------------------------------------------------
    subroutine matrxd(ierr)

    use eprmat
    use basis
    implicit none
    include 'limits.inc'
    include 'simparm.inc'
!    include 'basis.inc'
!    include 'eprmat.inc'
    include 'maxl.inc'
    include 'rndoff.inc'
    include 'physcn.inc'

    logical :: fnz,fld,flk0,diftrm,sumtrm,fldmd
    logical :: newlr,newjkr,newkr,newjqer,newmr,newpnr, &
    newlc,newjkc,newkc,newjqec,newmc,newpnc
    integer :: ioldlr,ioldjkr,ioldkr,ioldjqer,ioldmr,ioldpnr, &
    ioldlc,ioldjkc,ioldkc,ioldjqec,ioldmc,ioldpnc

    double precision :: c,ctemp,c1,c2,c3,c4,cdfd,cdfd1, &
    cdff,cdffu,cga0,cga2,cgam, &
    cliou,cliou0,clioua,clioug, &
    cnorm,cnk,cnl,cpmkr,cplkc,cplmc,ct, &
    d1,dsq2,dsq12,dsq18,dsq23, &
    fii,fki,ga,gg,ossc, &
    ra,ra1,ra2,rg,rg1,rg2,rpl,rmi,rz,rj,rp, &
    sa1,sg1,wliou,z,t1e,t2n,t1n

    double precision :: zero,one
    parameter (zero=0.0d0,one=1.0d0)

    double precision :: d2km(2,5,5)
    integer :: i,j,l,li,ierr,iper, &
    jqer,jqec,jqed,ipnr,ipnc,ipnd,ipndab,ipns, &
    iqnr,iqnc,iqnd,iqndab, &
    lr,lc,ld,ldabs,lrprod,jkr,jkc,jkd, &
    kr,kc,kd,ks,kdabs,ksabs,kdi,ksi,kip, &
    mr,mc,md,ms,mdabs, &
    nrow,ncol,nel,nelr,neli

    integer :: ipar
    double precision :: w3j
    external ipar,w3j

!######################################################################

!      if(abs(zeen).gt.rndoff) fnz=.true.

!----------------------------------------------------------------------
!     define constants
!----------------------------------------------------------------------

    dsq2=dsqrt(2.0D0)
    dsq12=one/dsq2
    dsq18=one/dsqrt(8.0D0)
    dsq23=dsqrt(2.0D0/3.0D0)

!----------------------------------------------------------------------
!     Scale diffusion constants
!----------------------------------------------------------------------

    ct=g0*betae/hbar
    rpl=0.5D0*(dx+dy)/ct
    rmi=0.25D0*(dx-dy)/ct
    rz=dz/ct
    rj=djf/ct
    rp=djfprp/ct
    ossc=oss/ct

    t2n=t2ndi/ct
    t1e=t1edi/(2.0D0*ct)
    t1n=t1ndi/(2.0D0*ct)

    fii=0.25D0*dble(in2*(in2+2))

!    Non-Brownian motional parameters

    if(ipdf /= 1) then
        pml=zero
        pmzz=zero
        pmxy=zero
    end if

!   Correction term for anisotropic viscosity

    if((ipdf == 2) .AND. (ist == 0)) then
        rpl=rpl+rp
        rz=rz+rp
    end if

!----------------------------------------------------------------------
!                                 2
!     Get Wigner rotation matrix D  (0,psi,0) for director tilt.
!                                 KM

!     Note: resultant matrix will be real-valued and d2km(2,i,j)=zero

!                              2
!           d2km(1,i+3,j+3) = d  (psi)
!                              ij
!----------------------------------------------------------------------

    call cd2km(d2km,zero,psi,zero)

!----------------------------------------------------------------------
!                               L
!     Calculate the quantities X  used in the potential-dependent
!                               K
!     portion of the diffusion operator
!----------------------------------------------------------------------

    call anxlk(rpl,rmi,rz)

!----------------------------------------------------------------------
!     Initialize counters
!----------------------------------------------------------------------

    nelr=0
    neli=0
    nel=0

    iper=0
    ioldjqer=0
    ioldlr=-1
    ioldjkr=0
    ioldkr=1000
    ioldmr=1000
    ioldpnr=-in2-1

!----------------------------------------------------------------------
!     **** Loop over rows ***
!----------------------------------------------------------------------

    do 200 nrow=1,ndimd
        jqer=djqe1(nrow)
        lr=dl1(nrow)
        jkr=djk1(nrow)
        kr=dk1(nrow)
        mr=dm1(nrow)
        ipnr=dpi1(nrow)
        iqnr=dqi1(nrow)
    
        newjqer=jqer.ne.ioldjqer
        newlr=newjqer.or.(lr.ne.ioldlr)
        newjkr=newlr.or.(jkr.ne.ioldjkr)
        newkr=newjkr.or.(kr.ne.ioldkr)
        newmr=newkr.or.(mr.ne.ioldmr)
        newpnr=newmr.or.(ipnr.ne.ioldpnr)
    
        if(newjqer) ioldjqer=jqer
    
        if(newlr) then
            ioldlr=lr
            lrprod=lr*(lr+1)
        end if
    
        if(newjkr) ioldjkr=jkr
    
        if(newkr) then
            ioldkr=kr
        
        !         ** Calculate isotropic part of diffusion operator
        
            c1=dble(lrprod)
            if(ipdf /= 0) then
                c2=one+pml*c1
                if(dabs(pl) > rndoff) c2=c2**pl
                c1=c1/c2
            end if
            cdfd=rpl*c1
        
            c1=dble(kr*kr)
            c2=rz
            c3=rpl
            if(ipdf /= 0) then
                c4=one+pmzz*c1
                if(dabs(pkzz) > rndoff) c4=c4**pkzz
                c2=c2/c4
                c4=one+pmxy*c1
                if(dabs(pkxy) > rndoff) c4=c4**pkxy
                c3=c3/c4
            end if
            cdfd=cdfd+c1*(c2-c3)
        
        !         ** Discrete jump motion
        
            if(ist /= 0) then
                i=kr/ist
                if((i*ist) /= kr) cdfd=cdfd+rj
            end if
        
            cdfd1=cdfd
        
        end if
    
        if(newmr) then
            ioldmr=mr
            cpmkr=ipar(mr+kr)
        
        !   Anisotropic viscosity term in diffusion operator
        
            if(ipdf == 2) cdfd=cdfd1+dble(mr*mr)*(rj-rp)
        
        end if
    
        if(newpnr) ioldpnr=ipnr
    
    !----------------------------------------------------------------------
    !   **** Loop over columns (note only upper diagonal) ****
    !----------------------------------------------------------------------
    
        ioldjqec=0
        ioldlc=-1
        ioldjkc=0
        ioldkc=1000
        ioldmc=1000
        ioldpnc=-in2-1
    
        jzmat(nrow)=neli+1
        kzmat(nrow)=nelr+1
    
        do 300 ncol=nrow,ndimd
            jqec=djqe1(ncol)
            lc=dl1(ncol)
            jkc=djk1(ncol)
            kc=dk1(ncol)
            mc=dm1(ncol)
            ipnc=dpi1(ncol)
            iqnc=dqi1(ncol)
        
            newjqec=jqec.ne.ioldjqec
            newlc=newjqec.or.(lc.ne.ioldlc)
            newjkc=newlc.or.(jkc.ne.ioldjkc)
            newkc=newjkc.or.(kc.ne.ioldkc)
            newmc=newkc.or.(mc.ne.ioldmc)
            newpnc=newmc.or.(ipnc.ne.ioldpnc)
        
            if(newjqec) then
                ioldjqec=jqec
                jqed=jqer-jqec
            end if
        
            if(newlc) then
                ioldlc=lc
                ld=lr-lc
                ldabs=iabs(ld)
                fld=ldabs.le.2
                cnl=dsqrt((2.D0*lr+1.D0)*(2.D0*lc+1.D0))
            end if
        
            if(newjkc) then
                ioldjkc=jkc
                jkd=jkr-jkc
            end if
        
            if(newkc) then
                ioldkc=kc
                kd=kr-kc
                ks=kr+kc
                kdabs=iabs(kd)
                ksabs=iabs(ks)
                cplkc=ipar(lc+kc)
            
                flk0=(ld.eq.0).and.(kd.eq.0)
            
                if((kr == 0) .AND. (kc == 0)) then
                    cnk=0.5D0
                else if((kr /= 0) .AND. (kc /= 0)) then
                    cnk=one
                else
                    cnk=dsq12
                end if
            
            !---------------------------------------------------------------------
            !       Calculate R(mu,l) as in Meirovitch, Igner, Igner, Moro, & Freed
            !       J. Phys. Chem. 77 (1982) p. 3935, Eqs. A42 & A44
            !---------------------------------------------------------------------
            
                ra=zero
                rg=zero
                if(fld) then
                    ra1=zero
                    rg1=zero
                    if(kdabs <= 2) then
                        if(jkd == 0) then
                            ga=fad(1,kd+3)
                            gg=fgd(1,kd+3)
                        else
                            ga=fad(2,kd+3)*jkr
                            gg=fgd(2,kd+3)*jkr
                        end if
                        z=w3j(lr,2,lc,kr,-kd,-kc)
                        ra1=ga*z
                        rg1=gg*z
                    end if
                
                    ra2=zero

                    rg2=zero
                    if(ksabs <= 2) then
                        if(jkd == 0) then
                            ga=fad(1,ks+3)
                            gg=fgd(1,ks+3)
                        else
                            ga=fad(2,ks+3)*jkr
                            gg=fgd(2,ks+3)*jkr
                        end if
                        z=w3j(lr,2,lc,kr,-ks,kc)
                        ra2=ga*z
                        rg2=gg*z
                    end if
                !                                      --- for g tensor
                    rg=rg1+cplkc*jkc*rg2

                !                                      --- for A tensor
                    if(in2 /= 0) ra=ra1+cplkc*jkc*ra2
                end if
            
            !----------------------------------------------------------------------
            !       Calculate off-diagonal terms of the diffusion operator,
            !       including the potential-dependent portion and the rhombic
            !       component of the isotropic diffusion operator.
            
            !       See Meirovitch, Igner, Igner, Moro,and Freed,
            !       J. Chem. Phys. 77 (1982) pp.3933-3934, and especially
            !       Eqns. A22-24 and A40 for more information.
            !----------------------------------------------------------------------
            
            !           -------- Rhombic part of isotropic diffusion operator
            
                if((ld == 0) .AND. (kr+2 == kc)) then
                    cdff=rmi*dsqrt(dble((lr+kc-1)*(lr+kc)* &
                    (lr-kc+1)*(lr-kc+2)))/cnk
                else
                    cdff=zero
                end if
            
            !           -------- Potential-dependent part of diffusion operator
            
                if((ipt /= 0) .AND. (ldabs <= lband) .AND. (ipar(ks) == 1) &
                 .AND. (jkd == 0) .AND. ((kdabs <= kband) &
                 .OR. (ksabs <= kband))) then
                
                    kdi=1+kdabs/2
                    ksi=1+ksabs/2
                    c1=cpmkr*cnl*cnk
                
                !             Loop over L and K indices in sum over terms in potential
                
                    do 444 l=0,lband,2
                        cdffu=0.0D0
                        li=1+l/2
                        if (ksi <= li .AND. abs(xlk(li,ksi)) >= rndoff) &
                        cdffu=cplkc*jkc*xlk(li,ksi)*w3j(lr,l,lc,kr,-ks,kc)
                        if (kdi <= li .AND. abs(xlk(li,kdi)) >= rndoff) &
                        cdffu=cdffu+xlk(li,kdi)*w3j(lr,l,lc,kr,-kd,-kc)
                    
                        if (abs(cdffu) > rndoff) cdff=cdff+ &
                        w3j(lr,l,lc,mr,0,-mr)*c1*cdffu
                    444 END DO
                
                end if
            end if
        
            if(newmc) then
                ioldmc=mc
                md=mr-mc
                ms=mr+mc
                mdabs=iabs(md)
                cplmc=ipar(lc+mc)
            
            !----------------------------------------------------------------------
            !     set flags and zero out matrix elements
            !     if no non-zero elements are possible
            !----------------------------------------------------------------------
            
                if(fld .OR. (md == 0)) then
                    fldmd=.true.
                else
                    cliou=0.0D0
                    cgam=0.0D0
                    fldmd=.false.
                end if
            
            !----------------------------------------------------------------------
            !     get appropriate 3-J symbols if non-zero Liouville
            !     operator matrix elements are possible
            !----------------------------------------------------------------------
            
                wliou=zero
                if(fld) wliou=w3j(lr,2,lc,mr,-md,-mc)
            
            end if
        
            if(newpnc) then
                ioldpnc=ipnc
                ipnd=ipnr-ipnc
                ipndab=iabs(ipnd)
                ipns=ipnr+ipnc
            
                cnorm=cnl*cnk*cpmkr
            
            !----------------------------------------------------------------------
            !     get director tilt dependence of Liouville
            !     operator matrix elements
            !----------------------------------------------------------------------
            
                d1=zero
                if((ipndab <= 2) .AND. (mdabs <= 2)) &
                d1=d2km(1,ipnd+3,md+3)*wliou
            
            end if
        
            if(fldmd) then
                iqnd=iqnr-iqnc
                iqndab=iabs(iqnd)
            
            !----------------------------------------------------------------------
            !    Liouville operator
            
            !    Note: see Meirovitch, Igner,Igner,Moro, and Freed
            !    J. Chem. Phys. 77 (1982) Appendix B p.3936-3937 for details.
            
            !          Sa and Sg in this program also contain the
            !          appropriate Clebsch-Gordon coefficient, and fki is
            !          Ki in text.
            
            !          The manner in which the isotropic contributions from the
            !          spin Hamiltonian are included in this program is not clear
            !          from the text, and  deserves some comments.
            !          When Lc=Lr, both the (script) l=0 and l=2 ISTO's must
            !          be included, whereas only l=2 terms appear off the diagonal.
            !          The Clebsch-Gordan coefficients for these terms are
            !          c(1,0,1,0;0,0) and c(1,+/-1,1,-/+1;0,0), or -/+ 1/sqrt(3),
            !          respectively. Here, cga0 is set to +/- 1 instead of -/+
            !          1/sqrt(3) because the extra factors are already included in g0
            !          and a0 (note sign change). The isotropic terms are not
            !          normalized by any factor.
            !	   The cnk factor is not needed because the l=0 ISTO components
            !	   have not been transformed by the K-symmetrization.
            !          The factor cnl*cpmkr is canceled by the product of
            !          w3j(lr,mr,0,0,lr,-mr) and w3j(lr,kr,0,0,lr,-kr), which
            !          therefore do not need to be calculated.
            
            !----------------------------------------------------------------------
            
                cliou=zero
            
                if(fld) then
                    clioua=zero
                    clioug=zero
                    cliou0=zero
                
                !        **   Hyperfine interaction   **
                
                    if((jqed == 0) .AND. (in2 /= 0) .AND. (mdabs <= 2) .AND. &
                    (ipndab == iqndab) .AND. (ipndab <= 1)) then
                    
                        cga0=zero
                        if(ipnd /= 0) then
                        !                                         *** pseudosecular term ***
                        !					      (delta pI .ne. 0)
                            kip=iqnr*iqnd+ipnr*ipnd
                            kip=kip*(kip-2)
                            fki=dsqrt(fii-0.25D0*kip)
                            sa1=-iqnd*fki*dsq18
                            cga2=dsq12
                        else
                        !                                         *** secular term ***
                        !					      (delta pS,pI .eq. 0)
                            sa1=ipnr*0.5D0
                            cga0=one
                            cga2=dsq23
                        end if
                        clioua=sa1*d1*ra*cga2
                        if(flk0 .AND. (jkd == 0) .AND. (md == 0)) &
                        cliou0=cliou0+sa1*cga0*a0
                    
                    end if
                
                !        **   Nuclear Zeeman interaction (isotropic)  **
                
                !              if(fnz.and.(jqed.ne.0).and.flk0.and.(jkd.eq.0).and.
                !     #          (md.eq.0).and.(ipnd.eq.0).and.(iqnd.eq.0)) then
                !                cliou0=cliou0+ipnr*zeen
                !              end if
                
                    cliou=cnorm*(clioua+clioug)+cliou0
                
                end if
            
            !----------------------------------------------------------------------
            !           Diffusion operator
            !----------------------------------------------------------------------
            
                cgam=zero
            
            !      **   Heisenberg exchange   **
            
                if((dabs(ossc) > rndoff) .AND. flk0 .AND. (jkd == 0) .AND. &
                (jqed == 0)) then
                    c=zero
                    ctemp=zero
                    if((md == 0) .AND. (ipnd == 0)) then
                        if(iqnd == 0) c=1.0D0
                        if((ipnr == 0) .AND. (lr == 0)) c=c-one/dble(in2+1)
                    end if
                    if((ms == 0) .AND. (ipns == 0)) then
                        if(iqnd == 0) ctemp=1.0D0
                        if((ipnr == 0) .AND. (lr == 0)) &
                        ctemp=ctemp-one/dble(in2+1)
                        ctemp=ctemp*jqec*cplmc
                    end if
                    cgam=ossc*(c+ctemp)*0.5d0
                end if
            
            !      **   Potential-dependent terms   **
            
                if((ipnd == 0) .AND. (iqnd == 0) .AND. (md == 0) .AND. &
                (jkd == 0) .AND. (jqed == 0))     cgam=cgam+cdff
            
            !      **   Isotropic electronic T1 relaxation (We)  **
            
                if((jqed == 0) .AND. (ld == 0) .AND. (jkd == 0) .AND. &
                (kd == 0) .AND. (iqnd == 0)) then
                
                    if((mr == mc) .AND. (ipnr == ipnc)) cgam=cgam+t1e
                    if((mr == -mc) .AND. (ipnr == -ipnc)) &
                    cgam=cgam+jqec*cplmc*t1e
                end if
            
            !      **   Isotropic nuclear T1 relaxation (Wn)  **
            !           (specifically for I=1/2 and I=1 case)
            
                if((jqed == 0) .AND. (ipnd == 0) .AND. (ld == 0) .AND. &
                (jkd == 0) .AND. (kd == 0) .AND. (md == 0)) then
                !                                   * I=1/2 case
                    if(in2 == 1) then
                        if(ipnr == 0) then
                            if(iqnd == 0) then
                                cgam=cgam+t1n
                            else if(iqndab == 2) then
                                cgam=cgam-t1n
                            end if
                        else if((ipnr == 1) .OR. (ipnr == -1)) then
                            cgam=cgam+t1n*7.0d0/6.0d0
                        end if
                    !                                   * I=1 case
                    else if(in2 == 2) then
                        if(ipnr == 0) then
                            if(iqnd == 0) then
                                if((iqnr == in2) .OR. (iqnr == -in2)) then
                                    cgam=cgam+t1n
                                else
                                    cgam=cgam+t1n*2.0d0
                                end if
                            else if(iqndab == 2) then
                                cgam=cgam-t1n
                            end if
                        else if((ipnr == 1) .OR. (ipnr == -1)) then
                            if(iqndab == 0) cgam=cgam+t1n*13.0d0/6.0d0
                            if(iqndab == 2) cgam=cgam-t1n
                        else if((ipnr == 2) .OR. (ipnr == -2)) then
                            cgam=cgam+t1n*11.0d0/3.0d0
                        end if
                    end if
                end if
            
            !      **   T2n relaxation HERE !!!   **
            
            end if
        
        
        !----------------------------------------------------------------------
        !         Store matrix elements for Rutishauser algorithm
        !         ***  stores upper triangular matrix only
        !----------------------------------------------------------------------
        
            if(nrow == ncol) then
                cgam=cgam+cdfd
                zdiag(1,nrow)=cgam
                zdiag(2,nrow)=-cliou
            else
            
                if(abs(cliou) > rndoff) then
                    nel=nel+1
                    if(nel > mxel) then
                        ierr=1
                        return
                    else
                        neli=neli+1
                        zmat(neli)=-cliou
                        izmat(neli)=ncol
                    end if
                end if
            
                if(abs(cgam) > rndoff) then
                    nel=nel+1
                    if(nel > mxel) then
                        ierr=1
                        return
                    else
                        nelr=nelr+1
                        zmat(mxel-nelr+1)=cgam
                        izmat(mxel-nelr+1)=ncol
                    end if
                end if
            end if
        
        !----------------------------------------------------------------------
        !     Increment column counter
        !----------------------------------------------------------------------
        
        !  ****test*****
        !         if (nrow.eq.1) then
        !         if (ncol.eq.1)       write(6,332)
        !         write(6,333) ncol,iqec,lc,jkc,kc,mc,ipnc,iqnc
        ! 332     format(2x,' ncol  iqec   lc    jkc   mc   ipnc  iqnc')
        ! 333     format(2x,8i6)
        !         end if
        
            if(ncol > mxdim) then
                ierr=2
                return
            end if
        !----------------------------------------------------------------------
        !         end loop over columns
        !----------------------------------------------------------------------
        300 END DO
    !----------------------------------------------------------------------
    !       end loop over rows
    !----------------------------------------------------------------------
    200 END DO

    ierr=0
    neltot=nel
    nelre=nelr
    nelim=neli

    jzmat(ndimd+1)=neli+1
    kzmat(ndimd+1)=nelr+1

    return
    end subroutine matrxd
