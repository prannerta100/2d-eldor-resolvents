!----------------------------------------------------------------------
!                         ====================
!                          common block DATAS
!                         ====================
!     Experimental data used for least-squares fitting and associated
!     parameters.  It also keeps the variable in series option.

!     NOTE: File 'limits.inc' should be included before this file
!           in the program source code.

!  ** User-specified parameters **

!    *Set by data command
!     data           : Experimental data array
!     idepnd         : Dependency flag for scale factor (sfac)
!     sishft         : shift flag

!    *Set by series command
!     siexp          : Experimental type (0=FID, 1=COSY, 2=SECSY, 3=ELDOR
!                                         4=Echo-ELDOR, 5=Stim. SECSY)
!     sicomb         : Spectral combination
!     snpt1,snpt2    : Number of data points in t1 & t2
!     sinit1,sstept1 : Initial and step values in t1 (ns)
!     sinit2,sstept2 : Initial and step values in t2 (ns)
!     stfix          : Fixed time (mixing time in ELDOR)
!     datwgt         : Weight for this data set (1=standard weight)

!  ** parameters calculated during fit **

!     cwdata         : cw-equivalent spectra extracted from 2D data array
!     cwsim          : cw-equivalent spectra simulated
!     cwsimtot       : sum simulated cw-equivalent spectra over sites
!     sfac           : Scale factor for each spectrum in series option
!     sshft          : Field offset (shift spectrum)
!     ndata          : Number of data points in each spectrum
!     ixsp           : Starting index of each spectrum in data array
!     ixspcw         : Starting index for cw-eq. spectrum in cwdata array
!     nspc           : Number of spectra
!     nptot          : Total number of points (all spectra)
!     nptotcw        : Total number of cw data points (all spectra)
!     dataid         : Identification string for each data file
!     sratio         : scaling ratio for each independent experimental
!                      spectrum to have same maximum intensity
!     stfft1	     : specify if fft is to be zero filled in t1
!     stfft2	     : same in t2.

!----------------------------------------------------------------------
    character dataid*30
    integer :: siexp,sicomb,snpt1,snpt2,ndata,ixsp,ixspcw,idepnd, &
    sishft,nspc,nptot,nptotcw
!    #        ,ssfft1,ssfft2

    double precision :: data,cwdata,cwsim,cwsimtot,sinit1,sstept1, &
    sinit2,sstept2,stfix,datwgt,sfac,sshft,sratio
    complex*16 ccwsim,ccwsimtot

    common /expdat/ data(mxpt),cwdata(mxnpt*MXSPEC), &
    cwsim(mxnpt*MXSPEC,MXSITE),cwsimtot(mxnpt*MXSPEC), &
    ccwsim(mxnpt*MXSPEC,MXSITE),ccwsimtot(mxnpt*MXSPEC), &
    sinit1(MXSPEC),sstept1(MXSPEC),sinit2(MXSPEC), &
    sstept2(MXSPEC),stfix(MXSPEC),datwgt(MXSPEC),sfac(MXSPEC), &
    sshft(MXSPEC),sratio(MXSPEC), &
    siexp(MXSPEC),sicomb(MXSPEC),snpt1(MXSPEC),snpt2(MXSPEC), &
    ndata(MXSPEC),ixsp(MXSPEC),ixspcw(MXSPEC),idepnd(MXSPEC+1), &
    sishft(MXSPEC),nspc,nptot,nptotcw,dataid(MXSPEC)

!     #	    ssfft1(MXSPEC),ssfft2(MXSPEC),


