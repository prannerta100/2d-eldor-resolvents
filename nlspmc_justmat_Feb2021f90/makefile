########################################################################
#
#	     ====================================================
#	     Make descriptor file for NLSPC and related programs
#	     ====================================================
#
#  VERSION FOR RUNNING UNDER windows 2000
#
########################################################################
#
.SUFFIXES: .f90 .inc .o .c

# See: https://computing.llnl.gov/tutorials/openMP/
F77=gfortran -std=legacy
#the one above is original

#F77=gfortran

#F77=ifort
#replaced gfortran with ifort

#FFLAGS = -c -O2 -g
FFLAGS = -c -O2 -g -fopenmp -std=gnu -m64 -mcmodel=medium -Wuninitialized
#FFLAGS = -c -O2 -g -fopenmp -std=gnu -fbacktrace -ffpe-trap=zero,overflow,underflow -fno-automatic -m64 -mcmodel=large
#FFLAGS = -c -O2 -g -std=gnu -fbacktrace -ffpe-trap=zero,overflow,underflow -m64 #-fPIC -mcmodel=large
#LIB = -L /usr/X11/lib64 -lX11 -L/usr/X11R6/lib64  -llapack -lblas
LIB = -L /usr/X11/lib64 -lX11 -L/usr/X11R6/lib64 

CC = cl
LIBS = dfor.lib
#LIBS2  = C:\PROGRA~1\MPIPro\LIB\MPIPro.lib sstatd.lib sstats.lib smathd.lib smaths.lib sf90mp.lib

LIBS2  = C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsl.lib C:\PROGRA~1\Micros~3\DF98\IMSL\LIB\imsls_err.lib

#	
#-----------------------------------------------------------------------
#			Definitions
#-----------------------------------------------------------------------
#
NLSPI = limits.inc names.inc parms.inc stdio.inc miscel.inc 
NLSPI2 = parmequ.inc simparm.inc
CMDI = $(NLSPI) $(NLSPI2) lpnam.inc names.inc basis.inc datas.inc 
CMDI2 = $(CMDI) lmcomm.inc parms.inc rndoff.inc basis.inc
VCH = vchange.o

LETI = limits.inc simparm.inc datas.inc parms.inc lpnam.inc
NLSI = limits.inc simparm.inc datas.inc lmcomm.inc stdio.inc parms.inc
FUNI = limits.inc simparm.inc datas.inc parms.inc 
MATI = limits.inc simparm.inc basis.inc eprmat.inc stvcom.inc
# CMDI = $(NLSI) lmcomm.inc parms.inc lpnam.inc
DATI = limits.inc datas.inc names.inc stdio.inc parms.inc 
CGO  = cscg.o zaypx.o zdotu2.o zaxpy2.o scmvm.o znormu.o 
EVCG = csval.o compev.o lump.o comgap.o isoev.o inverr.o cmtqli.o csvec.o
EPO1 = sim2d.o evcgf.o cd2km.o anxlk.o w3j.o ipar.o fbasis.o
EPO2 = matrxo.o matrxd.o stveco.o stvecd.o ccrints.o fz.o bessel.o plgndr.o
#EPO3 = spcalc.o convft.o fft.o switch.o zgemm.o zgemv.o lsame.o xerbla.o polint.o dftcor.o
EPO3 = convft.o fft.o switch.o zgemm.o zgemv.o lsame.o xerbla.o polint.o dftcor.o
NEWSTF = spectra.o comps.o vchange.o
NLSC = cmds.o datac.o letcmc.o addprm.o rmvprm.o srchc.o tensym.o convtc.o
NLSS = strutl.o lprmpt.o helpc.o ipfind.o
NLSD = setnm.o getdat.o wrfit.o
NLSH = xshft.o sscale.o mnbrak.o brent.o p1pfun.o
NLSF = fitp.o pfunnew.o pcheck.o $(NLSH) $(EPO1) $(EPO2) $(EPO3) $(EVCG) $(CGO)
NLSB = lmnls.o enorm.o dpmpar.o qrfac.o lmpar.o qrsolv.o covar.o mapxxx.o
NLSO = nlspmc.o $(NLSC) $(NLSS) $(NLSD) $(NLSF) $(NLSB) nlsinit.o ordrpr.o

#-----------------------------------------------------------------------
#		Object files
#-----------------------------------------------------------------------

addprm.o	: addprm.f90 $(CMDI)
anxlk.o		: anxlk.f90 limits.inc simparm.inc rndoff.inc
bessel.o	: bessel.f90 rndoff.inc pidef.inc
brent.o		: brent.f90 limits.inc stdio.inc parms.inc
cmds.o		: cmds.f90 $(CMDI2) 
cd2km.o		: cd2km.f90 rndoff.inc pidef.inc
cmtqli.o	: cmtqli.f90 limits.inc rndoff.inc
comgap.o	: comgap.f90 limits.inc
compev.o	: compev.f90 limits.inc simparm.inc
comps.o		: comps.f90 $(CMDI2)
convft.o	: convft.f90 limits.inc simparm.inc wkspcm.inc physcn.inc
convtc.o	: convtc.f90 limits.inc simparm.inc parms.inc stdio.inc lpnam.inc
cscg.o 		: cscg.f90 limits.inc parms.inc stdio.inc rndoff.inc
csval.o		: csval.f90 limits.inc stdio.inc simparm.inc rndoff.inc
csvec.o		: csvec.f90 limits.inc rndoff.inc
datac.o		: datac.f90 $(DATI) lmcomm.inc wkspcm.inc
evcgf.o		: evcgf.f90 limits.inc stdio.inc simparm.inc parms.inc basis.inc wkspcm.inc
fbasis.o	: fbasis.f90 limits.inc simparm.inc basis.inc stdio.inc
fft.o		: fft.f90
fitp.o  	: fitp.f90 $(NLSI) names.inc tdspec.inc lmtxt.inc parms.inc
ftest.o		: ftest.f90
fz.o		: fz.f90 limits.inc simparm.inc
getdat.o	: getdat.f90 limits.inc stdio.inc
helpc.o		: helpc.f90 stdio.inc
inverr.o	: inverr.f90 limits.inc simparm.inc
ipar.o		: ipar.f90
ipfind.o	: ipfind.f90 limits.inc simparm.inc parms.inc lpnam.inc
isoev.o		: isoev.f90 limits.inc simparm.inc
letcmc.o	: letcmc.f90 $(VCH) $(LETI) stdio.inc
lmnls.o		: lmnls.f90 mapxxx.f90 limits.inc parms.inc parms.inc
lump.o		: lump.f90 limits.inc simparm.inc
mapxxx.o	: mapxxx.f90 parms.inc
matrxo.o	: matrxo.f90 $(MATI) maxl.inc rndoff.inc physcn.inc
matrxd.o	: matrxd.f90 $(MATI) maxl.inc rndoff.inc physcn.inc
mnbrak.o	: mnbrak.f90 limits.inc parms.inc
nlsinit.o 	: nlsinit.f90 $(DATI) parms.inc simparm.inc parmequ.inc lmcomm.inc
# nlspmc.o	: nlspmc.f90 $(NLSI) parmequ.inc names.inc parms.inc
nlspmc.o	: nlspmc.f90 $(NLSPI) $(NLSPI2)
ordrpr.o	: ordrpr.f90 rndoff.inc
p1pfun.o	: p1pfun.f90 limits.inc parms.inc datas.inc lmcomm.inc
pcheck.o	: pcheck.f90 $(NLSI) names.inc maxl.inc rndoff.inc
pfunnew.o	: pfunnew.f90 $(VCH) $(FUNI) lpnam.inc tdspec.inc wkspcm.inc stdio.inc parms.inc
rmvprm.o	: rmvprm.f90 limits.inc parms.inc simparm.inc stdio.inc
scmvm.o 	: scmvm.f90 limits.inc rndoff.inc eprmat.inc
scspec.o	: scspec.f90 stdio.inc
setnm.o		: setnm.f90 limits.inc names.inc
#sim2d.o		: sim2d.f90 $(MATI) parmequ.inc parms.inc stdio.inc
spectra.o	: spectra.f90 $(NLSPI) $(NLSPI2) lpnam.inc
srchc.o		: srchc.f90 $(VCH) limits.inc simparm.inc parms.inc lpnam.inc stdio.inc
sscale.o	: sscale.f90 limits.inc datas.inc
strutl.o	: strutl.f90 stdio.inc
stveco.o 	: stveco.f90 limits.inc simparm.inc basis.inc stvcom.inc rndoff.inc
stvecd.o 	: stvecd.f90 limits.inc simparm.inc basis.inc stvcom.inc wkspcm.inc
tensym.o	: tensym.f90 limits.inc simparm.inc parms.inc lpnam.inc stdio.inc
tstjac.o	: tstjac.f90 stdio.inc
vchange.o	: vchange.f90 rndoff.inc
w3j.o		: w3j.f90 maxl.inc
whris.o		: whris.f90
wrfit.o		: wrfit.f90 limits.inc datas.inc stdio.inc
xshft.o		: xshft.f90
zaypx.o 	: zaypx.f90 limits.inc rndoff.inc
zaxpy2.o 	: zaxpy2.f90 limits.inc rndoff.inc
znormu.o	: znormu.f90 limits.inc rndoff.inc
zscsw.o		: zscsw.f90 limits.inc rndoff.inc
zdotu2.o		: zdotu2.f90 limits.inc rndoff.inc

lsame.o         : lsame.f90
zgemm.o         : zgemm.f90
zgemv.o         : zgemv.f90
xerbla.o        : xerbla.f90
#normal.o			 : normal.f9090
polint.o        : polint.f90 stdio.inc
dftcor.o        : dftcor.f90
#-----------------------------------------------------------------------
#		Executable files
#-----------------------------------------------------------------------


#nlspmc	: $(NLSO) $(NEWSTF) ftest.o
#	$(F77) /exe:nlspmc $(NLSO) $(NEWSTF) ftest.o $(LIBS2) $(LIBS)

nlspmc	: $(NLSO) $(NEWSTF) ftest.o
	$(F77) $(LOADFLG) $(NEWSTF) ftest.o -o $@ $(NLSO) $(LIB) $(LIB2) -lc

#-----------------------------------------------------------------------
#			Default actions
#-----------------------------------------------------------------------

.c.o   :
	cc $(CFLAGS) $*.c

.f90.o   :
	$(F77) $(FFLAGS) $*.f90

.inc.o :
	$(F77) $(FFLAGS) $*.f90
