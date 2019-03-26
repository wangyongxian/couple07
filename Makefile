#
# Unix Makefile for couple. The list of 57 subroutines is
# correct, but it has not been tested. The may be some missing
# tab characters.
#

FC = gfortran
OPT = -c

OBJ = couple.o abortc.o amat.o backsub.o bathy.o bdf.o \
      bessj0.o bessy0.o bmat.o botdeps.o bounds.o \
      cdfun.o cdj0.o cdsrt.o cdtan.o cexib.o cexibd.o \
      cexiw.o cexiwb.o cexiwd.o chdh.o \
      choldc.o cjdj.o \
      clineqs.o cpres.o \
      crscup.o csds.o csymeig.o cupb.o \
      cupbw.o cupw.o cupwb.o cvmult.o cwns.o dreal.o \
      file.o foresub.o egnval.o galrkn.o gama0.o \
      gama00.o genabc.o header.o hk01.o hk02.o \
      indep.o indexx.o input.o mgs.o mmult.o \
      propm.o prtmat.o prttl.o rcf.o rwprof.o \
      tmmult.o zcexib.o zcexiw.o

all: couple

couple: $(OBJ)
	$(FC) -o couple $(OBJ)

clean:
	rm *.o couple

couple.o: couple.f
	$(FC) $(OPT) couple.f
abortc.o: abortc.f
	$(FC) $(OPT) abortc.f
amat.o: amat.f 
	$(FC) $(OPT) amat.f
backsub.o: backsub.f 
	$(FC) $(OPT) backsub.f
bathy.o: bathy.f 
	$(FC) $(OPT) bathy.f
bdf.o: bdf.f 
	$(FC) $(OPT) bdf.f
bessj0.o: bessj0.f 
	$(FC) $(OPT) bessj0.f
bessy0.o: bessy0.f 
	$(FC) $(OPT) bessy0.f
bmat.o: bmat.f 
	$(FC) $(OPT) bmat.f
botdeps.o: botdeps.f
	$(FC) $(OPT) botdeps.f
bounds.o: bounds.f 
	$(FC) $(OPT) bounds.f
cdfun.o: cdfun.f
	$(FC) $(OPT) cdfun.f
cdj0.o: cdj0.f
	$(FC) $(OPT) cdj0.f
cdsrt.o: cdsrt.f
	$(FC) $(OPT) cdsrt.f
cdtan.o: cdtan.f
	$(FC) $(OPT) cdtan.f
cexib.o: cexib.f
	$(FC) $(OPT) cexib.f
cexibd.o: cexibd.f
	$(FC) $(OPT) cexibd.f
cexiw.o: cexiw.f
	$(FC) $(OPT) cexiw.f
cexiwb.o: cexiwb.f
	$(FC) $(OPT) cexiwb.f
cexiwd.o: cexiwd.f
	$(FC) $(OPT) cexiwd.f
chdh.o: chdh.f
	$(FC) $(OPT) chdh.f
choldc.o: choldc.f
	$(FC) $(OPT) choldc.f
cjdj.o: cjdj.f
	$(FC) $(OPT) cjdj.f
clineqs.o: clineqs.f
	$(FC) $(OPT) clineqs.f
cpres.o: cpres.f
	$(FC) $(OPT) cpres.f
crscup.o: crscup.f
	$(FC) $(OPT) crscup.f
csds.o: csds.f
	$(FC) $(OPT) csds.f
csymeig.o: csymeig.f
	$(FC) $(OPT) csymeig.f
cupb.o: cupb.f 
	$(FC) $(OPT) cupb.f
cupbw.o: cupbw.f
	$(FC) $(OPT) cupbw.f
cupw.o: cupw.f
	$(FC) $(OPT) cupw.f
cupwb.o: cupwb.f
	$(FC) $(OPT) cupwb.f
cvmult.o: cvmult.f
	$(FC) $(OPT) cvmult.f
cwns.o: cwns.f
	$(FC) $(OPT) cwns.f
dreal.o: dreal.f
	$(FC) $(OPT) dreal.f
file.o: file.f
	$(FC) $(OPT) file.f
foresub.o: foresub.f
	$(FC) $(OPT) foresub.f
egnval.o: egnval.f
	$(FC) $(OPT) egnval.f
galrkn.o: galrkn.f
	$(FC) $(OPT) galrkn.f
gama0.o: gama0.f
	$(FC) $(OPT) gama0.f
gama00.o: gama00.f
	$(FC) $(OPT) gama00.f
genabc.o: genabc.f
	$(FC) $(OPT) genabc.f
header.o: header.f
	$(FC) $(OPT) header.f
hk01.o: hk01.f
	$(FC) $(OPT) hk01.f
hk02.o: hk02.f
	$(FC) $(OPT) hk02.f
indep.o: indep.f
	$(FC) $(OPT) indep.f
indexx.o: indexx.f
	$(FC) $(OPT) indexx.f
input.o: input.f
	$(FC) $(OPT) input.f
mgs.o: mgs.f
	$(FC) $(OPT) mgs.f
mmult.o: mmult.f
	$(FC) $(OPT) mmult.f
propm.o: propm.f
	$(FC) $(OPT) propm.f
prtmat.o: prtmat.f
	$(FC) $(OPT) prtmat.f
prttl.o: prttl.f
	$(FC) $(OPT) prttl.f
rcf.o: rcf.f
	$(FC) $(OPT) rcf.f
rwprof.o: rwprof.f
	$(FC) $(OPT) rwprof.f
tmmult.o: tmmult.f
	$(FC) $(OPT) tmmult.f
zcexib.o: zcexib.f
	$(FC) $(OPT) zcexib.f
zcexiw.o: zcexiw.f
	$(FC) $(OPT) zcexiw.f
