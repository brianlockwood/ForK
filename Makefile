#	Kriging MAKEFILE

# Fortran Compiler
FC = gfortran
# Archive Command (needed to create .a file)
ARC = ar

# Flags for Fortran Compiler
FFLAGS = -ffree-line-length-none -cpp -DHAVELAPACK -L/usr/lib64/ -llapack -lblas -static

# Flags for Archive command
ARCFLAG = cr

#  File List (do not touch)
VARMOD = argument.f90 

MODS = covars.f90 opt.f90 choleskymod.f90

SUB = covarmatrix.f90 likelihood.f90 likelihood_mle.f90 magnitude.f90 covarmatrix_grad.f90 

WRAPL1=funcmod.f90

WRAPL2 = funcwrapper.f90 patternsearch.f90 simplexsearch.f90 hyperparameters_all.f90 hyperparameters_mle.f90

MAIN= buildkriging.f90 buildkrigingGEK.f90 krigingfuncpredict.f90 krigingfuncpredictGEK.f90 krigingfuncvariance.f90 krigingfuncvarianceGEK.f90 kriginggradpredict.f90 kriginggradpredictGEK.f90 kriginggradvariance.f90 kriginggradvarianceGEK.f90 krigingextremefuncpredict.f90 krigingmaxvariancepredict.f90 krigingmaxeipredict.f90 krigingextremefuncpredictGEK.f90 krigingmaxvariancepredictGEK.f90 krigingmaxeipredictGEK.f90 krigingfunccovar.f90 krigingfuncsample.f90 krigingeipredict.f90 krigingeipredictGEK.f90 

FILES = $(VARMOD) $(MODS) $(SUB) $(WRAPL1) $(WRAPL2) $(MAIN)

OBJS =  $(VARMOD:.f90=.o) $(MODS:.f90=.o) $(SUB:.f90=.o) $(WRAPL1:.f90=.o) $(WRAPL2:.f90=.o) $(MAIN:.f90=.o) $(PROGRAM:.f90=.o) 
MODULES = $(VARMOD:.f90=.mod) $(MODS:.f90=.mod) $(SUB:.f90=.mod) $(WRAPL1.f90=.mod)

LIB = kriginglib.a

all: $(LIB)
	chmod g+rw *.o
clean: 
	rm *.o
	rm *.mod
	rm *.a

$(LIB): $(OBJS)
	$(ARC) $(ARCFLAG) $(LIB) *.o 

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 
