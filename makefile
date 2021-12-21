#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                ifort (version: 14.0.2, 16.0.3, 17.0.1 linux)
#                gfortran (version: 6.3.0 linux and osx)
#                pgf90 (version: 17.10-0, linux)
#                nagfor (version 7.0, osx)
 F90 = gfortran
#F90 = nagfor
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 0
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK = 1
## Some compilers (like PGF90) do not have inverse hyperbolic functions: atanh, asinh, acosh
# NVHYP  = 1 : with intrinsic inverse hyperbolic functions
# NVHYP  = 0 : with external inverse hyperbolic functions (without intrinsic ones)
INVHYP  = 1
#=================================================================================

#=================================================================================
# If ExternalF90 is empty, F90 is unchanged
ifeq  ($(strip $(ExternalF90)),)
else
  F90 = $(ExternalF90)
endif
# If F90 is empty, F90=gfortran
ifeq  ($(strip $(F90)),)
  F90 = gfortran
endif
# If ExternalOPT is empty, OPT is unchanged
ifeq  ($(strip $(ExternalOPT)),)
else
  OPT = $(ExternalOPT)
endif
# If OPT is empty, OPT=0
ifeq  ($(strip $(OPT)),)
  OPT = 0
endif
# If ExternalOMP is empty, OMP is unchanged
ifeq  ($(strip $(ExternalOMP)),)
else
  OMP = $(ExternalOMP)
endif
# If OMP is empty, OMP=1
ifeq  ($(strip $(OMP)),)
  OMP = 1
endif
#=================================================================================

# Operating system, OS? automatic using uname:
OS=$(shell uname)


#=================================================================================
# for c++ preprocessing
#=================================================================================
 CPPpre    = -cpp
#=================================================================================
#=================================================================================

#=================================================================================
# nag compillation (nagfor)
#=================================================================================
ifeq ($(F90),nagfor)
   # for c++ preprocessing
   CPPpre = -fpp
   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      OMPFLAG = -openmp
   endif
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O4  $(OMPFLAG) -Ounroll=4 -v
   else
      #F90FLAGS = -O0 $(OMPFLAG) -g -C=all -mtrace=all
      #  -C=undefined is not compatible with -framework Accelerate (because the routines are not compilled with the -C=undefined option)
      #with -mtrace=all add information on the memmory allocation/deallocation.
      ifeq ($(OMP),0)
        F90FLAGS =  -O0 $(OMPFLAG) -g -gline -C=all -C=undefined
      else
        F90FLAGS = -O0 $(OMPFLAG) -g -C=all
      endif
   endif

   ifeq ($(LAPACK),1)
     F90LIB = -framework Accelerate
   else
     F90LIB =
   endif

   F90_VER = $(shell $(F90) -V 3>&1 1>&2 2>&3 | head -1 )

endif


#=================================================================================
#=================================================================================
# ifort compillation v17 v18 with mkl
#=================================================================================
ifeq ($(F90),ifort)
   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      OMPFLAG = -qopenmp
   endif
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O  $(OMPFLAG) -parallel -g -traceback
   else
      F90FLAGS = -O0 $(OMPFLAG) -check all -g -traceback
   endif

   ifeq ($(LAPACK),1)
     F90LIB = -mkl -lpthread
     #F90LIB = $(MKLROOT)/lib/libmkl_lapack95_ilp64.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_blas95_ilp64.a -lpthread
   else
     F90LIB = -lpthread
   endif

   F90_VER = $(shell $(F90) --version | head -1 )

endif
#=================================================================================
#=================================================================================

#=================================================================================
# pgf90 compillation v12 with mkl
#=================================================================================
ifeq ($(F90),pgf90)

   # With pgf90 invers hyperbolic functions are not present => INVHYP = 0
   INVHYP = 0
   # for c++ preprocessing
   CPPpre    = -Mpreprocess

   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
       OMPFLAG = -mp=allcores
       F90LIB = -lpthread
   endif
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O $(OMPFLAG) -fast -Mallocatable=03
   else
      F90FLAGS = -O0 $(OMPFLAG)      -Mallocatable=03 -Mbounds -Mchkstk -g
   endif

   ifeq ($(LAPACK),1)
     F90LIB += -lblas -llapack
   else
     F90LIB +=
   endif

   F90_VER = $(shell $(F90) --version | head -2 | tail -1 )

endif

#=================================================================================
#=================================================================================
# gfortran (osx and linux)
#=================================================================================
 ifeq ($(F90),gfortran)
   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      OMPFLAG = -fopenmp
   endif
   # OS management
   ifeq ($(LAPACK),1)
     ifeq ($(OS),Darwin)    # OSX
        # OSX libs (included lapack+blas)
        F90LIB = -framework Accelerate
     else                   # Linux
        # linux libs
        F90LIB = -llapack -lblas
        # linux libs with mkl and with openmp
        #F90LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
        # linux libs with mkl and without openmp
        #F90LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
     endif
   else
    F90LIB =
   endif
   #
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O5 -g -fbacktrace $(OMPFLAG) -funroll-loops -ftree-vectorize -falign-loops=16
   else
      #F90FLAGS = -Og -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
      #F90FLAGS = -O0 -fbounds-check -Wuninitialized
      F90FLAGS = -Og $(OMPFLAG) -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
   endif

   F90_VER = $(shell $(F90) --version | head -1 )

endif

QML_ver=$(shell awk '/QML/ {print $$3}' version-QML)
QML_path:=$(shell pwd)
#=================================================================================
#=================================================================================
$(info ***********************************************************************)
$(info ***********OS:           $(OS))
$(info ***********COMPILER:     $(F90))
$(info ***********COMPILER_VER: $(F90_VER))
$(info ***********OPTIMIZATION: $(OPT))
$(info ***********OpenMP:       $(OMPFLAG))
$(info ***********Arpack:       $(ARPACK))
$(info ***********F90FLAGS:     $(F90FLAGS))
$(info ***********F90LIB:       $(F90LIB))
$(info ***********INVHYP:       $(INVHYP))
$(info ***********QML_ver:      $(QML_ver))
$(info ***********QML_path:     $(QML_path))
$(info ***********************************************************************)


CPPSHELL_QML = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
               -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
               -D__COMPILER="'$(F90)'" \
               -D__COMPILER_VER="'$(F90_VER)'" \
               -D__COMPILER_OPT="'$(F90FLAGS)'" \
               -D__COMPILER_LIBS="'$(F90LIB)'" \
               -D__QMLPATH="'$(QML_path)'" \
               -D__QML_VER='"$(QML_ver)"'

CPPSHELL_INVHYP  = -D__INVHYP="$(INVHYP)"
CPPSHELL_MATRIX  = -D__LAPACK="$(LAPACK)"

F90_FLAGS = $(F90) $(F90FLAGS)
LYNK90 = $(F90_FLAGS)

 LIBS := $(PESLIB) $(ARPACKLIB) $(F90LIB)
 LYNKFLAGS = $(LIBS)


#
GRIDEXE   = Grid.x
ADIAEXE   = Adia.x
MODEXE    = ModLib.x
dnSEXE    = dnS.x
dnPolyEXE = dnPoly.x
TESTEXE   = testOOP.x
DriverEXE = Driver.x
ModLib    = libpot.a
QMLib     = libQMLib.a

DIR0      = $(QML_path)
DIROBJ    = $(DIR0)/OBJ
DIRSRC    = $(DIR0)/SRC
DIRLib    = $(DIRSRC)/QMLLib
DIRdnS    = $(DIRSRC)/QMLdnSVM
DIRdnPoly = $(DIRSRC)/QMLdnSVM
DIRdnMat  = $(DIRSRC)/QMLdnSVM
DIRModel  = $(DIRSRC)/QML
DIRAdia   = $(DIRSRC)/AdiaChannels
DIROpt    = $(DIRSRC)/Opt

#
OBJ_QML   = $(DIROBJ)/Empty_m.o \
                 $(DIROBJ)/Sigmoid_m.o $(DIROBJ)/Morse_m.o $(DIROBJ)/Poly1D_m.o \
                 $(DIROBJ)/H2_m.o $(DIROBJ)/Buck_m.o \
                 $(DIROBJ)/Template_m.o $(DIROBJ)/Test_m.o \
                 $(DIROBJ)/H2NSi_m.o $(DIROBJ)/H2SiN_m.o \
                 $(DIROBJ)/HNNHp_m.o $(DIROBJ)/HONO_m.o \
                 $(DIROBJ)/HNO3_m.o $(DIROBJ)/NO3_m.o \
                 $(DIROBJ)/CH5_m.o $(DIROBJ)/PH4_m.o \
                 $(DIROBJ)/HOO_DMBE_m.o \
                 $(DIROBJ)/H3_m.o \
                 $(DIROBJ)/HenonHeiles_m.o $(DIROBJ)/LinearHBond_m.o \
								 $(DIROBJ)/TwoD_MullerBrown_m.o \
                 $(DIROBJ)/Phenol_m.o $(DIROBJ)/TwoD_m.o \
                 $(DIROBJ)/PSB3_m.o $(DIROBJ)/Retinal_JPCB2000_m.o \
                 $(DIROBJ)/OneDSOC_1S1T_m.o $(DIROBJ)/OneDSOC_2S1T_m.o \
                 $(DIROBJ)/Tully_m.o

OBJ_AdiaLib    = $(DIROBJ)/MakeHinact_m.o

OBJ_Model      = $(DIROBJ)/Model_m.o $(DIROBJ)/Basis_m.o $(DIROBJ)/Opt_m.o $(DIROBJ)/IRC_m.o

OBJ_test       = $(DIROBJ)/TEST_OOP.o
OBJ_driver     = $(DIROBJ)/Model_driver.o
OBJ_grid       = $(DIROBJ)/TEST_grid.o
OBJ_adia       = $(DIROBJ)/TEST_Adia.o
OBJ_testmod    = $(DIROBJ)/TEST_model.o
OBJ_testdnS    = $(DIROBJ)/TEST_dnS.o
OBJ_testdnPoly = $(DIROBJ)/TEST_dnPoly.o
OBJ_testdriver = $(DIROBJ)/TEST_driver.o


OBJ_lib        = $(DIROBJ)/FiniteDiff_m.o \
                 $(DIROBJ)/dnMat_m.o $(DIROBJ)/dnPoly_m.o $(DIROBJ)/dnS_m.o \
                 $(DIROBJ)/UtilLib_m.o $(DIROBJ)/diago_m.o $(DIROBJ)/Matrix_m.o \
                 $(DIROBJ)/NumParameters_m.o

OBJ_all        = $(OBJ_lib) $(OBJ_Model) $(OBJ_QML)

#===============================================
#============= Main program ====================
#
.PHONY: all
all: dnS lib model grid driver readme
# model tests

# test_OOP
.PHONY: test
test:$(TESTEXE)

$(TESTEXE): $(OBJ_test) $(OBJ_all)
	$(LYNK90)   -o $(TESTEXE) $(OBJ_test) $(OBJ_all) $(LYNKFLAGS)

.PHONY: model QML_Test
model QML_Test:$(MODEXE)
	echo "model (QML) compilation: OK"
$(MODEXE): $(OBJ_testmod) $(OBJ_all)
	$(LYNK90)   -o $(MODEXE) $(OBJ_testmod) $(OBJ_all) $(LYNKFLAGS)
#
# grid
.PHONY: grid testgrid
grid testgrid:$(GRIDEXE)
$(GRIDEXE): $(OBJ_grid) $(OBJ_all)
	$(LYNK90)   -o $(GRIDEXE) $(OBJ_grid) $(OBJ_all) $(LYNKFLAGS)
#
# Adia
.PHONY: adia testadia
adia testadia:$(ADIAEXE)
$(ADIAEXE): $(OBJ_adia) $(OBJ_AdiaLib) $(OBJ_all)
	$(LYNK90)   -o $(ADIAEXE) $(OBJ_adia) $(OBJ_AdiaLib) $(OBJ_all) $(LYNKFLAGS)
#
# dnS
.PHONY: dns dnS testdns testdnS
dns dnS testdns testdnS:$(dnSEXE)
		echo "dnS compilation: OK"

$(dnSEXE): $(OBJ_testdnS) $(OBJ_lib)
	$(LYNK90)   -o $(dnSEXE) $(OBJ_testdnS) $(OBJ_lib) $(LYNKFLAGS)
	echo "dnS compilation: OK"

# dnS
.PHONY: dnpoly dnPoly testdnpoly testdnPoly
dnpoly dnPoly testdnpoly testdnPoly: $(dnPolyEXE)
		echo "OBJ_testdnPoly compilation: OK"

$(dnPolyEXE): $(OBJ_testdnPoly) $(OBJ_lib)
	$(LYNK90)   -o $(dnPolyEXE) $(OBJ_testdnPoly) $(OBJ_lib) $(LYNKFLAGS)
	echo "dnPoly compilation: OK"

#
#driver
.PHONY: driver
driver:$(DriverEXE)
$(DriverEXE): $(OBJ_testdriver) $(ModLib)
	$(LYNK90)   -o $(DriverEXE) $(OBJ_testdriver) -L$(DIR0) -lpot $(LYNKFLAGS)
#
#readme
.PHONY: readme
readme:
	bin/extractReadMe

#===============================================
#============= Model Lib =======================
#
.PHONY: lib
lib: $(ModLib) readme
	echo "create the library: ",$(ModLib)
$(ModLib): $(OBJ_driver) $(OBJ_all)
	ar -r $(ModLib) $(OBJ_driver) $(OBJ_all)
	rm -f $(QMLib)
	ln -s $(ModLib) $(QMLib)
#
#
#===============================================
#===============================================
.PHONY: clean
clean:
	rm -f  $(MODEXE) $(GRIDEXE) $(ADIAEXE) $(dnSEXE) $(DriverEXE) $(TESTEXE) $(dnPolyEXE)
	rm -f  $(ModLib) libQMLib.a
	rm -fr *.dSYM comp.log
	cd $(DIROBJ) ; rm -f *.o *.mod *.MOD
	@cd Tests && ./clean
	@echo "  done cleaning up the example directories"
	@cd DOC && ./clean
	@echo "  done cleaning up the documentation"
#===============================================
#===============================================
#
##################################################################################
### Model libraries
#
$(DIROBJ)/Empty_m.o:$(DIRModel)/Empty_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Empty_m.f90

$(DIROBJ)/Morse_m.o:$(DIRModel)/Morse_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Morse_m.f90

$(DIROBJ)/Poly1D_m.o:$(DIRModel)/Poly1D_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Poly1D_m.f90

$(DIROBJ)/H2_m.o:$(DIRModel)/H2_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/H2_m.f90

$(DIROBJ)/Template_m.o:$(DIRModel)/Template_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Template_m.f90

$(DIROBJ)/Test_m.o:$(DIRModel)/Test_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Test_m.f90

$(DIROBJ)/LinearHBond_m.o:$(DIRModel)/LinearHBond_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/LinearHBond_m.f90

$(DIROBJ)/TwoD_MullerBrown_m.o:$(DIRModel)/TwoD_MullerBrown_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/TwoD_MullerBrown_m.f90

$(DIROBJ)/Phenol_m.o:$(DIRModel)/Phenol_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Phenol_m.f90

$(DIROBJ)/PSB3_m.o:$(DIRModel)/PSB3_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/PSB3_m.f90

$(DIROBJ)/Retinal_JPCB2000_m.o:$(DIRModel)/Retinal_JPCB2000_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Retinal_JPCB2000_m.f90

$(DIROBJ)/HONO_m.o:$(DIRModel)/HONO_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/HONO_m.f90

$(DIROBJ)/HNO3_m.o:$(DIRModel)/HNO3_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/HNO3_m.f90

$(DIROBJ)/NO3_m.o:$(DIRModel)/NO3_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/NO3_m.f90

$(DIROBJ)/CH5_m.o:$(DIRModel)/CH5_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/CH5_m.f90

$(DIROBJ)/PH4_m.o:$(DIRModel)/PH4_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/PH4_m.f90

$(DIROBJ)/HOO_DMBE_m.o:$(DIRModel)/HOO_DMBE_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/HOO_DMBE_m.f90

$(DIROBJ)/H3_m.o:$(DIRModel)/H3_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/H3_m.f90

$(DIROBJ)/HNNHp_m.o:$(DIRModel)/HNNHp_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/HNNHp_m.f90

$(DIROBJ)/H2SiN_m.o:$(DIRModel)/H2SiN_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/H2SiN_m.f90

$(DIROBJ)/H2NSi_m.o:$(DIRModel)/H2NSi_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/H2NSi_m.f90

$(DIROBJ)/TwoD_m.o:$(DIRModel)/TwoD_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/TwoD_m.f90

$(DIROBJ)/Tully_m.o:$(DIRModel)/Tully_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Tully_m.f90

$(DIROBJ)/OneDSOC_1S1T_m.o:$(DIRModel)/OneDSOC_1S1T_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/OneDSOC_1S1T_m.f90
$(DIROBJ)/OneDSOC_2S1T_m.o:$(DIRModel)/OneDSOC_2S1T_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/OneDSOC_2S1T_m.f90

$(DIROBJ)/Buck_m.o:$(DIRModel)/Buck_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Buck_m.f90

$(DIROBJ)/Sigmoid_m.o:$(DIRModel)/Sigmoid_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/Sigmoid_m.f90

$(DIROBJ)/HenonHeiles_m.o:$(DIRModel)/HenonHeiles_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/HenonHeiles_m.f90
#
##################################################################################
### QModel
#
$(DIROBJ)/Model_m.o:$(DIRSRC)/Model_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_QML)  -c $(DIRSRC)/Model_m.f90
#
##################################################################################
#
#
##################################################################################
### AdiaChannels
#
$(DIROBJ)/MakeHinact_m.o:$(DIRAdia)/MakeHinact_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIRAdia)/MakeHinact_m.f90
$(DIROBJ)/Basis_m.o:$(DIRAdia)/Basis_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIRAdia)/Basis_m.f90
#
##################################################################################
#
##################################################################################
### Optimization + IRC
#
$(DIROBJ)/Opt_m.o:$(DIROpt)/Opt_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIROpt)/Opt_m.f90
$(DIROBJ)/IRC_m.o:$(DIROpt)/IRC_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIROpt)/IRC_m.f90
#
##################################################################################

#
##################################################################################
### Main + driver + tests
#
$(DIROBJ)/TEST_driver.o:$(DIRSRC)/TEST_driver.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_driver.f90
$(DIROBJ)/Model_driver.o:$(DIRSRC)/Model_driver.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/Model_driver.f90
$(DIROBJ)/TEST_model.o:$(DIRSRC)/TEST_model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_model.f90
$(DIROBJ)/TEST_grid.o:$(DIRSRC)/TEST_grid.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_grid.f90
$(DIROBJ)/TEST_OOP.o:$(DIRSRC)/TEST_OOP.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_OOP.f90
$(DIROBJ)/TEST_Adia.o:$(DIRSRC)/TEST_Adia.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_Adia.f90
#
##################################################################################
#
#
#
##################################################################################
### dnS libraries
#
$(DIROBJ)/dnS_m.o:$(DIRdnS)/dnS_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_INVHYP)  -c $(DIRdnS)/dnS_m.f90
$(DIROBJ)/TEST_dnS.o:$(DIRSRC)/TEST_dnS.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_INVHYP)  -c $(DIRSRC)/TEST_dnS.f90
#
##################################################################################
#
##################################################################################
### dnPoly libraries
#
$(DIROBJ)/dnPoly_m.o:$(DIRdnPoly)/dnPoly_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) -c $(DIRdnPoly)/dnPoly_m.f90
$(DIROBJ)/TEST_dnPoly.o:$(DIRSRC)/TEST_dnPoly.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIRSRC)/TEST_dnPoly.f90
#
##################################################################################
#
#
#
##################################################################################
### dnMat libraries
#
$(DIROBJ)/dnMat_m.o:$(DIRdnMat)/dnMat_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRdnMat)/dnMat_m.f90
#
##################################################################################
#
#
#
##################################################################################
### libraries
#
$(DIROBJ)/FiniteDiff_m.o:$(DIRLib)/FiniteDiff_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/FiniteDiff_m.f90
$(DIROBJ)/UtilLib_m.o:$(DIRLib)/UtilLib_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/UtilLib_m.f90
$(DIROBJ)/diago_m.o:$(DIRLib)/diago_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_MATRIX)  -c $(DIRLib)/diago_m.f90
$(DIROBJ)/Matrix_m.o:$(DIRLib)/Matrix_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_MATRIX)  -c $(DIRLib)/Matrix_m.f90
$(DIROBJ)/NumParameters_m.o:$(DIRLib)/NumParameters_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/NumParameters_m.f90
#
##################################################################################
#
#
##################################################################################
### dependencies
#
$(DIROBJ)/TEST_OOP.o:    $(OBJ_lib) $(OBJ_Model) $(OBJ_QML)
$(DIROBJ)/TEST_model.o:  $(OBJ_lib) $(OBJ_Model) $(OBJ_QML)
$(DIROBJ)/TEST_dnS.o:    $(OBJ_lib)
$(DIROBJ)/TEST_dnPoly.o: $(OBJ_lib)
$(DIROBJ)/TEST_driver.o: $(ModLib)
$(DIROBJ)/TEST_Adia.o:   $(ModLib) $(OBJ_AdiaLib)


$(DIROBJ)/Model_driver.o: $(OBJ_lib) $(OBJ_Model) $(OBJ_QML)

$(DIROBJ)/Opt_m.o: $(OBJ_lib) $(DIROBJ)/Model_m.o $(OBJ_QML)
$(DIROBJ)/IRC_m.o: $(OBJ_lib) $(DIROBJ)/Model_m.o $(DIROBJ)/Opt_m.o $(OBJ_QML)

$(DIROBJ)/Model_m.o: $(DIROBJ)/Basis_m.o $(OBJ_lib) $(OBJ_QML)

$(DIROBJ)/MakeHinact_m.o: $(DIROBJ)/Basis_m.o $(ModLib)


$(DIROBJ)/Empty_m.o:   $(OBJ_lib)
$(DIROBJ)/Morse_m.o:   $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/Poly1D_m.o:   $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/H2_m.o:      $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/Buck_m.o:    $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/Sigmoid_m.o: $(DIROBJ)/Empty_m.o $(OBJ_lib)

$(DIROBJ)/Template_m.o:          $(OBJ_lib) $(DIROBJ)/Empty_m.o $(DIROBJ)/Morse_m.o
$(DIROBJ)/Test_m.o:              $(OBJ_lib) $(DIROBJ)/Empty_m.o

$(DIROBJ)/HenonHeiles_m.o:       $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/Tully_m.o:             $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/OneDSOC_1S1T_m.o:      $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/OneDSOC_2S1T_m.o:      $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/TwoD_m.o:              $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/PSB3_m.o:              $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/Retinal_JPCB2000_m.o:  $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/HONO_m.o:              $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/HNO3_m.o:              $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/NO3_m.o:               $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/HOO_DMBE_m.o:          $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/H3_m.o:                $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/CH5_m.o:               $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/PH4_m.o:               $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/HNNHp_m.o:             $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/H2SiN_m.o:             $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/H2NSi_m.o:             $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/LinearHBond_m.o:       $(DIROBJ)/Empty_m.o $(OBJ_lib) \
                                          $(DIROBJ)/Morse_m.o $(DIROBJ)/Buck_m.o
$(DIROBJ)/TwoD_MullerBrown_m.o:  $(DIROBJ)/Empty_m.o $(OBJ_lib)
$(DIROBJ)/Phenol_m.o:            $(DIROBJ)/Empty_m.o $(OBJ_lib) \
                                       $(DIROBJ)/Morse_m.o $(DIROBJ)/Sigmoid_m.o
#
#
$(DIROBJ)/UtilLib_m.o:    $(DIROBJ)/NumParameters_m.o
$(DIROBJ)/dnS_m.o:        $(DIROBJ)/UtilLib_m.o $(DIROBJ)/NumParameters_m.o
$(DIROBJ)/dnPoly_m.o:     $(DIROBJ)/dnS_m.o $(DIROBJ)/UtilLib_m.o $(DIROBJ)/NumParameters_m.o
$(DIROBJ)/dnMat_m.o:      $(DIROBJ)/dnS_m.o $(DIROBJ)/diago_m.o $(DIROBJ)/UtilLib_m.o $(DIROBJ)/NumParameters_m.o
$(DIROBJ)/diago_m.o:      $(DIROBJ)/NumParameters_m.o
$(DIROBJ)/Matrix_m.o:     $(DIROBJ)/NumParameters_m.o
$(DIROBJ)/FiniteDiff_m.o: $(DIROBJ)/dnMat_m.o $(DIROBJ)/dnS_m.o $(DIROBJ)/NumParameters_m.o
#
############################################################################
### Documentation with doxygen
#
doxy:
	@echo "Installing documentation with doxygen"
	@cd DOC ; doxygen ModLib_doxygen_settings
#
############################################################################
