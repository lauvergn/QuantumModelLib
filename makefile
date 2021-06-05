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
   F90LIB = -mkl -lpthread

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
   ifeq ($(OS),Darwin)    # OSX
      # OSX libs (included lapack+blas)
      F90LIB = -framework Accelerate
   else                   # Linux
      # linux libs
      F90LIB = -llapack -lblas
      #
      # linux libs with mkl and with openmp
      #F90LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
      # linux libs with mkl and without openmp
      #F90LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
   endif
   #
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O5 -g -fbacktrace $(OMPFLAG) -funroll-loops -ftree-vectorize -falign-loops=16
   else
      F90FLAGS = -Og -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
      #F90FLAGS = -O0 -fbounds-check -Wuninitialized
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
CPPSHELL_DIAGO   = -D__LAPACK="$(LAPACK)"

F90_FLAGS = $(F90) $(F90FLAGS)
LYNK90 = $(F90_FLAGS)

 LIBS := $(PESLIB) $(F90LIB) $(ARPACKLIB)
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
DIRLib    = $(DIRSRC)/Lib
DIRdnS    = $(DIRSRC)/dnSLib
DIRdnPoly = $(DIRSRC)/dnPolyLib
DIRdnMat  = $(DIRSRC)/dnMatLib
DIRModel  = $(DIRSRC)/ModelLib
DIRAdia   = $(DIRSRC)/AdiaChannels
#
OBJ_ModelLib   = $(DIROBJ)/mod_EmptyModel.o \
                 $(DIROBJ)/mod_SigmoidModel.o $(DIROBJ)/mod_MorseModel.o $(DIROBJ)/mod_BuckModel.o \
                 $(DIROBJ)/mod_TemplateModel.o $(DIROBJ)/mod_TestModel.o \
                 $(DIROBJ)/mod_H2NSi_Model.o $(DIROBJ)/mod_H2SiN_Model.o \
                 $(DIROBJ)/mod_HNNHp_Model.o $(DIROBJ)/mod_HONO_Model.o \
                 $(DIROBJ)/mod_HNO3_Model.o \
                 $(DIROBJ)/mod_CH5_Model.o $(DIROBJ)/mod_PH4_Model.o \
								 $(DIROBJ)/mod_HOO_DMBE_Model.o \
                 $(DIROBJ)/mod_H3_Model.o \
                 $(DIROBJ)/mod_HenonHeilesModel.o $(DIROBJ)/mod_LinearHBondModel.o \
                 $(DIROBJ)/mod_PhenolModel.o $(DIROBJ)/mod_TwoD_Model.o \
                 $(DIROBJ)/mod_PSB3_Model.o $(DIROBJ)/mod_Retinal_JPCB2000_Model.o \
                 $(DIROBJ)/mod_OneDSOC_1S1T_Model.o $(DIROBJ)/mod_OneDSOC_2S1T_Model.o \
                 $(DIROBJ)/mod_TullyModel.o

OBJ_AdiaLib    = $(DIROBJ)/mod_AdiaChannels.o

OBJ_Model      = $(DIROBJ)/mod_Model.o $(DIROBJ)/QML_Basis_m.o

OBJ_test       = $(DIROBJ)/TEST_OOP.o
OBJ_driver     = $(DIROBJ)/Model_driver.o
OBJ_grid       = $(DIROBJ)/TEST_grid.o
OBJ_adia       = $(DIROBJ)/TEST_Adia.o
OBJ_testmod    = $(DIROBJ)/TEST_model.o
OBJ_testdnS    = $(DIROBJ)/TEST_dnS.o
OBJ_testdnPoly = $(DIROBJ)/TEST_dnPoly.o
OBJ_testdriver = $(DIROBJ)/TEST_driver.o


OBJ_lib        = $(DIROBJ)/QML_FiniteDiff_m.o \
                 $(DIROBJ)/QML_dnMat_m.o $(DIROBJ)/QML_dnPoly_m.o $(DIROBJ)/QML_dnS_m.o \
                 $(DIROBJ)/QML_UtilLib_m.o $(DIROBJ)/QML_diago_m.o \
                 $(DIROBJ)/QML_NumParameters_m.o

OBJ_all        = $(OBJ_lib) $(OBJ_Model) $(OBJ_ModelLib)

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

.PHONY: model testmodel
model testmodel:$(MODEXE)
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
	$(LYNK90)   -o $(DriverEXE) $(OBJ_testdriver) $(LYNKFLAGS) -L$(DIR0) -lpot
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
	rm -f  $(MODEXE) $(GRIDEXE) $(ADIAEXE) $(dnSEXE) $(DriverEXE) $(TESTEXE)
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
$(DIROBJ)/mod_EmptyModel.o:$(DIRModel)/mod_EmptyModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_EmptyModel.f90

$(DIROBJ)/mod_MorseModel.o:$(DIRModel)/mod_MorseModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_MorseModel.f90

$(DIROBJ)/mod_TemplateModel.o:$(DIRModel)/mod_TemplateModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_TemplateModel.f90

$(DIROBJ)/mod_TestModel.o:$(DIRModel)/mod_TestModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_TestModel.f90

$(DIROBJ)/mod_LinearHBondModel.o:$(DIRModel)/mod_LinearHBondModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_LinearHBondModel.f90

$(DIROBJ)/mod_PhenolModel.o:$(DIRModel)/mod_PhenolModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_PhenolModel.f90

$(DIROBJ)/mod_PSB3_Model.o:$(DIRModel)/mod_PSB3_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_PSB3_Model.f90

$(DIROBJ)/mod_Retinal_JPCB2000_Model.o:$(DIRModel)/mod_Retinal_JPCB2000_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_Retinal_JPCB2000_Model.f90

$(DIROBJ)/mod_HONO_Model.o:$(DIRModel)/mod_HONO_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_HONO_Model.f90

$(DIROBJ)/mod_HNO3_Model.o:$(DIRModel)/mod_HNO3_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_HNO3_Model.f90

$(DIROBJ)/mod_CH5_Model.o:$(DIRModel)/mod_CH5_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_CH5_Model.f90

$(DIROBJ)/mod_PH4_Model.o:$(DIRModel)/mod_PH4_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_PH4_Model.f90

$(DIROBJ)/mod_HOO_DMBE_Model.o:$(DIRModel)/mod_HOO_DMBE_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_HOO_DMBE_Model.f90

$(DIROBJ)/mod_H3_Model.o:$(DIRModel)/mod_H3_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_H3_Model.f90

$(DIROBJ)/mod_HNNHp_Model.o:$(DIRModel)/mod_HNNHp_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_HNNHp_Model.f90

$(DIROBJ)/mod_H2SiN_Model.o:$(DIRModel)/mod_H2SiN_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_H2SiN_Model.f90

$(DIROBJ)/mod_H2NSi_Model.o:$(DIRModel)/mod_H2NSi_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_H2NSi_Model.f90

$(DIROBJ)/mod_TwoD_Model.o:$(DIRModel)/mod_TwoD_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_TwoD_Model.f90

$(DIROBJ)/mod_TullyModel.o:$(DIRModel)/mod_TullyModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_TullyModel.f90

$(DIROBJ)/mod_OneDSOC_1S1T_Model.o:$(DIRModel)/mod_OneDSOC_1S1T_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_OneDSOC_1S1T_Model.f90
$(DIROBJ)/mod_OneDSOC_2S1T_Model.o:$(DIRModel)/mod_OneDSOC_2S1T_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_OneDSOC_2S1T_Model.f90

$(DIROBJ)/mod_BuckModel.o:$(DIRModel)/mod_BuckModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_BuckModel.f90

$(DIROBJ)/mod_SigmoidModel.o:$(DIRModel)/mod_SigmoidModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_SigmoidModel.f90

$(DIROBJ)/mod_HenonHeilesModel.o:$(DIRModel)/mod_HenonHeilesModel.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRModel)/mod_HenonHeilesModel.f90
#
##################################################################################
### QModel
#
$(DIROBJ)/mod_Model.o:$(DIRSRC)/mod_Model.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_QML)  -c $(DIRSRC)/mod_Model.f90
#
##################################################################################
#
#
##################################################################################
### AdiaChannels
#
$(DIROBJ)/mod_AdiaChannels.o:$(DIRAdia)/mod_AdiaChannels.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIRAdia)/mod_AdiaChannels.f90
$(DIROBJ)/QML_Basis_m.o:$(DIRAdia)/QML_Basis_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIRAdia)/QML_Basis_m.f90
#
##################################################################################
#
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
$(DIROBJ)/QML_dnS_m.o:$(DIRdnS)/QML_dnS_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_INVHYP)  -c $(DIRdnS)/QML_dnS_m.f90
$(DIROBJ)/TEST_dnS.o:$(DIRdnS)/TEST_dnS.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_INVHYP)  -c $(DIRdnS)/TEST_dnS.f90
#
##################################################################################
#
##################################################################################
### dnPoly libraries
#
$(DIROBJ)/QML_dnPoly_m.o:$(DIRdnPoly)/QML_dnPoly_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) -c $(DIRdnPoly)/QML_dnPoly_m.f90
$(DIROBJ)/TEST_dnPoly.o:$(DIRdnPoly)/TEST_dnPoly.f90
	cd $(DIROBJ) ; $(F90_FLAGS)  -c $(DIRdnPoly)/TEST_dnPoly.f90
#
##################################################################################
#
#
#
##################################################################################
### dnMat libraries
#
$(DIROBJ)/QML_dnMat_m.o:$(DIRdnMat)/QML_dnMat_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRdnMat)/QML_dnMat_m.f90
#
##################################################################################
#
#
#
##################################################################################
### libraries
#
$(DIROBJ)/QML_FiniteDiff_m.o:$(DIRLib)/QML_FiniteDiff_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/QML_FiniteDiff_m.f90
$(DIROBJ)/QML_UtilLib_m.o:$(DIRLib)/QML_UtilLib_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/QML_UtilLib_m.f90
$(DIROBJ)/QML_diago_m.o:$(DIRLib)/QML_diago_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_DIAGO)  -c $(DIRLib)/QML_diago_m.f90
$(DIROBJ)/QML_NumParameters_m.o:$(DIRLib)/QML_NumParameters_m.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/QML_NumParameters_m.f90
#
##################################################################################
#
#
##################################################################################
### dependencies
#
$(DIROBJ)/TEST_OOP.o:    $(OBJ_lib) $(OBJ_Model) $(OBJ_ModelLib)
$(DIROBJ)/TEST_model.o:  $(OBJ_lib) $(OBJ_Model) $(OBJ_ModelLib)
$(DIROBJ)/TEST_dnS.o:    $(OBJ_lib)
$(DIROBJ)/TEST_dnPoly.o: $(OBJ_lib)
$(DIROBJ)/TEST_driver.o: $(ModLib)
$(DIROBJ)/TEST_Adia.o:   $(ModLib) $(OBJ_AdiaLib)


$(DIROBJ)/Model_driver.o: $(OBJ_lib) $(OBJ_Model) $(OBJ_ModelLib)

$(DIROBJ)/mod_Model.o: $(DIROBJ)/QML_Basis_m.o $(OBJ_lib) $(OBJ_ModelLib)

$(DIROBJ)/mod_AdiaChannels.o: $(DIROBJ)/QML_Basis_m.o $(ModLib)


$(DIROBJ)/mod_EmptyModel.o: $(OBJ_lib)
$(DIROBJ)/mod_MorseModel.o: $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_BuckModel.o: $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_SigmoidModel.o: $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)

$(DIROBJ)/mod_TemplateModel.o: $(OBJ_lib) $(DIROBJ)/mod_EmptyModel.o $(DIROBJ)/mod_MorseModel.o
$(DIROBJ)/mod_TestModel.o: $(OBJ_lib) $(DIROBJ)/mod_EmptyModel.o

$(DIROBJ)/mod_HenonHeilesModel.o:        $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_TullyModel.o:              $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_OneDSOC_1S1T_Model.o:      $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_OneDSOC_2S1T_Model.o:      $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_TwoD_Model.o:              $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_PSB3_Model.o:              $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_Retinal_JPCB2000_Model.o:  $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_HONO_Model.o:              $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_HNO3_Model.o:              $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_HOO_DMBE_Model.o:          $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_H3_Model.o:                $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_CH5_Model.o:               $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_PH4_Model.o:               $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_HNNHp_Model.o:             $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_H2SiN_Model.o:             $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_H2NSi_Model.o:             $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib)
$(DIROBJ)/mod_LinearHBondModel.o:        $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib) \
                                         $(DIROBJ)/mod_MorseModel.o $(DIROBJ)/mod_BuckModel.o
$(DIROBJ)/mod_PhenolModel.o:             $(DIROBJ)/mod_EmptyModel.o $(OBJ_lib) \
                                         $(DIROBJ)/mod_MorseModel.o $(DIROBJ)/mod_SigmoidModel.o
#
#
$(DIROBJ)/QML_UtilLib_m.o:    $(DIROBJ)/QML_NumParameters_m.o
$(DIROBJ)/QML_dnS_m.o:        $(DIROBJ)/QML_UtilLib_m.o $(DIROBJ)/QML_NumParameters_m.o
$(DIROBJ)/QML_dnPoly_m.o:     $(DIROBJ)/QML_dnS_m.o $(DIROBJ)/QML_UtilLib_m.o $(DIROBJ)/QML_NumParameters_m.o
$(DIROBJ)/QML_dnMat_m.o:      $(DIROBJ)/QML_dnS_m.o $(DIROBJ)/QML_diago_m.o $(DIROBJ)/QML_UtilLib_m.o $(DIROBJ)/QML_NumParameters_m.o
$(DIROBJ)/QML_diago_m.o:      $(DIROBJ)/QML_NumParameters_m.o
$(DIROBJ)/QML_FiniteDiff_m.o: $(DIROBJ)/QML_dnMat_m.o $(DIROBJ)/QML_dnS_m.o $(DIROBJ)/QML_NumParameters_m.o
#
############################################################################
### Documentation with doxygen
#
doxy:
	@echo "Installing documentation with doxygen"
	@cd DOC ; doxygen ModLib_doxygen_settings
#
############################################################################
