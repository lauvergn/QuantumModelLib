#=================================================================================
#=================================================================================
# Compiler? 
#Possible values: (Empty: gfortran)
#                ifort (version: 14.0.2, 16.0.3, 17.0.1 linux)
#                gfortran (version: 6.3.0 linux and osx)
#                pgf90 (version: 17.10-0, linux)
#                nagfor (version 7.0, osx)
F90 = gfortran
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 0
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
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
 CPP    = -cpp
#=================================================================================
#=================================================================================

#=================================================================================
# nag compillation (nagfor)
#=================================================================================
ifeq ($(F90),nagfor)
   # for c++ preprocessing
   CPP = -fpp
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
      #  -C=undefined is not compatible with -framework Accelerate
      # −kind=byte and −dcfuns is not working
      #with -mtrace=all add information on the memmory allocation/deallocation.
      ifeq ($(OMP),0)
        F90FLAGS = -O0 $(OMPFLAG) -g -gline -C -C=alias -C=intovf
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
# ifort compillation v12 with mkl
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
   CPP    = -Mpreprocess

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
endif

QML_ver=$(shell awk '/QML/ {print $$3}' version-QML)
QML_path=$(shell pwd)
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


CPPSHELL_QML_ver_Path  = -D__QML_VER='"$(QML_ver)"' \
                         -D__QMLPATH="'$(QML_path)'"


CPPSHELL_INVHYP  = -D__INVHYP="$(INVHYP)"

F90_FLAGS = $(F90) $(F90FLAGS)
LYNK90 = $(F90_FLAGS)

 LIBS := $(PESLIB) $(F90LIB) $(ARPACKLIB)
 LYNKFLAGS = $(LIBS)


#
GRIDEXE   = Grid.x
MODEXE    = ModLib.x
dnSEXE    = dnS.x
DriverEXE = Driver.x
ModLib    = libpot.a
QMLib     = libQMLib.a

DIR0      = $(shell pwd)
DIROBJ    = $(DIR0)/OBJ
DIRSRC    = $(DIR0)/SRC
DIRLib    = $(DIRSRC)/Lib
DIRdnS    = $(DIRSRC)/dnSLib
DIRdnMat  = $(DIRSRC)/dnMatLib
DIRPot    = $(DIRSRC)/PotLib

#
OBJ_lib        = $(DIROBJ)/dnMatPot_Module.o $(DIROBJ)/dnS_Module.o $(DIROBJ)/Lib_module.o $(DIROBJ)/sub_diago.o $(DIROBJ)/sub_module_NumParameters.o
OBJ_Pot        = $(DIROBJ)/LinearHBondPotential_Module.o $(DIROBJ)/TullyPotential_Module.o \
                 $(DIROBJ)/PhenolPotential_Module.o $(DIROBJ)/PSB3Potential_Module.o \
                 $(DIROBJ)/SOC_1S1T_1DModel_Module.o $(DIROBJ)/SOC_2S1T_1DModel_Module.o \
                 $(DIROBJ)/HONOPotential_Module.o $(DIROBJ)/HNNHp_Module.o \
                 $(DIROBJ)/H2SiN_Module.o $(DIROBJ)/H2NSi_Module.o \
                 $(DIROBJ)/TemplatePotential_Module.o \
                 $(DIROBJ)/HenonHeilesPotential_Module.o \
                 $(DIROBJ)/TwoD_Potential_Module.o \
                 $(DIROBJ)/BuckinghamPotential_Module.o $(DIROBJ)/MorsePotential_Module.o $(DIROBJ)/SigmoidPotential_Module.o

OBJ_Model      = $(DIROBJ)/Model_Module.o

OBJ_driver     = $(DIROBJ)/Model_driver.o
OBJ_grid       = $(DIROBJ)/TEST_grid.o
OBJ_testmod    = $(DIROBJ)/TEST_model.o
OBJ_testdnS    = $(DIROBJ)/TEST_dnS.o
OBJ_testdriver = $(DIROBJ)/TEST_driver.o


OBJ_all        = $(OBJ_lib) $(OBJ_Pot) $(OBJ_Model)

#===============================================
#============= Main program ====================
#
all: dnS lib model grid driver readme
# model tests
model:$(MODEXE)
testmodel:$(MODEXE)
$(MODEXE): $(OBJ_testmod) $(OBJ_all)
	$(LYNK90)   -o $(MODEXE) $(OBJ_testmod) $(OBJ_all) $(LYNKFLAGS)
#
# grid
grid:$(GRIDEXE)
testgrid:$(GRIDEXE)
$(GRIDEXE): $(OBJ_grid) $(OBJ_all)
	$(LYNK90)   -o $(GRIDEXE) $(OBJ_grid) $(OBJ_all) $(LYNKFLAGS)
#
# dnS
dns:$(dnSEXE)
dnS:$(dnSEXE)
testdns:$(dnSEXE)
testdnS:$(dnSEXE)
$(dnSEXE): $(OBJ_testdnS) $(OBJ_lib)
	$(LYNK90)   -o $(dnSEXE) $(OBJ_testdnS) $(OBJ_lib) $(LYNKFLAGS)
#
#driver
driver:$(DriverEXE)
$(DriverEXE): $(OBJ_testdriver) $(ModLib)
	$(LYNK90)   -o $(DriverEXE) $(OBJ_testdriver) $(LYNKFLAGS) -L$(DIR0) -lpot
#
#readme
readme:
	bin/extractReadMe

#===============================================
#============= Model Lib =======================
#
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
clean: 
	rm -f  $(MODEXE) $(GRIDEXE) $(dnSEXE) $(DriverEXE) $(ModLib) libQMLib.a
	rm -fr *.dSYM
	cd $(DIROBJ) ; rm -f *.o *.mod *.MOD
	@cd Tests && ./clean
	@echo "  done cleaning up the example directories"
	@cd DOC && ./clean
	@echo "  done cleaning up the documentation"
#===============================================
#===============================================
#
$(DIROBJ)/TEST_grid.o:   $(OBJ_lib) $(OBJ_Pot) $(OBJ_Model)
$(DIROBJ)/TEST_model.o:  $(OBJ_lib) $(OBJ_Pot) $(OBJ_Model)
$(DIROBJ)/TEST_dnS.o:    $(OBJ_lib)
$(DIROBJ)/TEST_driver.o: $(ModLib)

$(DIROBJ)/Model_driver.o: $(OBJ_lib) $(OBJ_Pot) $(OBJ_Model)

$(DIROBJ)/Model_Module.o: $(OBJ_lib) $(OBJ_Pot)

$(DIROBJ)/BuckinghamPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/MorsePotential_Module.o: $(OBJ_lib)
$(DIROBJ)/SigmoidPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/HenonHeilesPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/TullyPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/SOC_1S1T_1DModel_Module.o: $(OBJ_lib)
$(DIROBJ)/SOC_2S1T_1DModel_Module.o: $(OBJ_lib)
$(DIROBJ)/TwoD_Potential_Module.o: $(OBJ_lib)
$(DIROBJ)/PSB3Potential_Module.o: $(OBJ_lib)
$(DIROBJ)/HONOPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/HNNHp_Module.o: $(OBJ_lib)
$(DIROBJ)/H2SiN_Module.o: $(OBJ_lib)
$(DIROBJ)/H2NSi_Module.o: $(OBJ_lib)
$(DIROBJ)/LinearHBondPotential_Module.o: $(OBJ_lib) $(DIROBJ)/MorsePotential_Module.o $(DIROBJ)/BuckinghamPotential_Module.o
$(DIROBJ)/PhenolPotential_Module.o: $(OBJ_lib) $(DIROBJ)/MorsePotential_Module.o $(DIROBJ)/SigmoidPotential_Module.o
$(DIROBJ)/TemplatePotential_Module.o: $(OBJ_lib) $(DIROBJ)/MorsePotential_Module.o


$(DIROBJ)/Lib_module.o: $(DIROBJ)/sub_module_NumParameters.o
$(DIROBJ)/dnS_Module.o: $(DIROBJ)/Lib_module.o $(DIROBJ)/sub_module_NumParameters.o
$(DIROBJ)/dnMatPot_Module.o: $(DIROBJ)/dnS_Module.o $(DIROBJ)/Lib_module.o $(DIROBJ)/sub_module_NumParameters.o
$(DIROBJ)/sub_diago.o: $(DIROBJ)/sub_module_NumParameters.o

##################################################################################
### Potential libraries
#
$(DIROBJ)/TemplatePotential_Module.o:$(DIRPot)/TemplatePotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/TemplatePotential_Module.f90

$(DIROBJ)/LinearHBondPotential_Module.o:$(DIRPot)/LinearHBondPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/LinearHBondPotential_Module.f90

$(DIROBJ)/PhenolPotential_Module.o:$(DIRPot)/PhenolPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/PhenolPotential_Module.f90

$(DIROBJ)/PSB3Potential_Module.o:$(DIRPot)/PSB3Potential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/PSB3Potential_Module.f90

$(DIROBJ)/HONOPotential_Module.o:$(DIRPot)/HONOPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/HONOPotential_Module.f90

$(DIROBJ)/HNNHp_Module.o:$(DIRPot)/HNNHp_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/HNNHp_Module.f90

$(DIROBJ)/H2SiN_Module.o:$(DIRPot)/H2SiN_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/H2SiN_Module.f90

$(DIROBJ)/H2NSi_Module.o:$(DIRPot)/H2NSi_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/H2NSi_Module.f90

$(DIROBJ)/TwoD_Potential_Module.o:$(DIRPot)/TwoD_Potential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/TwoD_Potential_Module.f90

$(DIROBJ)/TullyPotential_Module.o:$(DIRPot)/TullyPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/TullyPotential_Module.f90

$(DIROBJ)/SOC_1S1T_1DModel_Module.o:$(DIRPot)/SOC_1S1T_1DModel_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/SOC_1S1T_1DModel_Module.f90
$(DIROBJ)/SOC_2S1T_1DModel_Module.o:$(DIRPot)/SOC_2S1T_1DModel_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/SOC_2S1T_1DModel_Module.f90
$(DIROBJ)/MorsePotential_Module.o:$(DIRPot)/MorsePotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/MorsePotential_Module.f90

$(DIROBJ)/BuckinghamPotential_Module.o:$(DIRPot)/BuckinghamPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/BuckinghamPotential_Module.f90

$(DIROBJ)/SigmoidPotential_Module.o:$(DIRPot)/SigmoidPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/SigmoidPotential_Module.f90

$(DIROBJ)/HenonHeilesPotential_Module.o:$(DIRPot)/HenonHeilesPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/HenonHeilesPotential_Module.f90
#
##################################################################################
#
#
#
##################################################################################
### Model libraries
#
$(DIROBJ)/Model_Module.o:$(DIRSRC)/Model_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPP) $(CPPSHELL_QML_ver_Path)  -c $(DIRSRC)/Model_Module.f90

$(DIROBJ)/TEST_driver.o:$(DIRSRC)/TEST_driver.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_driver.f90
$(DIROBJ)/Model_driver.o:$(DIRSRC)/Model_driver.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/Model_driver.f90
$(DIROBJ)/TEST_model.o:$(DIRSRC)/TEST_model.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_model.f90
$(DIROBJ)/TEST_grid.o:$(DIRSRC)/TEST_grid.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/TEST_grid.f90
#
##################################################################################
#
#
#
##################################################################################
### dnS libraries
#
$(DIROBJ)/dnS_Module.o:$(DIRdnS)/dnS_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPP) $(CPPSHELL_INVHYP)  -c $(DIRdnS)/dnS_Module.f90
$(DIROBJ)/TEST_dnS.o:$(DIRdnS)/TEST_dnS.f90
	cd $(DIROBJ) ; $(F90_FLAGS) $(CPP) $(CPPSHELL_INVHYP)  -c $(DIRdnS)/TEST_dnS.f90
#
##################################################################################
#
#
#
##################################################################################
### dnMat libraries
#
$(DIROBJ)/dnMatPot_Module.o:$(DIRdnMat)/dnMatPot_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRdnMat)/dnMatPot_Module.f90
#
##################################################################################
#
#
#
##################################################################################
### libraries
#
$(DIROBJ)/Lib_module.o:$(DIRLib)/Lib_module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/Lib_module.f90
$(DIROBJ)/sub_diago.o:$(DIRLib)/sub_diago.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/sub_diago.f90
$(DIROBJ)/sub_module_NumParameters.o:$(DIRLib)/sub_module_NumParameters.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRLib)/sub_module_NumParameters.f90
#
##################################################################################
#
#
#
############################################################################
### Documentation with doxygen
#
doxy:
	@echo "Installing documentation with doxygen"
	@cd DOC ; doxygen ModLib_doxygen_settings
#
############################################################################
