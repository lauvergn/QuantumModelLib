#=================================================================================
#=================================================================================
# Compiler? 
#Possible values: 
#                ifort (version: 14.0.2, 16.0.3, 17.0.1 linux)
#                gfortran (version: 6.3.0 linux and osx)
#                pgf90 (version: 17.10-0, linux): problems because datanh, dasinh, dacosh are not available with this compiler !!!
F90 = gfortran
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 0
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
#=================================================================================


# Operating system, OS? automatic using uname: 
OS=$(shell uname)
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
      F90FLAGS = -O0 -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized
      #F90FLAGS = -O0 -fbounds-check -Wuninitialized
   endif
endif
#=================================================================================
#=================================================================================
$(info ***********************************************************************)
$(info ***********OS:           $(OS))
$(info ***********COMPILER:     $(F90))
$(info ***********OPTIMIZATION: $(OPT))
$(info ***********OpenMP:       $(OMPFLAG))
$(info ***********Arpack:       $(ARPACK))
$(info ***********F90FLAGS:     $(F90FLAGS))
$(info ***********F90LIB:       $(F90LIB))
$(info ***********************************************************************)


F90_FLAGS = $(F90) $(F90FLAGS)
LYNK90 = $(F90_FLAGS)

 LIBS := $(PESLIB) $(F90LIB) $(ARPACKLIB)
 LYNKFLAGS = $(LIBS)



#=================================================================================
# for c++ preprocessing
#=================================================================================
 CPP    = -cpp
#=================================================================================
#=================================================================================
#
GRIDEXE   = Grid.x
MODEXE    = ModLib.x
dnSEXE    = dnS.x
DriverEXE = Driver.x
ModLib    = libpot.a

DIR0      = $(shell pwd)
DIROBJ    = $(DIR0)/OBJ
DIRSRC    = $(DIR0)/SRC
DIRLib    = $(DIRSRC)/Lib
DIRdnS    = $(DIRSRC)/dnSLib
DIRdnMat  = $(DIRSRC)/dnMatLib
DIRPot    = $(DIRSRC)/PotLib

#
OBJ_lib        = $(DIROBJ)/dnMatPot_Module.o $(DIROBJ)/dnS_Module.o $(DIROBJ)/Lib_module.o $(DIROBJ)/sub_diago.o $(DIROBJ)/sub_module_NumParameters.o
OBJ_Pot        = $(DIROBJ)/LinearHBondPotential_Module.o $(DIROBJ)/TullyPotential_Module.o $(DIROBJ)/PhenolPotential_Module.o \
                 $(DIROBJ)/TemplatePotential_Module.o \
                 $(DIROBJ)/HenonHeilesPotential_Module.o \
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
all: dnS lib model grid driver
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
$(dnSEXE): $(OBJ_testdnS) $(OBJ_all)
	$(LYNK90)   -o $(dnSEXE) $(OBJ_testdnS) $(OBJ_all) $(LYNKFLAGS)
#
#driver
driver:$(DriverEXE)
$(DriverEXE): $(OBJ_testdriver) $(ModLib)
	$(LYNK90)   -o $(DriverEXE) $(OBJ_testdriver) $(LYNKFLAGS) -L$(DIR0) -lpot


#===============================================
#============= Model Lib =======================
#
lib: $(ModLib)
	echo "create the library: ",$(ModLib)
$(ModLib): $(OBJ_driver) $(OBJ_all)
	ar -r $(ModLib) $(OBJ_driver) $(OBJ_all)
#
#
#===============================================
#===============================================
clean: 
	rm -f  $(MODEXE) $(GRIDEXE) $(dnSEXE) $(DriverEXE) $(ModLib)
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
$(DIROBJ)/TEST_dnS.o:    $(OBJ_lib) $(OBJ_Pot) $(OBJ_Model)
$(DIROBJ)/TEST_driver.o: $(ModLib)

$(DIROBJ)/Model_driver.o: $(OBJ_lib) $(OBJ_Pot) $(OBJ_Model)

$(DIROBJ)/Model_Module.o: $(OBJ_lib) $(OBJ_Pot)

$(DIROBJ)/BuckinghamPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/MorsePotential_Module.o: $(OBJ_lib)
$(DIROBJ)/SigmoidPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/HenonHeilesPotential_Module.o: $(OBJ_lib)
$(DIROBJ)/TullyPotential_Module.o: $(OBJ_lib)
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
$(DIROBJ)/TullyPotential_Module.o:$(DIRPot)/TullyPotential_Module.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRPot)/TullyPotential_Module.f90
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
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRSRC)/Model_Module.f90

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
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRdnS)/dnS_Module.f90
$(DIROBJ)/TEST_dnS.o:$(DIRdnS)/TEST_dnS.f90
	cd $(DIROBJ) ; $(F90_FLAGS)   -c $(DIRdnS)/TEST_dnS.f90
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
