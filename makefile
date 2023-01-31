#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx)
 FC = gfortran
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 1
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK = 1
## how to get external libraries;  "loc" (default): from local zip file, Empty or something else (v0.5): from github
EXTLIB_TYPE = loc
#=================================================================================
#=================================================================================
ifeq ($(FC),)
  FFC      := gfortran
else
  FFC      := $(FC)
endif
ifeq ($(OPT),)
  OOPT      := 1
else
  OOPT      := $(OPT)
endif
ifeq ($(OMP),)
  OOMP      := 1
else
  OOMP      := $(OMP)
endif
ifeq ($(LAPACK),)
  LLAPACK      := 1
else
  LLAPACK      := $(LAPACK)
endif
#=================================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)


QML_ver  := $(shell awk '/QML/ {print $$3}' version-QML)
QML_path := $(shell pwd)

# Extension for the object directory and the library
ext_obj=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)

# library name
QMLIBA=libQMLib$(ext_obj).a
#=================================================================================

OBJ_DIR=OBJ/obj$(ext_obj)
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))

MOD_DIR=$(OBJ_DIR)
SRC_DIR=SRC
MAIN_DIR=APP
TESTS_DIR=Tests

CPPSHELL_QML = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
               -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
               -D__COMPILER="'$(FFC)'" \
               -D__COMPILER_VER="'$(FC_VER)'" \
               -D__COMPILER_OPT="'$(FFLAGS0)'" \
               -D__QMLPATH="'$(QML_path)'" \
               -D__QML_VER='"$(QML_ver)"'
#=================================================================================
# External Libraries directory
#
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(QML_path)/Ext_Lib
endif

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA    = $(QD_DIR)/libQD$(ext_obj).a

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(ext_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(ext_obj).a

EXTMOD_DIR = $(QDMOD_DIR) $(ADMOD_DIR)
EXTLib     = $(ADLIBA) $(QDLIBA)
#===============================================================================
#=================================================================================
#=================================================================================
# gfortran (osx and linux)
#=================================================================================
ifeq ($(FFC),gfortran)

  ifeq ($(OOPT),1)
    FFLAGS = -O5 -g -fbacktrace -funroll-loops -ftree-vectorize -falign-loops=16
  else
    FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
    #FFLAGS = -O0 -fbounds-check -Wuninitialized
  endif

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -fopenmp
  endif
  FFLAGS0 := $(FFLAGS)

  # where to store the .mod files
  FFLAGS +=-J$(MOD_DIR)

  # where to look the .mod files
  FFLAGS += -I$(QDMOD_DIR) -I$(ADMOD_DIR)

  # some cpreprocessing
  FFLAGS += -cpp $(CPPSHELL_QML)

  FLIB   = $(EXTLib)
  # OS management
  ifeq ($(LLAPACK),1)
    ifeq ($(OS),Darwin)    # OSX
      # OSX libs (included lapack+blas)
      FLIB += -framework Accelerate
    else                   # Linux
      # linux libs
      FLIB += -llapack -lblas
      #
      # linux libs with mkl and with openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
      # linux libs with mkl and without openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
    endif
  endif

  FC_VER = $(shell $(FFC) --version | head -1 )

endif
#=================================================================================
#=================================================================================

$(info ***********************************************************************)
$(info ***********OS:           $(OS))
$(info ***********COMPILER:     $(FFC))
$(info ***********COMPILER_VER: $(FC_VER))
$(info ***********OPTIMIZATION: $(OOPT))
$(info ***********OpenMP:       $(OOMP))
$(info ***********LAPACK:       $(LLAPACK))
$(info ***********FFLAGS:       $(FFLAGS))
$(info ***********FLIB:         $(FLIB))
$(info ***********EXTMOD_DIR:   $(EXTMOD_DIR))
$(info ***********ext_obj:      $(ext_obj))
$(info ***********************************************************************)

VPATH = $(TESTS_DIR):$(MAIN_DIR):$(SRC_DIR): \
        $(SRC_DIR)/QMLLib:$(SRC_DIR)/QML:$(SRC_DIR)/Opt:$(SRC_DIR)/AdiaChannels

MAIN=TEST_driver
TESTS=TEST_model

MAINSRCFILES=TEST_Adia.f90 TEST_OOP.f90 TEST_grid.f90 TEST_OMPloop.f90 TEST_driver.f90
TESTSRCFILES=TEST_model.f90
LIBSRCFILES=Model_driver.f90 Model_m.f90 \
            IRC_m.f90    Opt_m.f90 \
            Basis_m.f90  MakeHinact_m.f90

QMLSRCFILES=Buck_m.f90 H3_m.f90 NO3_m.f90  Template_m.f90 \
            CH5_m.f90 HCN_Murrell_m.f90    OneDSOC_1S1T_m.f90 Test_m.f90 \
            ClH2p_Botschwina_m.f90         HNNHp_m.f90 OneDSOC_2S1T_m.f90 Tully_m.f90 \
            ClH2p_m.f90  HNO3_m.f90        PH4_m.f90   TwoD_MullerBrown_m.f90 \
            Empty_m.f90  HONO_m.f90        PSB3_m.f90  TwoD_RJDI2014_m.f90 \
            H2NSi_m.f90  HOO_DMBE_m.f90    Phenol_m.f90           TwoD_Valahu2022_m.f90 \
            H2O_m.f90    HenonHeiles_m.f90 Poly1D_m.f90           TwoD_m.f90 \
            H2SiN_m.f90  LinearHBond_m.f90 Retinal_JPCB2000_m.f90   \
            H2_m.f90     Morse_m.f90       Sigmoid_m.f90 \
            FiniteDiff_m.f90  UtilLib_m.f90

QMLOBJ0=${QMLSRCFILES:.f90=.o}
QMLOBJ=$(addprefix $(OBJ_DIR)/, $(QMLOBJ0))



SRCFILES= $(QMLSRCFILES) $(LIBSRCFILES)
OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))

#===============================================
#============= Tests ===========================
#===============================================
.PHONY: ut UT
UT ut: $(TESTS).x
	@echo "model (QML) compilation: OK"
#	./$(TESTS).x |	grep "Number of error(s)"
#	@echo "  done Tests"


#===============================================
#============= all: lib, tests ...  ============
#===============================================
.PHONY: all
all: $(QMLIBA) $(MAIN).x $(TESTS).x
#===============================================
#============= Main executable and tests  ======
#===============================================
$(MAIN).x: $(OBJ_DIR)/$(MAIN).o $(QMLIBA)
	$(FFC) $(FFLAGS) -o $(MAIN).x  $(OBJ_DIR)/$(MAIN).o $(QMLIBA) $(FLIB)

$(TESTS).x: $(OBJ_DIR)/$(TESTS).o $(QMLIBA)
	$(FFC) $(FFLAGS) -o $(TESTS).x  $(OBJ_DIR)/$(TESTS).o $(QMLIBA) $(FLIB)
#===============================================
#============= Library: libQD.a  ===============
#===============================================
.PHONY: lib
lib: $(QMLIBA)

$(QMLIBA): $(OBJ)
	ar -cr $(QMLIBA) $(OBJ)
	@echo "  done Library: "$(QMLIBA)
#===============================================
#===============================================
#
#===============================================
#============= compilation =====================
#===============================================
$(OBJ_DIR)/%.o: %.f90
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<

#===============================================
#================ cleaning =====================
.PHONY: clean cleanall
clean:
	cd $(TESTS_DIR) ; ./clean
	rm -f $(OBJ_DIR)/*/*.o $(OBJ_DIR)/*.o
	rm -f *.log 
	rm -f TEST*.x
	@echo "  done cleaning"

cleanall : clean
	cd $(ExtLibDIR) ; ./cleanlib
	rm -fr OBJ/obj* OBJ/*mod build
	rm -f lib*.a
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := QuantumModelLib
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	cd $(ExtLibSAVEDIR) ; rm -rf $(BaseName)_devloc
	mkdir $(ExtLibSAVEDIR)/$(BaseName)_devloc
	cp -r * $(ExtLibSAVEDIR)/$(BaseName)_devloc
	cd $(ExtLibSAVEDIR) ; zip -r Save_$(BaseName)_devloc.zip $(BaseName)_devloc
	cd $(ExtLibSAVEDIR) ; rm -rf $(BaseName)_devloc
	@echo "  done zip"
#=== external libraries ========================
# AD_dnSVM + QDUTIL Lib
#===============================================
#
$(QDLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QD_DIR) || (cd $(ExtLibDIR) ; ./get_QDUtilLib.sh $(EXTLIB_TYPE))
	@test -d $(QD_DIR) || (echo $(QD_DIR) "does not exist" ; exit 1)
	cd $(QD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(AD_DIR) || (cd $(ExtLibDIR) ; ./get_AD_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(AD_DIR) || (echo $(AD_DIR) "does not exist" ; exit 1)
	cd $(AD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(AD_DIR) " in "$(BaseName)
#===============================================
#============= module dependencies =============
#===============================================
##################################################################################
### dependencies
#
$(OBJ_DIR)/UtilLib_m.o:    $(EXTLib)
$(OBJ_DIR)/FiniteDiff_m.o: $(EXTLib) $(OBJ_DIR)/UtilLib_m.o
$(OBJ_DIR)/Basis_m.o:      $(OBJ_DIR)/UtilLib_m.o

$(OBJ_DIR)/Empty_m.o:      $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/FiniteDiff_m.o

$(OBJ_DIR)/Morse_m.o:      $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/Poly1D_m.o:     $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/H2_m.o:         $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/Buck_m.o:       $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/Sigmoid_m.o:    $(OBJ_DIR)/Empty_m.o

$(OBJ_DIR)/Template_m.o:          $(OBJ_DIR)/Empty_m.o $(OBJ_DIR)/Morse_m.o
$(OBJ_DIR)/Test_m.o:              $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/HenonHeiles_m.o:       $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/Tully_m.o:             $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/OneDSOC_1S1T_m.o:      $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/OneDSOC_2S1T_m.o:      $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/TwoD_m.o:              $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/TwoD_RJDI2014_m.o:     $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/TwoD_Valahu2022_m.o:   $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/PSB3_m.o:              $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/Retinal_JPCB2000_m.o:  $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/HONO_m.o:              $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/HNO3_m.o:              $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/NO3_m.o:               $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/HOO_DMBE_m.o:          $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/H3_m.o:                $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/HCN_Murrell_m.o:       $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/ClH2p_m.o:             $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/ClH2p_Botschwina_m.o:  $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/CH5_m.o:               $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/PH4_m.o:               $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/HNNHp_m.o:             $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/H2SiN_m.o:             $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/H2NSi_m.o:             $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/LinearHBond_m.o:       $(OBJ_DIR)/Empty_m.o $(OBJ_DIR)/Morse_m.o $(OBJ_DIR)/Buck_m.o
$(OBJ_DIR)/TwoD_MullerBrown_m.o:  $(OBJ_DIR)/Empty_m.o
$(OBJ_DIR)/Phenol_m.o:            $(OBJ_DIR)/Empty_m.o $(OBJ_DIR)/Morse_m.o $(OBJ_DIR)/Sigmoid_m.o


$(OBJ_DIR)/Model_m.o:      $(QMLOBJ) $(OBJ_DIR)/Basis_m.o \
                           $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/FiniteDiff_m.o

$(OBJ_DIR)/Opt_m.o:        $(OBJ_DIR)/Model_m.o
$(OBJ_DIR)/IRC_m.o:        $(OBJ_DIR)/Opt_m.o $(OBJ_DIR)/Model_m.o
$(OBJ_DIR)/MakeHinact_m.o: $(OBJ_DIR)/Model_m.o

$(OBJ_DIR)/Model_driver.o: $(OBJ_DIR)/Model_m.o $(OBJ_DIR)/Opt_m.o $(OBJ_DIR)/IRC_m.o $(OBJ_DIR)/MakeHinact_m.o

$(OBJ_DIR)/TEST_model.o:   $(OBJ_DIR)/Model_m.o $(OBJ_DIR)/Opt_m.o $(OBJ_DIR)/IRC_m.o $(OBJ_DIR)/MakeHinact_m.o
$(OBJ_DIR)/TEST_driver.o:  $(OBJ_DIR)/Model_driver.o
#
############################################################################



#=================================================================================
#=================================================================================
# ifort compillation v17 v18 with mkl
#=================================================================================
ifeq ($(FFC),ifort)

  # opt management
  ifeq ($(OOPT),1)
      #F90FLAGS = -O -parallel -g -traceback
      FFLAGS = -O  -g -traceback
  else
      FFLAGS = -O0 -check all -g -traceback
  endif

  # where to store the modules
  FFLAGS +=-module $(MOD_DIR)

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -qopenmp
  endif

  # some cpreprocessing
  FFLAGS += -cpp $(CPPSHELL_QML)

  # where to look the .mod files
  FFLAGS += -I$(QDMOD_DIR) -I$(ADMOD_DIR)

  FLIB    = $(EXTLib)
  ifeq ($(LLAPACK),1)
    FLIB += -mkl -lpthread
  else
    FLIB += -lpthread
  endif

  FC_VER = $(shell $(F90) --version | head -1 )

endif
#=================================================================================
#=================================================================================