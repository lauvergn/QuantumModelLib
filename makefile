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
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT = 4
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
ifeq ($(INT),)
  IINT      := 4
else
  IINT      := $(INT)
endif
#=================================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)


QML_ver  := $(shell awk '/QML/ {print $$3}' version-QML)
QML_path := $(shell pwd)

# Extension for the object directory and the library
ext_obj=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(IINT)

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
    FFLAGS0 = -O5 -g .... 
  else
    FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
    FFLAGS0 = -Og -g
    #FFLAGS = -O0 -fbounds-check -Wuninitialized
  endif

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -fopenmp
    FFLAGS0 += -fopenmp
  endif

  # integer kind management
  ifeq ($(IINT),8)
    FFLAGS += -fdefault-integer-8
    FFLAGS0 += -fdefault-integer-8
  endif

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
      IntLIB = -framework Accelerate
    else                   # Linux
      # linux libs
      IntLIB = -llapack -lblas
      #
      # linux libs with mkl and with openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
      # linux libs with mkl and without openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
    endif
  endif
  FLIB += $(IntLIB)

  FC_VER = $(shell $(FFC) --version | head -1 )

endif
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# ifort compillation v17 v18 with mkl
#=================================================================================
ifeq ($(FFC),$(filter $(FFC),ifort ifx))

  # opt management
  ifeq ($(OOPT),1)
      FFLAGS = -O  -g -traceback -heap-arrays
  else
      FFLAGS = -O0 -check all -g -traceback -heap-arrays
  endif

  # integer kind management
  ifeq ($(IINT),8)
    FFLAGS += -i8
  endif

  # where to store the modules
  FFLAGS +=-module $(MOD_DIR)

  # omp management
  ifeq ($(OOMP),1)
    ifeq ($(FFC),ifort)
      FFLAGS += -qopenmp -parallel
    else # ifx
      FFLAGS += -qopenmp
    endif
  endif
  # some cpreprocessing

  FFLAGS += -cpp $(CPPSHELL_QML)

  # where to look the .mod files
  FFLAGS += -I$(QDMOD_DIR) -I$(ADMOD_DIR)

  FLIB    = $(EXTLib)

  ifneq ($(LLAPACK),1)
    ifeq ($(FFC),ifort)
      IntLIB = -mkl -lpthread
    else # ifx
      IntLIB = -qmkl -lpthread
    endif
  else
    IntLIB = -lpthread
  endif

  FLIB += $(IntLIB)

  FC_VER = $(shell $(FFC) --version | head -1 )

endif
#===============================================================================
# nag compillation (nagfor)
#===============================================================================
ifeq ($(FFC),nagfor)

  # opt management
  ifeq ($(OOPT),1)
      FFLAGS = -O4 -o -compatible -kind=byte -Ounroll=4 -s
  else
    ifeq ($(OOMP),0)
      ifeq ($(LLAPACK),0)
          FFLAGS = -O0 -g -gline -kind=byte -C -C=alias -C=intovf -C=undefined
      else
          FFLAGS = -O0 -g -gline -kind=byte -C -C=alias -C=intovf
      endif
    else
          FFLAGS = -O0 -g        -kind=byte -C -C=alias -C=intovf
    endif
  endif

  # integer kind management
  ifeq ($(IINT),8)
    FFLAGS += -i8
  endif

 # where to store the .mod files
  FFLAGS +=-mdir $(MOD_DIR)

# where to look the .mod files
  FFLAGS +=-I $(MOD_DIR)

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -openmp
  endif

  # lapack management with cpreprocessing
  FFLAGS += -fpp -D__LAPACK="$(LLAPACK)"

  # where to look .mod files
  FFLAGS += -I$(QDMOD_DIR) -I$(ADMOD_DIR)

  FLIB    = $(QDLIBA)

  # lapact management (default with openmp), with cpreprocessing
  ifeq ($(LLAPACK),1)
    ifeq ($(OS),Darwin)    # OSX
      # OSX libs (included lapack+blas)
      FLIB += -framework Accelerate
    else                   # Linux
      # linux libs
      FLIB += -llapack -lblas
    endif
  endif

  FC_VER = $(shell $(FFC) -V 3>&1 1>&2 2>&3 | head -1 )

endif
#=================================================================================
#=================================================================================
$(info ***********************************************************************)
$(info ***********OS:           $(OS))
$(info ***********COMPILER:     $(FFC))
$(info ***********COMPILER_VER: $(FC_VER))
$(info ***********OPTIMIZATION: $(OOPT))
$(info ***********OpenMP:       $(OOMP))
$(info ***********INT:          $(IINT))
$(info ***********LAPACK:       $(LLAPACK))
$(info ***********FFLAGS:       $(FFLAGS))
$(info ***********FLIB:         $(FLIB))
$(info ***********EXTMOD_DIR:   $(EXTMOD_DIR))
$(info ***********ext_obj:      $(ext_obj))
$(info ***********************************************************************)

VPATH = $(TESTS_DIR) $(MAIN_DIR) $(SRC_DIR)  \
        $(SRC_DIR)/QMLLib $(SRC_DIR)/QML $(SRC_DIR)/Opt $(SRC_DIR)/AdiaChannels

MAINS= TEST_driver TEST_VibAdia TEST_grid TEST_OMPloop
MAINSx=$(addsuffix .x,$(MAINS))
$(info ***********MAINSx:      $(MAINSx))

TESTS=TEST_model

MAINSRCFILES=TEST_VibAdia.f90 TEST_grid.f90 TEST_OMPloop.f90 TEST_driver.f90
TESTSRCFILES=TEST_model.f90

# liste of source files in SRCFILES
include ./fortranlist.mk


OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))

#===============================================
#============= Tests ===========================
#===============================================
.PHONY: test TEST
test TEST: $(TESTS).x
	@echo "model (QML) compilation: OK"
#
.PHONY: ut UT
UT ut: $(TESTS).x
	@echo "model (QML) compilation: OK"
	cd Tests ; ./run_test_QML $(FFC) $(OOPT) $(OOMP) $(LLAPACK) $(IINT) 1
#	cd Tests ; ../$(TESTS).x  < input.dat > res 2>error.log

#	./$(TESTS).x |	grep "Number of error(s)"
#	@echo "  done Tests"
#===============================================
#============= all: lib, tests ...  ============
#===============================================
.PHONY: all
all: $(QMLIBA) $(MAINSx) $(TESTS).x
#===============================================
#============= Main executable and tests  ======
#=============================================== libQMLibFull$(ext_obj).a
TEST_VibAdia.x: $(OBJ_DIR)/TEST_VibAdia.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_VibAdia.x  $(OBJ_DIR)/TEST_VibAdia.o libQMLibFull$(ext_obj).a $(IntLIB)
TEST_OOP.x: $(OBJ_DIR)/TEST_OOP.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_OOP.x  $(OBJ_DIR)/TEST_OOP.o libQMLibFull$(ext_obj).a $(IntLIB)
TEST_grid.x: $(OBJ_DIR)/TEST_grid.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_grid.x  $(OBJ_DIR)/TEST_grid.o libQMLibFull$(ext_obj).a $(IntLIB)
TEST_OMPloop.x: $(OBJ_DIR)/TEST_OMPloop.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_OMPloop.x  $(OBJ_DIR)/TEST_OMPloop.o libQMLibFull$(ext_obj).a $(IntLIB)
TEST_driver.x: $(OBJ_DIR)/TEST_driver.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_driver.x  $(OBJ_DIR)/TEST_driver.o libQMLibFull$(ext_obj).a $(IntLIB)
#
$(TESTS).x: $(OBJ_DIR)/$(TESTS).o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o $(TESTS).x  $(OBJ_DIR)/$(TESTS).o libQMLibFull$(ext_obj).a $(IntLIB)
#===============================================
#============= Library: libQD.a  ===============
#===============================================
.PHONY: lib
lib: $(QMLIBA)

$(QMLIBA): $(OBJ)
	ar -cr $(QMLIBA) $(OBJ)
	ar -cr libQMLibFull$(ext_obj).a $(OBJ) $(ADMOD_DIR)/*.o $(QDMOD_DIR)/*.o
	@echo "  done Library: "$(QMLIBA)
#===============================================
#===============================================
#
#===============================================
#============= compilation =====================
#===============================================
$(OBJ_DIR)/%.o: %.f90 | $(EXTLib)
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
#===============================================
#================ cleaning =====================
.PHONY: clean cleanall
clean:
	cd $(TESTS_DIR) ; ./clean
	rm -f $(OBJ_DIR)/*/*.o $(OBJ_DIR)/*.o
	rm -f *.log grid* res*
	rm -f TEST*.x
	@echo "  done cleaning"

cleanall : clean
	cd $(ExtLibDIR) ; ./cleanlib
	rm -fr OBJ/obj* OBJ/*mod build
	rm -f lib*.a
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#===============================================
#================ make the readme.txt ==========
.PHONY: readme
readme:
	./scripts/extractReadMe
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := QuantumModelLib
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	$(ExtLibSAVEDIR)/makezip.sh $(BaseName)
	cd $(ExtLibSAVEDIR) ; ./cp_QML.sh
	@echo "  done zip"
#=== external libraries ========================
# AD_dnSVM + QDUTIL Lib
#===============================================
#
$(QDLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QD_DIR) || (cd $(ExtLibDIR) ; ./get_QDUtilLib.sh $(EXTLIB_TYPE))
	@test -d $(QD_DIR) || (echo $(QD_DIR) "does not exist" ; exit 1)
	cd $(QD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(IINT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(AD_DIR) || (cd $(ExtLibDIR) ; ./get_AD_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(AD_DIR) || (echo $(AD_DIR) "does not exist" ; exit 1)
	cd $(AD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(IINT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(AD_DIR) " in "$(BaseName)
#
#===============================================
#============= make dependencies =============
#===============================================
.PHONY: dep
dependencies.mk fortranlist.mk dep:
	./scripts/dependency.sh
##################################################################################
### dependencies
#
$(OBJ_DIR)/TEST_model.o:    $(QMLIBA) | $(EXTLib)
$(OBJ_DIR)/TEST_driver.o:   $(QMLIBA) | $(EXTLib)
$(OBJ_DIR)/TEST_VibAdia.o:  $(QMLIBA) | $(EXTLib)
$(OBJ_DIR)/TEST_grid.o:     $(QMLIBA) | $(EXTLib)
$(OBJ_DIR)/TEST_OMPloop.o:  $(QMLIBA) | $(EXTLib)

include ./dependencies.mk
#
############################################################################
