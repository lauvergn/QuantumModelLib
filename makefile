#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx), ifx, ifort
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
#=================================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)


QML_ver  := $(shell awk '/QML/ {print $$3}' version-QML)

# Extension for the object directory and the library
ext_obj=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)

# library name
QMLIBA=libQMLib$(ext_obj).a
#=================================================================================
MAIN_path := $(shell pwd)

OBJ_DIR=OBJ/obj$(ext_obj)
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))

MOD_DIR=$(OBJ_DIR)
SRC_DIR=SRC
MAIN_DIR=APP
TESTS_DIR=Tests

CPPSHELL_QML = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
               -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
			   -D__QMLPATH="'$(MAIN_path)'" \
               -D__QML_VER='"$(QML_ver)"'
#=================================================================================
# External Libraries directory
#
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(MAIN_path)/Ext_Lib
endif
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA    = $(QD_DIR)/libQD$(ext_obj).a

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(ext_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(ext_obj).a

EXTMOD_DIR = $(QDMOD_DIR) $(ADMOD_DIR)

EXTMod     = -I$(QDMOD_DIR) -I$(ADMOD_DIR)
EXTLib     = $(ADLIBA) $(QDLIBA)
#===============================================================================
#=================================================================================
 #
#=================================================================================
# To deal with external compilers.mk file
CompilersDIR = $(MAIN_path)
ifeq ($(CompilersDIR),)
  include compilers.mk
else
  include $(CompilersDIR)/compilers.mk
endif
FFLAGS += $(CPPSHELL_QML)
#=================================================================================
#=================================================================================
#=================================================================================
$(info ***********************************************************************)
$(info ***********OS:           $(OS))
$(info ***********COMPILER:     $(FFC))
$(info ***********COMPILER_VER: $(FC_VER))
$(info ***********OPTIMIZATION: $(OOPT))
$(info ***********OpenMP:       $(OOMP))
$(info ***********INT:          $(INT))
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
#$(info ***********MAINSx:      $(MAINSx))

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
	cd Tests ; ./run_test_QML $(FFC) $(OOPT) $(OOMP) $(LLAPACK) $(INT) 1
#===============================================
#============= all: lib, tests ...  ============
#===============================================
.PHONY: all
all: $(QMLIBA) $(MAINSx) $(TESTS).x
#===============================================
#============= Main executable and tests  ======
#=============================================== libQMLibFull$(ext_obj).a
TEST_VibAdia.x: $(OBJ_DIR)/TEST_VibAdia.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_VibAdia.x  $(OBJ_DIR)/TEST_VibAdia.o libQMLibFull$(ext_obj).a $(FLIB)
TEST_OOP.x: $(OBJ_DIR)/TEST_OOP.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_OOP.x  $(OBJ_DIR)/TEST_OOP.o libQMLibFull$(ext_obj).a $(FLIB)
TEST_grid.x: $(OBJ_DIR)/TEST_grid.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_grid.x  $(OBJ_DIR)/TEST_grid.o libQMLibFull$(ext_obj).a $(FLIB)
TEST_OMPloop.x: $(OBJ_DIR)/TEST_OMPloop.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_OMPloop.x  $(OBJ_DIR)/TEST_OMPloop.o libQMLibFull$(ext_obj).a $(FLIB)
TEST_driver.x: $(OBJ_DIR)/TEST_driver.o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o TEST_driver.x  $(OBJ_DIR)/TEST_driver.o libQMLibFull$(ext_obj).a $(FLIB)
#
$(TESTS).x: $(OBJ_DIR)/$(TESTS).o $(QMLIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o $(TESTS).x  $(OBJ_DIR)/$(TESTS).o libQMLibFull$(ext_obj).a $(FLIB)
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
	rm -f RES_files DAT_files
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
#
#===============================================
#=== external libraries ========================
# AD_dnSVM + QDUTIL Lib
#===============================================
#
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib dev
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM dev
#
$(QDLIBA): getlib
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QDLIBA)
$(ADLIBA): getlib
	cd $(ExtLibDIR)/AD_dnSVM ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(ADLIBA) || (echo $(ADLIBA) "does not exist" ; exit 1)
	@echo "  done " $(ADLIBA)
#
#===============================================
#=== Add links to directories for fpm ==========
#===============================================
#
.PHONY: fpmlink
fpmlink: getlib
	ln -s Tests/RES_files RES_files
	ln -s Tests/DAT_files DAT_files

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
