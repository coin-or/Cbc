OptLevel := -O1
OptLevel := -g

DRIVER ?= driver

# Possible choices for driver:
#   sample3
#   sample4

# The location of the customized Makefiles
CoinDir := ../../
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

IncDir :=
IncDir += ${CoinIncDir}
IncDir += ${ClpIncDir}
IncDir += ${zlibIncDir}
IncDir += ${bzlibIncDir}

LibDir :=
LibDir += ${CoinLibDir}
LibDir += ${ClpLibDir}
LibDir += ${zlibLibDir}
LibDir += ${bzlibLibDir}
LibDir += $(GlpkLibDir)

LibName :=
LibName += ${CbcLibName}
LibName += ${CglLibName}
LibName += ${OsiClpLibName}
LibName += ${OsiCbcLibName}
LibName += ${OsiLibName}
LibName += ${ClpLibName}
LibName += ${CoinLibName}
LibName += ${zlibLibName}
LibName += ${bzlibLibName}
LibName += $(GlpkLibName)
ifeq ($(VolDefine),COIN_HAS_VOL)
LibName += $(VolLibName)
endif

Define :=
Define += ${ClpDefine}
Define += ${CoinDefine}
Define += ${zlibDefine}
Define += ${bzlibDefine}
#lapack
IncDir +=  $(lapackIncDir)
LibDir +=  $(lapackLibDir)
LibName += $(lapackLibName)
Define +=  $(lapackDefine)

CXXFLAGS += $(OPTFLAG)
CXXFLAGS += -DCOIN_USE_CLP 
###############################################################################

space:= $(empty) $(empty)
OptVersion := $(subst $(space),_,$(OptLevel))

TARGETDIR := $(UNAME)$(OptVersion)
DEPDIR := dep

VPATH := . : include : Junk : ${TARGETDIR} : ${DEPDIR}

#########################################################################

CXXFLAGS += $(addprefix -I,${IncDir})
CXXFLAGS += $(addprefix -D,${Define})

LIBDIRS := ${LibDir}
LIBS    := ${LibName}

LDFLAGS := $(addprefix -L,$(LIBDIRS))
LDFLAGS += $(call ifprefix,$(SHLINKPREFIX),$(LIBDIRS))
LDFLAGS += $(patsubst lib%,-l%,$(basename $(LIBS)))

ifeq ($(findstring Linux,${UNAME}),Linux)
	CXXFLAGS +=	-DREADLINE
#	LDFLAGS += -lhistory -lreadline -ltermcap
endif
CXXFLAGS  += -I. 
#CXXFLAGS += -DCHECK_CUT_COUNTS
#LDFLAGS += -lefence
###############################################################################

TESTSRC := $(DRIVER).cpp
TESTSRC += CbcCompareUser.cpp
TESTSRC += CbcSolver2.cpp
#TESTSRC += CbcSolver3.cpp
TESTSRC += CbcBranchUser.cpp
TESTSRC += CbcBranchLink.cpp
TESTSRC += CbcHeuristicUser.cpp
TESTSRC += CbcHeuristicGreedy.cpp
TESTOBJ := $(addprefix $(TARGETDIR)/, $(TESTSRC:.cpp=.o))
TESTDEP := $(addprefix $(DEPDIR)/, $(TESTSRC:.cpp=.d))

###############################################################################
# This one must be right before the rules

include ${MakefileDir}/Makefile.rules

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean doc unitTest libdepend

default: testit

###############################################################################


testit : libdepend $(TESTOBJ)
	@rm -rf Junk
	@echo ""
	@echo Creating testit
	@rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(TESTOBJ) $(LDFLAGS) $(SYSLD) -lm
#	${CP} $@ ..

###############################################################################

-include $(TESTDEP)
