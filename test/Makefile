include ../config/blasxprs_config.mak
include ../config/boost_config.mak
include ../config/armadillo_config.mak
include ../config/blaze_config.mak
include ../config/eigen_config.mak
include ../config/openblas/openblas_config.mak
include ../config/atlas/atlas_config.mak
include ../config/netlib/netlib_config.mak
include ../config/intel_mkl/intel_mkl_config.mak
include ../config/amd_acml/amd_acml_config.mak

WHAT=mtrx

DEFINES=-DBLASXPRS_USE_$(BLAS_INTERFACE) -DBLASXPRS_MTRX_$(MTRX_LAYOUT)

# OPENBLAS
ifeq ($(WHICH_BLAS), openblas) 
  BLAS_LIB=$(OPENBLAS_LIB)
  BLAS_LDLIBS = -lopenblas -lpthread 
  BLAS_INCLUDE=$(OPENBLAS_INCLUDE)
  DEFINES += -DBLASXPRS_WITH_OPENBLAS
endif

# NETLIB
ifeq ($(WHICH_BLAS), netlib)
  BLAS_LIB=$(NETLIB_BLAS_LIB)
  ifeq ($(BLAS_INTERFACE), CBLAS)
    BLAS_LDLIBS= -lcblas -lblas
  endif
  ifeq ($(BLAS_INTERFACE), F77BLAS)
    BLAS_LDLIBS= -lblas
  endif
  BLAS_INCLUDE=$(NETLIB_BLAS_INCLUDE)
  DEFINES += -DBLASXPRS_WITH_NETLIB_BLAS
endif

# ATLAS
ifeq ($(WHICH_BLAS), atlas)
  BLAS_LIB=$(ATLAS_LIB)
  BLAS_LDLIBS= -ltatlas
  BLAS_INCLUDE=$(ATLAS_INCLUDE)
  DEFINES += -DBLASXPRS_WITH_ATLAS
endif

# MKL
ifeq ($(WHICH_BLAS), mkl)
  BLAS_LIB=$(INTEL_MKL_LIB)
  BLAS_LDLIBS= -lmkl_rt
  BLAS_INCLUDE=$(INTEL_MKL_INCLUDE)
  DEFINES += -DBLASXPRS_WITH_INTEL_MKL
endif

# ACML
ifeq ($(WHICH_BLAS), acml)
  BLAS_LIB=$(AMD_ACML_LIB)
  BLAS_LDLIBS= -lacml_mp
  BLAS_INCLUDE=$(AMD_ACML_INCLUDE)
  DEFINES += -DBLASXPRS_WITH_AMD_ACML
endif

CXX=c++
CXXFLAGS=-std=c++14 -m64 -Wall -Wextra -pedantic -ggdb3 -fopenmp
INCLUDE=-I../include -I${BLAS_INCLUDE} -I${BOOST_INCLUDE}
LDFLAGS=-L${BLAS_LIB} ${BLAS_LDLIBS}
RUNPATH=${BLAS_LIB}

ifeq ($(WITH_ARMADILLO_UTIL), true)
  INCLUDE += -I${ARMADILLO_INCLUDE}
  LDFLAGS += -L${ARMADILLO_LIB} -larmadillo
  DEFINES += -DBLASXPRS_WITH_ARMADILLO_UTIL
  ifeq ($(TEST_AGAINST_ARMADILLO), true)
    DEFINES += -DBLASXPRS_TEST_AGAINST_ARMADILLO 
	RUNPATH := $(RUNPATH):${ARMADILLO_LIB}
  endif
endif

ifeq ($(CXX), pgc++)
  CXXFLAGS = -std=c++11 -m64 -pedantic 
endif

ifeq ($(CXX), sunCC)
  CXXFLAGS = -std=c++14 -m64 -Wall -Wextra -pedantic -fopenmp 
endif

$(WHAT) : $(WHAT).cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -o $@ $< $(INCLUDE) $(LDFLAGS) -Wl,-rpath=$(RUNPATH)
