## config

CXXFLAGS=-std=c++17
#HPCOPTS=-msse4.2 -axavx
HPCOPTS=-march=haswell -mtune=haswell
CPPFLAGS=-O3 ${HPCOPTS} -DNDEBUG -DSLURM #-DSLURMLOCK
#CPPFLAGS= -g -O0 -DDEBUG
CPPFLAGS += -fpic -Wall -Wextra -Wpedantic -Wvector-operation-performance
CXXFLAGS += -Wnon-virtual-dtor -Woverloaded-virtual #-Weffc++
LDFLAGS=${HPCOPTS}
#LDFLAGS=-g
LDLIBS=-lstdc++ -lm

# Mosek
CPPFLAGS += -I${MSKHOME}/h -I$(abspath ${MSKHOME}/../../examples/c)
LDLIBS += -L${MSKHOME}/bin -lmosek64

# Gurobi
CPPFLAGS += -I${GUROBI_HOME}/include

ifneq (,$(findstring taurus.hrsk.tu-dresden.de,$(shell hostname)))
  LDLIBS += -L${HOME}/gurobi8/build
else
  LDLIBS += -L${GUROBI_HOME}/src/build
endif

LDLIBS +=-lgurobi_c++ -L${GUROBI_HOME}/lib -lgurobi80

# HDF5
CXX = h5c++ -shlib
CC = h5cc -shlib

## default targets

.PHONY:
all: afsnd afsnd_dinkelbach afsndEE aftinEE afsnd2 hk leakage PA_afsnd hkEE

.PHONY:
clean:
	rm -f *.o *.a
	rm -rf .d

# Mosek
vpath scopt-ext.c $(abspath ${MSKHOME}/../../examples/c)
scopt-ext.o: CPPFLAGS:=$(filter-out -Wall -DDEBUG,$(CPPFLAGS))

scopt-ext.o: CPPFLAGS+=-Wno-unused-result


# Dependencies
OBJ_SIT = util.o
OBJ_RR = $(OBJ_SIT) Mosek.o Gurobi.o scopt-ext.o
OBJ_SND = $(OBJ_RR) HDF.o util.o
OBJ_HK = $(OBJ_RR) HDF.o HK.o util.o
OBJ_HK2 = $(OBJ_RR) HDF.o HK2.o util.o
OBJ_PA = PA_Solver.o HDF.o util.o Gurobi.o

## rules

afsnd: $(OBJ_SND)

afsnd_dinkelbach: $(OBJ_SND)

afsndEE: $(OBJ_SND)

aftinEE: $(OBJ_SND)

afsnd2: $(OBJ_SND)

hk: $(OBJ_HK)

hkEE: $(OBJ_HK)

hk2: $(OBJ_HK2)

leakage: $(OBJ_SIT)

PA_afsnd: $(OBJ_PA)

tinEE_dm: util.o HDF.o

tinEE_dm: CPPFLAGS+=-Itin_dinkelbach

tinEE_dm: CPPFLAGS:=$(filter-out -Wpedantic, $(CPPFLAGS))

.PRECIOUS: bench_%.o
bench_%.o: bench.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DNUMUE=$(firstword $(subst _, ,$*)) -DNUMNC=$(word 2,$(subst _, ,$*)) -USLURM -c bench.cpp -o $@

bench_%: bench_%.o $(OBJ_SND)
	$(CC) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

%.a: %.o 
	$(AR) $(ARFLAGS) $@ $?



### automatic dependency tracking

SRCS = $(wildcard *.cpp)

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.c = $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o : %.c
%.o : %.c $(DEPDIR)/%.d
	$(COMPILE.c) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

%.o : %.cc
%.o : %.cc $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

%.o : %.cxx
%.o : %.cxx $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d


include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS))))

