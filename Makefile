PROGRAM = blast

CXX_DBG = c++
CXX_CPU = c++
CXX_OMP = c++
CXX_GPU = nvcc

FLAG_ANY = -Wall -std=c++17 -MMD -MP -Ivapor/include
FLAG_DBG = -O0 -g
FLAG_CPU = -Ofast
FLAG_OMP = -Ofast -Xpreprocessor -fopenmp
FLAG_GPU = -O3 -x cu --extended-lambda

LIBS_ANY = -lhdf5
LIBS_DBG =
LIBS_CPU =
LIBS_OMP = -lomp
LIBS_GPU =

DEP_DBG = build/$(PROGRAM)_dbg.d
DEP_CPU = build/$(PROGRAM)_cpu.d
DEP_OMP = build/$(PROGRAM)_omp.d
DEP_GPU = build/$(PROGRAM)_gpu.d
DEP = $(DEP_DBG) $(DEP_CPU) $(DEP_OMP) $(DEP_GPU)

OBJ_DBG = build/$(PROGRAM)_dbg.o
OBJ_CPU = build/$(PROGRAM)_cpu.o
OBJ_OMP = build/$(PROGRAM)_omp.o
OBJ_GPU = build/$(PROGRAM)_gpu.o
OBJ = $(OBJ_DBG) $(OBJ_CPU) $(OBJ_OMP) $(OBJ_GPU)

EXE_DBG = bin/$(PROGRAM)_dbg
EXE_CPU = bin/$(PROGRAM)_cpu
EXE_OMP = bin/$(PROGRAM)_omp
EXE_GPU = bin/$(PROGRAM)_gpu
EXE = $(EXE_DBG) $(EXE_CPU) $(EXE_OMP) $(EXE_GPU)

-include Makefile.in

dbg: $(EXE_DBG)
cpu: $(EXE_CPU)
omp: $(EXE_OMP)
gpu: $(EXE_GPU)

$(EXE_DBG): $(OBJ_DBG) bin
	$(CXX_DBG) -o $@ $< $(LIBS_ANY) $(LIBS_DBG)
$(EXE_CPU): $(OBJ_CPU) bin
	$(CXX_CPU) -o $@ $< $(LIBS_ANY) $(LIBS_CPU)
$(EXE_OMP): $(OBJ_OMP) bin
	$(CXX_OMP) -o $@ $< $(LIBS_ANY) $(LIBS_OMP)
$(EXE_GPU): $(OBJ_GPU) bin
	$(CXX_GPU) -o $@ $< $(LIBS_ANY) $(LIBS_GPU)

$(OBJ_DBG): src/$(PROGRAM).cpp build
	$(CXX_DBG) -o $@ $< $(FLAG_ANY) $(FLAG_DBG) -c
$(OBJ_CPU): src/$(PROGRAM).cpp build
	$(CXX_CPU) -o $@ $< $(FLAG_ANY) $(FLAG_CPU) -c
$(OBJ_OMP): src/$(PROGRAM).cpp build
	$(CXX_OMP) -o $@ $< $(FLAG_ANY) $(FLAG_OMP) -c
$(OBJ_GPU): src/$(PROGRAM).cpp build
	$(CXX_GPU) -o $@ $< $(FLAG_ANY) $(FLAG_GPU) -c

all: $(EXE)

build:
	mkdir -p $@
bin:
	mkdir -p $@
clean:
	$(RM) $(DEP) $(OBJ) $(EXE)

-include $(DEP)
