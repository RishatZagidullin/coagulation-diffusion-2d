GCC = g++

LIB  = -lgomp -lm -ldl -lfftw3 -lblas -llapacke

GCC_COPT = -fopenmp -DFFTW -O3 #-g

VPATH = ../:../utils:../solvers:../ppm:../geometry:../solvers/coagulation

SOURCES_SOL = space2d_reg.cpp vector3d.cpp coagulation.cpp ppm.cpp blas_double.cpp blas_int.cpp cross.cpp matrix.cpp maxvol.cpp parallel_cross_omp.cpp tensor.cpp tensor_train.cpp convolution.cpp
SOURCES_MAIN = main.cpp

OBJECTS_SOL = $(SOURCES_SOL:%.cpp=$(OBJDIR)/%.o)
OBJECTS_MAIN = $(SOURCES_MAIN:%.cpp=$(OBJDIR)/%.o)

OBJDIR = build

all: tet.exe

tet.exe: $(OBJECTS_MAIN) $(OBJECTS_SOL)
	$(GCC) $^ -o $@ $(LIB) $(GCC_COPT)

$(OBJDIR)/%.o: %.cpp
	$(GCC) -c $(GCC_COPT) $< -o $@

$(OBJDIR)/%.o: %.c
	$(GCC) -c $(GCC_COPT) $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o *.exe
