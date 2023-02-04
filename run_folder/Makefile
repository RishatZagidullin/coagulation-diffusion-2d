GCC = g++

LIB  = -lgomp -lm -ldl

GCC_COPT = -fopenmp

VPATH = ../:../utils:../solvers:../ppm:../geometry

SOURCES_SOL = diffusion_reg.cpp vector3d.cpp coagulation.cpp ppm.cpp
SOURCES_MAIN = main.cpp

OBJECTS_SOL = $(SOURCES_SOL:%.cpp=$(OBJDIR)/%.o)
OBJECTS_MAIN = $(SOURCES_MAIN:%.cpp=$(OBJDIR)/%.o)

OBJDIR = build

all: tet.exe

tet.exe: $(OBJECTS_MAIN) $(OBJECTS_SOL)
	$(GCC) $^ -o $@ $(LIB)

$(OBJDIR)/%.o: %.cpp
	$(GCC) -c $(GCC_COPT) $< -o $@

$(OBJDIR)/%.o: %.c
	$(GCC) -c $(GCC_COPT) $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o *.exe
