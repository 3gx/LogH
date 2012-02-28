CXX = g++
CC  = gcc
LD  = g++

.SUFFIXES: .o .f90 .cpp .cu .s

OFLAGS = -O4 -g -m64 -Wall -msse3 -funroll-all-loops  
# OFLAGS += -ffast-math
# OFLAGS = -O0 -g -m64 -Wall -msse3 
# OFLAGS += -fopenmp
# OFLAGS = -O4 -m64 -Wall -msse3 -funroll-all-loops
CFLAGS = $(OFLAGS) -I$(CUDA_TK)/include  -D__$(MACOSX)__  -I/opt/local/include
# CFLAGS += -D_CH2SIG_
CXXFLAGS = $(CFLAGS) 


OFLAGS += -D_FPESIG_ENABLE_ 

PROG = nbody

OBJS = nbody.o

# LIBS = -L/opt/local/lib -lCGAL  -lgmp
# LIBS += -lefence


all: $(PROG)

$(PROG): $(OBJS)
	$(LD) $(CXXFLAGS) -o $@ $(OBJS)




.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

.cpp.s: 
	$(CXX) $(CXXFLAGS) -S $< -o $@

clean:
	/bin/rm -rf *.o version.h
	/bin/rm -rf $(PROG) lib$(VLIB).a

$(OBJS) : LogH.h LogH+TTL.h particle.h vector3.h kepler.h 
