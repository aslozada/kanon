PROG =	kanon

SRCS =	kanon.f90 kanon_commands_mod.f90 kanon_constants_mod.f90 \
	kanon_delaunay_mod.f90 kanon_eigen_mod.f90 kanon_error_mod.f90 \
	kanon_initialize_mod.f90 kanon_kinds_mod.f90 kanon_lexical_mod.f90 \
	kanon_molecule_mod.f90 kanon_tools_mod.f90 kanon_vector_mod.f90

OBJS =	kanon.o kanon_commands_mod.o kanon_constants_mod.o \
	kanon_delaunay_mod.o kanon_eigen_mod.o kanon_error_mod.o \
	kanon_initialize_mod.o kanon_kinds_mod.o kanon_lexical_mod.o \
	kanon_molecule_mod.o kanon_tools_mod.o kanon_vector_mod.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = gfortran
F90FLAGS = -Wall -Wno-tabs -g -O0 -fopenmp -ffpe-trap=zero -fbacktrace -fcheck=all -fbounds-check 
LDFLAGS = -Wall -Wno-tabs -g -O0 -fopenmp -ffpe-trap=zero -fbacktrace -fcheck=all -fbounds-check 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

kanon.o: kanon_delaunay_mod.o kanon_error_mod.o kanon_initialize_mod.o \
	kanon_kinds_mod.o kanon_molecule_mod.o kanon_tools_mod.o
kanon_commands_mod.o: kanon_error_mod.o kanon_lexical_mod.o
kanon_constants_mod.o: kanon_kinds_mod.o
kanon_delaunay_mod.o: kanon_error_mod.o kanon_kinds_mod.o \
	kanon_molecule_mod.o kanon_tools_mod.o
kanon_eigen_mod.o: kanon_error_mod.o kanon_kinds_mod.o
kanon_initialize_mod.o: kanon_commands_mod.o kanon_error_mod.o \
	kanon_lexical_mod.o kanon_tools_mod.o
kanon_molecule_mod.o: kanon_constants_mod.o kanon_eigen_mod.o \
	kanon_error_mod.o kanon_initialize_mod.o kanon_kinds_mod.o \
	kanon_lexical_mod.o kanon_tools_mod.o kanon_vector_mod.o
kanon_tools_mod.o: kanon_error_mod.o kanon_kinds_mod.o kanon_lexical_mod.o
kanon_vector_mod.o: kanon_error_mod.o kanon_kinds_mod.o
