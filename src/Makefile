# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = brightness climatology formod interpolate kernel limb nadir planck raytrace retrieval tab2bin

# Library directories...
LIBDIR = -L /usr/lib64 

# Include directories...
INCDIR = -I /usr/local/include

# OpenMP...
OMP = 1

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler with/without openMP and without MPI...
# -----------------------------------------------------------------------------

# Compiler...
CC = gcc

# CFLAGS...
CFLAGS = $(INCDIR) -pedantic-errors -Wall -Wextra -Werror -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Wno-long-long -Winline -fno-common -DHAVE_INLINE

# Flags for OpenMP...
ifdef OMP
CFLAGS += -fopenmp
endif

# -----------------------------------------------------------------------------
# Set flags for open MPI C wrapper compiler with/without openMP...
# -----------------------------------------------------------------------------
# # Compiler...
# CC = mpicc 

# # CFLAGS...
# CFLAGS = $(INCDIR) -pedantic-errors -Wall -Wextra -Werror -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Wno-long-long -Winline -fno-common -DHAVE_INLINE

# # Flags for OpenMP...
# ifdef OMP
# CFLAGS += -fopenmp
# endif
# -----------------------------------------------------------------------------

CFLAGS += -D WORK_QUEUE


# Flags for profiling...
ifdef PROF
CFLAGS += -O2 -g -pg
else
CFLAGS += -O3
endif

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC)
	rm -f *~

$(EXC): %: %.c  atmosphere.o continua.o control.o forwardmodel.o lineofsight.o misc.o retrievalmodel.o workqueue.o scatter.o Makefile
	$(CC) $(CFLAGS) -o $@ $< atmosphere.o continua.o control.o forwardmodel.o lineofsight.o misc.o retrievalmodel.o workqueue.o scatter.o  $(LDFLAGS)

atmosphere.o: atmosphere.c 
	$(CC) $(CFLAGS) -c atmosphere.c

continua.o: continua.c 
	$(CC) $(CFLAGS) -c continua.c

control.o: control.c
	$(CC) $(CFLAGS) -c control.c

forwardmodel.o: forwardmodel.c 
	$(CC) $(CFLAGS) -c forwardmodel.c

lineofsight.o: lineofsight.c 
	$(CC) $(CFLAGS) -c lineofsight.c

misc.o: misc.c 
	$(CC) $(CFLAGS) -c misc.c

retrievalmodel.o: retrievalmodel.c 
	$(CC) $(CFLAGS) -c retrievalmodel.c

scatter.o: scatter.c 
	$(CC) $(CFLAGS) -c scatter.c

workqueue.o: workqueue.c 
	$(CC) $(CFLAGS) -c workqueue.c



bak:
	tar czf ../bak/jurassic_`date +"%y%m%d%H%M"`.tgz Makefile *.c *.h

clean:
	rm -f $(EXC) *.o *~
