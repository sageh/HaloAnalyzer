SHELL = /bin/bash
# decide on compiler and flags depending on host
ifeq ($(HOSTNAME),kustaanheimo)
CC = gcc
CFLAGS = -march=pentium4m
else
ifeq ($(HOSTNAME),phobos2)
CC = gcc
CFLAGS = -march=i686
else
CC = gcc
CFLAGS =
endif
endif

# Possible debug and additional information output options
#DOPTS = -DDEBUG -DPROGRESS
DOPTS =

ALL_CFLAGS = -Wall -O3 -s $(CFLAGS) $(DOPTS)
LIBS = -lm

# Change the MPI libraries according to the system installed. The defaults
# are for MPICH2
MPILIBS = -lpmpich -lmpich

mpi_corr_BIN = mpi_correlator
mpi_corr_SRC = mpi_correlator.c

mpi_compcorr_BIN = mpi_compcorr
mpi_compcorr_SRC = mpi_compcorr.c

BIN = $(mpi_corr_BIN) $(mpi_compcorr_BIN)
SRC = $(mpi_corr_SRC) $(mpi_compcorr_SRC)

ifeq ($(HOSTNAME),kustaanheimo)
else
$(BIN) = $(BIN:.exe=)
endif

all: $(BIN) Makefile

$(mpi_corr_BIN): $(mpi_corr_SRC)
	$(CC) $(ALL_CFLAGS) -o $(mpi_corr_BIN) $(mpi_corr_SRC) $(LIBS) $(MPILIBS) 

$(mpi_compcorr_BIN): $(mpi_compcorr_SRC)
	$(CC) $(ALL_CFLAGS) -o $(mpi_compcorr_BIN) $(mpi_compcorr_SRC) $(LIBS) $(MPILIBS) 

# We have to make sure here that we get also the .exe suffixed executables
# that cygwin gcc seems to want to create. Stupid Windows.
.PHONY : clean
clean:
	rm -f $(BIN) $(BIN:=.exe) $(fg_OBJ)
