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

ALL_CFLAGS = -Wall -O3 -s $(CFLAGS)
LIBS = -lm

corr_BIN = correlator
corr_SRC = correlator.c

compcorr_BIN = compcorr
compcorr_SRC = compcorr.c

fg_BIN = fieldgen
fg_SRC = fieldgen.c mersenne/mt19937ar.c
fg_OBJ = fieldgen.o mt19937ar.o

BIN = $(corr_BIN) $(compcorr_BIN) $(fg_BIN)
SRC = $(corr_SRC) $(compcorr_SRC) $(fg_SRC)

ifeq ($(HOSTNAME),kustaanheimo)
else
$(BIN) = $(BIN:.exe=)
endif

all: $(BIN) Makefile

$(corr_BIN): $(corr_SRC)
	$(CC) $(ALL_CFLAGS) $(LIBS) -o $(corr_BIN) $(corr_SRC)

$(compcorr_BIN): $(compcorr_SRC)
	$(CC) $(ALL_CFLAGS) $(LIBS) -o $(compcorr_BIN) $(compcorr_SRC)

$(fg_BIN): $(fg_SRC)
	$(CC) $(ALL_CFLAGS) -c $(fg_SRC)
	$(CC) $(ALL_CFLAGS) $(LIBS) -o $(fg_BIN) $(fg_OBJ)

# We have to make sure here that we get also the .exe suffixed executables
# that cygwin gcc seems to want to create. Stupid Windows.
.PHONY : clean
clean:
	rm -f $(BIN) $(BIN:=.exe) $(fg_OBJ)
