SHELL = /bin/bash

# Decide on compiler and flags depending on host
ifeq ($(HOSTNAME),kustaanheimo)
CC = g++
CFLAGS = -march=pentium4m
# For cygwin/winXP combo
GLLIBS = -lglut32 -lglu32 -lopengl32
else
CC = g++
CFLAGS =
# Seems to work on Debian/Mesa
GLLIBS = -lglut
endif

# Possible #definable options:
# ROTATING_COLORS 
# -- Instead of single color intensity variation, use RGB cycling
DEFINED_OPTS = -DROTATING_COLORS
ALL_CFLAGS = -Wall -O6 `sdl-config --cflags` $(DEFINED_OPTS) $(CFLAGS)

# Libraries
LIBS = -lm $(GLLIBS) `sdl-config --libs` 

# Binaries
halodump_EXE = visualize_halodump 
halos_EXE = visualize_halos
EXES = $(halodump_EXE) $(halos_EXE)

# Object files
halodump_OBJS = visualize_halodump.o common.o
halos_OBJS = visualize_halos.o common.o
OBJS = $(halos_OBJS) $(halodump_OBJS)

# pattern rule to compile object files from C files
# might not work with make programs other than GNU make
%.o : %.cpp Makefile
	$(CC) $(ALL_CFLAGS) -c $< -o $@

all: $(EXES)

$(halodump_EXE): $(halodump_OBJS) Makefile
	$(CC) $(ALL_CFLAGS) $(halodump_OBJS) -o $(halodump_EXE) $(LIBS)

$(halos_EXE): $(halos_OBJS) Makefile
	$(CC) $(ALL_CFLAGS) $(halos_OBJS) -o $(halos_EXE) $(LIBS)

# Make sure that we also get possible .exe suffixed cygwin executables
.PHONY : clean
clean:
	rm -f $(OBJS) $(EXES) $(EXES:=.exe)
