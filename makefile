CC = g++
STD = --std=c++11
LGSL =

FLAGS = $(STD)
BINDIR = 
LDD = 

# OS-specific flags:
OSXFLAGS = -L/opt/local/lib
LINUXFLAGS = -O3 -Wall
OTHERFLAGS = -O2 -g -funroll-loops
LINUXLDD = -pthread
# Detect OS:
OS = $(shell uname)

# Set compile flags
ifeq ($(OS), Darwin)
	FLAGS = $(STD) -pg $(OSXFLAGS)
	LGSL = -lgsl
	LDD = $(OSXFLAGS)
endif
ifeq ($(OS), Linux)
	FLAGS = $(STD) $(LINUXFLAGS)
	LGSL = -lgsl -lgslcblas -lm
	LDD = $(LINUXLDD)
endif

# Set install directory
ifeq ($(OS), Darwin)
	BINDIR = $(HOME)/.local/bin/
endif
ifeq ($(OS), Linux)
	BINDIR = $(HOME)/.local/bin/
endif

# Recipes:
   
hopscotch : hopscotch.cpp params.h
	$(CC) hopscotch.cpp $(FLAGS) $(LGSL) -o hopscotch $(LDD)

install : hopscotch
	cp exec/hopscotch $(BINDIR)
 
