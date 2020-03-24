SOURCES := covid19.c

OBJECTS := $(SOURCES:%.c=%.o)

PROGRAM = covid19

CC := gcc
LD := $(CC)

#OMP = -fopenmp
OPT = -O3 -march=native -mtune=native
#DEBUG = -O0 -g
CPPFLAGS = 
CFLAGS = $(OPT) $(DEBUG) $(OMP)
LIBS := -lm

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS) $(LIBS)

clean:
	$(RM) $(PROGRAM) $(OBJECTS)

# vim:sw=8:sts=8:noexpandtab:
