SOURCES := covid19.c interventions.c COV_rand.c prop_distribution.c distance.c locale.c infections.c

OBJECTS := $(SOURCES:%.c=%.o)

PROGRAM = covid19

CC := gcc
LD := $(CC)

OMP = -fopenmp
OPT = -O3 -march=native -mtune=native
#DEBUG = -O0 -g
CPPFLAGS = -DUSE_LOCALE_DISTANCE
CFLAGS = $(OPT) $(DEBUG) $(OMP)
LDFLAGS = $(OPT) $(DEBUG) $(OMP)
LIBS := -lm

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(OMP) $(OBJECTS) $(LIBS)

$(OBJECTS): common.h

clean:
	$(RM) $(PROGRAM) $(OBJECTS)

# vim:sw=8:sts=8:noexpandtab:
