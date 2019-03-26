#
#  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
#  DATE WRITTEN:     Feb 2017
#  LAST MODIFIED:    Feb 2017
#  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
#  DESCRIPTION:      File to compile adaptive_picard_chebyshev code
#

IDIR = ../include
CC = gcc -std=c99 -g
# CC = gcc
CFLAGS = -I$(IDIR)

ODIR = obj
LDIR = ../lib

LIBS = -lm

_DEPS1 = adaptive_picard_chebyshev.h polydegree_segments.h picard_chebyshev_propagator.h prepare_propagator.h \
picard_iteration.h clenshaw_curtis_ivpII.h rv2elm.h c_functions.h FandG.h perigee_approx.h EGM2008.h \
lsq_chebyshev_fit.h chebyshev.h radial_gravity.h matrix_loader.h eci2ecef.h ecef2eci.h perturbed_gravity.h \
picard_error_feedback.h interpolate.h reosc_perigee.h
DEPS1 = $(patsubst %,$(IDIR)/%,$(_DEPS1))

_DEPS2 = clenshaw_curtis_ivpII.h lsq_chebyshev_fit.h chebyshev.h c_functions.h
DEPS2 = $(patsubst %,$(IDIR)/%,$(_DEPS2))

_OBJ1 = test.o adaptive_picard_chebyshev.o polydegree_segments.o picard_chebyshev_propagator.o prepare_propagator.o \
picard_iteration.o clenshaw_curtis_ivpII.o rv2elm.o c_functions.o FandG.o perigee_approx.o EGM2008.o \
lsq_chebyshev_fit.o chebyshev.o radial_gravity.o matrix_loader.o eci2ecef.o ecef2eci.o perturbed_gravity.o \
picard_error_feedback.o interpolate.o reosc_perigee.o
OBJ1 = $(patsubst %,$(ODIR)/%,$(_OBJ1))

_OBJ2 = matrix_builder.o clenshaw_curtis_ivpII.o lsq_chebyshev_fit.o chebyshev.o c_functions.o
OBJ2 = $(patsubst %,$(ODIR)/%,$(_OBJ2))

$(ODIR)/%.o: %.c $(DEPS1) $(DEPS2)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ1)
	gcc -g -o $@ $^ $(CFLAGS) $(LIBS)

matrix_builder: $(OBJ2)
	gcc -g -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
