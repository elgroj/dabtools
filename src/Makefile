CC=gcc
PROFILEFLAGS=-pg

OPTARGS=-ffast-math -mfpmath=sse
# -ffast-math: -fno-math-errno, -funsafe-math-optimizations, -ffinite-math-only, -fno-rounding-math,
#            -fno-signaling-nans and -fcx-limited-range
# -funsafe-math-optimizations: -fno-signed-zeros, -fno-trapping-math, -fassociative-math and
#            -freciprocal-math
# -fno-rounding-math: on by default; allows constant folding
# -mfpmath=sse
# -msse -msse2 -msse3 : enable instruction sets


CFLAGS=-O3 -W -Wall -g $(OPTARGS) $(PROFILEFLAGS)
LDFLAGS=-lm -lfftw3 -lpthread -g $(PROFILEFLAGS)

CFLAGS+=$(shell pkg-config --cflags librtlsdr)
LDFLAGS+=$(shell pkg-config --libs librtlsdr)

# Uncomment the following line to use the fast, but x86-only Spiral Viterbi function.
ENABLE_SPIRAL_VITERBI=yes

ifdef ENABLE_SPIRAL_VITERBI
  CFLAGS+=-DENABLE_SPIRAL_VITERBI -msse2 -msse3 -std=gnu99 
  VITERBI_OBJS=viterbi_spiral.o viterbi_spiral_sse16.o
else
  VITERBI_OBJS=viterbi.o
endif


OBJS=wf_sync.o wf_prstables.o wf_maths.o fic.o misc.o dab_tables.o input_wf.o depuncture.o input_sdr.o sdr_fifo.o sdr_sync.o dab.o $(VITERBI_OBJS)

all: dab2eti eti2mpa

dab2eti: $(OBJS) dab2eti.o
	$(CC) -o dab2eti dab2eti.o $(OBJS) $(LDFLAGS)

eti2mpa: eti2mpa.o
	$(CC) -o eti2mpa eti2mpa.o

clean:
	rm -f dab2eti eti2mpa *.o *~
