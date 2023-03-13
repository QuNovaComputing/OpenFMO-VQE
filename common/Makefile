#INC_DIR = -I.. -I../basis -I../integ

OBJS = ofmo-string.o ofmo-compress.o ofmo-input.o ofmo-data.o \
       ofmo-init.o ofmo-misc.o ofmo-prof.o \
       ofmo-debug.o ofmo-fragment-basis.o ofmo-fragment.o \
       log_mpe.o

CFLAGS = ${CFLAGS0}
#CFLAGS = ${CFLAGS}
#CFLAGS += ${INC_DIR}

all: ${OBJS}

clean:
	-rm -f ${OBJS}

