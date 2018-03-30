#INC_DIR = -I.. -I../common -I../scf -I../integ -I../matope -I../basis

OBJS = ofmo-worker-main.o ofmo-mserv-cont.o ofmo-init-dens.o \
       ofmo-inter-frag.o ofmo-monomer-data.o ofmo-projection.o \
       ofmo-calc-frag.o ofmo-counter.o

#CFLAGS += ${INC_DIR}

all: ${OBJS}

clean:
	-rm -f ${OBJS}

