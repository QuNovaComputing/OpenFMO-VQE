OBJS = ofmo-cutoff-core-dddd.o ofmo-cutoff-core-ssss-dpdp.o ofmo-cutoff.o \
    ofmo-ifc2c-core.o ofmo-ifc2c.o ofmo-ifc3c-os.o ofmo-ifc4c-os.o \
    ofmo-integ.o ofmo-oneint-core.o ofmo-oneint.o fmt.o \
    fmt-m.o fmt-method1.o fmt-method2.o fmt-method3.o \
    ofmo-twoint-buffer.o ofmo-twoint-core-dddd.o ofmo-twoint-core-dddp.o \
    ofmo-twoint-core-ddds.o ofmo-twoint-core-ddss-ddpp.o \
    ofmo-twoint-core-dpds-dpdp.o ofmo-twoint-core-dsds-dppp.o \
    ofmo-twoint-core-ssss-dspp.o ofmo-twoint-direct.o ofmo-twoint.o \
    ofmo-index.o ofmo-os-xxxx.o ofmo-oneint-gen.o \
    ofmo-rys-xxxx.o ofmo-root.o \
    ofmo-rys-psss.o ofmo-rys-psps.o ofmo-rys-ppss.o ofmo-rys-ppps.o \
    ofmo-rys-pppp.o ofmo-rys-dsss.o ofmo-rys-dsps.o ofmo-rys-dspp.o \
    ofmo-rys-dsds.o ofmo-rys-dpss.o ofmo-rys-dpps.o ofmo-rys-dppp.o \
    ofmo-rys-dpds.o ofmo-rys-dpdp.o ofmo-rys-ddss.o ofmo-rys-ddps.o \
    ofmo-rys-ddpp.o ofmo-rys-ddds.o ofmo-rys-dddp.o ofmo-rys-dddd.o \
    ofmo-os-psss.o ofmo-os-psps.o ofmo-os-ppss.o ofmo-os-ppps.o \
    ofmo-os-pppp.o ofmo-os-dsss.o ofmo-os-dsps.o ofmo-os-dspp.o \
    ofmo-os-dsds.o ofmo-os-dpss.o ofmo-os-dpps.o ofmo-os-dppp.o \
    ofmo-os-dpds.o ofmo-os-dpdp.o ofmo-os-ddss.o ofmo-os-ddps.o \
    ofmo-os-ddpp.o ofmo-os-ddds.o ofmo-os-dddp.o ofmo-os-dddd.o \
    ofmo-ifc4c-rys.o ofmo-ifc3c-rys.o ofmo-ifc4c.o


INC_DIR = -I..

CFLAGS += ${INC_DIR}
#CFLAGS2 += ${INC_DIR}
CFLAGS2 = $(patsubst -O3,-O2, $(CFLAGS))
CFLAGS1 = $(patsubst -O3,-O1, $(CFLAGS))

all :: $(OBJS)


OBJS1 = ofmo-os-xxxx.o ofmo-oneint-gen.o
$(OBJS1): CFLAGS := $(CFLAGS1)

OBJS2 = ofmo-twoint-core-dddd.o ofmo-rys-dddd.o ofmo-os-dddd.o
$(OBJS2): CFLAGS := $(CFLAGS2)

clean:
	-rm -f ${OBJS} *.f90 *.mod
