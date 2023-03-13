#TEST_RUN=YES
TEST_RUN=NO

OBJS_4 = cuda-drv.o cuda-fmt-drv.o cuda.o
OBJS_5 = cudalib.o cuda-fmt.o cuda-utils.o \
       cuda-twoint-core-ssss-dspp.o \
       cuda-twoint-core-dsds-dppp.o \
       cuda-twoint-direct.o \
       cuda-integ.o 
SRCS = cudalib.cu cuda-fmt.cu cuda-utils.cu \
       cuda-fmt-m.cu \
       cuda-root.cu \
       cuda-twoint-core-ssss-dspp.cu \
       cuda-os-dsss.cu cuda-os-dsps.cu cuda-os-dspp.cu cuda-os-dsds.cu \
       cuda-os-dpss.cu cuda-os-dpps.cu cuda-os-dppp.cu cuda-os-dpds.cu cuda-os-dpdp.cu \
       cuda-os-ddss.cu \
       cuda-rys-ppps.cu cuda-rys-dsps.cu \
       cuda-twoint-direct.cu \
       cuda-integ.cu \
       cuda-ifc4c.cu cuda-ifc4c-calc.cu cuda-ifc4c-os.cu
SRC_FMT_M = cuda-fmt-method3-8-8-12.cu \
       cuda-fmt-method2-8-10-12.cu \
       cuda-fmt-method2-8-8-12.cu \
       cuda-fmt-method1-8-10-12.cu 

INC_DIR = -I.. -I. -I../include
all :: $(OBJS_4)

cuda.cu: $(SRCS) cudalib.h cuda-fmt.h cuda-utils.h cuda-twoint-core.h cuda-twoint-direct.h cuda-twoint-core-os.h cuda-twoint-core-rys.h cuda-root.h cuda-fmt-m.h cuda-integ.h cuda-ifc4c.h cuda-ifc4c-calc.h $(SRC_FMT_M)
	cat $(SRCS) > $@
%.o : %.cu cudalib.h
	$(NVCC) -c $(NVCC_CFLAGS) $(NVCC_CPPFLAGS) $< -o $@


NVCC_CFLAGS += ${INC_DIR}
CC = ${MPICC}
CFLAGS += ${INC_DIR}


clean:
	-rm *.o cuda.cu
