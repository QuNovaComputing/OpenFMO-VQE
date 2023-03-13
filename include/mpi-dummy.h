#ifndef _MPI_DUMMY_DEF_H_
#define _MPI_DUMMY_DEF_H_

#include <stdio.h>
#include <sys/time.h>
#include <string.h>

#define MPI_COMM_WORLD	1
#define MPI_COMM_NULL   0

#define MPI_SUCCESS 0

//copy from mpi.h in impi
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)0x4c000101)
#define MPI_SIGNED_CHAR    ((MPI_Datatype)0x4c000118)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)0x4c000102)
#define MPI_BYTE           ((MPI_Datatype)0x4c00010d)
#define MPI_WCHAR          ((MPI_Datatype)0x4c00040e)
#define MPI_SHORT          ((MPI_Datatype)0x4c000203)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)0x4c000204)
#define MPI_INT            ((MPI_Datatype)0x4c000405)
#define MPI_UNSIGNED       ((MPI_Datatype)0x4c000406)
#define MPI_LONG           ((MPI_Datatype)0x4c000807)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)0x4c000808)
#define MPI_FLOAT          ((MPI_Datatype)0x4c00040a)
#define MPI_DOUBLE         ((MPI_Datatype)0x4c00080b)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)0x4c00100c)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)0x4c000809)
#define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_LONG_LONG      MPI_LONG_LONG_INT

//copy from mpi.h in impi
typedef int MPI_Op;
#define MPI_MAX     (MPI_Op)(0x58000001)
#define MPI_MIN     (MPI_Op)(0x58000002)
#define MPI_SUM     (MPI_Op)(0x58000003)
#define MPI_PROD    (MPI_Op)(0x58000004)
#define MPI_LAND    (MPI_Op)(0x58000005)
#define MPI_BAND    (MPI_Op)(0x58000006)
#define MPI_LOR     (MPI_Op)(0x58000007)
#define MPI_BOR     (MPI_Op)(0x58000008)
#define MPI_LXOR    (MPI_Op)(0x58000009)
#define MPI_BXOR    (MPI_Op)(0x5800000a)
#define MPI_MINLOC  (MPI_Op)(0x5800000b)
#define MPI_MAXLOC  (MPI_Op)(0x5800000c)
#define MPI_REPLACE (MPI_Op)(0x5800000d)


#define MPI_THREAD_SINGLE 0
#define MPI_THREAD_FUNNELED 1
#define MPI_THREAD_SERIALIZED 2
#define MPI_THREAD_MULTIPLE 3

#define MPI_IN_PLACE  (void *) -1

/* Results of the compare operations. */
#define MPI_IDENT     0
#define MPI_CONGRUENT 1
#define MPI_SIMILAR   2
#define MPI_UNEQUAL   3

typedef int MPI_Comm;
typedef int MPI_Status;

static int MPI_Init_thread( int *argc, char*** argv, int require,
	int* provided ) {
    *provided = require;
    return 0; }

static int MPI_Init( int *argc, char*** argv ) {
    return 0; }

static int MPI_Comm_size( MPI_Comm comm, int *nprocs ) {
    *nprocs = 1;
    return 0;
}

static int MPI_Comm_rank( MPI_Comm comm, int *myrank ) {
    *myrank = 0;
    return 0;
}

static void MPI_Abort( MPI_Comm comm, int errcode ) { exit( errcode ); }

static int MPI_Bcast( void* buf, int count, int type, int root,
	MPI_Comm comm ) { return 0; }

static int MPI_Finalize() { return 0; }

static double MPI_Wtime() {
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return (double)tv.tv_sec+(double)tv.tv_usec*1.e-6;
}

static int MPI_Recv( void* dest, int count, int dtype,
	int root, int tag, MPI_Comm comm, MPI_Status *status ) {
    return 0;
}

static int MPI_Send( void* src, int count, int dtype,
	int root, int tag, MPI_Comm comm ) {
    return 0;
}

static int MPI_Allreduce( void* src, void* buf, int count,
	int type, int op, MPI_Comm comm ) {
    int tsize;
    if ( type == MPI_DOUBLE ) {
	tsize = sizeof(double);
    } else if ( type == MPI_INT ) {
	tsize = sizeof(int);
    } else if ( type == MPI_CHAR ) {
	tsize = sizeof(char);
    } else if ( type == MPI_LONG_LONG_INT ) {
	tsize = sizeof(long long int);
    } else {
	printf("ERROR: Illegal type (%d)\n", type );
	return -1;
    }
    if (src != MPI_IN_PLACE) memcpy( buf, src, tsize*count );
    return 0;
}

static int MPI_Reduce( void* src, void* buf, int count,
	int type, int op, int root, MPI_Comm comm ) {
    int tsize;
    if ( type == MPI_DOUBLE ) {
	tsize = sizeof(double);
    } else if ( type == MPI_INT ) {
	tsize = sizeof(int);
    } else if ( type == MPI_CHAR ) {
	tsize = sizeof(char);
    } else if ( type == MPI_LONG_LONG_INT ) {
	tsize = sizeof(long long int);
    } else {
	printf("ERROR: Illegal type (%d)\n", type );
	return -1;
    }
    if (src != MPI_IN_PLACE) memcpy( buf, src, tsize*count );
    return 0;
}

static int MPI_Gather( void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                       void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                       int root, MPI_Comm comm)
{
    size_t tsize;
    if ( recvtype != sendtype ) {
	printf("ABORT: MPI_Gather() in mpi-dummy.h \n" );
	return -1;
    } else if ( recvtype == MPI_DOUBLE ) {
	tsize = sizeof(double);
    } else if ( recvtype == MPI_INT ) {
	tsize = sizeof(int);
    } else if ( recvtype == MPI_CHAR ) {
	tsize = sizeof(char);
    } else if ( recvtype == MPI_LONG_LONG_INT ) {
	tsize = sizeof(long long int);
    } else {
	printf("ERROR: Illegal type (%d)\n", recvtype );
	return -1;
    }
    memcpy( recvbuf, sendbuf, tsize*recvcnt );
    return 0;
}

static int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result)
{
  *result = MPI_UNEQUAL;
  if (comm1 == MPI_COMM_NULL || comm2 == MPI_COMM_NULL) return -1;
  *result = (comm1==comm2)? MPI_IDENT: MPI_UNEQUAL;
  return 0;
}

#endif
