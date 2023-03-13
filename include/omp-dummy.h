#ifndef _OMP_DUMMY_H_
#define _OMP_DUMMY_H_

#include <sys/time.h>

static int omp_get_max_threads() { return 1; }
static int omp_get_thread_num() { return 0; }
static int omp_get_num_threads() { return 1; }

static double omp_get_wtime() {
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return ( (double)tv.tv_sec + 1.e-6*(double)tv.tv_usec );
}

#endif
