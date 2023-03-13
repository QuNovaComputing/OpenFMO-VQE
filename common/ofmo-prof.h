#ifndef _OFMO_PROF_H_
#define _OFMO_PROF_H_

#include <stdio.h>
#include <stdlib.h>

//#include <mpi.h>
#include "ofmo-parallel.h"

extern FILE* fp_prof;

extern int ofmo_prof_init( const char *prof_name, MPI_Comm comm );
extern int ofmo_prof_reinit( const char *prof_name, MPI_Comm comm );


extern int ofmo_create_worker_prof_name( const char *input,
	const char *header,
	const int group_id, char *prof_name, size_t maxstrlen );

extern int ofmo_create_mserv_prof_name( const char *header,
	const int mserv_id, char *prof_name, size_t maxstrlen );

// スレッド毎のタイマー
extern int ofmo_create_thread_timer( const char *str, const int attr );
extern void ofmo_reset_all_thread_timer();
extern void ofmo_start_thread_timer( const int id, const int mythread );
extern void ofmo_acc_thread_timer( const int id, const int mythread );
extern void ofmo_set_thread_timer( const int id, const int mythread,
	const double val );
extern void ofmo_show_thread_timer_all();
extern void ofmo_show_thread_timers( int n, int list[] );

// プロセス毎のタイマー
extern int ofmo_create_proc_timer( const char *str, const int attr );
extern void ofmo_reset_all_proc_timer();
extern void ofmo_start_proc_timer( const int id );
extern void ofmo_acc_proc_timer( const int id );
extern void ofmo_show_proc_timer_all();
extern void ofmo_show_proc_timers( int n, int list[] );

#endif
