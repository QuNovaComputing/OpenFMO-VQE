#ifndef _LOG_MPE_H_
#define _LOG_MPE_H_

#define MAX_MPE_EVENT 10

//typedef enum {e_init, e_scc, e_scf, e_ifc4c, e_dimer, e_dimerscf, e_esdimer} mpe_event;
#define mpe_event int

#ifdef USE_MPE

void mpe_setup(void);
void mpe_def_event(mpe_event id, char *name, char *color);
void mpe_log_event_begin(mpe_event id);
void mpe_log_event_end(mpe_event id);
void mpe_log_event(mpe_event id, int stat);

#else /* !USE_MPE */

#define mpe_setup(a)
#define mpe_log_event_begin(a)
#define mpe_log_event_end(a)
#define mpe_log_event(a,b)

#endif /* USE_MPE */

#endif /* _LOG_MPE_H_ */
