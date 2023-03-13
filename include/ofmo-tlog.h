#ifndef _OFMO_TLOG_H_
#define _OFMO_TLOG_H_

#ifdef USE_TLOG

#ifndef USE_MPE // tlog
#include "tlog.h"
#define TLOG_INITIALIZE()       tlog_initialize()
#define TLOG_FINALIZE()         tlog_finalize()
//#define TLOG_LOG_EVENT(a)       tlog_log(TLOG_EVENT_##a)
//#define TLOG_LOG_IN(a)          tlog_log(TLOG_EVENT_##a##_IN)
//#define TLOG_LOG_OUT(a)         tlog_log(TLOG_EVENT_##a##_OUT)
#define TLOG_LOG_EVENT(a)       tlog_log(TLOG_EVENT_1+((a)-1))
#define TLOG_LOG_IN(a)          tlog_log(TLOG_EVENT_1_IN+((a)-1)*2)
#define TLOG_LOG_OUT(a)         tlog_log(TLOG_EVENT_1_OUT+((a)-1)*2)
#else /* USE_MPE */ // MPE

#include "log_mpe.h"

#define TLOG_INITIALIZE()       mpe_setup()
#define TLOG_FINALIZE()
#define TLOG_LOG_EVENT(a)
#define TLOG_LOG_IN(a)          mpe_log_event_begin(a)
#define TLOG_LOG_OUT(a)         mpe_log_event_end(a)

#endif

#else /* !USE_TLOG */

#define TLOG_INITIALIZE()
#define TLOG_FINALIZE()
#define TLOG_LOG_EVENT(a)
#define TLOG_LOG_IN(a)
#define TLOG_LOG_OUT(a)

#endif /* USE_TLOG */
#endif /* _OFMO_TLOG_H_ */
