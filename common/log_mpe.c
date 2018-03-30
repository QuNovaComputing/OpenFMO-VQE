#ifdef USE_MPE
#include "stdio.h"
#include "mpe.h"
#include "log_mpe.h"

static int mpe_eventID[MAX_MPE_EVENT][2];

void mpe_setup(void)
{
  int i;

  for (i=0; i<MAX_MPE_EVENT; i++) {
    mpe_eventID[i][0] = MPE_Log_get_event_number();
    mpe_eventID[i][1] = MPE_Log_get_event_number();
  }

  mpe_def_event(0, "E0", "white");
  mpe_def_event(1, "E1", "dark gray");
  mpe_def_event(2, "E2", "aquamarine");
  mpe_def_event(3, "E3", "green");
  mpe_def_event(4, "E4", "red");
  mpe_def_event(5, "E5", "blue");
  mpe_def_event(6, "E6", "yellow");
  mpe_def_event(7, "E7", "purple");
  mpe_def_event(8, "E8", "magenta");
  mpe_def_event(9, "E9", "pink");
}

void mpe_def_event(mpe_event id, char *name, char *color)
{
  MPE_Describe_state(mpe_eventID[id][0], mpe_eventID[id][1], name, color);
}

void mpe_log_event(mpe_event id, int stat)
{
  MPE_Log_event(mpe_eventID[id][stat], 0, NULL);
}

void mpe_log_event_begin(mpe_event id)
{
  MPE_Log_event(mpe_eventID[id][0], 0, NULL);
}

void mpe_log_event_end(mpe_event id)
{
  MPE_Log_event(mpe_eventID[id][1], 0, NULL);
}


#endif
