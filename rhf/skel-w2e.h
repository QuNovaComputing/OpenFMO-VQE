#ifndef _SKEL_W2E_H_
#define _SKEL_W2E_H_

#include "ofmo-def.h"

void setup_w2e(const MPI_Comm wcomm, const int optsync);
void start_w2e(void);
void set_w2e(const int Labcd);
void print_w2e(void);

#endif /* _SKEL_W2E_H_ */

