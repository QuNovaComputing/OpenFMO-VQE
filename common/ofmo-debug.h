#ifndef _OFMO_DEBUG_H_
#define _OFMO_DEBUG_H_

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-prof.h"
#include "ofmo-data.h"

extern int ofmo_show_input_data();
extern int ofmo_show_entire_molecule_data();
extern int ofmo_show_monomer_data();

#endif
