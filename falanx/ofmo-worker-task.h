
#ifndef _OFMO_WORKER_TASK_H_
#define _OFMO_WORKER_TASK_H_

#include <falanx_worker.h>

/* falanx task declarations */

/*
 * 初期密度行列計算
 */
void ofmo_init_dens_task(falanx_context_t* context, MPI_Comm worker);

/*
 * モノマー間距離計算
 */
void ofmo_distance_task(falanx_context_t* context, MPI_Comm worker);

/*
 * SCF計算
 */
void ofmo_scf_task(falanx_context_t* context, MPI_Comm worker);

/*
 * APPROXタスク
 */
void ofmo_approx_task(falanx_context_t* context, MPI_Comm worker);

#define OFMO_CMD_INIT_DENS      "OFMO_INIT_DENS"
#define OFMO_CMD_DISTANCE       "OFMO_DISTANCE"
#define OFMO_CMD_SCF            "OFMO_SCF"
#define OFMO_CMD_APPROX         "OFMO_APPROX"


#endif

