/**
 * @file ofmo-vqe.c
 * A file that defines the functions for vqe calling and data acquisition.
 * 
 * */
/**
 * @defgroup ofmo-vqe vqe related funtions
 * 
 * @ingroup ofmo
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _amp_carrier_ {
    int namp;
    double * amp;
    int * fock_vec;
};

static struct _amp_carrier_ **old_amp = NULL;
static struct _amp_carrier_ **amp = NULL;
static int _nfrag;

int ofmo_amp_alloc( const int nfrag ){
    _nfrag = nfrag;
    old_amp = (struct _amp_carrier_**) malloc( nfrag * sizeof(struct _amp_carrier_) );
    amp     = (struct _amp_carrier_**) malloc( nfrag * sizeof(struct _amp_carrier_) );
    for(int ifrag=0; ifrag < nfrag; ifrag++){
        old_amp[ifrag] = (struct _amp_carrier_ * ) malloc(sizeof(struct _amp_carrier_));
        amp[ifrag]     = (struct _amp_carrier_ * ) malloc(sizeof(struct _amp_carrier_));
        old_amp[ifrag]->namp = 0;
        old_amp[ifrag]->amp = NULL;
        old_amp[ifrag]->fock_vec = NULL;
        amp[ifrag]->namp = 0;
        amp[ifrag]->amp = NULL;
        amp[ifrag]->fock_vec = NULL;
    }
    return 0;
}

int ofmo_amp_dealloc(){
    for(int ifrag=0; ifrag < _nfrag; ifrag++){
        if(old_amp[ifrag]->amp != NULL)      free(old_amp[ifrag]->amp);
        if(old_amp[ifrag]->fock_vec != NULL) free(old_amp[ifrag]->fock_vec);
        if(amp[ifrag]->amp != NULL)          free(old_amp[ifrag]->amp);
        if(amp[ifrag]->fock_vec != NULL)     free(old_amp[ifrag]->fock_vec);
        free(old_amp[ifrag]);
        free(amp[ifrag]);
    }
    free(old_amp);
    free(amp);
    return 0;
}

int ofmo_get_amps( const int ifrag, int *namp, double **alpha, int **fock_vec){
    *namp = amp[ifrag]->namp;
    *alpha = amp[ifrag]->amp;
    *fock_vec = amp[ifrag]->fock_vec;
    return 0;
}

int ofmo_get_oldamps( const int ifrag, int *namp, double **old_alpha, int **old_fock_vec){
    *namp = old_amp[ifrag]->namp;
    *old_alpha = old_amp[ifrag]->amp;
    *old_fock_vec = old_amp[ifrag]->fock_vec;
    return 0;
}

int ofmo_update_amps(){

}

int ofmo_vqe_call( const int ifrag, const int nao, const double H[],
    const double ao_eri_val[], const short int ao_eri_idx4[], const size_t nstored_eri,
    const double S[], const double C[], const int nelec, double *energy){

    /* Generate integral file */
    char fpath[256];
    sprintf(fpath, "./integ_temp/temp_int_%d.dat", ifrag);
    ofmo_export_integ(fpath, nao, H, ao_eri_val, ao_eri_idx4, nstored_eri, S, C, nelec);

    /* Call VQE */
    char ofpath[256];
    sprintf(ofpath, "./result_temp/temp_res_%d.dat");
    char *args[] = {"python", "py_script.py", fpath, ofpath};
    exec_prog(args);

    /* Data acquisition */


    /* Update to the memory */

    /* Outputs */

    return 0;
}

int ofmo_export_integ(const* fpath, const int nao, const double H[],
    const double ao_eri_val[], const short int ao_eri_idx4[], const size_t nstored_eri,
    const double S[], const double C[], const int nelec){
    

    FILE *fp = fopen(fpath, "w");

    fprintf("NELEC\t%d\n", nelec);

    int nao2 = ( nao + 1 ) * nao / 2;
    int i,j,k,l;
    fprintf(fp, "OEI\n");
    int ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<i; j++, ij++){
            fprintf(fp, "%f\t", H[ij]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nOVERLAP\n");
    ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<i; j++, ij++){
            fprintf(fp, "%f\t", S[ij]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nMOCOEFF\n");
    ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<nao; j++){
            fprintf(fp, "%f\t", C[ij]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nERI\n");
    int ix, ix4;
    for(ix=0, ix4=0; ix<nstored_eri; ix++, ix4+=4){
        i = (int) ao_eri_idx4[ix4+0];
        j = (int) ao_eri_idx4[ix4+1];
        k = (int) ao_eri_idx4[ix4+2];
        l = (int) ao_eri_idx4[ix4+3];
        fprintf(fp, "%d\t%d\t%d\t%d\t%f\n", i, j, k, l, ao_eri_val[ix]);
    }

    fclose(fp);

    return 0;
}

static int exec_prog(const char **argv)
{
    pid_t   my_pid;
    int     status, timeout /* unused ifdef WAIT_FOR_COMPLETION */;

    if (0 == (my_pid = fork())) {
            if (-1 == execve(argv[0], (char **)argv , NULL)) {
                    perror("child process execve failed [%m]");
                    return -1;
            }
    }

    while (0 == waitpid(my_pid , &status , WNOHANG)) {
            sleep(1);
    }

    printf("%s WEXITSTATUS %d WIFEXITED %d [status %d]\n",
            argv[0], WEXITSTATUS(status), WIFEXITED(status), status);

    if (1 != WIFEXITED(status) || 0 != WEXITSTATUS(status)) {
            perror("%s failed, halt system");
            return -1;
    }

    return 0;
}

int ofmo_parse_result(const char *ofpath, const int ifrag){
    FILE * fp = fopen(ofpath, "r");
    
}