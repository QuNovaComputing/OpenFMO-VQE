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
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>


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
    int ifrag;
    for(ifrag = 0; ifrag < _nfrag; ifrag++){
        int namp = amp[ifrag]->namp;
        old_amp[ifrag]->namp = namp;
        
        if(old_amp[ifrag]->amp) free(old_amp[ifrag]->amp);
        if(old_amp[ifrag]->fock_vec) free(old_amp[ifrag]->fock_vec);
        old_amp[ifrag]->amp = (double *) malloc(namp * sizeof(double));
        old_amp[ifrag]->fock_vec = (int *) malloc(namp * sizeof(int));

        memcpy(old_amp[ifrag]->amp, amp[ifrag]->amp, namp * sizeof(double));
        memcpy(old_amp[ifrag]->fock_vec, amp[ifrag]->fock_vec, namp * sizeof(int));
    }
    return 0;
}

int ofmo_export_integ(const char* fpath, const int nao, const double H[],
    const double ao_eri_val[], const short int ao_eri_idx4[], const size_t nstored_eri,
    const double S[], const double C[], const double Enuc, const int nelec){
    

    FILE *fp = fopen(fpath, "w");

    fprintf(fp, "NELEC\t%d\n", nelec);
    fprintf(fp, "ENUC\t%f\n", Enuc);

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
    int     status;

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

int ofmo_parse_result(const char *ofpath, const int ifrag, double *energy){
    char line[256];
    ssize_t read;
    FILE * fp = fopen(ofpath, "r");
    double val;
    char idx_str[100];
    char *stop_idx_str;
    int namp, iamp;

    // Read size
    if(fgets(line, 256, fp) == NULL){
        return -1;
    }
    sscanf(line, "%d", &namp);
    amp[ifrag]->namp = namp;
    if(amp[ifrag] -> amp) free(amp[ifrag] -> amp);
    if(amp[ifrag] -> fock_vec) free(amp[ifrag] -> fock_vec);
    amp[ifrag]->amp = (double *) malloc(namp * sizeof(double));
    amp[ifrag]->fock_vec = (int *) malloc(namp * sizeof(int));

    // Read energy
    if(fgets(line, 256, fp) == NULL){
        return -1;
    }
    sscanf(line, "%lf", energy);

    // Read amplitudes
    for(iamp=0; iamp<namp; iamp++){
        if(fgets(line, 256, fp) == NULL){
            return -1;
        }
        sscanf(line, "%[0-1]\t%lf\n", idx_str, &val);
        amp[ifrag]->amp[iamp] = val;
        amp[ifrag]->fock_vec[iamp] = (int) strtol(idx_str, &stop_idx_str, 2);
    }

    fclose(fp);
    return 0;
}


int ofmo_vqe_call( const int ifrag, const int nao, const double H[],
    const double ao_eri_val[], const short int ao_eri_idx4[], const size_t nstored_eri,
    const double S[], const double C[], const int nelec, const double Enuc, double *energy){

    /* Generate integral file */
    char fpath[256];
    sprintf(fpath, "./integ_temp/temp_int_%d.dat", ifrag);
    ofmo_export_integ(fpath, nao, H, ao_eri_val, ao_eri_idx4, nstored_eri, S, C, Enuc, nelec);

#ifdef DEBUG
    return 0;
#endif

    /* Call VQE */
    char ofpath[256];
    sprintf(ofpath, "./result_temp/temp_res_%d.dat", ifrag);
    const char *args[64] = {"python", "py_script.py", fpath, ofpath, NULL};
    exec_prog(args);

    /* Data acquisition */
    ofmo_parse_result(ofpath, ifrag, energy);

    return 0;
}


int ofmo_posthf_density( const int na, const double * A, const int * fock_vec,
    const double C[], const int nao, const int max_nocc, double D[]){
    int i, j, ij, ia, ifk, ifk_s;
    int vec, noe;
    double alpha, dt;
    const int nsao = nao * 2;
    double * Ct, *ci, *cj;

    Ct  = (double*)malloc(sizeof(double) * nao * nao );

    ofmo_transpose_matrix(nao, C, Ct);

    ij = 0;
    for (i=0, ci=Ct; i<nao; i++, ci+=nao){
    for (j=0, cj=Ct; j<=i;  j++, cj+=nao){
        D[ij] = 0;
        for (ia=0; ia<na; ia++){
            vec = fock_vec[ia];
            alpha = A[ia];
            dt = 0;
            for(ifk=0; ifk<nsao; ifk+=2){
                noe = (vec>>ifk) & 1 + (vec>>(ifk + 1)) & 1;
                if(noe > 0){
                    ifk_s = ifk/2;
                    dt += ci[ifk_s] * cj[ifk_s] * ((double)noe) / 2 ;
                }
            }
            D[ij] += alpha * dt;
        }
        ij++;
    }
    }
    free(Ct);
    return 0;
}
