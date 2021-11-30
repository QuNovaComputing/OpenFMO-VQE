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
#include <assert.h>

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

int ofmo_export_integ(const int nmonomer, const char* fpath, const int nao, const double H[],
    const double mo_tei[], const double S[], const double C[], const double Enuc, const int nelec,
    const double ev[], const double energy){
    

    FILE *fp = fopen(fpath, "w");

    fprintf(fp, "NMONOMER\t%d\n", nmonomer);
    fprintf(fp, "NELEC\t%d\n", nelec);
    fprintf(fp, "ENUC\t%f\n", Enuc);
    fprintf(fp, "NBASIS\t%d\n", nao);
    fprintf(fp, "HF_ENERGY\t%f\n", energy);

    int nao2 = ( nao + 1 ) * nao / 2;
    int i,j,k,l;

    fprintf(fp, "ENERGIES\n");
    for(i=0; i<nao; i++){
        fprintf(fp, "%f\n", ev[i]);
    }

    fprintf(fp, "\nOEI\n");
    int ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<=i; j++, ij++){
            fprintf(fp, "%f\t", H[ij]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nOVERLAP\n");
    ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<=i; j++, ij++){
            fprintf(fp, "%f\t", S[ij]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nMOCOEFF\n");
    ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<nao; j++, ij++){
            fprintf(fp, "%f\t", C[ij]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nERI\n");
    const int nao_2 = nao * nao;
    const int nao_3 = nao_2 * nao;
    int imo, imo4, jmo, jmo3, kmo, kmo2, lmo, mo_idx;
    int count_mo = 0;
    for(imo=0, imo4=0; imo<nao;   imo++, imo4+=nao_3){
    for(jmo=0, jmo3=0; jmo<imo+1; jmo++, jmo3+=nao_2){
    for(kmo=0, kmo2=0; kmo<imo+1; kmo++, kmo2+=nao){
    for(lmo=0; lmo<kmo+1; lmo++, count_mo++){
        /* mo_idx = imo * nao * nao * nao
                +jmo * nao * nao
                +kmo * nao
                +lmo; */
        mo_idx = imo4 + jmo3 + kmo2 + lmo;
        fprintf(fp, "%d %d %d %d %.7f\n", imo, jmo, kmo, lmo, mo_tei[mo_idx]);
    }}}}
    fflush(stdout);
    printf("%s -> %d eris\n", fpath, count_mo);

    fclose(fp);

    return 0;
}

static int exec_prog(const char **argv)
{
    pid_t   my_pid;
    int     status;

    if (0 == (my_pid = fork())) {
            if (-1 == execvp(argv[0], (char **)argv)) {
                    perror("child process execve failed [%m]");
                    return -1;
            }
    }

    while (0 == waitpid(my_pid , &status , WNOHANG)) {
            sleep(1);
    }

    //printf("%s WEXITSTATUS %d WIFEXITED %d [status %d]\n",
    //        argv[0], WEXITSTATUS(status), WIFEXITED(status), status);

    if (1 != WIFEXITED(status) || 0 != WEXITSTATUS(status)) {
            perror("%s failed, halt system");
            return -1;
    }

    return 0;
}

int ofmo_parse_result(const char *ofpath, const int nmonomer, const int monomer_list[], double *energy){
    char line[256];
    ssize_t read;
    FILE * fp = fopen(ofpath, "r");
    double val;
    char idx_str[256];
    char *stop_idx_str;
    int namp, iamp;

    // Read energy
    if(fgets(line, 256, fp) == NULL){
        return -1;
    }
    sscanf(line, "%lf\n", &val);
    
    //printf("%s : energy=%lf, %s\n", ofpath, val, line);

    if (nmonomer==1){
        // Read size
        if(fgets(line, 256, fp) == NULL){
            return -1;
        }
        sscanf(line, "%d", &namp);
        // init amps
        /*int ifrag = monomer_list[0];
        amp[ifrag]->namp = namp;
        if(amp[ifrag] -> amp) free(amp[ifrag] -> amp);
        if(amp[ifrag] -> fock_vec) free(amp[ifrag] -> fock_vec);
        amp[ifrag]->amp = (double *) malloc(namp * sizeof(double));
        amp[ifrag]->fock_vec = (int *) malloc(namp * sizeof(int));
        */
       // Read amplitudes
        for(iamp=0; iamp<namp; iamp++){
            if(fgets(line, 256, fp) == NULL){
                return -1;
            }
            // sscanf(line, "%[0-1]\t%lf\n", idx_str, &val);
                int ifrag = monomer_list[0];
                //printf("%s : %s, %f\n", ofpath, idx_str, val);

                //amp[ifrag]->amp[iamp] = val;
                //amp[ifrag]->fock_vec[iamp] = (int) strtol(idx_str, &stop_idx_str, 2);
            }
    }

    fclose(fp);
    return 0;
}


int ofmo_vqe_ofpath( const int nmonomer, const int monomer_list[], const int iscc, const int mythread, char* ofpath, char* desc){
    if(mythread != 0){
        printf("thread is not 0.\n");
        return -1;
    }
    if(desc == NULL){
        desc = "\0";
    }
    if(nmonomer == 1){
        sprintf(ofpath, "./result_temp/res_%s_mono_%d_%d_%d.dat", desc, monomer_list[0], iscc, mythread);
    }
    else if(nmonomer == 2){
        sprintf(ofpath, "./result_temp/res_%s_dim_%d-%d_%d_%d.dat", desc, monomer_list[0], monomer_list[1], iscc, mythread);
    }
    else{
        printf("nmonomer is neither 1 nor 2. ( nmonomer = %d )\n", nmonomer);
        fflush(stdout);
        return -1;
    }
    return 0;
}

int ofmo_vqe_ifpath( const int nmonomer, const int monomer_list[], const int iscc, const int mythread, char* fpath, char* desc){
    if(mythread != 0){
        printf("thread is not 0.\n");
        return -1;
    }
    if(desc == NULL){
        desc = "\0";
    }
    if(nmonomer == 1){
        sprintf(fpath, "./integ_temp/int_%s_mono_%d_%d_%d.dat", desc, monomer_list[0], iscc, mythread);
    }
    else if(nmonomer == 2){
        sprintf(fpath, "./integ_temp/int_%s_dim_%d-%d_%d_%d.dat", desc, monomer_list[0], monomer_list[1], iscc, mythread);
    }
    else{
        printf("nmonomer is neither 1 nor 2. ( nmonomer = %d )\n", nmonomer);
        fflush(stdout);
        return -1;
    }
    return 0;
}

int ofmo_vqe_get_energy( const int nmonomer, const int monomer_list[], const int iscc, double * energy, char * desc ){
    char ofpath[256];
    const int mythread = 0;
    int ierr;
    FILE * fp;
    double val;
    char line[256];

    ierr = ofmo_vqe_ofpath(nmonomer, monomer_list, iscc, mythread, ofpath, desc);
    if(ierr < 0){
        return -1;
    }
    fp = fopen(ofpath, "r");
    // Read energy
    if(fgets(line, 256, fp) == NULL){
        return -1;
    }
    sscanf(line, "%lf\n", &val);

    *energy = val;

    fclose(fp);
    return 0;
}

int ofmo_vqe_get_amplitudes( const int ifrag, const int iscc, const int nso, int * namps, double ** amp, char *** fock, char * desc){

    char ofpath[256];
    const int mythread = 0;
    int ierr, monomer_list[1]={ifrag};
    int ival, iamp;
    FILE * fp;
    double val;
    char line[256], idx_str[256];
    
    ierr = ofmo_vqe_ofpath(1, monomer_list, iscc, mythread, ofpath, desc);
    if(ierr < 0){
        return -1;
    }
    fp = fopen(ofpath, "r");
    // Read energy
    if(fgets(line, 256, fp) == NULL){
        return -1;
    }
    sscanf(line, "%lf\n", &val);

    // Read size
    if(fgets(line, 256, fp) == NULL){
        return -1;
    }
    sscanf(line, "%d", &ival);

    // Allocate
    *namps = ival;
    *amp = (double *) malloc (sizeof(double) * (*namps));
    *fock = (char **) malloc (sizeof(char*) * (*namps));

    // Read amplitudes
    for(iamp=0; iamp<*namps; iamp++){
        if(fgets(line, 256, fp) == NULL){
            return -1;
        }
        sscanf(line, "%[0-1]\t%lf\n", idx_str, &val);
        (*fock)[iamp] = (char *) malloc (sizeof(char) * (nso + 1));
        memcpy((*fock)[iamp], idx_str, sizeof(char) * (nso + 1));
        (*amp)[iamp] = val;
    }

    fclose(fp);
    return 0;

}

int ofmo_vqe_call( const int mythread, const int nmonomer, const int monomer_list[], const int nao, const double H[],
    const double mo_tei[], const double S[], const double C[], const int nelec, const double Enuc, const double energy,
    const int iscc, const double ev[], char *desc){

    char * vqescr = NULL;
    ofmo_data_get_vals("vqescr", &vqescr);

    //printf("vqescr = %s, desc = %s\n", vqescr, desc);

    /* Generate integral file */
    char fpath[256];
    char ofpath[256];
    int ierr;
    ierr = ofmo_vqe_ifpath(nmonomer, monomer_list, iscc, mythread, fpath, desc);
    if(ierr < 0){
        return -1;
    }
    ierr = ofmo_vqe_ofpath(nmonomer, monomer_list, iscc, mythread, ofpath, desc);
    if(ierr < 0){
        return -1;
    }
    ofmo_export_integ(nmonomer, fpath, nao, H, mo_tei, S, C, Enuc, nelec, ev, energy);

    /* Call VQE */
    const char *args[64] = {"python", vqescr, fpath, ofpath, NULL};
    ierr = exec_prog(args);
    if(ierr < 0){
        return -1;
    }

    return 0;
}


int ofmo_vqe_posthf_density( const int na, const double * A, char ** fock_vec,
    const double C[], const int nao, double D[]){
    int i, j, ij, ia, ifk, ifk_s;
    int noe;
    char * vec;
    double alpha, dt;
    const int nsao = nao * 2;
    double * Ct, *ci, *cj;

    Ct  = (double*)malloc(sizeof(double) * nao * nao );

    ofmo_transpose_matrix(nao, C, Ct);

    ij = 0;
    for (i=0, ci=Ct; i<nao; i++, ci+=nao){
    for (j=0, cj=Ct; j<=i;  j++, cj+=nao){
        D[ij] = 0.0;
        for (ia=0; ia<na; ia++){
            vec = fock_vec[ia];
            alpha = A[ia];
            dt = 0;
            for(ifk=0; ifk<nsao; ifk+=2){
                noe = ((vec[ifk] == '1') & 1) + ((vec[ifk + 1] == '1') & 1);
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
