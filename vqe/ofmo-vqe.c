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
#include <sys/stat.h>
#include <unistd.h>
#include <unistd.h>
#include <sys/wait.h>
#include <assert.h>


int ofmo_export_integ(const int nmonomer, const char* fpath, const int nao, const double H[],
    const double U[], const double mo_tei[], const double S[], const double C[], const double Enuc, const int nelec,
    const double ev[], const double energy, const int homo, const int lumo, const int ent, const char *desc,
    double D[], char* ansatz, char* threshold){

    char dir_path[1024];
    sprintf(dir_path, "./integ_temp/%s/", desc);
    struct stat st = {0};

#pragma omp critical
{
    if (stat(dir_path, &st) == -1){
        mkdir(dir_path, 0700);
    }
}

    FILE *fp = fopen(fpath, "w");

    fprintf(fp, "ANSATZ\t%s\n", ansatz);
    fprintf(fp, "THRESHOLD\t%s\n", threshold);

    fprintf(fp, "NMONOMER\t%d\n", nmonomer);
    fprintf(fp, "NELEC\t%d\n", nelec);
    fprintf(fp, "ENUC\t%f\n", Enuc);
    fprintf(fp, "NBASIS\t%d\n", nao);
    fprintf(fp, "HOMO\t%d\n", homo);
    fprintf(fp, "LUMO\t%d\n", lumo);
    fprintf(fp, "ENT\t%d\n", ent);
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

    fprintf(fp, "\nENV_OEI\n");
    ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<=i; j++, ij++){
            fprintf(fp, "%f\t", U[ij]);
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

    fprintf(fp, "\nMOONS\n");
    int diagonal_index = 0;
    for (int interval = 2; interval != nao+2; ++interval) {
        fprintf(fp, "%lf\n", D[diagonal_index]);
        diagonal_index += interval;
    }

    fprintf(fp, "\nDENSITY\n");
    ij=0;
    for(i=0; i<nao; i++){
        for(j=0; j<=i; j++, ij++){
            fprintf(fp, "%f\t", D[ij]);
        }
        fprintf(fp, "\n");
    }
    /*
    for(i=0; i<nao; i++){
        for(j=0; j<nao; j++, ij++){
            fprintf(fp, "%lf\t", D[ij]);
        }
        fprintf(fp, "\n");
    }
    */

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
    // fflush(stdout);
    // printf("%s -> %d eris\n", fpath, count_mo);

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


int ofmo_vqe_ofpath( const int nmonomer, const int monomer_list[], const int iscc, const int mythread, char* ofpath, char* desc){
    if(mythread != 0){
        printf("thread is not 0.\n");
        return -1;
    }
    if(nmonomer > 2){
        printf("nmonomer is neither 1 nor 2. ( nmonomer = %d )\n", nmonomer);
        fflush(stdout);
        return -1;
    }
    if(desc == NULL){
        if(nmonomer == 1) sprintf(ofpath, "./result_temp/res_mono_%d_%d_%d.dat", monomer_list[0], iscc, mythread);
        else              sprintf(ofpath, "./result_temp/res_dim_%d-%d_%d_%d.dat", monomer_list[0], monomer_list[1], iscc, mythread);
    }
    else{
        if(nmonomer == 1) sprintf(ofpath, "./result_temp/%s/res_mono_%s_%d_%d_%d.dat", desc, desc, monomer_list[0], iscc, mythread);
        else              sprintf(ofpath, "./result_temp/%s/res_dim_%s_%d-%d_%d_%d.dat", desc, desc, monomer_list[0], monomer_list[1], iscc, mythread);
    }
    return 0;
}

int ofmo_vqe_ifpath( const int nmonomer, const int monomer_list[], const int iscc, const int mythread, char* fpath, char* desc){
    if(mythread != 0){
        printf("thread is not 0.\n");
        return -1;
    }
    if(nmonomer > 2){
        printf("nmonomer is neither 1 nor 2. ( nmonomer = %d )\n", nmonomer);
        fflush(stdout);
        return -1;
    }
    if(desc == NULL){
        if(nmonomer == 1) sprintf(fpath, "./integ_temp/int_mono_%d_%d_%d.dat", monomer_list[0], iscc, mythread);
        else              sprintf(fpath, "./integ_temp/int_dim_%d-%d_%d_%d.dat", monomer_list[0], monomer_list[1], iscc, mythread);
    }
    else{
        if(nmonomer == 1) sprintf(fpath, "./integ_temp/%s/int_mono_%s_%d_%d_%d.dat", desc, desc, monomer_list[0], iscc, mythread);
        else              sprintf(fpath, "./integ_temp/%s/int_dim_%s_%d-%d_%d_%d.dat", desc, desc, monomer_list[0], monomer_list[1], iscc, mythread);
    }
    return 0;
}

int ofmo_vqe_get_energy( const int nmonomer, const int monomer_list[], const int iscc, double * energy, double * dv, char * desc ){
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


    // Read dv
    if(fgets(line, 256, fp) == NULL){
        return -1;
    }
    sscanf(line, "%lf\n", &val);

    *dv = val;

    fclose(fp);
    return 0;
}

int ofmo_vqe_get_amplitudes( const int ifrag, const int iscc, const int nao, int * namps, double ** amp, char *** fock,
    double ** corr_mat, char * desc){

    char ofpath[256];
    const int mythread = 0;
    const int nso = nao * 2;
    int ierr, monomer_list[1]={ifrag};
    int ival, iamp, m, n, nm;
    FILE * fp;
    double val;
    char line[256], idx_str[256];
    const int size_corr_mat = ((nao) * (nao + 1))/2;
    
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

    // Read dv
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
    *corr_mat = (double *) malloc (sizeof(double) * (size_corr_mat));

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

    char * token;
    for(n=0, nm=0; n<nao; n++){
        if(fgets(line, 2048, fp) == NULL){
            return -1;
        }
        token = strtok(line, "\t");
        for(m=0; m<=n; m++, nm++){
            (*corr_mat)[nm] = atof(token);
            token = strtok(NULL, "\t");
        }
    }

    fclose(fp);
    return 0;

}

int ofmo_vqe_call( const int mythread, const int nmonomer, const int monomer_list[], const int nao, const double H[],
    const double U[], const double mo_tei[], const double S[], const double C[], const int nelec, const double Enuc, const double energy,
    const int iscc, const double ev[], char *desc, const int homo, const int lumo, const int ent, double D[] ){

    char *vqescr, *ansatz, *threshold = NULL;
    ofmo_data_get_vals("vqescr ansatz threshold",
                        &vqescr, &ansatz, &threshold);

    //printf("vqescr = %s, desc = %s\n", vqescr, desc);

    /* Generate integral file */
    char fpath[256];
    char ofpath[256];
    int ierr;
    ierr = ofmo_vqe_ifpath(nmonomer, monomer_list, iscc, mythread, fpath, desc);
    if(ierr < 0) return -1;
    ierr = ofmo_vqe_ofpath(nmonomer, monomer_list, iscc, mythread, ofpath, desc);
    if(ierr < 0) return -1;
    ofmo_export_integ(nmonomer, fpath, nao, H, U, mo_tei, S, C, Enuc, nelec, ev, energy,
            homo, lumo, ent, desc, D, ansatz, threshold);

    /* Call VQE */
    const char *args[64] = {"python", vqescr, fpath, ofpath, NULL};
    ierr = exec_prog(args);
    if(ierr < 0){
        return -1;
    }

    return 0;
}


int ofmo_vqe_posthf_density( const int na, const double * A, char ** fock_vec, const double corr_mat[],
    const double C[], const int nao, double D[]){
    int i, j, ij, ia, ifk, ifk_s, m, n, nm;
    int noe;
    char * vec;
    double alpha, dt;
    const int nsao = nao * 2;
    double * Ct, *ci, *cj;
    const int n_cmat = (nao / 2) * (nao - 1);

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
        for (n=0, nm=0; n< nao; n++){
        for (m=0; m<=n; m++, nm++){
            D[ij] += corr_mat[nm] * ci[m] * cj[n];
        }
        }
        ij++;
    }
    }
    free(Ct);
    return 0;
}
