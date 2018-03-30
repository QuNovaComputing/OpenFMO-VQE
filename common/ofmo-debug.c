#include <stdio.h>
#include <stdlib.h>
#include "ofmo-parallel.h"

#include "ofmo-def.h"
#include "ofmo-prof.h"
#include "ofmo-data.h"

int ofmo_show_input_data() {
    int iat, ifrag, ierr, imo, ifn, ibd;
    if ( fp_prof == NULL ) return 0;
    // $fmoxyz section
    int natom, *atomic_number;
    double *atom_x, *atom_y, *atom_z;
    ierr = ofmo_data_get_vals("natom atn atx aty atz",
	    &natom, &atomic_number, &atom_x, &atom_y, &atom_z );
    if ( ierr != 0 ) {
	fdbg(fp_prof, "no data in $fmoxyz section");
	return -1;
    }
    fprintf(fp_prof, "============= $fmoxyz section data =============\n");
    fprintf(fp_prof, "# of atom in entire molecule = %d\n", natom );
    fprintf(fp_prof, "------------------------------------------------\n");
    //                123456---12--123456789012-123456789012-123456789012
    fprintf(fp_prof, "   SN  : atn        x            y            z\n");
    fprintf(fp_prof, "------------------------------------------------\n");
    for ( iat=0; iat<natom; iat++ ) {
	fprintf(fp_prof, "%6d : %2d  %12.7f %12.7f %12.7f\n",
		(iat+1), atomic_number[iat],
		atom_x[iat], atom_y[iat], atom_z[iat] );
    }
    fprintf(fp_prof, "------------------------------------------------\n");
    // $fmo section
    int nfrag, *icharge, *at2frg, nbody;
    double lptc, laop, ldim;
    ierr = ofmo_data_get_vals("nfrag icharg at2frg nbody lptc laop ldim",
	    &nfrag, &icharge, &at2frg, &nbody, &lptc, &laop, &ldim );
    if ( ierr != 0 ) {
	fdbg(fp_prof, "no data in $fmo section\n");
	return -1;
    }
    fprintf(fp_prof, "\n");
    fprintf(fp_prof, "================ $fmo section data =============\n");
    fprintf(fp_prof, "# of fragment = %d\n", nfrag );
    fprintf(fp_prof, "nbody         = %d\n", nbody );
    fprintf(fp_prof, "lptc          = %4.1f\n", lptc );
    fprintf(fp_prof, "laop          = %4.1f\n", laop );
    fprintf(fp_prof, "ldim          = %4.1f\n", ldim );
    //--
    fprintf(fp_prof, "icharge = ");
    for ( ifrag=0; ifrag<nfrag; ifrag++ )
	fprintf(fp_prof, "%s%2d", ((ifrag%10)==0?"\n    " : " "),
		icharge[ifrag] );
    if ( (nfrag%10) != 0 ) fprintf(fp_prof, "\n");
    //--
    fprintf(fp_prof, "at2frg = ");
    for ( iat=0; iat<natom; iat++ )
	fprintf(fp_prof, "%s%4d", ((iat%10)==0?"\n    " : " "),
		at2frg[iat] );
    if ( (natom%10) != 0 ) fprintf(fp_prof, "\n");
    fprintf(fp_prof, "------------------------------------------------\n");
    // $fmolmo section
    int nlmo, nfunc, *ban;
    double *clmo;
    char *lmobs;
    ierr = ofmo_data_get_vals("nlmo nfunc clmo ban lmobs",
	    &nlmo, &nfunc, &clmo, &ban, &lmobs );
    if ( ierr != 0 ) {
	fdbg(fp_prof, "no data in $fmolmo section\n");
	return -1;
    }
    fprintf(fp_prof, "\n");
    fprintf(fp_prof, "=============== $fmolmo section data ===========\n");
    fprintf(fp_prof, "# of localized MO = %d\n", nlmo );
    fprintf(fp_prof, "# of basis func.  = %d\n", nfunc );
    fprintf(fp_prof, "basis name        = %s\n", lmobs );
    for ( imo=0; imo<nlmo; imo++ ) {
	fprintf(fp_prof, " %d %d ", ban[imo], !ban[imo] );
	for ( ifn=0; ifn<nfunc; ifn++ )
	    fprintf(fp_prof, "%s%10.6f", ((ifn+1)%5==0 ? "\n      ": " "),
		    clmo[imo*nfunc + ifn] );
	fprintf(fp_prof, "\n");
    }
    fprintf(fp_prof, "------------------------------------------------\n");
    // $fmobnd section
    int nbond, *welec, *woelec;
    ierr = ofmo_data_get_vals("nbond welec woelec",
	    &nbond, &welec, &woelec );
    if ( ierr != 0 ) {
	fdbg(fp_prof, "no data in $fmobnd section\n");
	return -1;
    }
    fprintf(fp_prof, "\n");
    fprintf(fp_prof, "=============== $fmobnd section data ===========\n");
    fprintf(fp_prof, "# of bond atom = %d\n", nbond );
    for ( ibd=0; ibd<nbond; ibd++ )
	fprintf(fp_prof, "%6d %5d\n", woelec[ibd], welec[ibd] );
    return 0;
}

int ofmo_show_entire_molecule_data() {
    int ierr, ibs, ics, iat, iw;
    if ( fp_prof == NULL ) return 0;
    int natom;
    int nsbs, *ushel_lqn, *ushel_ini, ncs, nao, nps, maxlqn, *atm_lcs;
    char **bslst;
    ierr = ofmo_data_get_vals("natom nsbs bslst maxlqn ncs nao nps "
	    "ushel_lqn ushel_ini atm_lcs",
	    &natom, &nsbs, &bslst, &maxlqn, &ncs, &nao, &nps,
	    &ushel_lqn, &ushel_ini, &atm_lcs );
    if ( ierr != 0 ) {
	fdbg(fp_prof, "no data of entire molecule\n");
	return -1;
    }
    fprintf(fp_prof, "========= basis data for entire molecule =======\n");
    fprintf(fp_prof, "nat= %d, ncs= %d, nao= %d, nps= %d\n",
	    natom, ncs, nao, nps );
    fprintf(fp_prof, "# of species of basis set = %d\n", nsbs );
    for ( ibs=0; ibs<nsbs; ibs++ )
	fprintf(fp_prof, "%d %16s\n", (ibs+1), bslst[ibs] );
    fprintf(fp_prof, "------------------------------------------------\n");
    //-
    fprintf(fp_prof, "ushel_lqn =");
    for ( ics=0; ics<ncs; ics++ )
	fprintf(fp_prof, "%s%2d", (ics%20==0 ? "\n    " : " "),
		ushel_lqn[ics] );
    if ( (ncs%20) != 0 ) fprintf(fp_prof, "\n");
    //-
    if      ( nao < 10 ) iw=1;
    else if ( nao < 100 ) iw=2;
    else if ( nao < 1000 ) iw=3;
    else if ( nao < 10000 ) iw=4;
    else if ( nao < 100000 ) iw=5;
    else if ( nao < 1000000 ) iw=6;
    fprintf(fp_prof, "ushel_ini =");
    for ( ics=0; ics<ncs; ics++ )
	fprintf(fp_prof, "%s%*d", (ics%10==0 ? "\n    " : " "), iw,
		ushel_ini[ics] );
    if ( (ncs%10) != 0 ) fprintf(fp_prof, "\n");
    //-
    if      ( ncs < 10 ) iw=1;
    else if ( ncs < 100 ) iw=2;
    else if ( ncs < 1000 ) iw=3;
    else if ( ncs < 10000 ) iw=4;
    else if ( ncs < 100000 ) iw=5;
    else if ( ncs < 1000000 ) iw=6;
    fprintf(fp_prof, "leading CS of atom =");
    for ( iat=0; iat<natom; iat++ )
	fprintf(fp_prof, "%s%*d", (iat%10==0 ? "\n    " : " "), iw,
		atm_lcs[iat] );
    if ( (natom%10) != 0 ) fprintf(fp_prof, "\n");
    return 0;
}

int ofmo_show_monomer_data() {
    int ifrag, ics, ips, ips0, ips1, ics0, ics1, lqn;
    int nfrag, ierr, maxlqn;
    char *CS="SPDFGHIJ";
    if ( fp_prof == NULL ) return 0;
    ierr = ofmo_data_get_vals("nfrag", &nfrag );
    if ( ierr != 0 ) {
	fdbg(fp_prof, "nfrag values is not set\n");
	return -1;
    }
    //
    int maxnfcs, maxnfao, maxnfps, maxnpspair;
    int *nfcs,*nfao, *nfps, **mlcs, *nfatom;
    int **mshel_tem, **mshel_atm, **mshel_add, **mshel_ini;
    double **mprim_exp, **mprim_coe;
    ierr = ofmo_data_get_vals("nfrag maxlqn nfatom "
	    "maxnfcs maxnfao maxnfps maxnpspair "
	    "nfcs nfao nfps mlcs "
	    "mshel_tem mshel_atm mshel_add mshel_ini "
	    "mprim_exp mprim_coe",
	    &nfrag, &maxlqn, &nfatom,
	    &maxnfcs, &maxnfao, &maxnfps, &maxnpspair,
	    &nfcs, &nfao, &nfps, &mlcs,
	    &mshel_tem, &mshel_atm, &mshel_add, &mshel_ini,
	    &mprim_exp, &mprim_coe );
    if ( ierr != 0 ) {
	fdbg(fp_prof, "no data of monomer basis\n");
	return -1;
    }
    fprintf(fp_prof, "============= monomer basis data ===============\n");
    fprintf(fp_prof, "# of monomer = %d\n", nfrag );
    fprintf(fp_prof, "max. orbital quantum number = %d\n", maxlqn );
    fprintf(fp_prof, "max. # of CS in monomer = %d\n", maxnfcs );
    fprintf(fp_prof, "max. # of AO in monomer = %d\n", maxnfao );
    fprintf(fp_prof, "max. # of PS in monomer = %d\n", maxnfps );
    fprintf(fp_prof, "max. # of nps pair      = %d\n", maxnpspair );
    //-
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	fprintf(fp_prof, "frag=%d  nat=%d, ncs=%d, nao=%d, nps=%d\n",
		ifrag, nfatom[ifrag], nfcs[ifrag], nfao[ifrag],
		nfps[ifrag] );
	for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	    ics0 = mlcs[ifrag][lqn];
	    ics1 = mlcs[ifrag][lqn+1];
	    fprintf(fp_prof, "%c-shell\n", CS[lqn] );
	    for ( ics=ics0; ics<ics1; ics++ ) {
		ips0 = mshel_add[ifrag][ics];
		ips1 = ips0 + mshel_tem[ifrag][ics];
		fprintf(fp_prof, "%3d %3d :",
			mshel_atm[ifrag][ics], mshel_ini[ifrag][ics] );
		for ( ips=ips0; ips<ips1; ips++ )
		    fprintf(fp_prof, "   %8.2e %8.2e",
			    mprim_exp[ifrag][ips], mprim_coe[ifrag][ips] );
		fprintf(fp_prof, "\n");
	    }
	}
	fprintf(fp_prof, "\n");
    }
    return 0;
}
