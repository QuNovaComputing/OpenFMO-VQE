#ifndef _OFMO_OS_XXXX_H
#define _OFMO_OS_XXXX_H

// for partiall direct SCF

extern int ofmo_twoint_os_ssss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_psss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_psps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_ppss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_ppps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_pppp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dsss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dsps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dspp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dsds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dpss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dpps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dppp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dpds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dpdp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_ddss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_ddps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_ddpp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_ddds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dddp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );

extern int ofmo_twoint_os_dddd(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs );


// for direct SCF

extern int ofmo_twoint_direct_os_ssss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_psss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_psps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_ppss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_ppps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_pppp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dsss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dsps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dspp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dsds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dpss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dpps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dppp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dpds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dpdp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_ddss(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_ddps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_ddpp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_ddds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dddp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_os_dddd(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_xxxx(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] );

#endif
