#ifndef _OFMO_RYS_XXXX_H
#define _OFMO_RYS_XXXX_H

// for partiall direct SCF

extern int ofmo_twoint_rys_ssss(
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

extern int ofmo_twoint_rys_psss(
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

extern int ofmo_twoint_rys_psps(
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

extern int ofmo_twoint_rys_ppss(
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

extern int ofmo_twoint_rys_ppps(
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

extern int ofmo_twoint_rys_pppp(
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

extern int ofmo_twoint_rys_dsss(
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

extern int ofmo_twoint_rys_dsps(
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

extern int ofmo_twoint_rys_dspp(
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

extern int ofmo_twoint_rys_dsds(
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

extern int ofmo_twoint_rys_dpss(
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

extern int ofmo_twoint_rys_dpps(
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

extern int ofmo_twoint_rys_dppp(
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

extern int ofmo_twoint_rys_dpds(
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

extern int ofmo_twoint_rys_dpdp(
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

extern int ofmo_twoint_rys_ddss(
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

extern int ofmo_twoint_rys_ddps(
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

extern int ofmo_twoint_rys_ddpp(
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

extern int ofmo_twoint_rys_ddds(
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

extern int ofmo_twoint_rys_dddp(
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

extern int ofmo_twoint_rys_dddd(
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

extern int ofmo_twoint_direct_rys_ssss(
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

extern int ofmo_twoint_direct_rys_psss(
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

extern int ofmo_twoint_direct_rys_psps(
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

extern int ofmo_twoint_direct_rys_ppss(
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

extern int ofmo_twoint_direct_rys_ppps(
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

extern int ofmo_twoint_direct_rys_pppp(
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

extern int ofmo_twoint_direct_rys_dsss(
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

extern int ofmo_twoint_direct_rys_dsps(
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

extern int ofmo_twoint_direct_rys_dspp(
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

extern int ofmo_twoint_direct_rys_dsds(
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

extern int ofmo_twoint_direct_rys_dpss(
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

extern int ofmo_twoint_direct_rys_dpps(
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

extern int ofmo_twoint_direct_rys_dppp(
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

extern int ofmo_twoint_direct_rys_dpds(
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

extern int ofmo_twoint_direct_rys_dpdp(
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

extern int ofmo_twoint_direct_rys_ddss(
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

extern int ofmo_twoint_direct_rys_ddps(
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

extern int ofmo_twoint_direct_rys_ddpp(
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

extern int ofmo_twoint_direct_rys_ddds(
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

extern int ofmo_twoint_direct_rys_dddp(
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

extern int ofmo_twoint_direct_rys_dddd(
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

extern int ofmo_twoint_direct_rys_xxxx(
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
