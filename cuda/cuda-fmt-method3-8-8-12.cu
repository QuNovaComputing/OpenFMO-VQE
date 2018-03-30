extern __shared__ double shared[];

#define DINV2_0 (1.0)
#define DINV2_1 (1.0/(2*1+1))
#define DINV2_2 (1.0/(2*2+1))
#define DINV2_3 (1.0/(2*3+1))
#define DINV2_4 (1.0/(2*4+1))
#define DINV2_5 (1.0/(2*5+1))
#define DINV2_6 (1.0/(2*6+1))
#define DINV2_7 (1.0/(2*7+1))
#define DINV2_8 (1.0/(2*8+1))

#define DINV_0  (1.0)
#define DINV_1  (1.0)
#define DINV_2  (1.0/2)
#define DINV_3  (1.0/3)
#define DINV_4  (1.0/4)
#define DINV_5  (1.0/5)
#define DINV_6  (1.0/6)
#define DINV_7  (1.0/7)
#define DINV_8  (1.0/8)
#define DINV_9  (1.0/9)

// automatically generated by gen_fmt_method3 function
// m = 0, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt0_method3( const double t, const double coef, double fmt[] ) {
    if ( t >= 26 ) {
        fmt[0] = coef * LDG(FMT_m_sqrt_pi_2) * sqrt(1.e0/t);
    } else {
      int it;
      //double t0, dt, *f, expt0, t2, ff[8];
      double t0, dt, *f, expt0, t2;
      //double *tbl  = (double *)&shared[0];
      double delta0, dhalf0;
      dhalf0 = (delta0 = LDG(FMT_m_delta[0])) * 0.5;
        it = (int)(t*4);
        t0 = delta0 * (double)it + dhalf0;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt0_table[it*2];
        //f  = &FMT_m_table[0][it*2];
        //f  = &FMT_m_table0[it*2];
#ifndef CUDA_FMT_M_SM
        f  = &FMT_m_table0[it*2];
#else
        //f  = (double *)&shared[it*2];
        f  = &shared[it*2];
#endif
#if 1
        double ff[8];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // F(m+k)(T0) (k=0, n-2)
#if 1
        ff[6] = LDG(FMT_m_dinv2[0+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[0+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[0+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[0+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[0+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[0+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[0+0]) * ( expt0 + t2*ff[1] );
#else
        ff[6] = DINV2_6 * ( expt0 + t2*ff[7] );
        ff[5] = DINV2_5 * ( expt0 + t2*ff[6] );
        ff[4] = DINV2_4 * ( expt0 + t2*ff[5] );
        ff[3] = DINV2_3 * ( expt0 + t2*ff[4] );
        ff[2] = DINV2_2 * ( expt0 + t2*ff[3] );
        ff[1] = DINV2_1 * ( expt0 + t2*ff[2] );
        ff[0] = DINV2_0 * ( expt0 + t2*ff[1] );
#endif
        // F[0]
        fmt[0] = ff[0] + dt*( ff[1] + LDG(FMT_m_dinv[2])*dt*( ff[2] + LDG(FMT_m_dinv[3])*dt*
                ( ff[3] + LDG(FMT_m_dinv[4])*dt*( ff[4] + LDG(FMT_m_dinv[5])*dt*
                ( ff[5] + LDG(FMT_m_dinv[6])*dt*( ff[6] + LDG(FMT_m_dinv[7])*dt*ff[7]))))));
        fmt[0] *= coef;
#else
        double ff;
        ff = LDM(f[0]);
        expt0       = LDM(f[1]);
        t0 = LDG(FMT_m_dinv[7])*dt*ff;
        for (int k=6; k>=2; k--) {
          ff = LDG(FMT_m_dinv2[0+k]) * ( expt0 + t2*ff );
          t0 = LDG(FMT_m_dinv[k])*dt*(ff+t0);
        }
        ff = LDG(FMT_m_dinv2[0+1]) * ( expt0 + t2*ff );
        t0 = dt*(ff+t0);
        //ff = LDG(FMT_m_dinv2[0+0]) * ( expt0 + t2*ff );
        ff = ( expt0 + t2*ff );
        t0 += ff;
        fmt[0] = t0 * coef;
#endif
    }
}

// automatically generated by gen_fmt_method3 function
// m = 1, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt1_method3( const double t, const double coef, double fmt[] ) {
    //int it, k;
    //double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 30 ) {
      double t2, sqrt2;
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =        t2 * fmt[0];
    } else {
      //int it, k;
      int it;
      double t0, dt, *f, expt0, t2, expt;
      double delta1 = LDG(FMT_m_delta[1]);
      double dhalf1 = delta1 * 0.5;
      //dhalf1 = (delta1 = LDG(FMT_m_delta[1])) * 0.5;
        it = (int)(t*4);
        t0 = delta1 * (double)it + dhalf1;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt1_table[it*2];
        //f  = &FMT_m_table[1][it*2];
        //f  = &FMT_m_table1[it*2];
#ifndef CUDA_FMT_M_SM
        f  = &FMT_m_table1[it*2];
#else
        f  = (double *)&shared[it*2];
#endif
#ifdef FERMI
        double ff[8];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
#ifndef FERMI
        double dtk[10];
        // (-dt)^k/k!
        dtk[1] = dt;
        for (int k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
#endif
        // F(m+k)(T0) (k=0, n-2)
#if 1
        ff[6] = LDG(FMT_m_dinv2[1+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[1+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[1+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[1+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[1+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[1+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[1+0]) * ( expt0 + t2*ff[1] );
#else
        ff[6] = DINV2_7 * ( expt0 + t2*ff[7] );
        ff[5] = DINV2_6 * ( expt0 + t2*ff[6] );
        ff[4] = DINV2_5 * ( expt0 + t2*ff[5] );
        ff[3] = DINV2_4 * ( expt0 + t2*ff[4] );
        ff[2] = DINV2_3 * ( expt0 + t2*ff[3] );
        ff[1] = DINV2_2 * ( expt0 + t2*ff[2] );
        ff[0] = DINV2_1 * ( expt0 + t2*ff[1] );
#endif
        // F[1]
#ifndef FERMI
        fmt[1] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[1] *= coef;
        // F[0]-F[0]
        t2 = t + t;
        fmt[0] = coef*expt + t2*fmt[1];
#else
        t0 = ff[0];
        t2 = dt;
        t0 += t2 * ff[1];
        expt = 1.e0 + t2;
#if 1
#pragma unroll
        for (int k=2; k<8; k++ ) {
          t2 = dt*LDG(FMT_m_dinv[k])*t2;
          t0 += t2 * ff[k];
          expt += t2;
        }
        expt += dt*LDG(FMT_m_dinv[8])*t2 * (1+dt*LDG(FMT_m_dinv[9]));
#else
        t2 = dt*t2*DINV_2;
        t0 += t2 * ff[2];
        expt += t2;
        t2 = dt*t2*DINV_3;
        t0 += t2 * ff[3];
        expt += t2;
        t2 = dt*t2*DINV_4;
        t0 += t2 * ff[4];
        expt += t2;
        t2 = dt*t2*DINV_5;
        t0 += t2 * ff[5];
        expt += t2;
        t2 = dt*t2*DINV_6;
        t0 += t2 * ff[6];
        expt += t2;
        t2 = dt*t2*DINV_7;
        t0 += t2 * ff[7];
        expt += t2;
        expt += dt*t2*DINV_8 * (1+dt*DINV_9);
#endif
        expt *= expt0;
        t0 *= coef;
        fmt[1] = t0;
        // F[0]-F[0]
        t2 = t + t;
        //fmt[0] = coef*expt + t2*fmt[1];
        fmt[0] = coef*expt + t2*t0;
#endif
#else // ----------------------
        //double dtk[8];
        double dtk[7];
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        //dtk[1] = dt;
        dtk[0] = dt;
#pragma unroll
        //for (int k=2; k<=7; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        for (int k=1; k<=6; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k+1])*dtk[k-1];
        //expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6] + dtk[7];
        expt = 1.e0 + dtk[0] + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6];
        //t0 = dt*LDG(FMT_m_dinv[8])*dtk[7];
        t0 = dt*LDG(FMT_m_dinv[8])*dtk[7-1];
        expt += t0 * (1e0 + dt*LDG(FMT_m_dinv[9]));
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        dt = LDM(f[0]);
        //t0 = dtk[7]*dt;
        t0 = dtk[7-1]*dt;
#pragma unroll
        for (int k=6; k>0; k--) {
          dt = LDG(FMT_m_dinv2[1+k]) * ( expt0 + t2*dt );
          //t0 += dtk[k]*dt;
          t0 += dtk[k-1]*dt;
        }
        t0 += LDG(FMT_m_dinv2[1]) * ( expt0 + t2*dt );
        t0 *= coef;
        fmt[1] = t0;
        // F[0]-F[0]
        t2 = t + t;
        fmt[0] = coef*expt + t2*t0;
#endif
    }
}

// automatically generated by gen_fmt_method3 function
// m = 2, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt2_method3( const double t, const double coef, double fmt[] ) {
    //int it, k;
    //double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 33 ) {
      double t2, sqrt2;
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
    } else {
      //int it, k;
      int it;
      //double t0, dt, *f, expt0, t2, ff[8], dtk[10], expt;
      //double t0, dt, *f, expt0, t2, ff[8],  expt;
      double t0, dt, *f, expt0, t2,  expt;
      double delta2, dhalf2;
      dhalf2 = (delta2 = LDG(FMT_m_delta[2])) * 0.5;
        it = (int)(t*4);
        t0 = delta2 * (double)it + dhalf2;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt2_table[it*2];
        //f  = &FMT_m_table[2][it*2];
#ifndef CUDA_FMT_M_SM
        f  = &FMT_m_table2[it*2];
#else
        f  = (double *)&shared[it*2];
#endif
#if 0
        double ff[8];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[2+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[2+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[2+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[2+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[2+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[2+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[2+0]) * ( expt0 + t2*ff[1] );
        // (-dt)^k/k!
        // F[2]
        t2 = dt; // dtk[1]
        t0 = ff[0] + t2*ff[1];
        expt = 1.e0 + t2;
#pragma unroll
        for (int k=2; k<=7; k++) {
          t2 = dt*LDG(FMT_m_dinv[k])*t2;
          t0 += t2*ff[k];
          expt += t2;
        }
        t0 *= coef;
        fmt[2] = t0;
        expt += dt*LDG(FMT_m_dinv[8])*t2 * (1 + dt*LDG(FMT_m_dinv[9]));
        expt *= expt0;
        // F[1]-F[0]
        t2 = t + t;
        expt      *= coef;
        //fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        //fmt[0] =              expt + t2 * fmt[1];
        t0 = LDG(FMT_m_dinv2[1]) * ( expt + t2 * t0 );
        fmt[1] = t0;
        fmt[0] =              expt + t2 * t0;
#else // -----------------------------
        //double dtk[8];
        double dtk[7];
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        //dtk[1] = dt;
        dtk[0] = dt;
#pragma unroll
        //for (int k=2; k<8; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        for (int k=1; k<7; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k+1])*dtk[k-1];
        // exp(-T)
        //expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6] + dtk[7];
        expt = 1.e0 + dtk[0] + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6];
        //expt += dt*LDG(FMT_m_dinv[8])*dtk[7]*(1+dt*LDG(FMT_m_dinv[9]));
        expt += dt*LDG(FMT_m_dinv[8])*dtk[7-1]*(1+dt*LDG(FMT_m_dinv[9]));
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        dt = LDM(f[0]); // ff[7]
        //t0 = dtk[7]*dt;
        t0 = dtk[7-1]*dt;
#pragma unroll
        for (int k=6; k>=1; k--) {
          dt = LDG(FMT_m_dinv2[2+k]) * ( expt0 + t2*dt );
          //t0 += dtk[k]*dt;
          t0 += dtk[k-1]*dt;
        }
        dt = LDG(FMT_m_dinv2[2+0]) * ( expt0 + t2*dt );
        t0 = (t0 + dt) * coef;
        fmt[2] = t0;
        // F[1]-F[0]
        t2 = t + t;
        expt      *= coef;
        //fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        //fmt[0] =              expt + t2 * fmt[1];
        t0 = LDG(FMT_m_dinv2[1]) * ( expt + t2 * t0 );
        fmt[1] = t0;
        fmt[0] =              expt + t2 * t0;
#endif
        /* -------------- ORIG 
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[2+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[2+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[2+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[2+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[2+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[2+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[2+0]) * ( expt0 + t2*ff[1] );
        // F[2]
        fmt[2] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[2] *= coef;
        // F[1]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
        ----------------- ORIG */
    }
}

// automatically generated by gen_fmt_method3 function
// m = 3, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt3_method3( const double t, const double coef, double fmt[] ) {
    //int it, k;
    //double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 36 ) {
      double t2, sqrt2;
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
    } else {
      int it;
      double t0, dt, *f, expt0, t2, expt;
      double delta3, dhalf3;
      dhalf3 = (delta3 = LDG(FMT_m_delta[3])) * 0.5;
        it = (int)(t*4);
        t0 = delta3 * (double)it + dhalf3;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt3_table[it*2];
        //f  = &FMT_m_table[3][it*2];
#ifndef CUDA_FMT_M_SM
        f  = &FMT_m_table3[it*2];
#else
        f  = (double *)&shared[it*2];
#endif
#if 0
      double ff[8];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[3+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[3+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[3+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[3+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[3+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[3+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[3+0]) * ( expt0 + t2*ff[1] );
        // (-dt)^k/k!
        t2 = dt; // dtk[1]
        t0 = ff[0] + t2*ff[1];
        expt = 1.e0 + t2;
#pragma unroll
        for (int k=2; k<=7; k++ ) {
          t2 = dt*LDG(FMT_m_dinv[k])*t2;
          t0 += t2*ff[k];
          expt += t2;
        }
        // exp(-T)
        expt += dt*LDG(FMT_m_dinv[8])*t2*(1+ dt*LDG(FMT_m_dinv[9]));
        expt *= expt0;
        t0 *= coef;
        fmt[3] = t0;
        // F[3]
        // F[2]-F[0]
        t2 = t + t;
        expt      *= coef;
        //fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * t0 );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
#else
      //double dtk[8];
      double dtk[7];
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        //dtk[1] = dt;
        dtk[0] = dt;
#pragma unroll
        //for (int k=2; k<8; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        for (int k=1; k<7; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k+1])*dtk[k-1];
        // exp(-T)
        //expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6] + dtk[7];
        expt = 1.e0 + dtk[0] + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6];
        //expt += dt*LDG(FMT_m_dinv[8])*dtk[8-1]*(1+dt*LDG(FMT_m_dinv[9]));
        expt += dt*LDG(FMT_m_dinv[8])*dtk[6]*(1+dt*LDG(FMT_m_dinv[9]));
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        dt = LDM(f[0]); // ff[7]
        //t0 = dtk[7]*dt;
        t0 = dtk[6]*dt;
#pragma unroll
        for (int k=6; k>0; k--) {
          dt = LDG(FMT_m_dinv2[3+k]) * ( expt0 + t2*dt );
          //t0 += dtk[k]*dt;
          t0 += dtk[k-1]*dt;
        }
        t0 += LDG(FMT_m_dinv2[3+0]) * ( expt0 + t2*dt );
        t0 *= coef;
        // F[3]
        fmt[3] = t0;
        // F[2]-F[0]
        t2 = t + t;
        expt      *= coef;
        //fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * t0 );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
#endif
        /* -------------- ORIG 
      double ff[8], dtk[10];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        dtk[1] = dt;
        for (int k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[3+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[3+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[3+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[3+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[3+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[3+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[3+0]) * ( expt0 + t2*ff[1] );
        // F[3]
        fmt[3] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[3] *= coef;
        // F[2]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
        ----------------- ORIG */
    }
}

// automatically generated by gen_fmt_method3 function
// m = 4, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt4_method3( const double t, const double coef, double fmt[] ) {
    //int it, k;
    //double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 39 ) {
      double t2, sqrt2;
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
        fmt[4] = 7.e0 * t2 * fmt[3];
    } else {
      int it;
      double t0, dt, *f, expt0, t2, expt;
      double delta4, dhalf4;
      dhalf4 = (delta4 = LDG(FMT_m_delta[4])) * 0.5;
        it = (int)(t*4);
        t0 = delta4 * (double)it + dhalf4;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt4_table[it*2];
        //f  = &FMT_m_table[4][it*2];
#ifndef CUDA_FMT_M_SM
        f  = &FMT_m_table4[it*2];
#else
        f  = (double *)&shared[it*2];
#endif
#if 0
      double ff[8];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[4+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[4+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[4+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[4+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[4+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[4+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[4+0]) * ( expt0 + t2*ff[1] );
        // (-dt)^k/k!
        t2 = dt; // dtk[1]
        t0 = ff[0] + t2*ff[1];
        expt = 1.e0 + t2;
#pragma unroll
        for (int k=2; k<=7; k++ ) {
          t2 = dt*LDG(FMT_m_dinv[k])*t2;
          t0 += t2*ff[k];
          expt += t2;
        }
        // exp(-T)
        expt += dt*LDG(FMT_m_dinv[8])*t2*(1+ dt*LDG(FMT_m_dinv[9]));
        expt *= expt0;
        t0 *= coef;
        // F[4]
        fmt[4] = t0;
        // F[3]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[3] = LDG(FMT_m_dinv2[3]) * ( expt + t2 * fmt[4] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
#else
      //double dtk[8];
      double dtk[7];
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        //dtk[1] = dt;
        dtk[0] = dt;
#pragma unroll
        //for (int k=2; k<8; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        for (int k=1; k<7; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k+1])*dtk[k-1];
        // exp(-T)
        //expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6] + dtk[7];
        expt = 1.e0 + dtk[0] + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6];
        //expt += dt*LDG(FMT_m_dinv[8])*dtk[8-1]*(1+dt*LDG(FMT_m_dinv[9]));
        expt += dt*LDG(FMT_m_dinv[8])*dtk[6]*(1+dt*LDG(FMT_m_dinv[9]));
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        dt = LDM(f[0]); // ff[7]
        //t0 = dtk[7]*dt;
        t0 = dtk[6]*dt;
#pragma unroll
        for (int k=6; k>0; k--) {
          dt = LDG(FMT_m_dinv2[4+k]) * ( expt0 + t2*dt );
          //t0 += dtk[k]*dt;
          t0 += dtk[k-1]*dt;
        }
        t0 += LDG(FMT_m_dinv2[4+0]) * ( expt0 + t2*dt );
        t0 *= coef;
        fmt[4] = t0;
        // F[3]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[3] = LDG(FMT_m_dinv2[3]) * ( expt + t2 * fmt[4] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
#endif
        /* -------------- ORIG 
      double ff[8], dtk[10];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        dtk[1] = dt;
        for (int k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[4+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[4+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[4+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[4+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[4+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[4+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[4+0]) * ( expt0 + t2*ff[1] );
        // F[4]
        fmt[4] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[4] *= coef;
        // F[3]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[3] = LDG(FMT_m_dinv2[3]) * ( expt + t2 * fmt[4] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
        ----------------- ORIG */
    }
}

// automatically generated by gen_fmt_method3 function
// m = 5, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt5_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 41 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
    } else {
      double delta5, dhalf5;
      dhalf5 = (delta5 = LDG(FMT_m_delta[5])) * 0.5;
        it = (int)(t*4);
        t0 = delta5 * (double)it + dhalf5;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt5_table[it*2];
        f  = &FMT_m_table[5][it*2];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[5+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[5+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[5+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[5+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[5+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[5+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[5+0]) * ( expt0 + t2*ff[1] );
        // F[5]
        fmt[5] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[5] *= coef;
        // F[4]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[4] = LDG(FMT_m_dinv2[4]) * ( expt + t2 * fmt[5] );
        fmt[3] = LDG(FMT_m_dinv2[3]) * ( expt + t2 * fmt[4] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 6, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt6_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 43 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
    } else {
      double delta6, dhalf6;
      dhalf6 = (delta6 = LDG(FMT_m_delta[6])) * 0.5;
        it = (int)(t*4);
        t0 = delta6 * (double)it + dhalf6;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt6_table[it*2];
        f  = &FMT_m_table[6][it*2];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[6+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[6+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[6+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[6+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[6+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[6+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[6+0]) * ( expt0 + t2*ff[1] );
        // F[6]
        fmt[6] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[6] *= coef;
        // F[5]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[5] = LDG(FMT_m_dinv2[5]) * ( expt + t2 * fmt[6] );
        fmt[4] = LDG(FMT_m_dinv2[4]) * ( expt + t2 * fmt[5] );
        fmt[3] = LDG(FMT_m_dinv2[3]) * ( expt + t2 * fmt[4] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 7, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt7_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 45 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
        fmt[7] = 13.e0 * t2 * fmt[6];
    } else {
      double delta7, dhalf7;
      dhalf7 = (delta7 = LDG(FMT_m_delta[7])) * 0.5;
        it = (int)(t*4);
        t0 = delta7 * (double)it + dhalf7;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt7_table[it*2];
        f  = &FMT_m_table[7][it*2];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[7+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[7+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[7+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[7+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[7+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[7+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[7+0]) * ( expt0 + t2*ff[1] );
        // F[7]
        fmt[7] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[7] *= coef;
        // F[6]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[6] = LDG(FMT_m_dinv2[6]) * ( expt + t2 * fmt[7] );
        fmt[5] = LDG(FMT_m_dinv2[5]) * ( expt + t2 * fmt[6] );
        fmt[4] = LDG(FMT_m_dinv2[4]) * ( expt + t2 * fmt[5] );
        fmt[3] = LDG(FMT_m_dinv2[3]) * ( expt + t2 * fmt[4] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 8, nexp = 8, eps = 1.0e-12
__device__ void gpu_fmt8_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[8], sqrt2, dtk[10], expt;
    if ( t >= 48 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * LDG(FMT_m_sqrt_pi2) * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
        fmt[7] = 13.e0 * t2 * fmt[6];
        fmt[8] = 15.e0 * t2 * fmt[7];
    } else {
      double delta8, dhalf8;
      dhalf8 = (delta8 = LDG(FMT_m_delta[8])) * 0.5;
        it = (int)(t*4);
        t0 = delta8 * (double)it + dhalf8;
        t2 = t0 + t0;
        dt = t0 - t;
        //f  = &fmt8_table[it*2];
        f  = &FMT_m_table[8][it*2];
        ff[8-1] = LDM(f[0]);
        expt0       = LDM(f[1]);
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<10; k++ ) dtk[k] = dt*LDG(FMT_m_dinv[k])*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7] + dtk[8] + dtk[9];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[6] = LDG(FMT_m_dinv2[8+6]) * ( expt0 + t2*ff[7] );
        ff[5] = LDG(FMT_m_dinv2[8+5]) * ( expt0 + t2*ff[6] );
        ff[4] = LDG(FMT_m_dinv2[8+4]) * ( expt0 + t2*ff[5] );
        ff[3] = LDG(FMT_m_dinv2[8+3]) * ( expt0 + t2*ff[4] );
        ff[2] = LDG(FMT_m_dinv2[8+2]) * ( expt0 + t2*ff[3] );
        ff[1] = LDG(FMT_m_dinv2[8+1]) * ( expt0 + t2*ff[2] );
        ff[0] = LDG(FMT_m_dinv2[8+0]) * ( expt0 + t2*ff[1] );
        // F[8]
        fmt[8] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5] + dtk[6]*ff[6]
                 + dtk[7]*ff[7];
        fmt[8] *= coef;
        // F[7]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[7] = LDG(FMT_m_dinv2[7]) * ( expt + t2 * fmt[8] );
        fmt[6] = LDG(FMT_m_dinv2[6]) * ( expt + t2 * fmt[7] );
        fmt[5] = LDG(FMT_m_dinv2[5]) * ( expt + t2 * fmt[6] );
        fmt[4] = LDG(FMT_m_dinv2[4]) * ( expt + t2 * fmt[5] );
        fmt[3] = LDG(FMT_m_dinv2[3]) * ( expt + t2 * fmt[4] );
        fmt[2] = LDG(FMT_m_dinv2[2]) * ( expt + t2 * fmt[3] );
        fmt[1] = LDG(FMT_m_dinv2[1]) * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}
