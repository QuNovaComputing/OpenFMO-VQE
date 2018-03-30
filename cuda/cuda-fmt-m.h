#ifndef _CUDA_FMT_M_H_
#define _CUDA_FMT_M_H_

#ifdef __cplusplus
extern "C" {
#endif

size_t cuda_FMT_m_get_size(const int m);
int cuda_FMT_m_Init(double **dtmp, const size_t *mtmp,
        const int *ndivs, const int mmax);
int cuda_FMT_m_Finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* _CUDA_FMT_M_H_ */
