#define ASSERT(x) assert(x)

#define CU_CHECK(stmt)                          \
  do                                            \
    {                                           \
      CUresult status = (stmt);                 \
      if (CUDA_SUCCESS != status) {             \
        printf("status = %d\n", status);        \
      }                                         \
      ASSERT(CUDA_SUCCESS == status);           \
    } while (0)

#define checkCudaErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors(cudaError err, const char *file, const int line)
{
    if (cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",
                file, line, (int)err, cudaGetErrorString(err));
        exit(-1);
    }
}
