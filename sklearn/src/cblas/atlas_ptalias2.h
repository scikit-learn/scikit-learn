#ifndef ATLAS_PTALIAS2_H
#define ATLAS_PTALIAS2_H
/*
 * Real BLAS
 */
   #define ATL_sger    ATL_stger
   #define ATL_sgemv   ATL_stgemv

   #define ATL_dger    ATL_dtger
   #define ATL_dgemv   ATL_dtgemv

/*
 * Complex BLAS
 */
   #define ATL_cgemv     ATL_ctgemv
   #define ATL_cgerc     ATL_ctgerc
   #define ATL_cgeru     ATL_ctgeru

   #define ATL_zgemv     ATL_ztgemv
   #define ATL_zgerc     ATL_ztgerc
   #define ATL_zgeru     ATL_ztgeru

#endif
