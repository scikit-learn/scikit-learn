#define ATLAS_PTALIAS1_H /* no threaded routs for Level 1 yet */
#ifndef ATLAS_PTALIAS1_H
#define ATLAS_PTALIAS1_H
/*
 * Real BLAS
 */
   #define ATL_dsdot   ATL_dstdot
   #define ATL_sdsdot  ATL_sdstdot
   #define ATL_sasum   ATL_stasum
   #define ATL_snrm2   ATL_stnrm2
   #define ATL_sdot    ATL_stdot
   #define ATL_saxpy   ATL_staxpy
   #define ATL_scopy   ATL_stcopy
   #define ATL_sscal   ATL_stscal
   #define ATL_sswap   ATL_stswap
   #define ATL_srotm   ATL_strotm
   #define ATL_srot    ATL_strot
   #define ATL_srotmg  ATL_strotmg
   #define ATL_srotg   ATL_strotg
   #define ATL_isamax  ATL_istamax

   #define ATL_dasum   ATL_dtasum
   #define ATL_dnrm2   ATL_dtnrm2
   #define ATL_ddot    ATL_dtdot
   #define ATL_daxpy   ATL_dtaxpy
   #define ATL_dcopy   ATL_dtcopy
   #define ATL_dscal   ATL_dtscal
   #define ATL_dswap   ATL_dtswap
   #define ATL_drotm   ATL_dtrotm
   #define ATL_drot    ATL_dtrot
   #define ATL_drotmg  ATL_dtrotmg
   #define ATL_drotg   ATL_dtrotg
   #define ATL_idamax  ATL_idtamax

/*
 * Complex BLAS
 */
   #define ATL_cdotc_sub ATL_ctdotc_sub
   #define ATL_cdotu_sub ATL_ctdotu_sub
   #define ATL_caxpy     ATL_ctaxpy
   #define ATL_ccopy     ATL_ctcopy
   #define ATL_cscal     ATL_ctscal
   #define ATL_cswap     ATL_ctswap
   #define ATL_icamax    ATL_ictamax
   #define ATL_csscal    ATL_cstscal
   #define ATL_scnrm2    ATL_sctnrm2
   #define ATL_scasum    ATL_sctasum

   #define ATL_zdotc_sub ATL_ztdotc_sub
   #define ATL_zdotu_sub ATL_ztdotu_sub
   #define ATL_zaxpy     ATL_ztaxpy
   #define ATL_zcopy     ATL_ztcopy
   #define ATL_zscal     ATL_ztscal
   #define ATL_zswap     ATL_ztswap
   #define ATL_izamax    ATL_iztamax
   #define ATL_zdscal    ATL_zdtscal
   #define ATL_dznrm2    ATL_dztnrm2
   #define ATL_dzasum    ATL_dztasum

#endif
