#ifndef ATLAS_REFALIAS1_H
#define ATLAS_REFALIAS1_H
/*
 * Real BLAS
 */
   #define ATL_dsdot   ATL_dsrefdot
   #define ATL_sdsdot  ATL_sdsrefdot
   #define ATL_sasum   ATL_srefasum
   #define ATL_snrm2   ATL_srefnrm2
   #define ATL_sdot    ATL_srefdot
   #define ATL_saxpy   ATL_srefaxpy
   #define ATL_scopy   ATL_srefcopy
   #define ATL_sscal   ATL_srefscal
   #define ATL_sswap   ATL_srefswap
   #define ATL_srotm   ATL_srefrotm
   #define ATL_srot    ATL_srefrot
   #define ATL_srotmg  ATL_srefrotmg
   #define ATL_srotg   ATL_srefrotg
   #define ATL_isamax  ATL_isrefamax

   #define ATL_dasum   ATL_drefasum
   #define ATL_dnrm2   ATL_drefnrm2
   #define ATL_ddot    ATL_drefdot
   #define ATL_daxpy   ATL_drefaxpy
   #define ATL_dcopy   ATL_drefcopy
   #define ATL_dscal   ATL_drefscal
   #define ATL_dswap   ATL_drefswap
   #define ATL_drotm   ATL_drefrotm
   #define ATL_drot    ATL_drefrot
   #define ATL_drotmg  ATL_drefrotmg
   #define ATL_drotg   ATL_drefrotg
   #define ATL_idamax  ATL_idrefamax

/*
 * Complex BLAS
 */
   #define ATL_cdotc_sub ATL_crefdotc_sub
   #define ATL_cdotu_sub ATL_crefdotu_sub
   #define ATL_caxpy     ATL_crefaxpy
   #define ATL_ccopy     ATL_crefcopy
   #define ATL_cscal     ATL_crefscal
   #define ATL_cswap     ATL_crefswap
   #define ATL_icamax    ATL_icrefamax
   #define ATL_csscal    ATL_csrefscal
   #define ATL_scnrm2    ATL_screfnrm2
   #define ATL_scasum    ATL_screfasum

   #define ATL_zdotc_sub ATL_zrefdotc_sub
   #define ATL_zdotu_sub ATL_zrefdotu_sub
   #define ATL_zaxpy     ATL_zrefaxpy
   #define ATL_zcopy     ATL_zrefcopy
   #define ATL_zscal     ATL_zrefscal
   #define ATL_zswap     ATL_zrefswap
   #define ATL_izamax    ATL_izrefamax
   #define ATL_zdscal    ATL_zdrefscal
   #define ATL_dznrm2    ATL_dzrefnrm2
   #define ATL_dzasum    ATL_dzrefasum

#endif
