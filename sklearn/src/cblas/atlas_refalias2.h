#ifndef ATLAS_REFALIAS2_H
#define ATLAS_REFALIAS2_H
/*
 * Real BLAS
 */
   #define ATL_sger2   ATL_srefger2
   #define ATL_sspr2   ATL_srefspr2
   #define ATL_ssyr2   ATL_srefsyr2
   #define ATL_sspr    ATL_srefspr
   #define ATL_ssyr    ATL_srefsyr
   #define ATL_sger    ATL_srefger
   #define ATL_stpsv   ATL_sreftpsv
   #define ATL_stbsv   ATL_sreftbsv
   #define ATL_strsv   ATL_sreftrsv
   #define ATL_stpmv   ATL_sreftpmv
   #define ATL_stbmv   ATL_sreftbmv
   #define ATL_strmv   ATL_sreftrmv
   #define ATL_sspmv   ATL_srefspmv
   #define ATL_ssbmv   ATL_srefsbmv
   #define ATL_ssymv   ATL_srefsymv
   #define ATL_sgbmv   ATL_srefgbmv
   #define ATL_sgemv   ATL_srefgemv

   #define ATL_dger2   ATL_drefger2
   #define ATL_dspr2   ATL_drefspr2
   #define ATL_dsyr2   ATL_drefsyr2
   #define ATL_dspr    ATL_drefspr
   #define ATL_dsyr    ATL_drefsyr
   #define ATL_dger    ATL_drefger
   #define ATL_dtpsv   ATL_dreftpsv
   #define ATL_dtbsv   ATL_dreftbsv
   #define ATL_dtrsv   ATL_dreftrsv
   #define ATL_dtpmv   ATL_dreftpmv
   #define ATL_dtbmv   ATL_dreftbmv
   #define ATL_dtrmv   ATL_dreftrmv
   #define ATL_dspmv   ATL_drefspmv
   #define ATL_dsbmv   ATL_drefsbmv
   #define ATL_dsymv   ATL_drefsymv
   #define ATL_dgbmv   ATL_drefgbmv
   #define ATL_dgemv   ATL_drefgemv

/*
 * Complex BLAS
 */
   #define ATL_cger2c    ATL_crefger2c
   #define ATL_cger2u    ATL_crefger2u
   #define ATL_chpr2     ATL_crefhpr2
   #define ATL_cher2     ATL_crefher2
   #define ATL_chpr      ATL_crefhpr
   #define ATL_cher      ATL_crefher
   #define ATL_cgerc     ATL_crefgerc
   #define ATL_cgeru     ATL_crefgeru
   #define ATL_ctpsv     ATL_creftpsv
   #define ATL_ctbsv     ATL_creftbsv
   #define ATL_ctrsv     ATL_creftrsv
   #define ATL_ctpmv     ATL_creftpmv
   #define ATL_ctbmv     ATL_creftbmv
   #define ATL_ctrmv     ATL_creftrmv
   #define ATL_chpmv     ATL_crefhpmv
   #define ATL_chbmv     ATL_crefhbmv
   #define ATL_chemv     ATL_crefhemv
   #define ATL_cgbmv     ATL_crefgbmv
   #define ATL_cgemv     ATL_crefgemv

   #define ATL_zger2c    ATL_zrefger2c
   #define ATL_zger2u    ATL_zrefger2u
   #define ATL_zhpr2     ATL_zrefhpr2
   #define ATL_zher2     ATL_zrefher2
   #define ATL_zhpr      ATL_zrefhpr
   #define ATL_zher      ATL_zrefher
   #define ATL_zgerc     ATL_zrefgerc
   #define ATL_zgeru     ATL_zrefgeru
   #define ATL_ztpsv     ATL_zreftpsv
   #define ATL_ztbsv     ATL_zreftbsv
   #define ATL_ztrsv     ATL_zreftrsv
   #define ATL_ztpmv     ATL_zreftpmv
   #define ATL_ztbmv     ATL_zreftbmv
   #define ATL_ztrmv     ATL_zreftrmv
   #define ATL_zhpmv     ATL_zrefhpmv
   #define ATL_zhbmv     ATL_zrefhbmv
   #define ATL_zhemv     ATL_zrefhemv
   #define ATL_zgbmv     ATL_zrefgbmv
   #define ATL_zgemv     ATL_zrefgemv

#endif
