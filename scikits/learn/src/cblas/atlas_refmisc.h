/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.2 -- December 25, 2000
 *
 * -- Suggestions,  comments,  bugs reports should be sent to the follo-
 *    wing e-mail address: atlas@cs.utk.edu
 *
 * Author         : Antoine P. Petitet
 * University of Tennessee - Innovative Computing Laboratory
 * Knoxville TN, 37996-1301, USA.
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the ATLAS group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------
 */
#ifndef ATL_REFMISC_H
#define ATL_REFMISC_H
/*
 * =====================================================================
 * Include files
 * =====================================================================
 */
#include <math.h>
#include "atlas_enum.h"
/*
 * =====================================================================
 * #define macro constants
 * =====================================================================
 */
#define    ATL_sNONE                     (-1.0f)
#define    ATL_sNTWO                     (-2.0f)
#define    ATL_sONE                      ( 1.0f)
#define    ATL_sZERO                     ( 0.0f)

#define    ATL_dNONE                     (-1.0)
#define    ATL_dNTWO                     (-2.0)
#define    ATL_dONE                      ( 1.0)
#define    ATL_dZERO                     ( 0.0)
/*
 * =====================================================================
 * # macro functions
 * =====================================================================
 */
#define    Msabs( a_ ) ( ( (a_) < ATL_sZERO ) ? -(a_) : (a_) )

#define    Mszero( a_r_, a_i_ )                                        \
           ( ( (a_r_) == ATL_sZERO ) && ( (a_i_) == ATL_sZERO ) )

#define    Msone( a_r_,  a_i_ )                                        \
           ( ( (a_r_) == ATL_sONE  ) && ( (a_i_) == ATL_sZERO ) )

#define    Msscl( a_r_, a_i_, c_r_, c_i_ )                             \
           {                                                           \
              register float tmp_r_, tmp_i_;                           \
              tmp_r_ = (a_r_) * c_r_ - (a_i_) * c_i_;                  \
              tmp_i_ = (a_r_) * c_i_ + (a_i_) * c_r_;                  \
              c_r_   = tmp_r_;                                         \
              c_i_   = tmp_i_;                                         \
           }
/*
 * Msdiv performs complex division in real arithmetic
 *    a_r_ + i * a_i_ = ( a_r_ + i * a_i_ ) / ( b_r_ + i * b_i_ );
 * The algorithm is due to Robert L. Smith and can be found in D. Knuth,
 * The art of Computer Programming, Vol.2, p.195
 */
#define    Msdiv( b_r_, b_i_, a_r_, a_i_ )                             \
           {                                                           \
              register float c_i_, c_r_, tmp1_, tmp2_;                 \
              if( Msabs( b_i_ ) < Msabs( b_r_ ) )                      \
              {                                                        \
                 tmp1_ = (b_i_) / (b_r_);                              \
                 tmp2_ = (b_r_) + (b_i_) * tmp1_;                      \
                 c_r_  = ( (a_r_) + (a_i_) * tmp1_ ) / tmp2_;          \
                 c_i_  = ( (a_i_) - (a_r_) * tmp1_ ) / tmp2_;          \
              }                                                        \
              else                                                     \
              {                                                        \
                 tmp1_ = (b_r_) / (b_i_);                              \
                 tmp2_ = (b_i_) + (b_r_) * tmp1_;                      \
                 c_r_  = (  (a_i_) + (a_r_) * tmp1_ ) / tmp2_;         \
                 c_i_  = ( -(a_r_) + (a_i_) * tmp1_ ) / tmp2_;         \
              }                                                        \
              a_r_ = c_r_;                                             \
              a_i_ = c_i_;                                             \
           }

#define    Mdabs( a_ ) ( ( (a_) < ATL_dZERO ) ? -(a_) : (a_) )

#define    Mdzero( a_r_, a_i_ )                                        \
           ( ( (a_r_) == ATL_dZERO ) && ( (a_i_) == ATL_dZERO ) )

#define    Mdone( a_r_, a_i_ )                                         \
           ( ( (a_r_) == ATL_dONE  ) && ( (a_i_) == ATL_dZERO ) )

#define    Mdscl( a_r_, a_i_, c_r_, c_i_ )                             \
           {                                                           \
              register double tmp_r_, tmp_i_;                          \
              tmp_r_ = (a_r_) * c_r_ - (a_i_) * c_i_;                  \
              tmp_i_ = (a_r_) * c_i_ + (a_i_) * c_r_;                  \
              c_r_   = tmp_r_;                                         \
              c_i_   = tmp_i_;                                         \
           }
/*
 * Mddiv performs complex division in real arithmetic
 *    a_r_ + i * a_i_ = ( a_r_ + i * a_i_ ) / ( b_r_ + i * b_i_ );
 * The algorithm is due to Robert L. Smith and can be found in D. Knuth,
 * The art of Computer Programming, Vol.2, p.195
 */
#define    Mddiv( b_r_, b_i_, a_r_, a_i_ )                             \
           {                                                           \
              register double c_i_, c_r_, tmp1_, tmp2_;                \
              if( Mdabs( b_i_ ) < Mdabs( b_r_ ) )                      \
              {                                                        \
                 tmp1_ = (b_i_) / (b_r_);                              \
                 tmp2_ = (b_r_) + (b_i_) * tmp1_;                      \
                 c_r_  = ( (a_r_) + (a_i_) * tmp1_ ) / tmp2_;          \
                 c_i_  = ( (a_i_) - (a_r_) * tmp1_ ) / tmp2_;          \
              }                                                        \
              else                                                     \
              {                                                        \
                 tmp1_ = (b_r_) / (b_i_);                              \
                 tmp2_ = (b_i_) + (b_r_) * tmp1_;                      \
                 c_r_  = (  (a_i_) + (a_r_) * tmp1_ ) / tmp2_;         \
                 c_i_  = ( -(a_r_) + (a_i_) * tmp1_ ) / tmp2_;         \
              }                                                        \
              a_r_ = c_r_;                                             \
              a_i_ = c_i_;                                             \
           }

#define    Mmin( a_, b_ ) ( ( (a_) < (b_) ) ?  (a_) : (b_) )

#define    Mmax( a_, b_ ) ( ( (a_) > (b_) ) ?  (a_) : (b_) )

#define    Mmul( a_r_, a_i_, b_r_, b_i_, c_r_, c_i_ )                  \
           {                                                           \
              c_r_ = (a_r_) * (b_r_) - (a_i_) * (b_i_);                \
              c_i_ = (a_r_) * (b_i_) + (a_i_) * (b_r_);                \
           }

#define    Mmla( a_r_, a_i_, b_r_, b_i_, c_r_, c_i_ )                  \
           {                                                           \
              c_r_ += (a_r_) * (b_r_) - (a_i_) * (b_i_);               \
              c_i_ += (a_r_) * (b_i_) + (a_i_) * (b_r_);               \
           }

#define    Mmls( a_r_, a_i_, b_r_, b_i_, c_r_, c_i_ )                  \
           {                                                           \
              c_r_ -= (a_r_) * (b_r_) - (a_i_) * (b_i_);               \
              c_i_ -= (a_r_) * (b_i_) + (a_i_) * (b_r_);               \
           }

#define    Mset( a_r_, a_i_, b_r_, b_i_ )                              \
           {                                                           \
              b_r_ = (a_r_);                                           \
              b_i_ = (a_i_);                                           \
           }

#define    Mselscal( al_, a_ )                                         \
           {                                                           \
              if(      (al_) == ATL_sZERO ) { (a_)  = ATL_sZERO; }     \
              else if( (al_) != ATL_sONE  ) { (a_) *= (al_);     }     \
           }

#define    Mdelscal( al_, a_ )                                         \
           {                                                           \
              if(      (al_) == ATL_dZERO ) { (a_)  = ATL_dZERO; }     \
              else if( (al_) != ATL_dONE  ) { (a_) *= (al_);     }     \
           }

#define    Mcelscal( al_r_, al_i_, a_r_, a_i_ )                        \
           {                                                           \
              if( Mszero( (al_r_), (al_i_) ) )                         \
              { (a_r_) = (a_i_) = ATL_sZERO; }                         \
              else if( ! Msone( (al_r_), (al_i_) ) )                   \
              { Msscl( (al_r_), (al_i_), (a_r_), (a_i_) ); }           \
           }

#define    Mzelscal( al_r_, al_i_, a_r_, a_i_ )                        \
           {                                                           \
              if( Mdzero( (al_r_), (al_i_) ) )                         \
              { (a_r_) = (a_i_) = ATL_dZERO; }                         \
              else if( ! Mdone( (al_r_), (al_i_) ) )                   \
              { Mdscl( (al_r_), (al_i_), (a_r_), (a_i_) ); }           \
           }

#define    Msvscal( n_, al_, x_, incx_ )                               \
           {                                                           \
            int i_, ix_;                                               \
            if(      (al_) == ATL_sZERO )                              \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx_) )   \
             { (x_)[ix_] = ATL_sZERO; }                                \
            }                                                          \
            else if( (al_) != ATL_sONE )                               \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx_) )   \
             { (x_)[ix_] *= (al_); }                                   \
            }                                                          \
           }

#define    Mdvscal( n_, al_, x_, incx_ )                               \
           {                                                           \
            int i_, ix_;                                               \
            if(      (al_) == ATL_dZERO )                              \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx_) )   \
             { (x_)[ix_] = ATL_dZERO; }                                \
            }                                                          \
            else if( (al_) != ATL_dONE )                               \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx_) )   \
             { (x_)[ix_] *= (al_); }                                   \
            }                                                          \
           }

#define    Mcvscal( n_, al_, x_, incx_ )                               \
           {                                                           \
            int i_, ix_, incx2_ = ( 2 * (incx_) );                     \
            if( Mszero( (al_)[0], (al_)[1] ) )                         \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx2_) )  \
             { (x_)[ix_] = (x_)[ix_+1] = ATL_sZERO; }                  \
            }                                                          \
            else if( ! Msone( (al_)[0], (al_)[1] ) )                   \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx2_) )  \
             { Msscl( (al_)[0], (al_)[1], (x_)[ix_], (x_)[ix_+1] ); }  \
            }                                                          \
           }

#define    Mzvscal( n_, al_, x_, incx_ )                               \
           {                                                           \
            int i_, ix_, incx2_ = ( 2 * (incx_) );                     \
            if( Mdzero( (al_)[0], (al_)[1] ) )                         \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx2_) )  \
             { (x_)[ix_] = (x_)[ix_+1] = ATL_dZERO; }                  \
            }                                                          \
            else if( ! Mdone( (al_)[0], (al_)[1] ) )                   \
            {                                                          \
             for( i_ = 0, ix_ = 0; i_ < (n_); i_++, ix_ += (incx2_) )  \
             { Mdscl( (al_)[0], (al_)[1], (x_)[ix_], (x_)[ix_+1] ); }  \
            }                                                          \
           }

#define    Msgescal( m_, n_, al_, a_, lda_ )                           \
           {                                                           \
            int i_, iaij_, j_, jaj_;                                   \
            if(      (al_) == ATL_sZERO )                              \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += (lda_) )  \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 1 ) \
              { (a_)[iaij_] = ATL_sZERO; }                             \
             }                                                         \
            }                                                          \
            else if( (al_) != ATL_sONE )                               \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += (lda_) )  \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 1 ) \
              { (a_)[iaij_] *= (al_); }                                \
             }                                                         \
            }                                                          \
           }

#define    Mdgescal( m_, n_, al_, a_, lda_ )                           \
           {                                                           \
            int i_, iaij_, j_, jaj_;                                   \
            if(      (al_) == ATL_dZERO )                              \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += (lda_) )  \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 1 ) \
              { (a_)[iaij_] = ATL_dZERO; }                             \
             }                                                         \
            }                                                          \
            else if( (al_) != ATL_dONE )                               \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += (lda_) )  \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 1 ) \
              { (a_)[iaij_] *= (al_); }                                \
             }                                                         \
            }                                                          \
           }

#define    Mcgescal( m_, n_, al_, a_, lda_ )                           \
           {                                                           \
            int i_, iaij_, j_, jaj_, lda2_ = ( (lda_) << 1 );          \
            if( Mszero( (al_)[0], (al_)[1] ) )                         \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += lda2_ )   \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 2 ) \
              { (a_)[iaij_] = (a_)[iaij_+1] = ATL_sZERO; }             \
             }                                                         \
            }                                                          \
            else if( ! Msone( (al_)[0], (al_)[1] ) )                   \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += lda2_ )   \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 2 ) \
              {                                                        \
              Msscl( (al_)[0], (al_)[1], (a_)[iaij_], (a_)[iaij_+1] ); \
              }                                                        \
             }                                                         \
            }                                                          \
           }

#define    Mzgescal( m_, n_, al_, a_, lda_ )                           \
           {                                                           \
            int i_, iaij_, j_, jaj_, lda2_ = ( (lda_) << 1 );          \
            if( Mdzero( (al_)[0], (al_)[1] ) )                         \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += lda2_ )   \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 2 ) \
              { (a_)[iaij_] = (a_)[iaij_+1] = ATL_dZERO; }             \
             }                                                         \
            }                                                          \
            else if( ! Mdone( (al_)[0], (al_)[1] ) )                   \
            {                                                          \
             for( j_ = 0, jaj_ = 0; j_ < (n_); j_++, jaj_ += lda2_ )   \
             {                                                         \
              for( i_ = 0, iaij_ = jaj_; i_ < (m_); i_++, iaij_ += 2 ) \
              {                                                        \
              Mdscl( (al_)[0], (al_)[1], (a_)[iaij_], (a_)[iaij_+1] ); \
              }                                                        \
             }                                                         \
            }                                                          \
           }

#endif
/*
 * End of atlas_refmisc.h
 */
