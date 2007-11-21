/******************************************************************************

 mpn_extras.c
 
 Extra functions for manipulating mpn's and limbs.

 Copyright (C) 2006, William Hart

 mp_limb_t mpn_divmod_1_preinv was adapted from GMP, (C) Free Software Foundation

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "long_extras.h"
#include "memory-manager.h"
#include "mpn_extras.h"
#include "ZmodF_poly.h"
#include "ZmodF_mul.h"
#include "F_mpn_mul-tuning.h"

/*=======================================================================================

    Performs division by a limb d and places the quotient in qp and returns the 
    remainder. Requires a single limb approximation to 1/d as input. If the most
    significant bit of d is not 1 it expects d to be shifted left (by norm bits)
    until the most significant bit is 1 before the inverse is computed. However 
    the original d should be supplied to the function, not the shifted d. 
    
    This code has been adapted from code found in the GMP package version 4.2.1
    (divrem_1.c) (C) Free Software Foundation
*/

mp_limb_t F_mpn_divmod_1_preinv(mp_limb_t * qp, mp_limb_t * up, 
                                  unsigned long un, mp_limb_t d, mp_limb_t dinv, unsigned long norm)
{
  mp_size_t  n;
  mp_size_t  i;
  mp_limb_t  n1, n0;
  mp_limb_t  r = 0;

  n = un;
  if (n == 0)
    return 0;
  
  qp += (n - 1);   /* Make qp point at most significant quotient limb */

  if ((d & (1L<<(FLINT_BITS-1))) != 0)
  {
     if (un != 0)
     {
        /* High quotient limb is 0 or 1, skip a divide step. */
	    mp_limb_t q;
	    r = up[un - 1];
	    q = (r >= d);
	    *qp-- = q;
	    r -= (d & -q);
	    n--;
	    un--;
	 }

     /* Multiply-by-inverse, divisor already normalized. */
     for (i = un - 1; i >= 0; i--)
     {
        n0 = up[i];
        udiv_qrnnd_preinv (*qp, r, r, n0, d, dinv);
        qp--;
     }
     return r;
  } else
  {
     /* Most significant bit of divisor == 0.  */
     
     /* Skip a division if high < divisor (high quotient 0).  Testing here
	 before normalizing will still skip as often as possible.  */
     if (un != 0)
	 {
	    n1 = up[un - 1];
	    if (n1 < d)
        {
           r = n1;
	       *qp-- = 0;
	       n--;
	       if (n == 0) return r;
	       un--;
        }
	 }  

     d <<= norm;
     r <<= norm;

     if (un != 0)
     {
        n1 = up[un - 1];
        r |= (n1 >> (FLINT_BITS - norm));
        for (i = un - 2; i >= 0; i--)
		{
		  n0 = up[i];
		  udiv_qrnnd_preinv (*qp, r, r, 
				     ((n1 << norm) | (n0 >> (FLINT_BITS - norm))), d, dinv);
		  qp--;
		  n1 = n0;
		}
        udiv_qrnnd_preinv (*qp, r, r, n1 << norm, d, dinv);
        qp--;
     }
     
     return r >> norm;
  }
}

mp_limb_t F_mpn_addmul(mp_limb_t * rp, mp_limb_t * s1p, unsigned long s1n, 
                                      mp_limb_t * s2p, unsigned long s2n)
{
   if (s2n == 0) return 0;
   
   mp_limb_t carry;
   
   carry = mpn_addmul_1(rp, s1p, s1n, s2p[0]);
   for (unsigned long i = 1; i < s2n; i++)
   {
      carry = mpn_add_1(rp+i+s1n-1, rp+i+s1n-1, 1, carry); 
      if (s2p[i]) carry += mpn_addmul_1(rp+i, s1p, s1n, s2p[i]);
   }
   return carry;
}

/*=====================================================================================

   Fast Integer Multiplication Code
   
=====================================================================================*/

unsigned long MUL_TWK_VALS[MUL_TWK_COUNT][3] = 
{
   {2000, 2140, 1024},
   {2140, 2430, 64},
   {2430, 2580, 1024},
   {2580, 2700, 64},
   {2700, 2880, 4096},
   {2880, 3850, 16},
   {3850, 4220, 4},
   {4220, 4400, 1024},
   {4400, 4850, 16},
   {4850, 5700, 1024},
   {5700, 7900, 4},
   {7900, 8900, 1024},
   {8900, 97000, 4},
   {97000, 127000, 1},
   {127000, 262000, 4},
   {262000, 517000, 1},
   {517000, 1050000, 4},
   {1050000, 2060000, 1},
   {2060000, 4230000, 4},
   {4230000, 8350000, 1}
};

unsigned long SQR_TWK_VALS[SQR_TWK_COUNT][3] = 
{
   {1564, 1994, 16},
   {1994, 2952, 64},
   {2952, 5921, 16},
   {5921, 32575, 4},
   {32575, 40006, 16},
   {40006, 66526, 4},
   {66526, 127370, 1},
   {127370, 257473, 4},
   {257473, 520507, 1},
   {520507, 1050000, 4},
   {1050000, 2060000, 1},
   {2060000, 4230000, 4},
   {4230000, 8350000, 1}
};

/*
   Splits an mpn into segments of length coeff_limbs and stores in a ZmodF_poly
   in zero padded coefficients of length output_limbs, for use in FFT 
   convolution code. Assumes that the input is total_limbs in length. 
   Used by the large integer multiplication code 
   (F_mpn_mul and F_mpn_mul_precomp and F_mpn_mul_trunc)
*/

void F_mpn_FFT_split(ZmodF_poly_t poly, mp_limb_t * limbs, unsigned long total_limbs,
                               unsigned long coeff_limbs, unsigned long output_limbs)
{
   unsigned long length = (total_limbs-1)/coeff_limbs + 1;
   unsigned long i, j, skip;
   
   for (skip = 0, i = 0; skip+coeff_limbs <= total_limbs; skip+=coeff_limbs, i++)
   {
      for (j = 0; j < output_limbs; j += 8) FLINT_PREFETCH(poly->coeffs[i+1], j);
      
      F_mpn_clear(poly->coeffs[i], output_limbs+1);
      // convert a coefficient
      F_mpn_copy(poly->coeffs[i], limbs+skip, coeff_limbs);
   }
   if (i < length) F_mpn_clear(poly->coeffs[i], output_limbs+1);
   if (total_limbs > skip) F_mpn_copy(poly->coeffs[i], limbs+skip, total_limbs-skip);
   
   poly->length = length;
}

/*
   Recombines coefficients of a ZmodF_poly after doing a convolution. Assumes 
   each of the coefficients of the ZmodF_poly is output_limbs long, that each 
   of the coefficients is being shifted by a multiple of coeff_limbs and added
   to an mpn which is total_limbs long. It is assumed that the mpn has been 
   zeroed in advance.
   Used by the large integer multiplication code (F_mpn_mul and F_mpn_mul_precomp
   and F_mpn_mul_trunc)
*/

void F_mpn_FFT_combine(mp_limb_t * res, ZmodF_poly_t poly, unsigned long coeff_limbs, 
                             unsigned long output_limbs, unsigned long total_limbs)
{
   unsigned long skip, i, j;
   unsigned long length = poly->length;
   
   for (skip = 0, i = 0; (i < length) && (skip+output_limbs <= total_limbs); i++, skip+=coeff_limbs)
   { 
      for (j = 0; j < output_limbs; j += 8) FLINT_PREFETCH(poly->coeffs[i+1], j);
      mpn_add(res+skip, res+skip, output_limbs+1, poly->coeffs[i], output_limbs);      
   } 
   while ((skip < total_limbs) && (i < length))
   {
      mpn_add(res+skip, res+skip, total_limbs - skip, poly->coeffs[i], FLINT_MIN(total_limbs - skip, output_limbs));
      i++;
      skip+=coeff_limbs;
   }  

}

mp_limb_t __F_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2, unsigned long twk)
{
   unsigned long length = 1;
   unsigned long log_length = 0;
   
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long s1 = (FLINT_BIT_COUNT(data1[limbs1-1]) + FLINT_BIT_COUNT(data2[limbs2-1]) <= FLINT_BITS);
   unsigned long total_limbs = coeff_limbs - s1;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 0;
   
   if (twk > 64)
   {
      length = 2;
      log_length = 1;
      while ((1<<(log_length-1)) < output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
      while (twk > 64)
      {
         log_length--;
         length>>=1;
         twk>>=2;
      }
      coeff_limbs = (limbs1+limbs2-1)/length+1;
      while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
      output_bits = (2*coeff_limbs+1)*FLINT_BITS;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
      while ((output_bits%3) != 0) output_bits+=(1<<(log_length-1));
      coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
      if ((long) coeff_limbs < 1) coeff_limbs = 1;
      length1 = (limbs1-1)/coeff_limbs+1;
      length2 = (limbs2-1)/coeff_limbs+1;
      log_length = 1;
      while ((1<<log_length) < length1 + length2) log_length++;
      length = (1<<log_length);        
   }
   else
   {
      while (twk*length < 2*output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
   }
         
   n = output_bits/FLINT_BITS;
   //printf("n= %ld\n",n);
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld\n", length1, length2, output_bits, coeff_limbs, log_length);
#endif   
   ZmodF_poly_t poly1;
   ZmodF_poly_stack_init(poly1, log_length, n, 1);
   F_mpn_FFT_split(poly1, data1, limbs1, coeff_limbs, n);
   
   if ((data1 == data2) && (limbs1 == limbs2))
   {
      // identical operands case
      ZmodF_poly_convolution(poly1, poly1, poly1);
   }
   else
   {
      // distinct operands case
      ZmodF_poly_t poly2;
      ZmodF_poly_stack_init(poly2, log_length, n, 1);
      F_mpn_FFT_split(poly2, data2, limbs2, coeff_limbs, n);

      ZmodF_poly_convolution(poly1, poly1, poly2);

      ZmodF_poly_stack_clear(poly2);
   }
   
   ZmodF_poly_normalise(poly1);
   
   F_mpn_clear(res, total_limbs);
   
   F_mpn_FFT_combine(res, poly1, coeff_limbs, 2*coeff_limbs+1, total_limbs);
   ZmodF_poly_stack_clear(poly1);
   
   if (s1) return 0;
   else return res[limbs1+limbs2-1];
}

mp_limb_t __F_mpn_mul_trunc(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2, 
                                      unsigned long twk, unsigned long trunc)
{
   unsigned long length = 1;
   unsigned long log_length = 0;
   
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 0;
   
   if (twk > 64)
   {
      length = 2;
      log_length = 1;
      while ((1<<(log_length-1)) < output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1L) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
      while (twk > 64)
      {
         log_length--;
         length>>=1;
         twk>>=2;
      }
      coeff_limbs = (limbs1+limbs2-1)/length+1;
      while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
      output_bits = (2*coeff_limbs+1)*FLINT_BITS;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
      while ((output_bits%3) != 0) output_bits+=(1<<(log_length-1));
      coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
      if ((long) coeff_limbs < 1L) coeff_limbs = 1;
      length1 = (limbs1-1)/coeff_limbs+1;
      length2 = (limbs2-1)/coeff_limbs+1;
      log_length = 1;
      while ((1<<log_length) < length1 + length2) log_length++;
      length = (1<<log_length);        
   }
   else
   {
      while (twk*length < 2*output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1L) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
   }
         
   n = output_bits/FLINT_BITS;
   //printf("n= %ld\n",n);
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld\n", length1, length2, output_bits, coeff_limbs, log_length);
#endif   
   ZmodF_poly_t poly1;
   ZmodF_poly_stack_init(poly1, log_length, n, 1);
   F_mpn_FFT_split(poly1, data1, limbs1, coeff_limbs, n);
   
   if (data1 == data2 && limbs1 == limbs2)
   {
      // identical operands case
      ZmodF_poly_convolution_trunc(poly1, poly1, poly1, (trunc-1)/coeff_limbs + 1);
   }
   else
   {
      // distinct operands case
      ZmodF_poly_t poly2;
      ZmodF_poly_stack_init(poly2, log_length, n, 1);
      F_mpn_FFT_split(poly2, data2, limbs2, coeff_limbs, n);

      ZmodF_poly_convolution_trunc(poly1, poly1, poly2, (trunc-1)/coeff_limbs + 1);

      ZmodF_poly_stack_clear(poly2);
   }
   
   ZmodF_poly_normalise(poly1);
   
   F_mpn_clear(res, trunc);
   
   F_mpn_FFT_combine(res, poly1, coeff_limbs, 2*coeff_limbs+1, trunc);
   ZmodF_poly_stack_clear(poly1);
   
   return res[trunc - 1];
}

/*
   Multiply two integers in mpn format
   
   WARNING: This function requires limbs1+limbs2 output limbs when limbs1+limbs2
   < FLINT_FFT_LIMBS_CROSSOVER but may require one less limb otherwise. The function
   will return 0 if it did not require (and indeed did not zero) the extra limb,
   otherwise it returns the (non zero) value of this high limb after multiplication.
   
   Assumes neither of limbs1, limbs2 is zero. 
*/

mp_limb_t F_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2)
{
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long twk;
   
   if (coeff_limbs/2 < FLINT_FFT_LIMBS_CROSSOVER) 
   {
      return mpn_mul(res, data1, limbs1, data2, limbs2);
   } 
   
   if (data1 != data2)
   {
      if (coeff_limbs/2 < MUL_TWK_SMALL_CUTOFF) twk = MUL_TWK_SMALL_DEFAULT;
      else if (coeff_limbs/2 > MUL_TWK_LARGE_CUTOFF) twk = MUL_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < MUL_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs/2 < MUL_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < MUL_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = MUL_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   } else
   {
      if (coeff_limbs/2 < SQR_TWK_SMALL_CUTOFF) twk = SQR_TWK_SMALL_DEFAULT;
      else if (coeff_limbs/2 > SQR_TWK_LARGE_CUTOFF) twk = SQR_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < SQR_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs/2 < SQR_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < SQR_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = SQR_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   }

   return __F_mpn_mul(res, data1, limbs1, data2, limbs2, twk);
}

/*
   Multiply two integers in mpn format truncating to _trunc_ output limbs
   Assumes none of limbs1, limbs2 and trunc is zero. 
*/

mp_limb_t F_mpn_mul_trunc(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                        mp_limb_t * data2, unsigned long limbs2, unsigned long trunc)
{
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long twk;
   
   if (coeff_limbs/2 < FLINT_FFT_LIMBS_CROSSOVER) 
   {
      return mpn_mul(res, data1, limbs1, data2, limbs2);
   } 
   
   if (data1 != data2)
   {
      if (coeff_limbs/2 < MUL_TWK_SMALL_CUTOFF) twk = MUL_TWK_SMALL_DEFAULT;
      else if (coeff_limbs/2 > MUL_TWK_LARGE_CUTOFF) twk = MUL_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < MUL_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs/2 < MUL_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < MUL_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = MUL_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   } else
   {
      if (coeff_limbs/2 < SQR_TWK_SMALL_CUTOFF) twk = SQR_TWK_SMALL_DEFAULT;
      else if (coeff_limbs/2 > SQR_TWK_LARGE_CUTOFF) twk = SQR_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < SQR_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs/2 < SQR_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < SQR_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = SQR_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   }

   return __F_mpn_mul_trunc(res, data1, limbs1, data2, limbs2, twk, trunc);
}

/*   
   Precompute an FFT for integer multiplication.
   Assumes neither of limbs1, limbs2 is zero. 
*/

void F_mpn_mul_precomp_init(F_mpn_precomp_t precomp, mp_limb_t * data1, unsigned long limbs1, unsigned long limbs2)
{
   unsigned long length = 1;
   unsigned long log_length = 0;
   
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 0;
   
   unsigned long twk;
   
   if (coeff_limbs/2 < MUL_TWK_SMALL_CUTOFF) twk = MUL_TWK_SMALL_DEFAULT;
   else if (coeff_limbs/2 > MUL_TWK_LARGE_CUTOFF) twk = MUL_TWK_LARGE_DEFAULT;
   else 
   {  
      for (unsigned long twk_count = 0; twk_count < MUL_TWK_COUNT; twk_count++)
      {
         if ((coeff_limbs/2 < MUL_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < MUL_TWK_VALS[twk_count][1])) continue;
         else 
         {
            twk = MUL_TWK_VALS[twk_count][2];
            break;
         }
      }
   }
   if (twk > 64)
   {
      length = 2;
      log_length = 1;
      while ((1<<(log_length-1)) < output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1L) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
      while (twk > 64)
      {
         log_length--;
         length>>=1;
         twk>>=2;
      }
      coeff_limbs = (limbs1+limbs2-1)/length+1;
      while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
      output_bits = (2*coeff_limbs+1)*FLINT_BITS;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
      while ((output_bits%3) != 0) output_bits+=(1<<(log_length-1));
      coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
      if ((long) coeff_limbs < 1L) coeff_limbs = 1;
      length1 = (limbs1-1)/coeff_limbs+1;
      length2 = (limbs2-1)/coeff_limbs+1;
      log_length = 1;
      while ((1<<log_length) < length1 + length2) log_length++;
      length = (1<<log_length);        
   }
   else
   {
      while (twk*length < 2*output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1L) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
   }
      
   n = output_bits/FLINT_BITS;
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld\n", length1, length2, output_bits, coeff_limbs, log_length);
#endif   
   ZmodF_poly_p poly1;
   poly1 = (ZmodF_poly_p) malloc(sizeof(ZmodF_poly_struct));
   ZmodF_poly_stack_init(poly1, log_length, n, 1);
   F_mpn_FFT_split(poly1, data1, limbs1, coeff_limbs, n);
   
   ZmodF_poly_FFT(poly1, length1 + length2 - 1);
   precomp->type = FFT_PRE;
   precomp->length = length1;
   precomp->length2 = length2;
   precomp->coeff_limbs = coeff_limbs;
   precomp->limbs1 = limbs1;
   precomp->limbs2 = limbs2;
   precomp->poly = poly1; 
   precomp->msl_bits = FLINT_BIT_COUNT(data1[limbs1-1]);
}

void F_mpn_mul_precomp_clear(F_mpn_precomp_t precomp)
{
   if (precomp->type == FFT_PRE) 
   {
      ZmodF_poly_stack_clear(precomp->poly);
      free(precomp->poly);
   }   
}

/*   
   Compute an integer multiplication given a precomputed FFT for 
   one of the integers.
   Assumes neither of limbs1, limbs2 is zero. 
*/

mp_limb_t F_mpn_mul_precomp(mp_limb_t * res, mp_limb_t * data2, unsigned long limbs2, F_mpn_precomp_t precomp)
{
   ZmodF_poly_t poly2;
   ZmodF_poly_stack_init(poly2, precomp->poly->depth, precomp->poly->n, 1);
   int s1 = (FLINT_BIT_COUNT(data2[limbs2-1]) + precomp->msl_bits <= FLINT_BITS);
   
   F_mpn_FFT_split(poly2, data2, limbs2, precomp->coeff_limbs, precomp->poly->n);
   
   /*for (unsigned long i = poly2->length; i < precomp->length2; i++)
   {
      F_mpn_clear(poly2->coeffs[i], poly2->n+1);
   }
   poly2->length = precomp->length2;*/
   
   ZmodF_poly_FFT(poly2, precomp->length+poly2->length-1);
   ZmodF_poly_pointwise_mul(poly2, poly2, precomp->poly);
   ZmodF_poly_IFFT(poly2);
   ZmodF_poly_rescale(poly2);
   
   ZmodF_poly_normalise(poly2);
   F_mpn_clear(res, precomp->limbs1 + limbs2 - s1);
   
   F_mpn_FFT_combine(res, poly2, precomp->coeff_limbs, 2*precomp->coeff_limbs+1, precomp->limbs1 + limbs2 - s1);
   
   ZmodF_poly_stack_clear(poly2);
   
   if (s1) return 0;
   else return res[precomp->limbs1+limbs2-1];
}


