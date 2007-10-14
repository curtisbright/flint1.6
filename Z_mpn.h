/****************************************************************************

Z_mpn.h: Z arithmetic

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_Z_MPN_H
#define FLINT_Z_MPN_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "ZmodF_poly.h"

/****************************************************************************/

typedef enum {FFT_PRE, KAR_PRE} precomp_type;
 
typedef struct
{
   precomp_type type;
   ZmodF_poly_p poly;
   unsigned long length;
   unsigned long length2;
   unsigned long coeff_limbs;
   unsigned long limbs1;
   unsigned long limbs2;
   unsigned long msl_bits;
} Z_mpn_precomp_s;

typedef Z_mpn_precomp_s Z_mpn_precomp_t[1]; 

typedef Z_mpn_precomp_s * Z_mpn_precomp_p;

mp_limb_t Z_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2);
                                      
void Z_mul(mpz_t res, mpz_t a, mpz_t b);

void __Z_mul(mpz_t res, mpz_t a, mpz_t b, unsigned long twk);

mp_limb_t Z_mpn_mul_trunc(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                        mp_limb_t * data2, unsigned long limbs2, unsigned long trunc);

void Z_mpn_mul_precomp_init(Z_mpn_precomp_t precomp, mp_limb_t * data1, unsigned long limbs1, unsigned long limbs2);

void Z_mpn_mul_precomp_clear(Z_mpn_precomp_t precomp);

mp_limb_t Z_mpn_mul_precomp(mp_limb_t * res, mp_limb_t * data2, unsigned long limbs2, Z_mpn_precomp_t precomp);
      
#ifdef __cplusplus
 }
#endif
                                 
#endif
