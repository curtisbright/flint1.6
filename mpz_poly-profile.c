/****************************************************************************

mpz_poly-profile.c

Profiling for mpz_poly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler-main.h"
#include "mpz_poly.h"
#include "flint.h"
#include "test-support.h"
#include <string.h>
#include <math.h>


// ============================================================================


/*
this function samples multiplying polynomials of lengths len1 and len2
using mpz_poly_mul_karatsuba

arg should point to an unsigned long, giving the coefficient bitlengths
*/
void sample_mpz_poly_mul_karatsuba_mixlengths(
     unsigned long len1, unsigned long len2, void* arg, unsigned long count)
{
   unsigned long bits = *(unsigned long*) arg;
   
   mpz_poly_t poly1, poly2, poly3;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_poly_init(poly3);

   mpz_t x;
   mpz_init(x);
   for (unsigned long i = 0; i < len1; i++)
   {
      mpz_urandomb(x, randstate, bits);
      if (random_ulong(2)) mpz_neg(x, x);
      mpz_poly_set_coeff(poly1, i, x);
   }
   for (unsigned long i = 0; i < len2; i++)
   {
      mpz_urandomb(x, randstate, bits);
      if (random_ulong(2)) mpz_neg(x, x);
      mpz_poly_set_coeff(poly2, i, x);
   }
   mpz_clear(x);
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
      mpz_poly_mul_karatsuba(poly3, poly1, poly2);

   prof_stop();
   
   mpz_poly_clear(poly3);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
}


char* profDriverString_mpz_poly_mul_karatsuba_mixlengths(char* params)
{
   return "mpz_poly_mul_karatubsa for distinct input lengths and fixed\n"
   "coefficient size. Parameters are: max length; length skip; coefficient size (in bits)\n";
}

char* profDriverDefaultParams_mpz_poly_mul_karatsuba_mixlengths()
{
   return "50 1 100";
}


void profDriver_mpz_poly_mul_karatsuba_mixlengths(char* params)
{
   unsigned long max_length, skip, bits;

   sscanf(params, "%ld %ld %ld", &max_length, &skip, &bits);

   prof2d_set_sampler(sample_mpz_poly_mul_karatsuba_mixlengths);

   test_support_init();

   for (unsigned long len1 = skip; len1 <= max_length; len1 += skip)
      for (unsigned long len2 = skip; len2 <= len1; len2 += skip)
         prof2d_sample(len1, len2, &bits);

   test_support_cleanup();
}


// end of file ****************************************************************
