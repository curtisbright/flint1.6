/****************************************************************************

test-support.c: Support code for test modules

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <string.h>
#include "flint.h"
#include "test-support.h"


gmp_randstate_t randstate;


// a bunch of global mpz's, guaranteed to be init'd
mpz_t test_mpz[TEST_MPZ_COUNT];


void test_support_init()
{
   gmp_randinit_default(randstate);
   
   for (unsigned long i = 0; i < TEST_MPZ_COUNT; i++)
      mpz_init(test_mpz[i]);
}


void test_support_cleanup()
{
   gmp_randclear(randstate);

   for (unsigned long i = 0; i < TEST_MPZ_COUNT; i++)
      mpz_clear(test_mpz[i]);
}


unsigned long random_ulong(unsigned long max)
{
   return gmp_urandomm_ui(randstate, max);
}


mp_limb_t random_limb()
{
   return gmp_urandomb_ui(randstate, FLINT_BITS);
}


void random_limbs(mp_limb_t* dest, unsigned long limbs)
{
   mpz_rrandomb(test_mpz[0], randstate, limbs*FLINT_BITS);

   if (random_ulong(2))
   {
      // GMP always sets the high bit equal to 1,
      // so with probability 1/2 we flip all the bits
      mpz_set_ui(test_mpz[1], 1);
      mpz_mul_2exp(test_mpz[1], test_mpz[1], limbs*FLINT_BITS);
      mpz_sub_ui(test_mpz[1], test_mpz[1], 1);
      mpz_sub(test_mpz[0], test_mpz[1], test_mpz[0]);
   }

   memset(dest, 0, limbs * sizeof(mp_limb_t));
   mpz_export(dest, NULL, -1, sizeof(mp_limb_t), 0, 0, test_mpz[0]);
}


void mpz_to_mpn(mp_limb_t* dest, unsigned long limbs, mpz_t src)
{
   memset(dest, 0, limbs * sizeof(mp_limb_t));
   mpz_export(dest, NULL, -1, sizeof(mp_limb_t), 0, 0, src);
}


void mpn_to_mpz(mpz_t dest, mp_limb_t* src, unsigned long limbs)
{
   mpz_import(dest, limbs, -1, sizeof(mp_limb_t), 0, 0, src);
}


// end of file ****************************************************************
