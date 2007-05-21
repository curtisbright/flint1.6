/****************************************************************************

ZmodFpoly-profile.c

Profiling for ZmodFpoly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler-main.h"
#include "ZmodFpoly.h"
#include "flint.h"


// ============================================================================


void sample_ZmodFpoly_FFT(unsigned long length, unsigned long n,
                          unsigned long count)
{
   unsigned long m = ceil_log2(2*length);
   
   ZmodFpoly_t poly;
   ZmodFpoly_init(poly, m, n, 1);
   poly->length = length;
   
   // todo: need to generate random data here
   
   for (unsigned long i = 0; i < count; i++)
      ZmodFpoly_FFT(poly, 2*length);
   
   ZmodFpoly_clear(poly);
}



char* prof2dDriverString_ZmodFpoly_FFT(int argc, char* argv[])
{
   return "ZmodFpoly_FFT over various transform lengths and coefficient sizes";
}


void prof2dDriver_ZmodFpoly_FFT(int argc, char* argv[])
{
   prof2d_set_sampler(sample_ZmodFpoly_FFT);

   for (unsigned long length = 200; length < 250; length++)
   {
      unsigned long m = ceil_log2(2*length);

      unsigned long n_skip = (1 << m) / (4*FLINT_BITS_PER_LIMB);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n < 4*n_skip; n += n_skip)
      {
         prof2d_sample(length, n);
      }
   }
}


// ============================================================================



// end of file ****************************************************************
