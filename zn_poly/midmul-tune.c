/*
   midmul-tune.c:  tuning program for middle product algorithms
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <math.h>
#include "support.h"
#include "zn_poly_internal.h"
#include "profiler.h"


/*
   For each modulus size, finds approximate crossover points between fallback
   and fft middle product algorithms.
   
   Store these in the global crossover table, and writes some logging
   information to flog.
*/
void tune_midmul(FILE* flog, int verbose)
{
   unsigned bits;
   
   fprintf(flog, "midmul FFT: ");
   fflush(flog);

   // how long we are willing to wait for each profile run
   const double speed = 0.0001;

   // run tuning process for each modulus size
   for (bits = 2; bits <= ULONG_BITS; bits++)
   {
      unsigned long crossover;

      double result[2];

      profile_midmul_info_t info[2];
      info[0]->n = info[1]->n = (1UL << (bits - 1)) + 1;

      info[0]->algo = ALGO_MIDMUL_FALLBACK;
      info[1]->algo = ALGO_MIDMUL_FFT;

      // find an upper bound, where 2nd algorithm appears to be safely
      // ahead of 1st algorithm
      size_t upper;
      int found = 0;
      for (upper = 45; upper <= 16384 && !found; upper = 2 * upper)
      {
         info[0]->len = info[1]->len = upper;
         
         result[0] = profile(NULL, NULL, profile_midmul, info[0], speed);
         result[1] = profile(NULL, NULL, profile_midmul, info[1], speed);
         
         if (result[1] < 0.95 * result[0])
            found = 1;
      }

      if (!found)
      {
         // couldn't find a reasonable upper bound
         crossover = SIZE_MAX;
         goto done;
      }
      
      // find a lower bound, where 1st algorithm appears to be safely
      // ahead of 2nd algorithm
      size_t lower;
      found = 0;
      for (lower = upper/2; lower >= 2 && !found; lower = lower / 2)
      {
         info[0]->len = info[1]->len = lower;
         
         result[0] = profile(NULL, NULL, profile_midmul, info[0], speed);
         result[1] = profile(NULL, NULL, profile_midmul, info[1], speed);
         
         if (result[1] > 1.05 * result[0])
            found = 1;
      }

      if (!found)
      {
         // couldn't find a reasonable lower bound
         crossover = 0;
         goto done;
      }

      // subdivide [lower, upper] into intervals and sample at each endpoint
      double ratio = (double) upper / (double) lower;
      const int max_intervals = 30;
      size_t points[max_intervals + 1];
      double score[max_intervals + 1];
      unsigned i;
      for (i = 0; i <= max_intervals; i++)
      {
         points[i] = ceil(lower * pow(ratio, (double) i / max_intervals));
         info[0]->len = info[1]->len = points[i];
         result[0] = profile(NULL, NULL, profile_midmul, info[0], speed);
         result[1] = profile(NULL, NULL, profile_midmul, info[1], speed);
         score[i] = result[1] / result[0];
      }
      
      // estimate crossover point
      unsigned count = 0;
      for (i = 0; i <= max_intervals; i++)
         if (score[i] > 1.0)
            count++;
      crossover = (size_t)
           ceil(lower * pow(ratio, (double) count / (max_intervals + 1)));

      done:

      if (verbose)
      {
         fprintf(flog, "\nbits = %u, cross to FFT at ", bits);
         
         if (crossover == SIZE_MAX)
            fprintf(flog, "infinity");
         else
            fprintf(flog, "%lu", crossover);
      }
      else
         fprintf(flog, ".");

      fflush(flog);

      tuning_info[bits].midmul_fft_crossover = crossover;
   }

   fprintf(flog, "\n");
}


// end of file ****************************************************************