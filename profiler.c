/******************************************************************************

 Timing/profiling

 (C) 2007 William Hart and David Harvey

******************************************************************************/

#include "profiler.h"
#include <math.h>


/*
clock_last[i] is the last read clock value for clock #i.
clock_accum[i] is the total time attributed to clock #i so far.
These should not be read directly; use get_clock(i) instead.
 */
double clock_last[FLINT_NUM_CLOCKS];
double clock_accum[FLINT_NUM_CLOCKS];


void prof_repeat(double* min, double* max, profile_target_t target, void* arg)
{
   // the number of timings that were at least DURATION_THRESHOLD microseconds:
   unsigned long good_count = 0;
   double max_time, min_time;

   // first try one loop
   unsigned long num_trials = 1;
   init_clock(0);
   target(arg, num_trials);
   double last_time = get_clock(0);

   // loop until we have enough good times
   while (1)
   {
      double per_trial = last_time / num_trials;
      
      // if the last recorded time was long enough, record it
      if (last_time > DURATION_THRESHOLD)
      {
         if (good_count)
         {
            if (per_trial > max_time)
               max_time = per_trial;
            if (per_trial < min_time)
               min_time = per_trial;
         }
         else
            max_time = min_time = per_trial;

         if (++good_count == 5)
         {
            // we've got enough data
            break;
         }
      }

      // adjust num_trials so that the elapsed time gravitates towards
      // DURATION_TARGET; num_trials can be changed by a factor of
      // at most 25%, and must be at least 1
      if (last_time < 0.0001)
         last_time = 0.0001;
      double adjust_ratio = DURATION_TARGET / last_time;
      if (adjust_ratio > 1.25)
         adjust_ratio = 1.25;
      if (adjust_ratio < 0.75)
         adjust_ratio = 0.75;
      num_trials = (unsigned long) ceil(adjust_ratio * num_trials);
      // just to be safe:
      if (num_trials == 0)
         num_trials = 1;

      // run another trial
      init_clock(0);
      target(arg, num_trials);
      last_time = get_clock(0);
   }
   
   // store results
   if (min)
      *min = min_time;
   if (max)
      *max = max_time;
}



// end of file ****************************************************************
