/*
   Copyright 2005, 2006 Damien Stehl�.
   Copyright 2009, 2010 William Hart
   Copyright 2009, 2010 Andy Novocin

   This program implements ideas from the paper "Floating-point LLL Revisited", 
   by Phong Nguyen and Damien Stehl�, in the Proceedings of Eurocrypt'2005, 
   Springer-Verlag; and was partly inspired by Shoup's NTL library: 
   http://www.shoup.net/ntl/
 
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2 of the License, or (at your
   option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.
*/
/****************************************************************************

   F_mpz_LLL.c: LLL reduction on the rows of an F_mpz_mat_t

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctype.h>
#include <limits.h>
#include "gmp.h"
#include "flint.h"
#include "profiler.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL.h"
#include "mpfr.h"
#include "mpfr_mat.h"
#include "d_mat.h"
#include <signal.h>

/****************************************************************************

   The various Babai's: check_Babai, check_Babai_heuristic_d, check_Babai_heuristic

****************************************************************************/

#define PROFILE 0

#if PROFILE
   // LLL profiles are on, checking babai total, time spent updating B-- 
   // but not including the time in advanced babai which is totaled
   
   double babai_start, babai_stop;
   double babai_total = 0;

   double adv_babai_start, adv_babai_stop;
   double adv_babai_total = 0;

   double update_start, update_stop;
   double update_total = 0.0;

   double convert_start, convert_stop;
   double convert_total = 0.0;

   double inner_start, inner_stop;
   double inner_total = 0.0;

   double ldexp_start, ldexp_stop;
   double ldexp_total = 0.0;

   double hbabai_start, hbabai_stop;
   double hbabai_total = 0;

   double hadv_babai_start, hadv_babai_stop;
   double hadv_babai_total = 0;

   double hupdate_start, hupdate_stop;
   double hupdate_total = 0.0;

   double hconvert_start, hconvert_stop;
   double hconvert_total = 0.0;

   double hinner_start, hinner_stop;
   double hinner_total = 0.0;

   double hldexp_start, hldexp_stop;
   double hldexp_total = 0.0;

#endif

int global_flag = 0;

typedef void (*sighandler_t)(int);

void handler(int varn){

    printf("I've been handled\n");
    global_flag++;

    return;
}

/*
   Computes the scalar product of two vectors of doubles vec1 and vec2, which are 
   respectively double approximations (up to scaling by a power of 2) to rows k and
   j in the exact integer matrix B. If massive cancellation is detected an exact
   computation is made.

   The exact computation is scaled by 2^-exp_adj, where exp_adj = r2 + r1 where
   r2 is the exponent for row j and r1 is the exponent for row k (i.e. row j is 
   notionally thought of as being multiplied by 2^r2, etc).

   The final scalar product computed by this function is then notionally the return
   value times 2^exp_adj.
*/
double heuristic_scalar_product(double * vec1, double * vec2, ulong n, 
								F_mpz_mat_t B, ulong k, ulong j, long exp_adj)
{
   double sum = _d_vec_scalar_product(vec1, vec2, n);
   double tmp = _d_vec_norm(vec1, n);
   double tmp2 = _d_vec_norm(vec2, n);

   tmp = ldexp(tmp*tmp2, -70);
   tmp2 = sum*sum;

   if (tmp2 <= tmp)
   {
      F_mpz_t sp;
      F_mpz_init(sp);
      ulong exp;
      _F_mpz_vec_scalar_product(sp, B->rows[k], B->rows[j], n);
      sum = F_mpz_get_d_2exp(&exp, sp);
      sum = ldexp(sum, exp - exp_adj);
      F_mpz_clear(sp);
   }

   return sum;
} 

/*
   Performs floating point size reductions of the kappa-th row of B by all of 
   the previous rows, uses d_mats mu and r for storing the GSO data.
   While the double array s will contain the size of the kappa-th row if it 
   were moved into position i.  The d_mat appB is an approximation of 
   B with each row receiving an exponent stored in expo which gets populated 
   only when needed.  The d_mat appSP is an approximation of the Gram matrix 
   whose entries are scalar products of the rows of B.  The index a is the 
   smallest row index which will be reduced from the kappa-th row.  Index zeros 
   is the number of zero rows in the matrix.  Kappamax is the highest index 
   which has been size-reduced so far, and n is the number of columns you want to 
   consider.  The output is the value -1 if the process fails (usually do to 
   insufficient precision) or 0 if everything was successful.  These descriptions 
   will be true for the future Babai procedures as well.
*/

int check_Babai (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;

   long loops = 0;

   do
   {
      test = 0;

      loops++;
      if (loops > 200){
         return -1;
      }
            
      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
#if PROFILE
             inner_start = get_cycle_counter();
#endif
	         appSP[kappa][j] = _d_vec_scalar_product(appB[kappa], appB[j], n);
#if PROFILE
             inner_stop = get_cycle_counter();
             inner_total = inner_total + inner_stop - inner_start;
#endif
	      }
	  	  
          if (j > zeros + 2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         rtmp = appSP[kappa][j] - tmp;
	         
			 for (k = zeros + 2; k < j - 1; k++)
		      {
		         tmp = mu[j][k] * r[kappa][k];
		         rtmp = rtmp - tmp;
		      }
	         
			 tmp = mu[j][j-1] * r[kappa][j-1];
	         r[kappa][j] = rtmp - tmp;
          } else if (j == zeros+2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         r[kappa][j] = appSP[kappa][j] - tmp;
	      } else r[kappa][j] = appSP[kappa][j];

	      mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	  {

	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
#if PROFILE
          ldexp_start = get_cycle_counter();
#endif
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
#if PROFILE
          ldexp_stop = get_cycle_counter();
          ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif	  
	      if (tmp > halfplus) 
	      {
	         test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		      {		  
		         if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
                 {
		             for (k = zeros + 1; k < j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
#if PROFILE
                     update_start = get_cycle_counter();
#endif
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
#if PROFILE
                    update_stop = get_cycle_counter();
                    update_total = update_total + update_stop - update_start;
#endif
		  
		         } else          /* otherwise X is -1 */ 
                 {
                     for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
#if PROFILE
                     update_start = get_cycle_counter();
#endif
                     _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
#if PROFILE
                     update_stop = get_cycle_counter();
                     update_total = update_total + update_stop - update_start;
#endif
                 }
		      } else   /* we must have |X| >= 2 */
		      {
#if PROFILE
                 ldexp_start = get_cycle_counter();
#endif
		         tmp = ldexp (mu[kappa][j] , -exponent);
#if PROFILE
                 ldexp_stop = get_cycle_counter();
                 ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
	             if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		            tmp = rint (tmp); 
		      
		            for (k = zeros + 1; k < j; k++)
			        {
			           rtmp = tmp * mu[j][k];
#if PROFILE
                       ldexp_start = get_cycle_counter();
#endif
			           rtmp = ldexp (rtmp, exponent);
#if PROFILE
                       ldexp_stop = get_cycle_counter();
                       ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
			           mu[kappa][k] = mu[kappa][k] - rtmp;
			        }

		            xx = (long) tmp;
		      
                    if (xx > 0L)
                    { 
#if PROFILE
                       update_start = get_cycle_counter();
#endif
                       _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
#if PROFILE
                       update_stop = get_cycle_counter();
                       update_total = update_total + update_stop - update_start;
#endif
                    } else
                    {
#if PROFILE
                       update_start = get_cycle_counter();
#endif
                       _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
#if PROFILE
                       update_stop = get_cycle_counter();
                       update_total = update_total + update_stop - update_start;
#endif
                    } 
                } else
		        {
#if PROFILE
                   ldexp_start = get_cycle_counter();
#endif
		           tmp = frexp(mu[kappa][j], &exponent); 
#if PROFILE
                   ldexp_stop = get_cycle_counter();
                   ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
		           tmp = tmp * MAX_LONG;
		           xx = (long) tmp;
		           exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		           /* This case is extremely rare: never happened for me */
		           if (exponent <= 0) 
			       {
//                     printf("rare case kappa = %d, j = %d ******\n", kappa, j);
			           xx = xx << -exponent;
			           exponent = 0;
			  
			           if (xx > 0)
                       {
#if PROFILE
                          update_start = get_cycle_counter();
#endif
                          _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
#if PROFILE
                          update_stop = get_cycle_counter();
                          update_total = update_total + update_stop - update_start;
#endif
                       } else
                       {
#if PROFILE
                          update_start = get_cycle_counter();
#endif
                          _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
#if PROFILE
                          update_stop = get_cycle_counter();
                          update_total = update_total + update_stop - update_start;
#endif
                       }
              			    
			           for (k = zeros + 1; k < j; k++)
			           {
                          rtmp = ((double) xx) * mu[j][k];
#if PROFILE
                          ldexp_start = get_cycle_counter();
#endif
			              rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
#if PROFILE
                          ldexp_stop = get_cycle_counter();
                          ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
			              mu[kappa][k] = mu[kappa][k] - rtmp;
			           }
			           } else
			           {
			              if (xx > 0)
                          {
#if PROFILE
                             update_start = get_cycle_counter();
#endif
                             _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
#if PROFILE
                             update_stop = get_cycle_counter();
                             update_total = update_total + update_stop - update_start;
#endif
                          } else
                          {
#if PROFILE
                              update_start = get_cycle_counter();
#endif
                              _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
#if PROFILE
                              update_stop = get_cycle_counter();
                              update_total = update_total + update_stop - update_start;
#endif
                           }
			               
						   for (k = zeros + 1; k < j; k++)
			               {
			                  rtmp = ((double) xx) * mu[j][k];
#if PROFILE
                              ldexp_start = get_cycle_counter();
#endif
			                  rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
#if PROFILE
                              ldexp_stop = get_cycle_counter();
                              ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
			                  mu[kappa][k] = mu[kappa][k] - rtmp;
					       }
				       }	    
			       }
		       }
		   }
	   }

       if (test)   /* Anything happened? */
	   {
#if PROFILE
          convert_start = get_cycle_counter();
#endif
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
#if PROFILE
          convert_stop = get_cycle_counter();
          convert_total = convert_total + convert_stop - convert_start;
#endif
	      aa = zeros + 1;
	      
		  for (i = zeros + 1; i <= kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      
		  for (i = kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   }
   } while (test);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
#if PROFILE
      inner_start = get_cycle_counter();
#endif
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
#if PROFILE
      inner_stop = get_cycle_counter();
      inner_total = inner_total + inner_stop - inner_start;
#endif
   }
   
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }

   return 0;
}

/* 
   Same as check_Babai but using the heuristic inner product rather than a purely 
   floating point inner product.  The heuristic will compute at full precision 
   when there is cancellation.
*/

int check_Babai_heuristic_d (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;

   int loops = 0;

   do
   {
      test = 0;
            
      loops++;
      if (loops > 200)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	  {	  
	     if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	     {
#if PROFILE
            hinner_start = get_cycle_counter();
#endif
            //### This is different -----
            appSP[kappa][j] = heuristic_scalar_product(appB[kappa], appB[j], 
				                             n, B, kappa, j, expo[kappa]+expo[j]);
#if PROFILE
            hinner_stop = get_cycle_counter();
            hinner_total = hinner_total + hinner_stop - hinner_start;
#endif
         }
	  	  
         if (j > zeros + 2)
	     {
	        tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	        rtmp = appSP[kappa][j] - tmp;
	        for (k = zeros + 2; k < j - 1; k++)
		    {
		       tmp = mu[j][k] * r[kappa][k];
		       rtmp = rtmp - tmp;
		    }
	        
			tmp = mu[j][j-1] * r[kappa][j-1];
	        r[kappa][j] = rtmp - tmp;
         } else if (j == zeros+2)
	     {
	        tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	        r[kappa][j] = appSP[kappa][j] - tmp;
	     } else r[kappa][j] = appSP[kappa][j];

	     mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	  {

	     /* test of the relaxed size-reduction condition */
	     tmp = fabs(mu[kappa][j]);
#if PROFILE
         hldexp_start = get_cycle_counter();
#endif
	     tmp = ldexp(tmp, expo[kappa] - expo[j]);
#if PROFILE
         hldexp_stop = get_cycle_counter();
         hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
	  
	     if (tmp > halfplus) 
	     {
	        test = 1; 
	        exponent = expo[j] - expo[kappa];
	      
	        /* we consider separately the cases X = +-1 */     
	        if (tmp <= onedothalfplus)   
		    {		  
		       if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
               {
		          for (k = zeros + 1; k < j; k++)
			      {
#if PROFILE
                     hldexp_start = get_cycle_counter();
#endif
			         tmp = ldexp (mu[j][k], exponent);
#if PROFILE
                     hldexp_stop = get_cycle_counter();
                     hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			         mu[kappa][k] =  mu[kappa][k] - tmp; 
			      }
#if PROFILE
                  hupdate_start = get_cycle_counter();
#endif		      
		          _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
#if PROFILE
                  hupdate_stop = get_cycle_counter();
                  hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif		  
		       } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			      {
#if PROFILE
                     hldexp_start = get_cycle_counter();
#endif
			         tmp = ldexp (mu[j][k], exponent);
#if PROFILE
                     hldexp_stop = get_cycle_counter();
                     hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			         mu[kappa][k] = mu[kappa][k] + tmp;
			      }

#if PROFILE
                  hupdate_start = get_cycle_counter();
#endif		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
#if PROFILE
                  hupdate_stop = get_cycle_counter();
                  hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
               }
		       } else   /* we must have |X| >= 2 */
		       {
#if PROFILE
                  hldexp_start = get_cycle_counter();
#endif
		          tmp = ldexp (mu[kappa][j] , -exponent);
#if PROFILE
                  hldexp_stop = get_cycle_counter();
                  hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif

	              if ((tmp < (double) MAX_LONG) && (tmp > (double) -MAX_LONG))  
		          {
		             tmp = rint (tmp); 
		      
		             for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
#if PROFILE
                        hldexp_start = get_cycle_counter();
#endif
			            rtmp = ldexp (rtmp, exponent);
#if PROFILE
                        hldexp_stop = get_cycle_counter();
                        hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		             xx = (long) tmp;
		      
                     if (xx > 0L)
                     { 
#if PROFILE
                        hupdate_start = get_cycle_counter();
#endif
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
#if PROFILE
                        hupdate_stop = get_cycle_counter();
                        hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                     } else
                     {
#if PROFILE
                        hupdate_start = get_cycle_counter();
#endif
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
#if PROFILE
                        hupdate_stop = get_cycle_counter();
                        hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                     } 
                  } else
		          {
#if PROFILE
                     hldexp_start = get_cycle_counter();
#endif
		             tmp = frexp(mu[kappa][j], &exponent); 
#if PROFILE
                     hldexp_stop = get_cycle_counter();
                     hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
		             tmp = tmp * MAX_LONG;
		             xx = (long) tmp;
		             exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		             /* This case is extremely rare: never happened for me */
		             if (exponent <= 0) 
			         {
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                        {
#if PROFILE
                           hupdate_start = get_cycle_counter();
#endif
                           _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
#if PROFILE
                           hupdate_stop = get_cycle_counter();
                           hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                        } else
                        {
#if PROFILE
                           hupdate_start = get_cycle_counter();
#endif
                           _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
#if PROFILE
                           hupdate_stop = get_cycle_counter();
                           hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                        }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                           rtmp = ((double) xx) * mu[j][k];
#if PROFILE
                           hldexp_start = get_cycle_counter();
#endif
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
#if PROFILE
                           hldexp_stop = get_cycle_counter();
                           hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                        {
#if PROFILE
                           hupdate_start = get_cycle_counter();
#endif
                           _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
#if PROFILE
                           hupdate_stop = get_cycle_counter();
                           hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                        } else
                        {
#if PROFILE
                           hupdate_start = get_cycle_counter();
#endif
                           _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
#if PROFILE
                           hupdate_stop = get_cycle_counter();
                           hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                        }
			            
						for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
#if PROFILE
                           hldexp_start = get_cycle_counter();
#endif
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
#if PROFILE
                           hldexp_stop = get_cycle_counter();
                           hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					    }
				     }	    
			      }
		       }
		    }
	    }

        if (test)   /* Anything happened? */
	    {
#if PROFILE
           hconvert_start = get_cycle_counter();
#endif
	       expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
#if PROFILE
           hconvert_stop = get_cycle_counter();
           hconvert_total = hconvert_total + hconvert_stop - hconvert_start;
#endif
	       aa = zeros + 1;
	        
		   for (i = zeros + 1; i <= kappa; i++) 
	          appSP[kappa][i] = NAN;//0.0/0.0;
	        
		   for (i = kappa + 1; i <= kappamax; i++) 
	       appSP[i][kappa] = NAN;//0.0/0.0;
       }
   } while (test);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
      // ### This is different -------
#if PROFILE
      hinner_start = get_cycle_counter();
#endif
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
#if PROFILE
      hinner_stop = get_cycle_counter();
      hinner_total = hinner_total + hinner_stop - hinner_start;
#endif
   }
   
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }

   return 0;
}

/* The mpfr version of check_Babai_heuristic_d, also inherits some temp variables.*/

int check_Babai_heuristic(int kappa, F_mpz_mat_t B, __mpfr_struct **mu, __mpfr_struct **r, __mpfr_struct *s, 
       __mpfr_struct **appB, __mpfr_struct **appSP, 
       int a, int zeros, int kappamax, int n, mpfr_t tmp, mpfr_t rtmp, mp_prec_t prec)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   F_mpz_t ztmp, X;
   F_mpz_init(ztmp);
   F_mpz_init(X);
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;

   long loops = 0;

   do
   {
      test = 0;
      loops++;
      if (loops > 200)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	  {	  
	      if ( mpfr_nan_p(appSP[kappa] + j) ) // if appSP[kappa] + j == NAN
	      {
              if (!(_mpfr_vec_scalar_product2(appSP[kappa] + j, appB[kappa], appB[j], n, prec) ) ){
//In this case a heuristic told us that some cancelation probably happened so we just compute the scalar product at full precision
               _F_mpz_vec_scalar_product(ztmp, B->rows[kappa], B->rows[j], n);
               F_mpz_get_mpfr(appSP[kappa] + j, ztmp);
          }
	  }
      
      if (j > zeros + 2)
	  {
	      mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	      mpfr_sub(rtmp, appSP[kappa] + j, tmp, GMP_RNDN);
	      
		  for (k = zeros + 2; k < j - 1; k++)
		  {
		       mpfr_mul(tmp, mu[j] + k, r[kappa] + k, GMP_RNDN);
		       mpfr_sub(rtmp, rtmp, tmp, GMP_RNDN);
		  }
	      
		  mpfr_mul(tmp, mu[j] + j - 1, r[kappa] + j - 1, GMP_RNDN);
	      mpfr_sub(r[kappa] + j, rtmp, tmp, GMP_RNDN);
       } else if (j == zeros+2)
	   {
	      mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	      mpfr_sub(r[kappa] + j, appSP[kappa] + j, tmp, GMP_RNDN);
	   } else mpfr_set(r[kappa] + j, appSP[kappa] + j, GMP_RNDN);

	   mpfr_div(mu[kappa] + j, r[kappa] + j, r[j] + j, GMP_RNDN);
   }
      
   /* **************************** */
   /* Step3--5: compute the X_j's  */
   /* **************************** */
      
   for (j = kappa - 1; j > zeros; j--)
   {
	   /* test of the relaxed size-reduction condition */
       mpfr_abs(tmp, mu[kappa] + j, GMP_RNDN); 
	  
	   if (mpfr_cmp_d(tmp, halfplus) > 0) 
	   {
	      test = 1; 
	      
	      /* we consider separately the cases X = +-1 */     
	      if (mpfr_cmp_d(tmp, onedothalfplus) <= 0)   
		  {
              int sgn = mpfr_sgn(mu[kappa] + j);		  
		      if (sgn >= 0)   /* in this case, X is 1 */
              {
		          for (k = zeros + 1; k < j; k++)
			      {
                     mpfr_sub(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			      }
		      
		          _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		      } else          /* otherwise X is -1 */ 
              {
                  for (k=zeros+1; k<j; k++)
			      {
			         mpfr_add(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			      }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
               }
		   } else   /* we must have |X| >= 2 */
		   {
               mpfr_round(tmp, mu[kappa] + j);
		       for (k = zeros + 1; k < j; k++)
			   {
			       mpfr_mul(rtmp, tmp, mu[j] + k, GMP_RNDN);
			       mpfr_sub(mu[kappa] + k, mu[kappa] + k, rtmp, GMP_RNDN);
			   }

	           if (mpfr_get_exp(tmp) < CPU_SIZE_1 - 2)  
		       {
                   /* X is stored in a long signed int */
                   xx = mpfr_get_si(tmp, GMP_RNDN);		      
                   if (xx > 0L)
                   { 
                      _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                   } else
                   {
                      _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
                   } 
               } else
		       {
                   exponent = F_mpz_set_mpfr_2exp(ztmp, tmp);
                   if (exponent <= 0)
				   {
                       F_mpz_div_2exp(ztmp, ztmp, -exponent);
                       _F_mpz_vec_submul_F_mpz(B->rows[kappa], B->rows[j], n, ztmp);
                   } else
				   {
                      _F_mpz_vec_submul_2exp_F_mpz(B->rows[kappa], B->rows[j], n, ztmp, exponent);
                   }
			    }
		     }
		  }
	   }

       if (test)   /* Anything happened? */
	   {
	      _F_mpz_vec_to_mpfr_vec(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      
		  for (i = zeros + 1; i <= kappa; i++) 
	         mpfr_set_nan(appSP[kappa] + i);//0.0/0.0;
	      
		  for (i = kappa + 1; i <= kappamax; i++) 
	         mpfr_set_nan(appSP[i] + kappa);//0.0/0.0;
	   }
   } while (test);

   if (mpfr_nan_p(appSP[kappa] + kappa)) 
   {
      _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
   }

   mpfr_set(s + zeros + 1, appSP[kappa] + kappa, GMP_RNDN);

   for (k = zeros + 1; k < kappa - 1; k++)
   {
      mpfr_mul( tmp, mu[kappa] + k, r[kappa] + k, GMP_RNDN);
      mpfr_sub( s + k + 1, s + k, tmp, GMP_RNDN);
   }

   mpfr_set(r[kappa] + kappa, s + kappa - 1, GMP_RNDN);

   F_mpz_clear(ztmp);
   F_mpz_clear(X);

   return 0;
}

/* 
   takes the scalar product of two dpe vectors with the same exponents by entry,
   mantissas in vec1 and vec2, n entries to be considered, while cexpo holds the
   exponents.
*/

double d_2exp_vec_scalar_product(double * vec1, double * vec2, int n, int *cexpo)
{
  double sum;

  sum = vec1[0] * vec2[0];
  for (long i = 1; i < n; i++)
     sum += ldexp(vec1[i] * vec2[i], 2*cexpo[i]);

  return sum;
} 

/*  same as above but returning the square of the l_2 norm of a single vector */

double d_2exp_vec_norm(double * vec, int n, int *cexpo)
{
  double sum;

  sum = vec[0] * vec[0];
  for (long i = 1 ; i < n ; i++)
     sum += ldexp(vec[i] * vec[i], 2*cexpo[i]);

  return sum;

} 

/* Computes the largest number of non-zero entries after the diagonal. */

ulong getShift(F_mpz_mat_t B)
{
   ulong n = B->c;
   ulong shift = 0;
   for (ulong i = 0; i < B->r; i++)
   {
      ulong j;
      for (j = n - 1; j >= i + shift + 1 && F_mpz_size(B->rows[i] + j) == 0L; j--);  
      
      if (shift < j - i) shift = j - i;
      
   }

   return shift;
}

/* 
   This is a Babai which is used when size reducing a vector beyond an index which 
   LLL has reached.  cur_kappa is the index behind which we can assume is LLL 
   reduced, while kappa is the vector to be reduced.  This procedure only size reduces
   kappa by vectors up to cur_kappa NOT kappa-1.

*/

int advance_check_Babai (int cur_kappa, int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   ulong temp_expo;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;

   long loops = 0;

   do
   {
      test = 0;

      loops++;
      if (loops > 200)
	  {
         return -1;
      }
            
      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < cur_kappa; j++)
	  {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
	         appSP[kappa][j] = _d_vec_scalar_product(appB[kappa], appB[j], n);
	      }
	  	  
          if (j > zeros + 2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         rtmp = appSP[kappa][j] - tmp;
	         
			 for (k = zeros + 2; k < j - 1; k++)
		     {
		         tmp = mu[j][k] * r[kappa][k];
		         rtmp = rtmp - tmp;
		     }
	         
			 tmp = mu[j][j-1] * r[kappa][j-1];
	         r[kappa][j] = rtmp - tmp;
          } else if (j == zeros+2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         r[kappa][j] = appSP[kappa][j] - tmp;
	      } else r[kappa][j] = appSP[kappa][j];

	      mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      for (j = cur_kappa - 1; j > zeros; j--)
	  {

	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
	  
	      if (tmp > halfplus) 
	      {
	         test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		     {		  
		         if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
                 {
		             for (k = zeros + 1; k < j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
		      
		             _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
                 {
                     for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                     _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
                 }
		     } else   /* we must have |X| >= 2 */
		     {
		         tmp = ldexp (mu[kappa][j] , -exponent);

	             if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		             tmp = rint (tmp); 
		      
		             for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
			            rtmp = ldexp (rtmp, exponent);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		             xx = (long) tmp;
		      
                     if (xx > 0L)
                     { 
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                     } else
                     {
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
                     } 
                  } else
		          {
		             tmp = frexp(mu[kappa][j], &exponent); 
		             tmp = tmp * MAX_LONG;
		             xx = (long) tmp;
		             exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		             /* This case is extremely rare: never happened for me */
		             if (exponent <= 0) 
			         {
//                      printf("rare case kappa = %d, j = %d ******************************************************\n", kappa, j);
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                        {
                           _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
                        } else
                        {
                           _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
                        }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                           rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                        {
                           _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
                        } else
                        {
                           _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
                        }
			            
						for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					    }
				     }	    
			      }
		      } 
		   }
	   }

       if (test == 1)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      
		  for (i = zeros + 1; i <= cur_kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      
		  for (i = cur_kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   } else
       {
          for (i = zeros + 1; i <= cur_kappa; i++)
             appSP[kappa][i] = NAN;
       }
   } while (test == 1);

   if (test == 0)
      return 0;
   else
      return -2;
}

/*
   Same as above but using the heuristic scalar product.
*/

int advance_check_Babai_heuristic_d (int cur_kappa, int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   ulong temp_expo;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;


   int loops = 0;

   do
   {
      test = 0;
            
      loops++;
      if (loops > 200)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < cur_kappa; j++)
	  {	  
	     if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	     {
//### This is different -----
            appSP[kappa][j] = heuristic_scalar_product(appB[kappa], appB[j], n, B, kappa, j, expo[kappa]+expo[j]);
//---------------------------
         }
	  	  
         if (j > zeros + 2)
	     {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         rtmp = appSP[kappa][j] - tmp;
	         
			 for (k = zeros + 2; k < j - 1; k++)
		     {
		        tmp = mu[j][k] * r[kappa][k];
		        rtmp = rtmp - tmp;
		     }
	         
			 tmp = mu[j][j-1] * r[kappa][j-1];
	         r[kappa][j] = rtmp - tmp;
         } else if (j == zeros+2)
	     {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         r[kappa][j] = appSP[kappa][j] - tmp;
	     } else r[kappa][j] = appSP[kappa][j];

	     mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = cur_kappa - 1; j > zeros; j--)
	  {
	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
	  
	      if (tmp > halfplus) 
	      {
	         test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		     {		  
		        if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
                {
		            for (k = zeros + 1; k < j; k++)
			        {
			           tmp = ldexp (mu[j][k], exponent);
			           mu[kappa][k] =  mu[kappa][k] - tmp; 
			        }
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
                 {
                     for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                     _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
                 }
		      } else   /* we must have |X| >= 2 */
		      {
		         tmp = ldexp (mu[kappa][j] , -exponent);

	             if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		             tmp = rint (tmp); 
		      
		             for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
			            rtmp = ldexp (rtmp, exponent);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		             xx = (long) tmp;
		      
                     if (xx > 0L)
                     { 
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                     } else
                     {
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
                     } 
                 } else
		         {
		             tmp = frexp(mu[kappa][j], &exponent); 
		             tmp = tmp * MAX_LONG;
		             xx = (long) tmp;
		             exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		             /* This case is extremely rare: never happened for me */
		             if (exponent <= 0) 
			         {
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                        {
                           _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
                        } else
                        {
                           _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
                        }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                           rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                        {
                           _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
                        } else
                        {
                           _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
                        }
			             
						for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					    }
				     }	    
			      }
		      }
		   }
	   }

       if (test == 1)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      
		  for (i = zeros + 1; i <= cur_kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      
		  for (i = cur_kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   } else
       {
          for (i = zeros + 1; i <= cur_kappa; i++)
             appSP[kappa][i] = NAN;
       }
   } while (test == 1);

   if (test == 0)
      return 0;
   else
      return -2;
}

/*
   Same as above but for mpfr and not doubles.
*/

int advance_check_Babai_heuristic(int cur_kappa, int kappa, F_mpz_mat_t B, __mpfr_struct **mu, __mpfr_struct **r, __mpfr_struct *s, 
       __mpfr_struct **appB, __mpfr_struct **appSP, 
       int a, int zeros, int kappamax, int n, mpfr_t tmp, mpfr_t rtmp, mp_prec_t prec)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   F_mpz_t ztmp, X;
   F_mpz_init(ztmp);
   F_mpz_init(X);
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;

   long loops = 0;

   do
   {
      test = 0;
      loops++;
      if (loops > 200)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < cur_kappa; j++)
	  {	  
	      if ( mpfr_nan_p(appSP[kappa] + j) ) // if appSP[kappa] + j == NAN
	      {
              if (!(_mpfr_vec_scalar_product2(appSP[kappa] + j, appB[kappa], appB[j], n, prec)))
			  {
                 // In this case a heuristic told us that some cancelation probably happened 
				 // so we just compute the scalar product at full precision
              
			     _F_mpz_vec_scalar_product(ztmp, B->rows[kappa], B->rows[j], n);
                 F_mpz_get_mpfr(appSP[kappa] + j, ztmp);
              }
	      }
          
		  if (j > zeros + 2)
	      {
	          mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	          mpfr_sub(rtmp, appSP[kappa] + j, tmp, GMP_RNDN);
	          for (k = zeros + 2; k < j - 1; k++)
		      {
		         mpfr_mul(tmp, mu[j] + k, r[kappa] + k, GMP_RNDN);
		         mpfr_sub(rtmp, rtmp, tmp, GMP_RNDN);
		      }
	          mpfr_mul(tmp, mu[j] + j - 1, r[kappa] + j - 1, GMP_RNDN);
	          mpfr_sub(r[kappa] + j, rtmp, tmp, GMP_RNDN);
          } else if (j == zeros+2)
	      {
	          mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	          mpfr_sub(r[kappa] + j, appSP[kappa] + j, tmp, GMP_RNDN);
	      } else mpfr_set(r[kappa] + j, appSP[kappa] + j, GMP_RNDN);

	      mpfr_div(mu[kappa] + j, r[kappa] + j, r[j] + j, GMP_RNDN);
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = cur_kappa - 1; j > zeros; j--)
	  {
	      /* test of the relaxed size-reduction condition */
          mpfr_abs(tmp, mu[kappa] + j, GMP_RNDN); 
	  
	      if ( mpfr_cmp_d(tmp, halfplus) > 0) 
	      {
	         test = 1; 
	      
	         /* we consider separately the cases X = +-1 */     
	         if (mpfr_cmp_d(tmp, onedothalfplus) <= 0)   
		     {
                 int sgn = mpfr_sgn(mu[kappa] + j);		  
		         if (sgn >= 0)   /* in this case, X is 1 */
                 {
		             for (k = zeros + 1; k < j; k++)
			         {
                        mpfr_sub(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			         }
		      
		             _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
                 {
                     for (k=zeros+1; k<j; k++)
			         {
			            mpfr_add(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			         }
		      
                     _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
                 }
		     } else   /* we must have |X| >= 2 */
		     {
                 mpfr_round(tmp, mu[kappa] + j);
		         for (k = zeros + 1; k < j; k++)
			     {
			         mpfr_mul(rtmp, tmp, mu[j] + k, GMP_RNDN);
			         mpfr_sub(mu[kappa] + k, mu[kappa] + k, rtmp, GMP_RNDN);
			     }

	             if (mpfr_get_exp(tmp) < CPU_SIZE_1 - 2)  
		         {
                     /* X is stored in a long signed int */
                     xx = mpfr_get_si(tmp, GMP_RNDN);		      
                     if (xx > 0L)
                     { 
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                     } else
                     {
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
                     } 
                 } else
		         {
                     exponent = F_mpz_set_mpfr_2exp(ztmp, tmp);
                     if (exponent <= 0)
					 {
                        F_mpz_div_2exp(ztmp, ztmp, -exponent);
                        _F_mpz_vec_submul_F_mpz(B->rows[kappa], B->rows[j], n, ztmp);
                     } else
					 {
                        _F_mpz_vec_submul_2exp_F_mpz(B->rows[kappa], B->rows[j], n, ztmp, exponent);
                     }
			     }
		      }
		   }
	   }

       if (test)   /* Anything happened? */
	   {
	      _F_mpz_vec_to_mpfr_vec(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      
		  for (i = zeros + 1; i <= kappa; i++) 
	         mpfr_set_nan(appSP[kappa] + i);//0.0/0.0;
	      
		  for (i = kappa + 1; i <= kappamax; i++) 
	         mpfr_set_nan(appSP[i] + kappa);//0.0/0.0;
	   }
   } while (test);

   if (mpfr_nan_p(appSP[kappa] + kappa)) 
   {
      _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
   }

   F_mpz_clear(ztmp);
   F_mpz_clear(X);

   return 0;
}

/* 
   Similar to above with a check for zero vectors after each size reduction.
*/

int check_Babai_heuristic_d_zero_vec (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;


   int loops = 0;

   do
   {
      test = 0;
            
      loops++;
      if (loops > 200)
      {
#if TRACE
         F_mpz_mat_print_pretty(B);
         d_mat_print(mu, expo, kappa+1, kappa+1);
#endif
#if TRACE
         printf("kappa = %d \n", kappa);
#endif
		 return -1;
      }

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	  {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
//### This is different -----
             appSP[kappa][j] = heuristic_scalar_product(appB[kappa], appB[j], n, B, kappa, j, expo[kappa]+expo[j]);
//---------------------------
          }
	  	  
          if (j > zeros + 2)
	      {
	          tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	          rtmp = appSP[kappa][j] - tmp;
	          
			  for (k = zeros + 2; k < j - 1; k++)
		      {
		         tmp = mu[j][k] * r[kappa][k];
		         rtmp = rtmp - tmp;
		      }
	          
			  tmp = mu[j][j-1] * r[kappa][j-1];
	          r[kappa][j] = rtmp - tmp;
          } else if (j == zeros+2)
	      {
	          tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	          r[kappa][j] = appSP[kappa][j] - tmp;
	      } else r[kappa][j] = appSP[kappa][j];

	      mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	  {

	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
	  
	      if (tmp > halfplus) 
	      {
	         test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		     {		  
		         if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
                 {
		             for (k = zeros + 1; k < j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
		      
		             _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
                 {
                     for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                     _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
                 }
		     } else   /* we must have |X| >= 2 */
		     {
		         tmp = ldexp (mu[kappa][j] , -exponent);

	             if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		             tmp = rint (tmp); 
		      
		             for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
			            rtmp = ldexp (rtmp, exponent);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		             xx = (long) tmp;
		      
                     if (xx > 0L)
                     { 
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                     } else
                     {
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
                     } 
                 } else
		         {
		             tmp = frexp(mu[kappa][j], &exponent); 
		             tmp = tmp * MAX_LONG;
		             xx = (long) tmp;
		             exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		             /* This case is extremely rare: never happened for me */
		             if (exponent <= 0) 
			         {
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                        {
                           _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
                        } else
                        {
                           _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
                        }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                           rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                        {
                           _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
                        } else
                        {
                           _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
                        }
			            
						for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					    }
				     }	    
			      }
		      }
		  }
	   }

       if (test == 1)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
          if (expo[kappa] != 0)
          {
   	         aa = zeros + 1;
	         
			 for (i = zeros + 1; i <= kappa; i++) 
	            appSP[kappa][i] = NAN;//0.0/0.0;
	         
			 for (i = kappa + 1; i <= kappamax; i++) 
	            appSP[i][kappa] = NAN;//0.0/0.0;
          } else
            test = 10;
	   }
   } while (test == 1);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
      // ### This is different -------
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
   }
   
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }

   if (test == 0)
      return 0;
   else
      return 10;
}

// This is a mildly greedy version, tries the fast version (doubles only) unless that fails then 
// switches to heuristic version for only one loop and right back to fast... reduces B in place!

int LLL_d(F_mpz_mat_t B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   ulong shift = getShift(B);

   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN; //0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int num_failed_fast = 0;
   int babai_ok = 0;
   int heuristic_fail = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   
      if (num_failed_fast < 20)
      {
         babai_ok = check_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n)); 
      }
      else{
         babai_ok = -1;
      }

      if (babai_ok == -1)
      {
         num_failed_fast++;
         heuristic_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n)); 
      }

      if (heuristic_fail == -1)
      {
         // Need to switch to mpfr
         return -1;
      }

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

      if (tmp <= s[kappa-1]) 
	   {
	      alpha[kappa] = kappa;
	      tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	      r[kappa][kappa] = s[kappa-1] - tmp;
	      kappa++;
	   } else
	   {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		     {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
          } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
	         r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   } 
  
   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);

   return 0;
}

/* 
   This LLL reduces B in place and only uses the heuristic inner products 
   which attempt to detect cancellations and otherwise use doubles.
*/

int LLL_d_heuristic(F_mpz_mat_t B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
   {
	   //### This is different -----
       appSP[i][i] = _d_vec_norm(appB[i], n); 
   } while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;

   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n);//FLINT_MIN(kappamax + 1 + shift, n)); 
      if (babai_fail == -1)
         return -1;

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

      if (tmp <= s[kappa-1]) 
	  {
	      alpha[kappa] = kappa;
	      tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	      r[kappa][kappa] = s[kappa-1] - tmp;
	      kappa++;
	  } else
	  {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		     {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
          } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)//fixme should be EPS
	      {
	         zeros++;
	         kappa++;
             //### This is different -----
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
             r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   }

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);

   return 0;
}

/*
   This is LLL using mpfr with the given precision for the underlying GSO,
   reduces B in place.  The mpfr2 refers to the way the mpfr_t's are 
   initialized.
*/
int LLL_mpfr2(F_mpz_mat_t B, mp_prec_t prec)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   __mpfr_struct ** mu, ** r, ** appB, ** appSP;
   __mpfr_struct * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init2(rtmp, prec);
   mpfr_init2(tmp, prec);

   mu = mpfr_mat_init2(d, d, prec);
   r = mpfr_mat_init2(d, d, prec);
   appB = mpfr_mat_init2(d, n, prec);
   appSP = mpfr_mat_init2(d, d, prec);

   s = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));
   appSPtmp = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));

   for (i = 0; i < d+1; i++){
      mpfr_init2(s + i, prec);
      mpfr_init2(appSPtmp + i, prec);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i] + j);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_to_mpfr_vec(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      _mpfr_vec_norm2(appSP[i] + i, appB[i], n, prec); 
   while ( (mpfr_sgn(appSP[i] + i) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) mpfr_set(r[i] + i, appSP[i] + i, GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int babai_fail = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      babai_fail = check_Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, n,  tmp, rtmp, prec);

      if (babai_fail == -1)
         return -1;
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1] + kappa - 1, ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s + kappa - 1) <= 0) 
	  {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa] + kappa - 1, r[kappa] + kappa - 1, GMP_RNDN);
         mpfr_sub(r[kappa] + kappa, s + kappa - 1, tmp, GMP_RNDN);
	      kappa++;
	  } else
	  {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		     {
               mpfr_mul_d(tmp, r[kappa-1] + kappa - 1, ctt, GMP_RNDN);
	         }
          } while ((kappa >= zeros + 2) && (mpfr_cmp(s + kappa - 1,tmp) <= 0));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      mpfr_set(r[kappa] + kappa, s + kappa, GMP_RNDN);
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp + i, appSP[kappa2] + i, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp + i, appSP[i] + kappa2, GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j, GMP_RNDN);	      
	         mpfr_set(appSP[i] + kappa, appSPtmp + i - 1, GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j - 1, GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j] + i, appSP[j] + i - 1, GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa] + i, appSPtmp + i, GMP_RNDN);
	      mpfr_set(appSP[kappa] + kappa, appSPtmp + kappa2, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i] + kappa, appSPtmp + i, GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa] + kappa) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
	         mpfr_set(r[kappa] + kappa, appSP[kappa] + kappa, GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < d+1; i++){
      mpfr_clear(s + i);
      mpfr_clear(appSPtmp + i);
   }


   mpfr_mat_clear(mu, d, d);
   mpfr_mat_clear(r, d, d);
   mpfr_mat_clear(appB, d, n);
   mpfr_mat_clear(appSP, d, d);
   free(s);
   free(appSPtmp);

   return 0;
}

/* 
   A wrapper of LLL_mpfr2.  This currently begins with prec == 53, 
   then for the first 20 loops increases the precision one limb at a time.
   After 20 loops it doubles the precision each time.  There is a proof that
   this will eventually work.
*/

int LLL_mpfr(F_mpz_mat_t B)
{

   mp_prec_t prec;
   prec = 53;

   int result = -1;
   int num_loops = 1;
   while ((result == -1) && (prec < MPFR_PREC_MAX))
   {
      result = LLL_mpfr2(B, prec);

#if TRACE
      printf("called LLL_mpfr with prec = %ld\n", prec);
#endif

      if (result == -1)
      {
         if (num_loops < 20)
            prec = prec + 53;
         else
            prec = prec*2;
         num_loops++;
      }
   }
   
   if (result >= 0)
      return result;
   else
      return -1;
}

/* 
   A wrapper of the above procedures.  Begins with the greediest version, then 
   adapts to heuristic inner products only, then finally to mpfr if needed.
*/

int LLL_wrapper(F_mpz_mat_t B){

   int res = LLL_d(B);
   if (res >= 0)
   { 
	  //hooray worked first time
      return res;
   } else if (res == -1)
   { 
	  //just in case the fast/heuristic switch has any impact
      res = LLL_d_heuristic(B);
   }

   if (res == -1)
   { 
	  //Now try the mpfr version
#if PROFILE
      printf("mpfr called\n");
#endif
      res = LLL_mpfr(B);
   }

   if (res >= 0)
   { 
	  //finally worked
      return res;
   } else // we've got big problems if this is the exit...
      return -1;
}

/**************************

    LLL_with_removals

    The following procedures are similar to the above but they accept an a priori bound for the 
    square length of any vector of interest.  That is to say, if the user is NOT interested in vectors
    with squared l_2 norm > B then use these procedures instead of full LLL.

***************************/



// Same as LLL_d but with the removal bound.  Output is the new dimension of 
// B if removals are desired.

int LLL_d_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
	   expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int num_failed_fast = 0;
   int babai_ok = 0;
   int heuristic_fail = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   
      if (num_failed_fast < 500)
      {
         babai_ok = check_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); //FLINT_MIN(kappamax + 1 + shift, n)); 
      }
      else{
         babai_ok = -1;
      }

      if (babai_ok == -1)
      {
         num_failed_fast++;
         heuristic_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, kappamax, n); 
      }

      if (heuristic_fail == -1)
      {
         // need to switch to mpfr...
         return -1;
      }

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

      if (tmp <= s[kappa-1]) 
	  {
	      alpha[kappa] = kappa;
	      tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	      r[kappa][kappa] = s[kappa-1] - tmp;
	      kappa++;
	  } else
	  {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		     {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
          } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
	         r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   } 

   // Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;

   // ldexp might not be the right choice as we move on... 
   // should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      // d_rii is the G-S length of ith vector divided by 2 (we shouldn't 
	  // make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i], 2*expo[i] - 1);

	  if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }
  
   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);

   return newd;
}

/*
   Same as LLL_d_heuristic but with the removal bound
*/

int LLL_d_heuristic_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN; //0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
   {
	  //### This is different -----
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   } while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); //FLINT_MIN(kappamax + 1 + shift, n)); 
      if (babai_fail == -1)
         return -1;

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

      if (tmp <= s[kappa-1]) 
	  {
	      alpha[kappa] = kappa;
	      tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	      r[kappa][kappa] = s[kappa-1] - tmp;
	      kappa++;
	  } else
	  {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		     {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
          } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)//fixme should be EPS
	      {
	         zeros++;
	         kappa++; 
             //### This is different -----
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
             r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   }

   // newd stuff here... 
   // Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
   
   // ldexp might not be the right choice as we move on... 
   // should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      // d_rii is the G-S length of ith vector divided by 2 (we shouldn't 
	  // make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i], 2*expo[i] - 1);

	  if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);
   
   return newd;
}

/*
   same as LLL_mpfr2 but with the removal bound
   mpfr_init2 not mpfr_init and set_prec not set_default_prec
*/

int LLL_mpfr2_with_removal(F_mpz_mat_t B, mp_prec_t prec, F_mpz_t gs_B)
{
   int kappa, kappa2, d, D, n, i, j, zeros, kappamax;
   __mpfr_struct ** mu, ** r, ** appB, ** appSP;
   __mpfr_struct * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;
   D = d;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init2(rtmp, prec);
   mpfr_init2(tmp, prec);

   mu = mpfr_mat_init2(d, d, prec);
   r = mpfr_mat_init2(d, d, prec);
   appB = mpfr_mat_init2(d, n, prec);
   appSP = mpfr_mat_init2(d, d, prec);

   s = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));
   appSPtmp = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));

   for (i = 0; i < d+1; i++){
      mpfr_init2(s + i, prec);
      mpfr_init2(appSPtmp + i, prec);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i] + j);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_to_mpfr_vec(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      _mpfr_vec_norm2(appSP[i] + i, appB[i], n, prec); 
   while ( (mpfr_sgn(appSP[i] + i) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;

   if (zeros < d - 1) mpfr_set(r[i] + i, appSP[i] + i, GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int babai_fail = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      babai_fail = check_Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, n,  tmp, rtmp, prec);

      if (babai_fail == -1)
         return -1;
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1] + kappa - 1, ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s + kappa - 1) <= 0) 
	  {
	      alpha[kappa] = kappa;
          mpfr_mul(tmp, mu[kappa] + kappa - 1, r[kappa] + kappa - 1, GMP_RNDN);
          mpfr_sub(r[kappa] + kappa, s + kappa - 1, tmp, GMP_RNDN);
	      kappa++;
	  } else
	  {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		     {
               mpfr_mul_d(tmp, r[kappa-1] + kappa - 1, ctt, GMP_RNDN);
	         }
          } while ( (kappa >= zeros + 2) && (mpfr_cmp(s + kappa - 1,tmp) <= 0) );

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      mpfr_set(r[kappa] + kappa, s + kappa, GMP_RNDN);
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp + i, appSP[kappa2] + i, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp + i, appSP[i] + kappa2, GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j, GMP_RNDN);	      
	         mpfr_set(appSP[i] + kappa, appSPtmp + i - 1, GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j - 1, GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j] + i, appSP[j] + i - 1, GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa] + i, appSPtmp + i, GMP_RNDN);
	      mpfr_set(appSP[kappa] + kappa, appSPtmp + kappa2, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i] + kappa, appSPtmp + i, GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa] + kappa) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
	         mpfr_set(r[kappa] + kappa, appSP[kappa] + kappa, GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 

   // newd stuff goes here...
   int ok = 1;
   long newd = d;

   F_mpz_get_mpfr(tmp, gs_B);

   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      // tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      mpfr_set(rtmp, r[i] + i, GMP_RNDN);
      mpfr_div_d(rtmp, rtmp, 8.0, GMP_RNDN);

	  ok = mpfr_cmp(rtmp, tmp);
      if (ok > 0)
	  {
         newd--;
      }
   }
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < D+1; i++){
      mpfr_clear(s + i);
      mpfr_clear(appSPtmp + i);
   }

   mpfr_mat_clear(mu, D, D);
   mpfr_mat_clear(r, D, D);
   mpfr_mat_clear(appB, D, n);
   mpfr_mat_clear(appSP, D, D);
   free(s);
   free(appSPtmp);

   return newd;
}

/* 
   wrapper of the mpfr-based LLL with removal of vectors activated
*/

int LLL_mpfr_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{

   mp_prec_t prec;
   prec = 53;

   int result = -1;
   int num_loops = 1;
   while ((result == -1) && (prec < MPFR_PREC_MAX))
   {
#if PROFILE
      printf("mpfr LLL with prec = %ld\n", prec);
#endif
      result = LLL_mpfr2_with_removal(B, prec, gs_B);
      if (result == -1){
         if (num_loops < 20)
            prec = prec + 53;
         else
            prec = prec*2;
         num_loops++;
      }
   }
   
   if (result >= 0)
      return result;
   else
      return -1;
}

/* 
   Wrapper of the base case LLL with the addition of the removal boundary.
*/

int LLL_wrapper_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int res = LLL_d_with_removal(B, gs_B);
   if (res >= 0)
   { 
	  //hooray worked first time
      return res;
   } else if (res == -1) //just in case the fast/heuristic switch has any impact
      res = LLL_d_heuristic_with_removal(B, gs_B);

   if (res == -1)
   { 
	  //Now try the mpfr version
#if PROFILE
      printf("using mpfr\n");
#endif
      res = LLL_mpfr_with_removal(B, gs_B);
   }

   if (res >= 0) //finally worked
      return res;
   else //we've got big problems if this is the exit...
      return -1;
}

/* 
   LLL specialized to knapsack-type lattices, it performs early size reductions 
   occasionally which makes things faster in the knapsack case, also uses the 
   removal bounds.
*/

int knapsack_LLL_wrapper_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int res = knapsack_LLL_d_with_removal(B, gs_B);
   if (res >= 0)
   { 
	  // hooray worked first time
      return res;
   } else if (res == -1) //just in case the fast/heuristic switch has any impact
      res = LLL_d_heuristic_with_removal(B, gs_B);

   if (res == -1)
   { 
	  // Now try the mpfr version
#if PROFILE
      printf("called mpfr!!\n");
#endif
      res = LLL_mpfr_with_removal(B, gs_B);
   }

   if (res >= 0) // finally worked
      return res;
   else // we've got big problems if this is the exit...
      return -1;
}

// Engine of the knapsack_LLL This is a mildly greedy version, 
// tries the fast version unless that fails then 
// switches to heuristic version for only one loop and right back to fast... 

int knapsack_LLL_d_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, D, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;

   int copy_kappa, copy_kappamax;
   double ** copy_mu, ** copy_r, ** copy_appB, ** copy_appSP;
   double * copy_s;
   int * copy_expo, * copy_alpha;

   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;

   n = B->c;
   d = B->r;
   D = d;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   copy_alpha = (int *) malloc(d * sizeof(int)); 
   copy_expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   copy_mu = d_mat_init(d, d);
   copy_r = d_mat_init(d, d);
   copy_appB = d_mat_init(d, n);
   copy_appSP = d_mat_init(d, d);


   s = (double *) malloc (d * sizeof(double));
   copy_s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN; //0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int num_failed_fast = 0;
   int babai_ok = 0;
   int heuristic_fail = 0;
   long new_kappa, newvec, newvec_max;

   int last_vec = 0;
      newvec = 0;
      newvec_max = 1;
   long num_loops = 0;

   while (kappa < d)
   {
      num_loops++;
      last_vec = 0;
      if (kappa == d - 1){
         last_vec = 1;
      }

      new_kappa = 0;
      if (kappa > kappamax)
      {
         // In the first time we hit a new kappa we're going to size-reduce in advance...
         kappamax = kappa; 
         newvec++;

         if (newvec > newvec_max){
            newvec_max = newvec_max * 2;
            newvec = 0;
            new_kappa = 1;
         }
      }

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */ 
  
      if (num_failed_fast < 150)
      {
#if PROFILE
         babai_start = get_cycle_counter();
#endif
         babai_ok = check_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); 
#if PROFILE
         babai_stop = get_cycle_counter();
         babai_total = babai_total + babai_stop - babai_start;
#endif
      }
      else
      {
         babai_ok = -1;
      }

      if (babai_ok == -1)
      {
#if PROFILE
         hbabai_start = get_cycle_counter();
#endif
         num_failed_fast++;
         heuristic_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, kappamax, n); 
#if PROFILE
         hbabai_stop = get_cycle_counter();
         hbabai_total = hbabai_total + hbabai_stop - hbabai_start;
#endif
      }


      if (heuristic_fail == -1)
      {
         free(alpha);
         free(expo);
         free(copy_alpha);
         free(copy_expo);
         d_mat_clear(mu);
         d_mat_clear(r);
         d_mat_clear(appB);
         d_mat_clear(appSP);
         d_mat_clear(copy_mu);
         d_mat_clear(copy_r);
         d_mat_clear(copy_appB);
         d_mat_clear(copy_appSP);
         free(s);
         free(copy_s);
         free(appSPtmp);
         
		 // need to switch to mpfr...
         return -1;
      }

      // End of the real Babai part...
      if (new_kappa == 1)
      {
         // Perhaps a bit naive but just making a copy of everything but B, running ahead
         // to kappa = d, without upsetting LLL...
         copy_kappa = kappa + 1;
         copy_kappamax = copy_kappa;

         for (copy_kappa = d-1; copy_kappa >  kappa; copy_kappa--)
         {
#if PROFILE
            adv_babai_start = get_cycle_counter();
#endif
            babai_ok = advance_check_Babai(kappa, copy_kappa, B, mu, r, copy_s, appB, expo, appSP, alpha[copy_kappa], zeros, copy_kappamax, n);
#if PROFILE
            adv_babai_stop = get_cycle_counter();
            adv_babai_total = adv_babai_total + adv_babai_stop - adv_babai_start;
#endif

            heuristic_fail = 0;
            if (babai_ok == -1)
            {
#if PROFILE
               hadv_babai_start = get_cycle_counter();
#endif
               heuristic_fail = advance_check_Babai_heuristic_d(kappa, copy_kappa, B, mu, r, copy_s, appB, expo, appSP, alpha[copy_kappa], zeros, copy_kappamax, n);
#if PROFILE
               hadv_babai_stop = get_cycle_counter();
               hadv_babai_total = hadv_babai_total + hadv_babai_stop - hadv_babai_start;
#endif
            }
         }
      }

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  

      if (last_vec == 12)
	  {
	     tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
         d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
         d_gs_B = ldexp( d_gs_B, exp);
         d_rii = ldexp(s[kappa-1] - tmp, 2*expo[kappa] - 1);
#if PROFILE
         printf("last vector reached\n kappa = d-1 kappa = %ld\n", kappa);
         printf("%5f r[%d] and gs_B = %5f\n", d_rii, d-1, d_gs_B);
#endif
         if (d_rii > d_gs_B)
		 {
#if PROFILE
            r[kappa][kappa] = s[kappa-1] - tmp;
            for (i = d-1; (i >= 0); i--)
            {
               // d_rii is the G-S length of ith vector divided by 2 
			   // (we shouldn't make a mistake and remove something valuable)
               d_rii = ldexp(r[i][i],  2*expo[i] - 1);
               printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
               if (d_rii < 1) printf("r[i][i] = %5f and expo = %ld\n", r[i][i], 2*expo[i] - 1);
               if (d_rii > d_gs_B) printf("big... \n"); 
            }
#endif
            break;
         }
      }

      tmp = r[kappa-1][kappa - 1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa - 1] - expo[kappa]));

      if (tmp <= s[kappa - 1]) 
      {
          alpha[kappa] = kappa;
	      tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	      r[kappa][kappa] = s[kappa-1] - tmp;
	      kappa++;
	  } else
	  {

          /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		      {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
          } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
	         r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   } 

#if PROFILE
   printf(" num loops = %ld \n", num_loops);
#endif

   // Use the newd stuff here...
   // ldexp might not be the right choice as we move on... 
   // should make a straight d_2exp comparison
   ok = 1;
   newd = d;
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      // d_rii is the G-S length of ith vector divided by 2 
	  // (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i], 2*expo[i] - 1);
#if PROFILE
      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
#endif
      if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }
  
   free(alpha);
   free(expo);
   free(copy_alpha);
   free(copy_expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   d_mat_clear(copy_mu);
   d_mat_clear(copy_r);
   d_mat_clear(copy_appB);
   d_mat_clear(copy_appSP);
   free(s);
   free(copy_s);
   free(appSPtmp);

   return newd;
}

//  Same as previous LLL_d_heuristic_with_removal but with advanced size-reduction

int knapsack_LLL_d_heuristic_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN; //0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
   {
	  //### This is different -----
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   } while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; 

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, 
		              alpha[kappa], zeros, kappamax, n); 
      
	  if (babai_fail == -1)
         return -1;

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

      if (tmp <= s[kappa-1]) 
	  {
	      alpha[kappa] = kappa;
	      tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	      r[kappa][kappa] = s[kappa-1] - tmp;
	      kappa++;
	  } else
	  {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		      {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
          } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0) //fixme should be EPS
	      {
	         zeros++;
	         kappa++;
             //### This is different -----
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
             r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   }

   // newd stuff here... 
   // Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
   
   // ldexp might not be the right choice as we move on... 
   // should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      // d_rii is the G-S length of ith vector divided by 2 
	  // (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i], 2*expo[i] - 1);
#if PROFILE
      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
#endif
      if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);

   return newd;
}

/*
   Same as LLL_mpfr2_with_removal but with advanced size-reduction
*/
int knapsack_LLL_mpfr2_with_removal(F_mpz_mat_t B, mp_prec_t prec, F_mpz_t gs_B)
{
   int kappa, kappa2, d, D, n, i, j, zeros, kappamax;
   __mpfr_struct ** mu, ** r, ** appB, ** appSP;
   __mpfr_struct * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;
   D = d;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5 ;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init2(rtmp, prec);
   mpfr_init2(tmp, prec);

   mu = mpfr_mat_init2(d, d, prec);
   r = mpfr_mat_init2(d, d, prec);
   appB = mpfr_mat_init2(d, n, prec);
   appSP = mpfr_mat_init2(d, d, prec);

   s = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));
   appSPtmp = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));

   for (i = 0; i < d+1; i++){
      mpfr_init2(s + i, prec);
      mpfr_init2(appSPtmp + i, prec);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i] + j); //0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
	   _F_mpz_vec_to_mpfr_vec(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      _mpfr_vec_norm2(appSP[i] + i, appB[i], n, prec); 
   while ( (mpfr_sgn(appSP[i] + i) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) mpfr_set(r[i] + i, appSP[i] + i, GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int babai_fail = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax = kappa; 

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      babai_fail = check_Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, n,  tmp, rtmp, prec);

      if (babai_fail == -1)
         return -1;
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1] + kappa - 1, ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s + kappa - 1) <= 0) 
	   {
	      alpha[kappa] = kappa;
          mpfr_mul(tmp, mu[kappa] + kappa - 1, r[kappa] + kappa - 1, GMP_RNDN);
          mpfr_sub(r[kappa] + kappa, s + kappa - 1, tmp, GMP_RNDN);
	      kappa++;
	   } else
	   {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		      {
               mpfr_mul_d(tmp, r[kappa-1] + kappa - 1, ctt, GMP_RNDN);
	         }
          } while ( (kappa >= zeros + 2) && (mpfr_cmp(s + kappa - 1,tmp) <= 0) );

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      mpfr_set(r[kappa] + kappa, s + kappa, GMP_RNDN);
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;


	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp + i, appSP[kappa2] + i, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp + i, appSP[i] + kappa2, GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j, GMP_RNDN);	      
	         mpfr_set(appSP[i] + kappa, appSPtmp + i - 1, GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j - 1, GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j] + i, appSP[j] + i - 1, GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa] + i, appSPtmp + i, GMP_RNDN);
	      mpfr_set(appSP[kappa] + kappa, appSPtmp + kappa2, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i] + kappa, appSPtmp + i, GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa] + kappa) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
	         mpfr_set(r[kappa] + kappa, appSP[kappa] + kappa, GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 

   // newd stuff goes here...
   int ok = 1;
   long newd = d;

   F_mpz_get_mpfr(tmp, gs_B);

   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      // tmp_gs is the G-S length of ith vector divided by 2 
	  // (we shouldn't make a mistake and remove something valuable)
      mpfr_set(rtmp, r[i] + i, GMP_RNDN);
      mpfr_div_d(rtmp, rtmp, 8.0, GMP_RNDN);

	  ok = mpfr_cmp(rtmp, tmp);
      if (ok > 0)
	  {
         newd--;
      }
   }
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < D+1; i++)
   {
      mpfr_clear(s + i);
      mpfr_clear(appSPtmp + i);
   }

   mpfr_mat_clear(mu, D, D);
   mpfr_mat_clear(r, D, D);
   mpfr_mat_clear(appB, D, n);
   mpfr_mat_clear(appSP, D, D);
   free(s);
   free(appSPtmp);

   return newd;
}

/*
   A wrapper of LLL with mpfr GSO and removal and advanced size-reduction
*/

int knapsack_LLL_mpfr_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
#if TRACE
   printf("******************warning mpfr**************\n");
#endif

   mp_prec_t prec;
   prec = 53;

   int result = -1;
   int num_loops = 1;
   while ((result == -1) && (prec < MPFR_PREC_MAX))
   {
      result = knapsack_LLL_mpfr2_with_removal(B, prec, gs_B);
      if (result == -1){
         if (num_loops < 20)
            prec = prec + 53;
         else
            prec = prec*2;
         num_loops++;
      }
   }
   
   if (result >= 0)
      return result;
   else
      return -1;
}

/*
   A wrapper of the various knapsack_LLL's starts with fastest moves to heuristic and finally mpfr
*/

int knapsack_LLL_with_removal(F_mpz_mat_t B, F_mpz_t gs_B){

   int res = knapsack_LLL_d_with_removal(B, gs_B);
   if (res >= 0) //hooray worked first time
      return res;
   else if (res == -1) //just in case the fast/heuristic switch has any impact
      res = knapsack_LLL_d_heuristic_with_removal(B, gs_B);

   if (res == -1) //Now try the mpfr version
      res = knapsack_LLL_mpfr_with_removal(B, gs_B);

   if (res >= 0) //finally worked
      return res;
   else //we've got big problems if this is the exit...
      return -1;
}

/*
   An LLL adapted to searching for zero vectors when not full rank.
*/


int LLL_d_zero_vec_heuristic_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = (4*DELTA + 1)/5;
   halfplus = (4*ETA + .5)/5;
   onedothalfplus = 1.0+halfplus;
	
   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN; //0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
   {
	  //### This is different -----
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   } while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;

   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; 

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d_zero_vec(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); 
      if (babai_fail == -1)
         return -1;
      
#if TRACE
	  if (babai_fail == 10)
      {
         // This means that the kappa^th vector is a 0-vector... 
         printf(" found a 0 vector!\n");

      }
#endif

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

      if (tmp <= s[kappa-1]) 
	  {
	      alpha[kappa] = kappa;
	      tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	      r[kappa][kappa] = s[kappa-1] - tmp;
	      kappa++;
	  } else
	  {

	      /* ******************************************* */
	      /* Step5: Find the right insertion index kappa */
          /* kappa2 remains the initial kappa            */
	      /* ******************************************* */  

	      kappa2 = kappa;
	      do
	      {
	         kappa--;
	         if (kappa > zeros + 1) 
		      {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
          } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

          for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
          for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
          B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0) //fixme should be EPS
	      {
	         zeros++;
	         kappa++;
             //### This is different -----
 	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
             r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   }

   // newd stuff here... 
   // Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
   
   // ldexp might not be the right choice as we move on... 
   // should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      //d_rii is the G-S length of ith vector divided by 2 
	  //(we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i], 2*expo[i] - 1);

	  if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);

   return newd;
}

/*
   A wrapper of the zero vector hunting LLL for B of not full rank
*/

int LLL_wrapper_zero_vec_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int res;

   res = LLL_d_zero_vec_heuristic_with_removal(B, gs_B);

   if (res == -1)
   { 
	  // Now try the mpfr version
      abort();
      res = LLL_mpfr_with_removal(B, gs_B);
   }

   if (res >= 0) // finally worked
      return res;
   else // we've got big problems if this is the exit...
      return -1;
}

/**

   Stripped down LLL_d to just calculate G-S lengths.  This only computes GSOs.

**/

void gs_Babai(int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
      for (j = aa; j < kappa; j++)
	   {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
            //### This is different -----
            appSP[kappa][j] = heuristic_scalar_product(appB[kappa], appB[j], n, B, kappa, j, expo[kappa]+expo[j]);
          }
	  	  
          if (j > zeros + 2)
	      {
	         tmp = mu[j][zeros + 1] * r[kappa][zeros + 1];
	         rtmp = appSP[kappa][j] - tmp;
	         
			 for (k = zeros + 2; k < j - 1; k++)
		     {
		         tmp = mu[j][k] * r[kappa][k];
		         rtmp = rtmp - tmp;
		     }
	         
			 tmp = mu[j][j-1] * r[kappa][j-1];
	         r[kappa][j] = rtmp - tmp;
          } else if (j == zeros+2)
	      {
	         tmp = mu[j][zeros + 1] * r[kappa][zeros + 1];
	         r[kappa][j] = appSP[kappa][j] - tmp;
	      } else r[kappa][j] = appSP[kappa][j];

	      mu[kappa][j] = r[kappa][j] / r[j][j];
      }

      for (j = kappa - 1; j > zeros; j--)
	   {
	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
      }
      
   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
      // ### This is different -------
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
   }

   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k + 1] = s[k] - tmp;
   }
}

//A fast and dirty approximation algorithm for the gs vectors of B.

ulong F_mpz_mat_gs_d( F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;
	
   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN; //0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  

   i = 0; 
  
   do
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;   

   while (kappa < d)
   {
       if (kappa > kappamax) kappamax = kappa; 

       gs_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); 
	   alpha[kappa] = kappa;
	   tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	   r[kappa][kappa] = s[kappa-1] - tmp;
	   kappa++;
   }

   ulong ok, newd, exp;
   double d_gs_B, d_rii; 
   ok = 1;
   newd = d;
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      // d_rii is the G-S length of ith vector divided by 2 
	  // (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i], 2*expo[i] - 1);
      
	  if (d_rii != NAN)
	  {
         if (d_rii > d_gs_B) newd--;
         else (ok = 0);
      } else 
		 (ok = 0);
   }

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);

   return newd;
}

/***************************************

   U_LLL is a new style of LLL which does adjoins an identity matrix to the input 
   lattice then scales the lattice down to new_size bits and reduces this augmented
   lattice.  This tends to be more stable numerically than traditional LLL which means
   higher dimensions can be attacked using doubles.  In each iteration a new identity matrix is adjoined
   is used.  Optimized for factoring polynomials.

****************************************/

int U_LLL_with_removal(F_mpz_mat_t FM, long new_size, F_mpz_t gs_B)
{

   long r, c, bits, i, j;
   int full_prec = 1;
   int done = 0;

    clock_t timer1, timer2;
    timer1 = clock();

#if PROFILE
   clock_t lll_start, lll_stop, lll_total, sum_start, sum_stop;
#endif

   int is_U_I;

#if PROFILE
   lll_total = 0;
   sum_start = clock();
#endif

   r = FM->r;
   c = FM->c;
   bits = FLINT_ABS(F_mpz_mat_max_bits(FM));

   F_mpz_mat_t U;

   F_mpz_mat_t I;
   F_mpz_mat_init_identity(I, r);

   F_mpz_mat_t full_U;
   F_mpz_mat_init_identity(full_U, r);

   F_mpz_mat_t big_FM;
   F_mpz_mat_init(big_FM, r, c + r);

   F_mpz_mat_t full_data;
   F_mpz_mat_init(full_data, r, c);

   F_mpz_mat_t trunc_data;
   F_mpz_mat_init(trunc_data, r, c);

   long mbits;

   int k = 1;

   int newd;
   long prev_mbits = bits;

   if (bits > new_size)
   {
      full_prec = 0;
      
	  //do some truncating
      for ( i = 0; i < r; i++)
         for ( j = 0; j < c; j++)
            F_mpz_set(full_data->rows[i]+j, FM->rows[i]+j);

      mbits = FLINT_ABS(F_mpz_mat_max_bits(full_data));

      if (mbits - new_size > 0)
	  {
         F_mpz_mat_resize(trunc_data, full_data->r, full_data->c);
		 F_mpz_mat_div_2exp(trunc_data, full_data, (ulong) (mbits - new_size));
         
		 // Make this iterate over i and j, make a LARGE lattice which has identity 
		 // in one corner and FM in the other
         for ( i = 0; i < r; i++)
		 {
            for (j = 0; j < i; j++)
               F_mpz_set_ui(big_FM->rows[i] + j, 0L);
            F_mpz_set_ui(big_FM->rows[i] + i, 1L);
            
			for (j = i + 1; j < r; j++)
               F_mpz_set_ui(big_FM->rows[i] + j, 0L);
            
			for (j = r; j < r + c; j++)
               F_mpz_set(big_FM->rows[i] + j, trunc_data->rows[i] + j - r);
         }
      } else
	  {
         full_prec = 1;
      }
   }

   while (done == 0)
   {
      k++;
      if (full_prec == 0)
	  {
#if PROFILE
         lll_start = clock();
#endif
		 knapsack_LLL_wrapper_with_removal(big_FM, gs_B);
#if PROFILE
         lll_stop = clock();
#endif
      } else
	  {
#if PROFILE
         lll_start = clock();
#endif
		 newd = knapsack_LLL_wrapper_with_removal(FM, gs_B);
#if PROFILE
         lll_stop = clock();
#endif
      }

#if PROFILE
      lll_total = lll_total + lll_stop - lll_start;
#endif

      if (full_prec == 1)
         done = 1;
      else 
	  {
         // add more bits
         F_mpz_mat_window_init(U, big_FM, 0, 0, big_FM->r, r);

#if PROFILE
         printf("U bits == %ld\n", FLINT_ABS(F_mpz_mat_max_bits(U)));
#endif

		 is_U_I = F_mpz_mat_equal(U, I);

         // do some truncating
         F_mpz_mat_mul_classical(full_data, U, full_data);

         mbits = FLINT_ABS(F_mpz_mat_max_bits(full_data));
         fprintf(stderr, "mbits is %ld\n", mbits);
         timer2 = clock();
         fprintf(stderr, " spent a total of %f seconds in ULLL\n", (double) (timer2 - timer1) / (double) CLOCKS_PER_SEC);

         if (global_flag > 0){

            F_mpz_mat_print_pretty(full_data);
            fflush(stdout);
            global_flag = 0;

         }
         
         // make this condition better?
         if ((mbits - new_size > 0) &&  (mbits <= prev_mbits - (long)(new_size/4)) && is_U_I == 0)
		 {
            F_mpz_mat_div_2exp(trunc_data, full_data, (ulong) (mbits - new_size));
         } else
		 {
            full_prec = 1;
         }

         prev_mbits = mbits;

         if (full_prec == 1)
		 {
            // can switch to FM, no need for a new identity
            for ( i = 0; i < r; i++){
               for (j = 0; j < c; j++)
                  F_mpz_set(FM->rows[i]+j, full_data->rows[i] + j);
            }
         } else
		 {
            // keep with the big_FM concept
            for ( i = 0; i < r; i++)
			{
               for (j = 0; j < i; j++)
                  F_mpz_set_ui(big_FM->rows[i]+j, 0L);
               F_mpz_set_ui(big_FM->rows[i]+i, 1L);
               
			   for (j = i+1; j < r; j++)
                  F_mpz_set_ui(big_FM->rows[i]+j, 0L);
               
			   for (j = r; j < r+c; j++)
                  F_mpz_set(big_FM->rows[i]+j, trunc_data->rows[i] + j-r);
            }
         }
         
		 F_mpz_mat_window_clear(U);
      }

   }

#if PROFILE
   sum_stop = clock();

   printf(" spent a total of %3f seconds on regular Babai\n", (double) babai_total / 2.4E9);
   printf(" of which %3f cycles spent doing ldexps\n", (double) ldexp_total / 2.4E9);
   printf(" of which %3f cycles spent updating full precision B\n", (double) update_total / 2.4E9);
   printf(" of which %3f cycles spent converting full precision B\n", (double) convert_total / 2.4E9);
   printf(" of which %3f cycles spent computing inner products\n", (double) inner_total /2.4E9);
   printf(" spent a total of %3f seconds on advanced Babai\n", (double) adv_babai_total /2.4E9);

   printf(" spent a total of %3f seconds on h regular Babai\n", (double) hbabai_total / 2.4E9);
   printf(" of which %3f cycles spent doing h ldexps\n", (double) hldexp_total / 2.4E9);
   printf(" of which %3f cycles spent updating h full precision B\n", (double) hupdate_total / 2.4E9);
   printf(" of which %3f cycles spent converting h full precision B\n", (double) hconvert_total / 2.4E9);
   printf(" of which %3f cycles spent computing h inner products\n", (double) hinner_total /2.4E9);
   printf(" spent a total of %3f seconds on h advanced Babai\n", (double) hadv_babai_total /2.4E9);

   printf(" spent a total of %f seconds on inner LLL\n", (double) lll_total / (double)CLOCKS_PER_SEC);

   printf(" spent a total of %f seconds in ULLL\n", (double) (sum_stop - sum_start) / (double) CLOCKS_PER_SEC);
#endif

   F_mpz_mat_clear(full_data);
   F_mpz_mat_clear(trunc_data);
   F_mpz_mat_clear(big_FM);
   F_mpz_mat_clear(I);
   F_mpz_mat_clear(full_U);

   return newd;
}

/* 
   A wrapper of U_LLL_with_removal the most numerically stable LLL.  This is the default LLL
   which reduced B in place.
*/


void LLL(F_mpz_mat_t B)
{
   F_mpz_t temp;
   F_mpz_init(temp);

   U_LLL_with_removal(B, 250L, temp);

   F_mpz_clear(temp);

}
