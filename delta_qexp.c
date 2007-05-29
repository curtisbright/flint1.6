/*
   Demo FLINT program for computing the q-expansion of the delta function.
   
   (C) 2007 David Harvey
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"


void print_poly(Zpoly_mpn_t x_mpn)
{
   Zpoly_t x;
   Zpoly_init2(x, x_mpn->length);
   _Zpoly_mpn_convert_out(x, x_mpn);
   Zpoly_print(stdout, x);
   Zpoly_clear(x);
   /*for (int i = 0; i < x->length; i++)
   {
      if ((signed long) x->coeffs[i*(x->limbs+1)] > 0)
         printf("+%d ", x->coeffs[i*(x->limbs+1)+1]);
      else if ((signed long) x->coeffs[i*(x->limbs+1)] < 0)
         printf("-%d ", x->coeffs[i*(x->limbs+1)+1]);
      else
         printf("0 ");
   }
   printf("\n");*/
}


int main(int argc, char* argv[])
{
   if (argc != 2)
   {
      printf("Syntax: delta_qexp <integer>\n");
      printf("where <integer> is the number of terms to compute\n");
      return 0;
   }
   
   // number of terms to compute
   long n = atoi(argv[1]);

   // initialise polynomial objects
   Zpoly_mpn_t F, F2, F4, F8;
   Zpoly_mpn_init(F, 2*n, 1);
   Zpoly_mpn_init(F2, 2*n, 1);
   Zpoly_mpn_init(F4, 2*n, 2);
   Zpoly_mpn_init(F8, 2*n, 2);
   
   F->length = n;
   for (long i = 0; i < n; i++)
      _Zpoly_mpn_set_coeff_ui(F, i, 0);
      
   // set F := whatever it's supposed to be
   for (long i = 0; 1; i++)
   {
      long index = i * (i+1) / 2;
      if (index >= n)
         break;
      _Zpoly_mpn_set_coeff_si(F, index, (i & 1) ? -(2*i+1) : (2*i+1));
   }

   print_poly(F); printf("\n");

   // compute F^8, truncated to length n
   _Zpoly_mpn_mul_KS(F2, F, F);
   F2->length = n;
   print_poly(F2); printf("\n");

   _Zpoly_mpn_mul_KS(F4, F2, F2);
   F4->length = n;
   print_poly(F4);  printf("\n");  

   _Zpoly_mpn_mul_KS(F8, F4, F4);
   F8->length = n;
   print_poly(F8);  printf("\n");  
   
   // print out last coefficient (or at least one word of it)
//   printf("coefficient of q^%d is %d\n", n-1, _Zpoly_mpn_get_coeff_ui(Fpow, n-1));
   
   // clean up
   Zpoly_mpn_clear(F8);
   Zpoly_mpn_clear(F4);
   Zpoly_mpn_clear(F2);
   Zpoly_mpn_clear(F);
   return 0;
}


// Here's what it *should* look like:
/*
int main()
{
   if (argc != 2)
   {
      printf("Syntax: delta_qexp <integer>\n");
      printf("where <integer> is the number of terms to compute\n");
      return 0;
   }
   
   // number of terms to compute
   long n = atoi(argv[1]);

   // initialise polynomial objects
   Zpoly_mpn_t F, Fpow;
   Zpoly_mpn_init(F);
   Zpoly_mpn_init(Fpow);

   // set F := whatever it's supposed to be
   for (long i = 0; 1; i++)
   {
      long index = i * (i+1) / 2;
      if (index >= n)
         break;
      Zpoly_mpn_set_coeff_si(F, index, (i & 1) ? -(2*i+1) : (2*i+1));
   }

   // compute F^8, truncated to length n
   Zpoly_mpn_mul(Fpow, F, F);
   Zpoly_mpn_truncate(Fpow, n);

   Zpoly_mpn_mul(Fpow, Fpow, Fpow);
   Zpoly_mpn_truncate(Fpow, n);

   Zpoly_mpn_mul(Fpow, Fpow, Fpow);
   Zpoly_mpn_truncate(Fpow, n);
   
   // print out last coefficient
   mpz_t x;
   mpz_init(x);
   Zpoly_mpn_get_coeff_mpz(Fpow, x, n-1);
   
   gmp_printf("coefficient of q^%d is %Z\n", n-1, x);
   
   // clean up
   mpz_clear(x);
   Zpoly_mpn_clear(Fpow);
   Zpoly_mpn_clear(F);
   return 0;
}
*/
