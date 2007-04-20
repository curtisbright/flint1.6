/****************************************************************************

Zpoly_mpn-test.c: Test code for Zpoly_mpn.c and Zpoly_mpn.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"
#include "flint-manager.h"

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

gmp_randstate_t Zpoly_test_randstate;

unsigned long randint(unsigned long randsup) 
{
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

void randpoly(Zpoly_t pol, unsigned long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   if (pol->coeffs) Zpoly_clear(pol);
   Zpoly_init3(pol, length, maxbits);
   for (unsigned long i = 0; i < length; i++)
   {
       bits = randint(maxbits);
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, Zpoly_test_randstate, bits);
          if (randint(2)) mpz_neg(temp,temp);
       }
       Zpoly_set_coeff(pol, i, temp);
       
   }
   
   mpz_clear(temp);
} 

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
int test_Zpoly_mpn_convert()
{
   Zpoly_t test_poly, test_poly2;
   Zpoly_mpn_t test_mpn_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          Zpoly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits); 
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          _Zpoly_mpn_convert_out(test_poly2, test_mpn_poly);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = Zpoly_equal(test_poly, test_poly2);
      }   
      Zpoly_mpn_clear(test_mpn_poly);
   }
   
   mpz_clear(temp);
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   
   return result;
}

int test_Zpoly_mpn_getset_ui()
{
   Zpoly_t test_poly;
   Zpoly_mpn_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long coeff, coeff_bits, coeff_num;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS_PER_LIMB);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(Zpoly_test_randstate, coeff_bits);
              coeff_num = randint(length);
#if DEBUG
              printf("Index = %ld, bits = %ld, coeff = %ld\n", coeff_num, coeff_bits, coeff);
#endif
              _Zpoly_mpn_set_coeff_ui(test_mpn_poly, coeff_num, coeff);
              result = (_Zpoly_mpn_get_coeff_ui(test_mpn_poly, coeff_num) == coeff);
          }
      }
      Zpoly_mpn_clear(test_mpn_poly);
   }
   
   return result; 
}

int test_Zpoly_mpn_getset_si()
{
   Zpoly_t test_poly;
   Zpoly_mpn_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS_PER_LIMB-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(Zpoly_test_randstate, coeff_bits);
              coeff_num = randint(length);
#if DEBUG
              printf("Index = %ld, bits = %ld, coeff = %ld\n", coeff_num, coeff_bits, coeff);
#endif
              if (randint(2)) sign = -1L; else sign = 1L;
              coeff = sign*coeff;
              _Zpoly_mpn_set_coeff_si(test_mpn_poly, coeff_num, coeff);
              result = ((_Zpoly_mpn_get_coeff_si(test_mpn_poly, coeff_num) == coeff) && (_Zpoly_mpn_get_coeff_ui(test_mpn_poly, coeff_num) == sign*coeff));
          }
      }
      Zpoly_mpn_clear(test_mpn_poly);
   }
   
   return result; 
}

int test_Zpoly_mpn_get_coeff_ptr()
{
   Zpoly_t test_poly;
   Zpoly_mpn_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS_PER_LIMB-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(Zpoly_test_randstate, coeff_bits);
              coeff_num = randint(length);
#if DEBUG
              printf("Index = %ld, bits = %ld, coeff = %ld\n", coeff_num, coeff_bits, coeff);
#endif
              if (randint(2)) sign = -1L; else sign = 1L;
              coeff = sign*coeff;
              _Zpoly_mpn_set_coeff_si(test_mpn_poly, coeff_num, coeff);
              if (coeff == 0) sign = 0;
              result = (_Zpoly_mpn_get_coeff_ptr(test_mpn_poly, coeff_num)[0] == sign);
          }
      }
      Zpoly_mpn_clear(test_mpn_poly);
   }
   
   return result; 
}

int test_Zpoly_mpn_normalise()
{
   Zpoly_t test_poly;
   Zpoly_mpn_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long nz_coeff;
   long sign;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          
          nz_coeff = randint(length+1)-1;
          if (randint(2)) sign = -1L; else sign = 1;
          if (nz_coeff != -1L) _Zpoly_mpn_set_coeff_si(test_mpn_poly, nz_coeff, sign*1000);
          for (unsigned long i = nz_coeff+1; i < length; i++)
          _Zpoly_mpn_set_coeff_ui(test_mpn_poly, i, 0);
              
          _Zpoly_mpn_normalise(test_mpn_poly);
#if DEBUG
          printf("length = %ld, nonzero coefficient = %ld\n",_Zpoly_mpn_length(test_mpn_poly), nz_coeff);
#endif              
          result = (_Zpoly_mpn_length(test_mpn_poly) == nz_coeff+1);
      }
      Zpoly_mpn_clear(test_mpn_poly);
   }
   
   return result; 
}

int test_Zpoly_mpn_getset_coeff()
{
   Zpoly_t test_poly;
   Zpoly_mpn_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          
          mp_limb_t * coeff1 = (mp_limb_t *) calloc(test_mpn_poly->limbs,sizeof(mp_limb_t));
          mp_limb_t * coeff2 = (mp_limb_t *) calloc(test_mpn_poly->limbs,sizeof(mp_limb_t));
          
          sign = _Zpoly_mpn_get_coeff(coeff1, test_mpn_poly, randint(test_mpn_poly->length));
          rand_coeff = randint(test_mpn_poly->length);
          _Zpoly_mpn_set_coeff(test_mpn_poly, rand_coeff, coeff1, sign, test_mpn_poly->limbs);
          sign2 = _Zpoly_mpn_get_coeff(coeff2, test_mpn_poly, rand_coeff);
          
          for (unsigned long i = 0; (i < test_mpn_poly->limbs) && (result == 1); i++)
          {
              result = (coeff1[i] == coeff2[i]);
          }
          
          free(coeff1);
          free(coeff2);
                    
          if (sign != sign2) result = 0;
      }
      Zpoly_mpn_clear(test_mpn_poly);
   }
   
   return result; 
}

int test_Zpoly_mpn_setequal()
{
   Zpoly_t test_poly;
   Zpoly_mpn_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long altered_coeff, extra_zeroes;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      Zpoly_mpn_init(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          Zpoly_mpn_realloc(test_mpn_poly, length+extra_zeroes);
          Zpoly_mpn_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          for (unsigned long i = 0; i < extra_zeroes; i++)
          {
             _Zpoly_mpn_set_coeff_ui(test_mpn_poly, length+i, 0);
          }
          test_mpn_poly->length = length;
          _Zpoly_mpn_set(test_mpn_poly2, test_mpn_poly);
          test_mpn_poly->length = length+extra_zeroes;
          result = _Zpoly_mpn_equal(test_mpn_poly2, test_mpn_poly); 
      }
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          Zpoly_mpn_realloc(test_mpn_poly, length+extra_zeroes);
          Zpoly_mpn_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          for (unsigned long i = 0; i < extra_zeroes; i++)
          {
             _Zpoly_mpn_set_coeff_ui(test_mpn_poly, length+i, 0);
          }
          test_mpn_poly->length = length;
          _Zpoly_mpn_set(test_mpn_poly2, test_mpn_poly);
          test_mpn_poly->length = length+extra_zeroes;
          altered_coeff = randint(length);
          test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)+1]++;
          if (test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] == 0)
             test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] = 1;
          result = !_Zpoly_mpn_equal(test_mpn_poly2, test_mpn_poly); 
      }
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          Zpoly_mpn_realloc(test_mpn_poly, length+extra_zeroes);
          Zpoly_mpn_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          for (unsigned long i = 0; i < extra_zeroes; i++)
          {
             _Zpoly_mpn_set_coeff_ui(test_mpn_poly, length+i, 0);
          }
          test_mpn_poly->length = length;
          _Zpoly_mpn_set(test_mpn_poly2, test_mpn_poly);
          test_mpn_poly->length = length+extra_zeroes;
          altered_coeff = randint(length);
          test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)]*=-1L;
          if (test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] == 0)
             test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] = 1;
          
          result = !_Zpoly_mpn_equal(test_mpn_poly2, test_mpn_poly); 
      }
      
      Zpoly_mpn_clear(test_mpn_poly);
      Zpoly_mpn_clear(test_mpn_poly2);         
   }
   
   return result; 
}



void Zpoly_mpn_test_all()
{
   int success, all_success = 1;

   RUN_TEST(Zpoly_mpn_convert);
   RUN_TEST(Zpoly_mpn_getset_ui);
   RUN_TEST(Zpoly_mpn_getset_si);
   RUN_TEST(Zpoly_mpn_get_coeff_ptr);
   RUN_TEST(Zpoly_mpn_normalise);
   RUN_TEST(Zpoly_mpn_getset_coeff);
   RUN_TEST(Zpoly_mpn_setequal);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   gmp_randinit_default(Zpoly_test_randstate);
   Zpoly_mpn_test_all();

   return 0;
}


