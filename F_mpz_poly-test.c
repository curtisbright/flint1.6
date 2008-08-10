/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

F_mpz_poly-test.c: Test code for F_mpz_poly.c and F_mpz_poly.h

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include "flint.h"
#include "mpz_poly.h"
#include "F_mpz_poly.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "test-support.h"
#include "zmod_poly.h"

#define VARY_BITS 1
#define SIGNS 1
#define SPARSE 1 
#define ITER 1

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

void mpz_randpoly(mpz_poly_t pol, long length, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_ensure_alloc(pol, length);
	mpz_poly_zero(pol);
   
   for (long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = z_randint(maxbits+1);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
#if SPARSE
          if (z_randint(10) == 1) mpz_rrandomb(temp, randstate, bits);
          else mpz_set_ui(temp, 0);
#else
          mpz_rrandomb(temp, randstate, bits);
#endif
#if SIGNS
          if (z_randint(2)) mpz_neg(temp,temp);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

void F_mpz_randpoly(F_mpz_poly_t poly, ulong length, ulong bits)
{
	mpz_poly_t m_poly;
	mpz_poly_init(m_poly);
	mpz_randpoly(m_poly, length, bits);
	mpz_poly_to_F_mpz_poly(poly, m_poly);
	mpz_poly_clear(m_poly);
}

int test_F_mpz_poly_convert()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_add()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_add(res, F_poly1, F_poly2);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_add(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_add(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_add(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_add(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_add(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_sub()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_sub(res, F_poly1, F_poly2);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_sub(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_sub(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_sub(res1, m_poly2, m_poly1);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_sub(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_sub(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_getset_si()
{
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   long coeff, coeff2;
   ulong coeff_bits, coeff_num;
   
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      for (ulong count2 = 0; (count2 < 500) && result == 1; count2++)
      {
          coeff_bits = z_randint(FLINT_BITS);
          coeff = z_randbits(coeff_bits);
          coeff_num = z_randint(F_poly->length);
              
			 if (z_randint(2)) coeff = -coeff;
              
			 if (F_poly->length)
          {
              F_mpz_poly_set_coeff_si(F_poly, coeff_num, coeff);
              coeff2 = F_mpz_poly_get_coeff_si(F_poly, coeff_num);
				  
				  result = (coeff2 == coeff);
				  if (!result)
				  {
					  printf("Error: length = %ld, coeff_num = %ld, coeff = %ld, coeff2 = %ld\n", length, coeff_num, coeff, coeff2);
				  }
          }
      }

      F_mpz_poly_clear(F_poly);
   }
   
   return result; 
}

int test_F_mpz_poly_getset_ui()
{
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   ulong coeff, coeff2;
   ulong coeff_bits, coeff_num;
   
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      for (ulong count2 = 0; (count2 < 500) && result == 1; count2++)
      {
          coeff_bits = z_randint(FLINT_BITS+1);
          coeff = z_randbits(coeff_bits);
          coeff_num = z_randint(F_poly->length);
                  
			 if (F_poly->length)
          {
              F_mpz_poly_set_coeff_ui(F_poly, coeff_num, coeff);
              coeff2 = F_mpz_poly_get_coeff_ui(F_poly, coeff_num);
				  
				  result = (coeff2 == coeff);
				  if (!result)
				  {
					  printf("Error: length = %ld, coeff_num = %ld, coeff = %ld, coeff2 = %ld\n", length, coeff_num, coeff, coeff2);
				  }
          }
      }

      F_mpz_poly_clear(F_poly);
   }
   
   return result; 
}

int test_F_mpz_poly_getset_mpz()
{
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   mpz_t coeff, coeff2;
   ulong coeff_bits, coeff_num;
	mpz_init(coeff);
	mpz_init(coeff2);
   
   for (ulong count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      for (ulong count2 = 0; (count2 < 300) && result == 1; count2++)
      {
          coeff_bits = z_randint(200);
          mpz_rrandomb(coeff, randstate, coeff_bits);
          coeff_num = z_randint(F_poly->length);
              
			 if (z_randint(2)) mpz_neg(coeff, coeff);
              
			 if (F_poly->length)
          {
              F_mpz_poly_set_coeff_mpz(F_poly, coeff_num, coeff);
              F_mpz_poly_get_coeff_mpz(coeff2, F_poly, coeff_num);
				  
				  result = (mpz_cmp(coeff2, coeff) == 0);
				  if (!result)
				  {
					  gmp_printf("Error: length = %ld, coeff_num = %ld, coeff = %Zd, coeff2 = %Zd\n", length, coeff_num, coeff, coeff2);
				  }
          }
      }

      F_mpz_poly_clear(F_poly);
   }
   
   mpz_clear(coeff);
	mpz_clear(coeff2);
	return result; 
}

int test_F_mpz_poly_set()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      F_mpz_poly_set(F_poly2, F_poly1);
		F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
      F_mpz_poly_clear(F_poly2);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

void F_mpz_poly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(F_mpz_poly_convert); 
   RUN_TEST(F_mpz_poly_add); 
   RUN_TEST(F_mpz_poly_sub); 
   RUN_TEST(F_mpz_poly_getset_si); 
   RUN_TEST(F_mpz_poly_getset_ui); 
   RUN_TEST(F_mpz_poly_getset_mpz); 
   RUN_TEST(F_mpz_poly_set); 
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}


