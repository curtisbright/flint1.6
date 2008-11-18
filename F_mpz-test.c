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

F_mpz-test.c: Test code for F_mpz.c and F_mpz.h

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include "flint.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "test-support.h"
#include "F_mpz.h"

#define VARY_BITS 1 // random coefficients have random number of bits up to the limit given
#define SIGNS 1 // random coefficients will be randomly signed
#define ITER 1 // if you want all tests to run longer, increase this

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // allows easy switching of debugging code on and off when debugging (if inserted)
#define DEBUG2 1 

/*
   Generate a random F_mpz_t with the given number of bits.
	If SIGNS is 1 then the value will be randomly signed, else it will be positive.
	Warning : do not use this function to test F_mpz_set_mpz!!
*/

void F_mpz_random(F_mpz_t f, ulong bits)
{
	if (bits = 0)
	{
		F_mpz_zero(f);
      return;
	}
	
	mpz_t temp;
	mpz_init(temp);
	
	mpz_rrandomb(temp, randstate, bits);
#if SIGNS
	if (z_randint(2)) mpz_neg(temp, temp);
#endif
   
	F_mpz_set_mpz(f, temp);

   mpz_clear(temp);
}

int test_F_mpz_getset_si()
{
   F_mpz_t f;
   int result = 1;
   ulong bits, val_bits;
	long val, val2;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      
		for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_random(f, bits); 

		   // Generate a random long
         val_bits = z_randint(FLINT_BITS - 1);
         val = z_randbits(val_bits);
         if (z_randint(2)) val = -val;
              
	      for (ulong i = 0; i < 100; i++) 
			{
				F_mpz_set_si(f, val);
            val2 = F_mpz_get_si(f);
			}

		   result = (val2 == val);
		   if (!result)
	      {
			   printf("Error: val = %ld, val2 = %ld\n", val, val2);
		   }
		}

      F_mpz_clear(f);
   }
   
   return result; 
}

int test_F_mpz_getset_ui()
{
   F_mpz_t f;
   int result = 1;
   ulong bits, val_bits, val, val2;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_random(f, bits); 

		   // Generate a random long
         val_bits = z_randint(FLINT_BITS);
         val = z_randbits(val_bits);
              
	      for (ulong i = 0; i < 100; i++) 
			{
				F_mpz_set_ui(f, val);
            val2 = F_mpz_get_ui(f);
			}

		   result = (val2 == val);
		   if (!result)
	      {
			   printf("Error: val = %u, val2 = %u\n", val, val2);
		   }
		}

      F_mpz_clear(f);
   }
   
   return result; 
}

int test_F_mpz_getset_mpz()
{
   F_mpz_t f;
   int result = 1;
   ulong bits;
   mpz_t val, val2;
   ulong val_bits;
	mpz_init(val);
	mpz_init(val2);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_init2(f, z_randint(10));
      
      F_mpz_random(f, bits); 
		    
      // set random coeffs in the poly
		for (ulong count2 = 0; (count2 < 100) && result == 1; count2++)
      {
         val_bits = z_randint(200);
         mpz_rrandomb(val, randstate, val_bits);
              
			if (z_randint(2)) mpz_neg(val, val);
              
			for (ulong i = 0; i < 100; i++) 
			{
				F_mpz_set_mpz(f, val);
            F_mpz_get_mpz(val2, f);
			}
				  
			result = (mpz_cmp(val2, val) == 0);
			if (!result)
			{
			   gmp_printf("Error: val = %Zd, val2 = %Zd\n", val, val2);
			}
      }

      F_mpz_clear(f);
   }
   
   mpz_clear(val);
	mpz_clear(val2);
	return result; 
}

int test_F_mpz_set()
{
   mpz_t m1, m2;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits, length;
   
   mpz_init(m1); 
   mpz_init(m2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits = z_randint(200) + 1;
      F_mpz_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      for (ulong i = 0; i < 100; i++) F_mpz_set(f2, f1);
		F_mpz_get_mpz(m2, f2);
          
      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
	// check aliasing
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      
      bits = z_randint(200) + 1;
      F_mpz_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      for (ulong i = 0; i < 100; i++) F_mpz_set(f1, f1);
		F_mpz_get_mpz(m2, f1);
          
      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   
   return result;
}

/*int test_F_mpz_add()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
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
   
	// test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
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
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
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

int test_F_mpz_sub()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
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
   
   // test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
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
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
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

int test_F_mpz_bits()
{
   mpz_poly_t m_poly1;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length, next_bits;
	long sign, mpz_bits, test_bits;
   
   mpz_poly_init(m_poly1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_bits = 0;
      sign = 1L;
      for (ulong i = 0; i < m_poly1->length; i++)
      {
         next_bits = mpz_sizeinbase(m_poly1->coeffs[i], 2);
         if (next_bits > mpz_bits) mpz_bits = next_bits;
         if (mpz_sgn(m_poly1->coeffs[i]) < 0L) sign = -1L;
      }
      mpz_bits = sign*mpz_bits;

		test_bits = F_mpz_poly_max_bits(F_poly);
          
      result = (mpz_bits == test_bits); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, m_poly1->length = %ld\n", length, bits, m_poly1->length);
         printf("mpz_bits = %ld, test_bits = %ld\n", mpz_bits, test_bits);
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   
   return result;
}

int test_F_mpz_size()
{
   mpz_poly_t m_poly1;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length, next_limbs;
	long sign, mpz_limbs, test_limbs;
   
   mpz_poly_init(m_poly1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(500) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_limbs = 0;
      for (ulong i = 0; i < m_poly1->length; i++)
      {
         next_limbs = mpz_size(m_poly1->coeffs[i]);
         if (next_limbs > mpz_limbs) mpz_limbs = next_limbs;
      }
      
		test_limbs = F_mpz_poly_max_limbs(F_poly);
          
      result = (mpz_limbs == test_limbs); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, m_poly1->length = %ld\n", length, bits, m_poly1->length);
         printf("mpz_limbs = %ld, test_limbs = %ld\n", mpz_limbs, test_limbs);
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   
   return result;
}

int test_F_mpz_equal()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   
   // random polys unlikely to be equal, test against mpz_poly_equal
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      
		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
          
      result = (mpz_poly_equal(m_poly1, m_poly2) == F_mpz_poly_equal(F_poly1, F_poly2)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
   }

	// polys are equal
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      F_mpz_poly_set(F_poly2, F_poly1);
          
      result = (F_mpz_poly_equal(F_poly1, F_poly2)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
   }

	// polys are aliased
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
          
      result = (F_mpz_poly_equal(F_poly1, F_poly1)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
   }

	mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_neg()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, F_poly4;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // negate and negate back
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      F_mpz_randpoly(F_poly1, length1, bits1);
      
		F_mpz_poly_neg(F_poly2, F_poly1);
      F_mpz_poly_neg(F_poly2, F_poly2);
          
      result = (F_mpz_poly_equal(F_poly1, F_poly2)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
   }

	// sub equals negate and add, included aliased negation
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(F_poly4);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      bits2 = z_randint(200) + 1;
      length2 = z_randint(100);
      F_mpz_randpoly(F_poly1, length1, bits1);
      F_mpz_randpoly(F_poly2, length2, bits2);
      
		F_mpz_poly_sub(F_poly3, F_poly1, F_poly2);
      F_mpz_poly_neg(F_poly2, F_poly2);
      F_mpz_poly_add(F_poly4, F_poly1, F_poly2);
          
      result = (F_mpz_poly_equal(F_poly3, F_poly4)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(F_poly4);
   }

   return result;
}

int test_F_mpz_swap()
{
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   unsigned long bits, length, length2;
   
   // reverse and reverse back, second time aliased
	for (unsigned long count1 = 0; (count1 < 20000) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);   
      
      length = z_randint(100);
      length2 = length + z_randint(200);
      
      F_mpz_randpoly(F_poly, length, bits); 
                
      F_mpz_poly_reverse(F_poly2, F_poly, length2);
		F_mpz_poly_reverse(F_poly2, F_poly2, length2);
           
      result = F_mpz_poly_equal(F_poly2, F_poly);
      if (!result)
		{
			printf("Error: length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
			printf("F_poly->length = %ld, F_poly2->length = %ld\n", F_poly->length, F_poly2->length);
		}

      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }
   
   // completely aliased reverse and reverse back
	for (unsigned long count1 = 0; (count1 < 20000) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);   
      
      length = z_randint(100);
      length2 = length + z_randint(200);
       
      F_mpz_randpoly(F_poly, length, bits); 
         
      F_mpz_poly_set(F_poly2, F_poly);
      F_mpz_poly_reverse(F_poly, F_poly, length2);
      F_mpz_poly_reverse(F_poly, F_poly, length2);
           
      result = F_mpz_poly_equal(F_poly2, F_poly);
		if (!result)
		{
			printf("Error: length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
			printf("F_poly->length = %ld, F_poly2->length = %ld\n", F_poly->length, F_poly2->length);
		}
      
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }

   return result; 
}

int test_F_mpz_mul_2exp()
{
   F_mpz_poly_t F_poly, F_poly2, F_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long shift;
   
   F_mpz_poly_init(F_poly);
   F_mpz_poly_init(F_poly2);
   F_mpz_poly_init(F_poly3);
      
	// left shift followed by right shift
	for (unsigned long count1 = 0; (count1 < 700) && (result == 1); count1++)
   {
      bits = z_randint(500)+ 1;
      length = z_randint(500);        
      shift = z_randint(100);
		    
		F_mpz_randpoly(F_poly, length, bits); 
      F_mpz_poly_set(F_poly3, F_poly);
      
		F_mpz_poly_left_shift(F_poly2, F_poly, shift); 
      F_mpz_poly_right_shift(F_poly, F_poly2, shift);
      
      result = F_mpz_poly_equal(F_poly3, F_poly);
		if (!result) printf("Error: bits = %ld, length = %ld, shift = %ld\n", bits, length, shift);
   }

	// F_poly3 is used uninitialised after this point so clear it
	F_mpz_poly_clear(F_poly3);

   // left shift followed by right shift, completely aliased
	for (unsigned long count1 = 0; (count1 < 700) && (result == 1); count1++)
   {
      bits = z_randint(500)+ 1;
      length = z_randint(500);        
      shift = z_randint(100);
		    
		F_mpz_randpoly(F_poly, length, bits); 
      F_mpz_poly_set(F_poly2, F_poly);
      
		F_mpz_poly_left_shift(F_poly, F_poly, shift); 
      F_mpz_poly_right_shift(F_poly, F_poly, shift);
      
      result = F_mpz_poly_equal(F_poly2, F_poly);
		if (!result) printf("Error: bits = %ld, length = %ld, shift = %ld\n", bits, length, shift);
   }
   
   // explicit check of right shift
	for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
   { 
       bits = z_randint(500)+ 1;
       length = z_randint(500);        
       if (length) shift = z_randint(length);
		 else shift = 0;
		
       do F_mpz_randpoly(F_poly, length, bits); 
       while (F_poly->length < length);

       F_poly3->length = F_poly->length - shift;
       F_poly3->coeffs = F_poly->coeffs + shift;
		 F_poly3->mpz_coeffs = F_poly->mpz_coeffs;
		 
		 F_mpz_poly_right_shift(F_poly2, F_poly, shift);      
          
		 result = F_mpz_poly_equal(F_poly3, F_poly2);
		 if (!result) printf("Error: bits = %ld, length = %ld, shift = %ld\n", bits, length, shift);
   }

   F_mpz_poly_clear(F_poly);
   F_mpz_poly_clear(F_poly2);
  
   return result; 
}

int test_F_mpz_mul_ui()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, bits2, length;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      
      mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_poly_scalar_mul_ui(F_poly2, F_poly, mult);
          
      mpz_poly_init(m_poly2);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mult == 0L) result = (F_poly2->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased multiply
	for (ulong count1 = 0; (count1 < 20000) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);        

      F_mpz_poly_init(F_poly);
      
		mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
          
		F_mpz_poly_scalar_mul_ui(F_poly, F_poly, mult);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly); 
          
      if (mult == 0L) result = (F_poly->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   mpz_clear(temp);
   
   return result; 
}

int test_F_mpz_mul_si()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, bits2, length;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      
      mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_poly_scalar_mul_si(F_poly2, F_poly, mult);
          
      mpz_poly_init(m_poly2);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mult == 0L) result = (F_poly2->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            if (mult < 0L) 
			   {
				   mpz_mul_ui(temp, m_poly->coeffs[i], -mult);
				   mpz_neg(temp, temp);
			   } else mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased multiply
	for (ulong count1 = 0; (count1 < 20000) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);        

      F_mpz_poly_init(F_poly);
      
		mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
          
		F_mpz_poly_scalar_mul_si(F_poly, F_poly, mult);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly); 
          
      if (mult == 0L) result = (F_poly->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            if (mult < 0L) 
			   {
				   mpz_mul_ui(temp, m_poly->coeffs[i], -mult);
				   mpz_neg(temp, temp);
			   } else mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }          
	   }
		
	   if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   mpz_clear(temp);
   
   return result; 
}

int test_F_mpz_mul_mpz()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, bits2, length;
   mpz_t temp, mult;
   mpz_init(temp);
   mpz_init(mult);
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 10000) && (result == 1); count1++)
   {
      bits = z_randint(300)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      
      mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(200);
		mpz_rrandomb(mult, randstate, bits2);
		if (z_randint(2)) mpz_neg(mult, mult);
      
		F_mpz_poly_scalar_mul_mpz(F_poly2, F_poly, mult);
          
      mpz_poly_init(m_poly2);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mpz_sgn(mult) == 0) result = (F_poly2->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			gmp_printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %Zd\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased multiply
	for (ulong count1 = 0; (count1 < 10000) && (result == 1); count1++)
   {
      bits = z_randint(300)+ 1;
      length = z_randint(200);        

      F_mpz_poly_init(F_poly);
      
		mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(200);
		mpz_rrandomb(mult, randstate, bits2);
		if (z_randint(2)) mpz_neg(mult, mult);
          
		F_mpz_poly_scalar_mul_mpz(F_poly, F_poly, mult);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly); 
          
      if (mpz_sgn(mult) == 0) result = (F_poly->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

	   if (!result) 
		{
			gmp_printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %Zd\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   mpz_clear(temp);
   mpz_clear(mult);
   
   return result; 
}

int test_F_mpz_addmul_ui()
{
   mpz_poly_t m_poly1, m_poly2, m_poly3, m_poly4;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, bits2, bits3, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(m_poly3); 
   mpz_poly_init(m_poly4); 

   for (ulong count1 = 0; (count1 < 4000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly);
		F_mpz_poly_init(F_poly2);
		
		length = z_randint(100) + 1;

      bits = z_randint(200) + 1;
      do mpz_randpoly(m_poly1, length, bits);
		while (m_poly1->length != length);
      
		bits2 = z_randint(200) + 1;
      do mpz_randpoly_dense(m_poly2, length, bits2);
		while (m_poly2->length != length);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);

		bits3 = z_randint(FLINT_BITS+1);
		ulong x = z_randbits(bits3);

		for (ulong i = 0; i < length - 1; i++) // don't touch last coefficient for normalisation
		{
			_F_mpz_addmul_ui(F_poly, i, F_poly2, i, x);
			mpz_addmul_ui(m_poly1->coeffs[i], m_poly2->coeffs[i], x);
		}

      F_mpz_poly_to_mpz_poly(m_poly2, F_poly);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, bits3 = %ld\n", length, bits, bits2, bits3);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }
   
   // sparse polynomial and aliasing
	for (ulong count1 = 0; (count1 < 4000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly);
		
		length = z_randint(100) + 1;

      bits = z_randint(200) + 1;
      do mpz_randpoly(m_poly1, length, bits);
		while (m_poly1->length != length);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		bits3 = z_randint(FLINT_BITS+1);
		ulong x = z_randbits(bits3);

		for (ulong i = 0; i < length - 1; i++) // don't touch last coefficient for normalisation
		{
			_F_mpz_addmul_ui(F_poly, i, F_poly, i, x);
			mpz_addmul_ui(m_poly1->coeffs[i], m_poly1->coeffs[i], x);
		}

      F_mpz_poly_to_mpz_poly(m_poly2, F_poly);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits3 = %ld\n", length, bits, bits3);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   mpz_poly_clear(m_poly3);
   mpz_poly_clear(m_poly4);
   
   return result;
}

int test_F_mpz_submul_ui()
{
   mpz_poly_t m_poly1, m_poly2, m_poly3, m_poly4;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, bits2, bits3, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(m_poly3); 
   mpz_poly_init(m_poly4); 

   for (ulong count1 = 0; (count1 < 4000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly);
		F_mpz_poly_init(F_poly2);
		
		length = z_randint(100) + 1;

      bits = z_randint(200) + 1;
      do mpz_randpoly(m_poly1, length, bits);
		while (m_poly1->length != length);
      
		bits2 = z_randint(200) + 1;
      do mpz_randpoly_dense(m_poly2, length, bits2);
		while (m_poly2->length != length);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);

		bits3 = z_randint(FLINT_BITS+1);
		ulong x = z_randbits(bits3);

		for (ulong i = 0; i < length - 1; i++) // don't touch last coefficient for normalisation
		{
			_F_mpz_submul_ui(F_poly, i, F_poly2, i, x);
			mpz_submul_ui(m_poly1->coeffs[i], m_poly2->coeffs[i], x);
		}

      F_mpz_poly_to_mpz_poly(m_poly2, F_poly);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, bits3 = %ld\n", length, bits, bits2, bits3);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }
   
   // sparse polynomial and aliasing
	for (ulong count1 = 0; (count1 < 4000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly);
		
		length = z_randint(100) + 1;

      bits = z_randint(200) + 1;
      do mpz_randpoly(m_poly1, length, bits);
		while (m_poly1->length != length);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		bits3 = z_randint(FLINT_BITS+1);
		ulong x = z_randbits(bits3);

		for (ulong i = 0; i < length - 1; i++) // don't touch last coefficient for normalisation
		{
			_F_mpz_submul_ui(F_poly, i, F_poly, i, x);
			mpz_submul_ui(m_poly1->coeffs[i], m_poly1->coeffs[i], x);
		}

      F_mpz_poly_to_mpz_poly(m_poly2, F_poly);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits3 = %ld\n", length, bits, bits3);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   mpz_poly_clear(m_poly3);
   mpz_poly_clear(m_poly4);
   
   return result;
}*/

void F_mpz_poly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(F_mpz_getset_ui); 
   RUN_TEST(F_mpz_getset_si); 
   RUN_TEST(F_mpz_getset_mpz); 
   RUN_TEST(F_mpz_set); 
   /*RUN_TEST(F_mpz_bits); 
   RUN_TEST(F_mpz_size); 
   RUN_TEST(F_mpz_equal); 
   RUN_TEST(F_mpz_swap); 
   RUN_TEST(F_mpz_neg); 
   RUN_TEST(F_mpz_add); 
   RUN_TEST(F_mpz_sub); 
   RUN_TEST(F_mpz_mul_2exp); 
   RUN_TEST(F_mpz_div_2exp); 
   RUN_TEST(F_mpz_mul_ui); 
   RUN_TEST(F_mpz_mul_si); 
   RUN_TEST(F_mpz_mul_mpz); 
   RUN_TEST(F_mpz_addmul_ui); 
   RUN_TEST(F_mpz_submul_ui);*/ 
	
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


