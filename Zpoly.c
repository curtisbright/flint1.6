/****************************************************************************

Zpoly.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <string.h>
#include "flint.h"
#include "Zpoly.h"
#include "flint-manager.h"

/****************************************************************************

   _Zpoly_* layer

****************************************************************************/

/* Normalise polynomial so that leading coefficient is nonzero */

void _Zpoly_normalise(Zpoly_t poly)
{
   while (poly->length && !mpz_sgn(poly->coeffs[poly->length-1]))
      poly->length--;
}

/* Set output poly to input poly */

void _Zpoly_set(Zpoly_t output, Zpoly_t input)
{
   FLINT_ASSERT(output->alloc >= input->length);
    
   for (unsigned long i = 0; i < input->length; i++)
      mpz_set(output->coeffs[i], input->coeffs[i]);
   
   output->length = input->length;
}

// todo: probably the "raw" version of this should assume something like:
// input1 is at least as long as input2. The hand-holding one would then
// call it with the parameters in the correct order.

/* Return 1 if polynomials are equal, 0 otherwise. 
   Polynomials do not need to be normalised. */

int _Zpoly_equal(Zpoly_t input1, Zpoly_t input2)
{
   int is_input1_longer = (input1->length > input2->length);

   unsigned long short_length =
           is_input1_longer ? input2->length : input1->length;

   unsigned long i;
   for (i = 0; i < short_length; i++)
      if (mpz_cmp(input1->coeffs[i], input2->coeffs[i]))
         return 0;

   if (is_input1_longer)
   {
      for (; i < input1->length; i++)
         if (mpz_sgn(input1->coeffs[i]))
            return 0;
   }
   else
   {
      for (; i < input2->length; i++)
         if (mpz_sgn(input2->coeffs[i]))
            return 0;
   }

   return 1;
}

/* Set output poly to sum of input polys */

// todo: we need several versions of this:
// 1) the lowest layer should assume that input1 is at least as long as input2
// 2) a middle layer which doesn't assume anything about the lengths, but assumes
//    the output has enough room for the result
// 3) the hand-holding layer does memory reallocation and everything.

void _Zpoly_add(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2)
{
   FLINT_ASSERT(output->alloc >= input1->length);
   FLINT_ASSERT(output->alloc >= input2->length);

   int is_input1_longer = (input1->length > input2->length);

   unsigned long short_length =
           is_input1_longer ? input2->length : input1->length;
   
   unsigned long i;
   for (i = 0; i < short_length; i++)
      mpz_add(output->coeffs[i], input1->coeffs[i], input2->coeffs[i]);
      
   if (is_input1_longer)
      for (; i < input1->length; i++)
         mpz_set(output->coeffs[i], input1->coeffs[i]);
   else
      for (; i < input2->length; i++)
         mpz_set(output->coeffs[i], input2->coeffs[i]);

   output->length = i;
}

/* Set output poly to input1 minus input2 */

// todo: see remarks above for _Zpoly_add, but for sub we need even an
// extra version, since things behave differently when the parameter order is
// switched.

void _Zpoly_sub(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2)
{
   FLINT_ASSERT(output->alloc >= input1->length);
   FLINT_ASSERT(output->alloc >= input2->length);

   int is_input1_longer = (input1->length > input2->length);

   unsigned long short_length =
           is_input1_longer ? input2->length : input1->length;
   
   unsigned long i;
   for (i = 0; i < short_length; i++)
      mpz_sub(output->coeffs[i], input1->coeffs[i], input2->coeffs[i]);
      
   if (is_input1_longer)
      for (; i < input1->length; i++)
         mpz_set(output->coeffs[i], input1->coeffs[i]);
   else
      for (; i < input2->length; i++)
         mpz_neg(output->coeffs[i], input2->coeffs[i]);

   output->length = i;
}

/* Negate input poly */

void _Zpoly_negate(Zpoly_t output, Zpoly_t input)
{
   FLINT_ASSERT(output->alloc >= input->length);
   
   for (unsigned long i = 0; i < input->length; i++)
      mpz_neg(output->coeffs[i], input->coeffs[i]);
   
   output->length = input->length;
}


void _Zpoly_scalar_mul(Zpoly_t poly, mpz_t x)
{
   abort();
}

void _Zpoly_scalar_mul_ui(Zpoly_t poly, unsigned long x)
{
   abort();
}

void _Zpoly_scalar_mul_si(Zpoly_t poly, long x)
{
   abort();
}

void _Zpoly_scalar_div(Zpoly_t poly, mpz_t x)
{
   abort();
}

void _Zpoly_scalar_div_ui(Zpoly_t poly, unsigned long x)
{
   abort();
}

/* Multiply two input polynomials */

void _Zpoly_mul(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2)
{
   FLINT_ASSERT(output != input1);
   FLINT_ASSERT(output != input2);

   // naive multiplication for now....
   // todo: plug in actual multiplication code :-)
   _Zpoly_mul_naive(output, input1, input2);
}
                           
/* Naieve schoolboy polynomial multiplication routine */

void _Zpoly_mul_naive(Zpoly_t output, Zpoly_t input1,
                             Zpoly_t input2)
{
   FLINT_ASSERT(output != input1);
   FLINT_ASSERT(output != input2);

   if (!input1->length || !input2->length)
   {
      // one of the inputs is zero
      output->length = 0;
      return;
   }
   
   output->length = input1->length + input2->length - 1;
   FLINT_ASSERT(output->alloc >= output->length);

   for (unsigned long i = 0; i < output->length; i++)
      mpz_set_ui(output->coeffs[i], 0);
   
   for (unsigned long i = 0; i < input1->length; i++)
      for (unsigned long j = 0; j < input2->length; j++)
         mpz_addmul(output->coeffs[i+j], input1->coeffs[i], input2->coeffs[j]);
}


void _Zpoly_mul_karatsuba(Zpoly_t output, Zpoly_t input1,
                                 Zpoly_t input2)
{
   abort();
}

void _Zpoly_sqr(Zpoly_t output, Zpoly_t input)
{
   abort();
}

void _Zpoly_sqr_naive(Zpoly_t output, Zpoly_t input)
{
   abort();
}

void _Zpoly_sqr_karatsuba(Zpoly_t output, Zpoly_t input)
{
   abort();
}

void _Zpoly_left_shift(Zpoly_t output, Zpoly_t input,
                              unsigned long n)
{
   abort();
}

void _Zpoly_right_shift(Zpoly_t output, Zpoly_t input,
                               unsigned long n)
{
   abort();
}


void _Zpoly_div(Zpoly_t quotient, Zpoly_t input1,
                       Zpoly_t input2)
{
   abort();
}

void _Zpoly_rem(Zpoly_t remainder, Zpoly_t input1,
                       Zpoly_t input2)
{
   abort();
}

void _Zpoly_div_rem(Zpoly_t quotient, Zpoly_t remainder,
                           Zpoly_t input1, Zpoly_t input2)
{
   abort();
}

void _Zpoly_gcd(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2)
{
   abort();
}

void _Zpoly_xgcd(Zpoly_t a, Zpoly_t b, Zpoly_t output,
                        Zpoly_t input1, Zpoly_t input2)
{
   abort();
}

void _Zpoly_content(mpz_t content, Zpoly_t a)
{
   abort();
}

/****************************************************************************

   Zpoly_* layer  ("non-raw" functions)

****************************************************************************/

/* Create polynomial with one allocated coefficient and length zero */

void Zpoly_init(Zpoly_t poly)
{
   poly->coeffs = (mpz_t*) flint_malloc(sizeof(mpz_t));
   mpz_init(poly->coeffs[0]);
   poly->alloc = 1;
   poly->length = 0;
}

/* Create polynomial with "alloc" allocated coefficients and length zero */

void Zpoly_init2(Zpoly_t poly, unsigned long alloc)
{
   poly->coeffs = (mpz_t*) flint_malloc(sizeof(mpz_t) * alloc);
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init(poly->coeffs[i]);
   poly->alloc = alloc;
   poly->length = 0;
}

/* Create polynomial with "alloc" allocated coefficients each with allocated
   space for at least the given number of bits, and length zero */

void Zpoly_init3(Zpoly_t poly, unsigned long alloc,
                     unsigned long coeff_bits)
{
   poly->coeffs = (mpz_t*) flint_malloc(sizeof(mpz_t) * alloc);
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init2(poly->coeffs[i], coeff_bits);
   poly->alloc = alloc;
   poly->length = 0;
}

/* Shrink or expand a polynomial to "alloc" coefficients */

void Zpoly_realloc(Zpoly_t poly, unsigned long alloc)
{
   if (alloc < poly->alloc)
   {
      // shrink available space; need to mpz_clear mpz's in the tail
      for (unsigned long i = alloc; i < poly->alloc; i++)
         mpz_clear(poly->coeffs[i]);
   }
   
   poly->coeffs = (mpz_t*) flint_realloc(poly->coeffs, sizeof(mpz_t) * alloc);

   // create new mpz's if necessary
   for (unsigned long i = poly->alloc; i < alloc; i++)
      mpz_init(poly->coeffs[i]);

   poly->alloc = alloc;
   
   // truncate actual data if necessary
   if (poly->length > alloc)
      poly->length = alloc;
}

/* Expand poly to "alloc" allocated coefficients, ensuring that it is made at
   least twice the current allocated size */

void Zpoly_ensure_space_IMPL(Zpoly_t poly, unsigned long alloc)
{
   if (alloc < 2*poly->alloc)
      alloc = 2*poly->alloc;
   Zpoly_realloc(poly, alloc);
}

/* Free all the poly coefficients */

void Zpoly_clear(Zpoly_t poly)
{
   for (unsigned long i = 0; i < poly->alloc; i++)
      mpz_clear(poly->coeffs[i]);
   flint_free(poly->coeffs);
}

/* Return a pointer to the given coefficient, or NULL if the poly isn't
   that long */

mpz_t* Zpoly_get_coeff_ptr(Zpoly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return NULL;
   return &poly->coeffs[n];
}

/* Set output to the given polynomial coefficient, or to zero if the poly
   isn't that long */

void Zpoly_get_coeff(mpz_t output, Zpoly_t poly,
                         unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(output, 0);
   else
      mpz_set(output, poly->coeffs[n]);
}

/* Return the given coefficient as an unsigned long (or zero) */

unsigned long Zpoly_get_coeff_ui(Zpoly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return 0;
   return mpz_get_ui(poly->coeffs[n]);
}

/* Return the given coefficient as a signed int (or zero) */

long Zpoly_get_coeff_si(Zpoly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return 0;
   return mpz_get_si(poly->coeffs[n]);
}

/* Make poly longer and set the coefficients to the given string of coeffs */

int Zpoly_set_from_string(Zpoly_t output, char* s)
{
   const char* whitespace = " \t\n\r";
   
   output->length = 0;
   
   while (1)
   {
      // skip whitespace
      s += strspn(s, whitespace);
      if (*s == 0)
         return 1;
         
      // ensure sufficient space, and grab coefficient
      Zpoly_ensure_space(output, output->length + 1);
      if (!gmp_sscanf(s, "%Zd", output->coeffs[output->length]))
         return 0;
      output->length++;
      
      // jump to next whitespace
      s += strcspn(s, whitespace);
   }
}

/* Return the length of a string to output a polynomial in base 10 */

unsigned long Zpoly_get_string_size(Zpoly_t poly)
{
   unsigned long size = 0;
   for (unsigned long i = 0; i < poly->length; i++)
      // (+2 is for the sign and a space)
      size += mpz_sizeinbase(poly->coeffs[i], 10) + 2;
   // (+1 is for the null terminator)
   return size + 1;
}

/* Set output to polynomial as a null terminated ascii string */

void Zpoly_get_as_string(char* output, Zpoly_t poly)
{
   if (!poly->length)
   {
      output[0] = 0;
      return;
   }
   
   for (unsigned long i = 0; i < poly->length; i++)
   {
      mpz_get_str(output, 10, poly->coeffs[i]);
      output += strlen(output);
      *output = ' ';
      output++;
   }
   
   output--;
   *output = 0;
}

/* Print polynomial as string to the given stream */

void Zpoly_print(FILE* output, Zpoly_t poly)
{
   unsigned long size = Zpoly_get_string_size(poly);
   char* buf = flint_malloc(size);
   Zpoly_get_as_string(buf, poly);
   fprintf(output, buf);
   flint_free(buf);
}

/* Set the given polynomial coefficient to the given value */

void Zpoly_set_coeff(Zpoly_t poly, unsigned long n, mpz_t x)
{
   Zpoly_ensure_space(poly, n+1);
   _Zpoly_set_coeff(poly, n, x);
   if (poly->length <= n)
      poly->length = n+1;
}

/* Set the given polynomial coefficient to the given unsigned long */

void Zpoly_set_coeff_ui(Zpoly_t poly, unsigned long n, unsigned long x)
{
   Zpoly_ensure_space(poly, n+1);
   _Zpoly_set_coeff_ui(poly, n, x);
   if (poly->length <= n)
      poly->length = n+1;
}

/* Set the given polynomial coefficient to the given signed int */

void Zpoly_set_coeff_si(Zpoly_t poly, unsigned long n, long x)
{
   Zpoly_ensure_space(poly, n+1);
   _Zpoly_set_coeff_si(poly, n, x);
   if (poly->length <= n)
      poly->length = n+1;
}

/* Set the output poly equal to the input poly */

void Zpoly_set(Zpoly_t output, Zpoly_t input)
{
   Zpoly_ensure_space(output, input->length);
   _Zpoly_set(output, input);
}

/* Set the output poly to the sum of the input polys */

void Zpoly_add(Zpoly_t output, Zpoly_t input1, Zpoly_t input2)
{
   Zpoly_ensure_space(output, FLINT_MAX(input1->length, input2->length));
   _Zpoly_add(output, input1, input2);
}

/* Set the output poly to input1 minus input2 */

void Zpoly_sub(Zpoly_t output, Zpoly_t input1, Zpoly_t input2)
{
   Zpoly_ensure_space(output, FLINT_MAX(input1->length, input2->length));
   _Zpoly_sub(output, input1, input2);
}

/* Negate the polynomial */

void Zpoly_negate(Zpoly_t output, Zpoly_t input)
{
   Zpoly_ensure_space(output, input->length);
   _Zpoly_negate(output, input);
}

void Zpoly_scalar_mul(Zpoly_t poly, mpz_t x)
{
   abort();
}

void Zpoly_scalar_mul_ui(Zpoly_t poly, unsigned long x)
{
   abort();
}

void Zpoly_scalar_mul_si(Zpoly_t poly, long x)
{
   abort();
}

void Zpoly_scalar_div(Zpoly_t poly, mpz_t x)
{
   abort();
}

void Zpoly_scalar_div_ui(Zpoly_t poly, unsigned long x)
{
   abort();
}

/* Set output poly to the product of the input polys */

void Zpoly_mul(Zpoly_t output, Zpoly_t input1, Zpoly_t input2)
{
   if (!input1->length || !input2->length)
   {
      // one of the inputs is zero
      output->length = 0;
      return;
   }
   
   unsigned long output_length = input1->length + input2->length - 1;

   if (output == input1 || output == input2)
   {
      // if output is inplace, need a temporary
      Zpoly_t temp;
      Zpoly_init2(temp, output_length);
      _Zpoly_mul(temp, input1, input2);
      _Zpoly_swap(output, temp);
      Zpoly_clear(temp);
   }
   else
   {
      // not inplace; just call directly
      Zpoly_ensure_space(output, output_length);
      _Zpoly_mul(output, input1, input2);
   }
}

/* Naieve schoolboy polynomial multiplication routine */

void Zpoly_mul_naive(Zpoly_t output, Zpoly_t input1,
                         Zpoly_t input2)
{
   if (!input1->length || !input2->length)
   {
      // one of the inputs is zero
      output->length = 0;
      return;
   }
   
   unsigned long output_length = input1->length + input2->length - 1;

   if (output == input1 || output == input2)
   {
      // if output is inplace, need a temporary
      Zpoly_t temp;
      Zpoly_init2(temp, output_length);
      _Zpoly_mul_naive(temp, input1, input2);
      _Zpoly_swap(output, temp);
      Zpoly_clear(temp);
   }
   else
   {
      // not inplace; just call directly
      Zpoly_ensure_space(output, output_length);
      _Zpoly_mul_naive(output, input1, input2);
   }
}


// ----------------------------------------------------------------------------
// A few support functions for naive KS multiplication.

/*
Sets y = \sum_{i=0}^{len-1} x[i] * 2^(ki)
Running time should be O(k*len*log(len))
*/

void Zpoly_mul_naive_KS_pack(mpz_t y, mpz_t* x, unsigned long len,
                                 unsigned long k)
{
   if (len == 1)
      mpz_set(y, x[0]);
   else
   {
      mpz_t temp;
      mpz_init(temp);
      unsigned long half = len/2;
      Zpoly_mul_naive_KS_pack(temp, x, half, k);
      Zpoly_mul_naive_KS_pack(y, x + half, len - half, k);
      mpz_mul_2exp(y, y, half*k);
      mpz_add(y, y, temp);
      mpz_clear(temp);
   }
}


/*
Inverse operation of Zpoly_mul_naive_KS_pack
(note: y is destroyed)
*/

void Zpoly_mul_naive_KS_unpack(mpz_t* x, unsigned long len, mpz_t y,
                                   unsigned long k)
{
   if (len == 1)
      mpz_set(x[0], y);
   else
   {
      mpz_t temp;
      mpz_init(temp);
      unsigned long half = len/2;
      if (mpz_tstbit(y, k*half - 1))
      {
         mpz_cdiv_q_2exp(temp, y, half*k);
         mpz_cdiv_r_2exp(y, y, half*k);
      }
      else
      {
         mpz_fdiv_q_2exp(temp, y, half*k);
         mpz_fdiv_r_2exp(y, y, half*k);
      }
      Zpoly_mul_naive_KS_unpack(x, half, y, k);
      Zpoly_mul_naive_KS_unpack(x + half, len - half, temp, k);
      mpz_clear(temp);
   }
}


/*
Counts maximum number of bits in abs(x->coeffs[i])
*/

unsigned long Zpoly_mul_naive_KS_get_max_bits(Zpoly_t x)
{
   unsigned long bits = 0, temp, i;
   for (i = 0; i < x->length; i++)
   {
      temp = mpz_sizeinbase(x->coeffs[i], 2);
      if (temp > bits)
         bits = temp;
   }
   return bits;
}

/* A simple Kronecker segmentation multiplication routine */

void Zpoly_mul_naive_KS(Zpoly_t output, Zpoly_t input1,
                            Zpoly_t input2)
{
   if (!input1->length || !input2->length)
   {
      // one of the inputs is zero
      output->length = 0;
      return;
   }
   
   mpz_t z1;
   mpz_t z2;
   mpz_init(z1);
   mpz_init(z2);

   unsigned long output_len = input1->length + input2->length - 1;
   unsigned long bits1 = Zpoly_mul_naive_KS_get_max_bits(input1);
   unsigned long bits2 = Zpoly_mul_naive_KS_get_max_bits(input2);
   unsigned long bits = bits1 + bits2 + 2 + ceil_log2(output_len);

   Zpoly_mul_naive_KS_pack(z1, input1->coeffs, input1->length, bits);
   Zpoly_mul_naive_KS_pack(z2, input2->coeffs, input2->length, bits);
   mpz_mul(z1, z1, z2);
   Zpoly_ensure_space(output, output_len);
   Zpoly_mul_naive_KS_unpack(output->coeffs, output_len, z1, bits);
   output->length = output_len;

   mpz_clear(z1);
   mpz_clear(z2);
}

/* A simple Kronecker substitution squaring routine */

void Zpoly_naive_KS_sqr(Zpoly_t output, Zpoly_t input)
{
   if (!input->length)
   {
      // input is zero
      output->length = 0;
      return;
   }
   
   mpz_t z;
   mpz_init(z);

   unsigned long output_len = 2*input->length - 1;
   unsigned long bits = 2 * Zpoly_mul_naive_KS_get_max_bits(input)
                          + 2 + ceil_log2(output_len);

   Zpoly_mul_naive_KS_pack(z, input->coeffs, input->length, bits);
   mpz_mul(z, z, z);
   Zpoly_ensure_space(output, output_len);
   Zpoly_mul_naive_KS_unpack(output->coeffs, output_len, z, bits);
   output->length = output_len;
   
   mpz_clear(z);
}


void Zpoly_mul_karatsuba(Zpoly_t output, Zpoly_t input1,
                             Zpoly_t input2)
{
   abort();
}

void Zpoly_sqr(Zpoly_t output, Zpoly_t input)
{
   abort();
}

void Zpoly_sqr_naive(Zpoly_t output, Zpoly_t input)
{
   abort();
}

void Zpoly_sqr_karatsuba(Zpoly_t output, Zpoly_t input)
{
   abort();
}

void Zpoly_left_shift(Zpoly_t output, Zpoly_t input,
                          unsigned long n)
{
   abort();
}

void Zpoly_right_shift(Zpoly_t output, Zpoly_t input,
                           unsigned long n)
{
   abort();
}

void Zpoly_div(Zpoly_t quotient, Zpoly_t input1,
                   Zpoly_t input2)
{
   abort();
}

void Zpoly_rem(Zpoly_t remainder, Zpoly_t input1,
                   Zpoly_t input2)
{
   abort();
}

void Zpoly_div_rem(Zpoly_t quotient, Zpoly_t remainder,
                       Zpoly_t input1, Zpoly_t input2)
{
   abort();
}

void Zpoly_gcd(Zpoly_t output, Zpoly_t input1,
                   Zpoly_t input2)
{
   abort();
}

void Zpoly_xgcd(Zpoly_t a, Zpoly_t b, Zpoly_t output,
                    Zpoly_t input1, Zpoly_t input2)
{
   abort();
}

void Zpoly_content(mpz_t content, Zpoly_t a)
{
   abort();
}

// *************** end of file
