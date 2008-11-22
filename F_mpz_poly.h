/*============================================================================

    F_mpz_poly.h: Polynomials over Z (FLINT 2.0 polynomials)

    Copyright (C) 2008, William Hart 

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

==============================================================================*/

#ifndef FLINT_F_MPZ_POLY_H
#define FLINT_F_MPZ_POLY_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "mpz_poly.h"
#include "zmod_poly.h"
#include "flint.h"
#include "F_mpz.h"

/*
   WARNING : 
	=======

	When implementing functions for F_mpz_poly, memory leaks can
	occur if the length of a polynomial is shortened without using
	_F_mpz_poly_set_length (which demotes thus releases any mpz_t's 
	that were being used beyong the new length). 

	In some cases it is desirable to use F_mpz_poly_truncate instead. 
	It won't increase the length of a polynomial, but if it is made
	shorter it also normalises it.
*/

/*==============================================================================

   F_mpz_poly_t
   -----------

F_mpz_poly_t represents a dense polynomial in Z[x] 

"coeffs" is an array of F_mpz's (longs), one for each coefficient

There are two things each entry in this array can represent:

1) If the most significant two bits are 01, then the entry represents
an index into an array F_mpz_arr of mpz_t's, defined in the F_mpz module, 
and the mpz_t in that array contains the coefficient.

2) Otherwise, the entry represents a signed coefficient
whose absolute value is no more than FLINT_BIT - 2 bits in length. The
coefficient is stored in twos complement format.

"alloc" is the number of allocated coefficients. Obviously we always have
alloc >= length.

"length" is the length of the polynomial. If length == 0, this is the zero 
polynomial. All functions normalise so that the top coefficient is non-zero.

================================================================================*/
 
typedef struct
{
   F_mpz * coeffs;
   ulong alloc;
   ulong length;
} F_mpz_poly_struct;

// F_mpz_poly_t allows reference-like semantics for F_mpz_poly_struct
typedef F_mpz_poly_struct F_mpz_poly_t[1];

/*===============================================================================

	Memory management

================================================================================*/

/** 
   \fn     void F_mpz_poly_init(F_mpz_poly_t poly)
   \brief  Initialise a polynomial of length zero with zero allocated coefficients
*/
void F_mpz_poly_init(F_mpz_poly_t poly);

/** 
   \fn     void F_mpz_poly_init2(F_mpz_poly_t poly, const ulong alloc)
   \brief  Initialise a polynomial of length zero with the given number of allocated 
	        coefficients
*/
void F_mpz_poly_init2(F_mpz_poly_t poly, const ulong alloc);

/** 
   \fn     void F_mpz_poly_realloc(F_mpz_poly_t poly, const ulong alloc)
   \brief  Reallocates poly to have space for precisely the given number of 
	        coefficients. If alloc = 0, then the polynomial is cleared. If
			  the alloc is smaller than the current length of the polynomial
			  the polynomial is truncated and normalised.
*/
void F_mpz_poly_realloc(F_mpz_poly_t poly, const ulong alloc);

/** 
   \fn     void F_mpz_poly_fit_length(F_mpz_poly_t poly, const ulong length)
   \brief  Expands poly, if necessary, so that it has space for the given number
	        of coefficients. This function never shrinks the polynomial, only 
			  expands it.
*/
void F_mpz_poly_fit_length(F_mpz_poly_t poly, const ulong length);

/** 
   \fn     void F_mpz_poly_clear(F_mpz_poly_t poly)
   \brief  Clear the polynomial, releasing any memory it was using.
*/
void F_mpz_poly_clear(F_mpz_poly_t poly);

/*===============================================================================

	Normalisation

================================================================================*/

/** 
   \fn     void _F_mpz_poly_normalise(F_mpz_poly_t poly)
   \brief  Normalise poly so that the leading coefficient is nonzero or the 
	        polynomial has length zero
*/
void _F_mpz_poly_normalise(F_mpz_poly_t poly);

/*===============================================================================

	Coefficient operations

================================================================================*/

/** 
   \fn     void F_mpz_poly_set_coeff_si(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Set coefficient n to the signed long value x. Coefficients are numbered
	        from the constant coefficient, starting at zero.
*/
void F_mpz_poly_set_coeff_si(F_mpz_poly_t poly, ulong n, const long x);

/** 
   \fn     void F_mpz_poly_get_coeff_si(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Return coefficient n of poly as a signed long. If n is greater than the degree
	        of poly, then zero is returned.
*/
long F_mpz_poly_get_coeff_si(const F_mpz_poly_t poly, const ulong n);

/** 
   \fn     void F_mpz_poly_set_coeff_ui(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Set coefficient n to the unsigned long value x. Coefficients are numbered
	        from the constant coefficient, starting at zero.
*/
void F_mpz_poly_set_coeff_ui(F_mpz_poly_t poly, ulong n, const ulong x);

/** 
   \fn     void F_mpz_poly_get_coeff_ui(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Return coefficient n of poly as an usigned long. If n is greater than the degree
	        of poly, then zero is returned.
*/
ulong F_mpz_poly_get_coeff_ui(const F_mpz_poly_t poly, const ulong n);

/** 
   \fn     void F_mpz_poly_set_coeff_mpz(F_mpz_poly_t poly, ulong n, const mpz_t x)
   \brief  Set coefficient n to the mpz_t value x. Coefficients are numbered
	        from the constant coefficient, starting at zero.
*/
void F_mpz_poly_set_coeff_mpz(F_mpz_poly_t poly, ulong n, const mpz_t x);

/** 
   \fn     void F_mpz_poly_get_coeff_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong n)
   \brief  Return coefficient n of poly as an mpz_t. If n is greater than the degree
	        of poly, then zero is returned.
*/
void F_mpz_poly_get_coeff_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong n);

/*===============================================================================

	Attributes

================================================================================*/

/** 
   \fn     long F_mpz_poly_degree(const fmpz_poly_t poly)
   \brief  Returns the degree of the polynomial. If the polynomial is zero, then
	        minus one is returned for the degree.
*/
static inline 
long F_mpz_poly_degree(const F_mpz_poly_t poly)
{
   return poly->length - 1;
}

/** 
   \fn     unsigned long F_mpz_poly_length(const fmpz_poly_t poly)
   \brief  Returns the length of the polynomial. The zero polynomial has length zero.
*/
static inline 
unsigned long F_mpz_poly_length(const F_mpz_poly_t poly)
{
   return poly->length;
}


/*===============================================================================

	Truncation

================================================================================*/

/** 
   \fn     void F_mpz_poly_set_length(F_mpz_poly_t poly, const ulong length)
   \brief  Set the length of the polynomial to the given length. Assumes that
	        all the coefficients are valid coefficients, and does not normalise.
			  If the poly is made shorter, mpz_t's past the end are demoted.
*/
static inline
void _F_mpz_poly_set_length(F_mpz_poly_t poly, const ulong length)
{
	if (poly->length > length) // demote coefficients beyond new length
   {
      for (ulong i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);
		
   } 

	poly->length = length;
}

/** 
   \fn     void F_mpz_poly_truncate(F_mpz_poly_t poly, const ulong length)
   \brief  Truncate the polynomial to the given length. It is permissible for
	        length to be greater than the current length of the polynomial, in 
			  which case nothing will happen.
*/
static inline
void F_mpz_poly_truncate(F_mpz_poly_t poly, const ulong length)
{
	if (poly->length > length) // only truncate if necessary
   {
      for (ulong i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);
		poly->length = length;
      _F_mpz_poly_normalise(poly);
   }  
}

/*===============================================================================

	Conversions

================================================================================*/

/** 
   \fn     void mpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const mpz_poly_t m_poly)
   \brief  Convert an mpz_poly_t to an F_mpz_poly_t
*/
void mpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const mpz_poly_t m_poly);

/** 
   \fn     void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly)
   \brief  Convert an F_mpz_poly_t to an mpz_poly_t
*/
void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly);


/*===============================================================================

	Assignment

================================================================================*/

/** 
   \fn     void F_mpz_poly_zero(F_mpz_poly_t poly)
   \brief  Sets poly to the zero polynomial
*/
static inline 
void F_mpz_poly_zero(F_mpz_poly_t poly)
{
   _F_mpz_poly_set_length(poly, 0);
}

/** 
   \fn     void F_mpz_poly_zero_coeffs(F_mpz_poly_t poly, const ulong n)
   \brief  Zero the first n coefficients of poly regardless of its length
*/
static inline
void F_mpz_poly_zero_coeffs(F_mpz_poly_t poly, const ulong n)
{
	if (n >= poly->length) _F_mpz_poly_set_length(poly, 0);
	else 
	   for (ulong i = 0; i < n; i++)
		   F_mpz_zero(poly->coeffs + i);
}

/** 
   \fn     void F_mpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2)
   \brief  Sets poly1 to equal poly2
*/
void F_mpz_poly_set(F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_swap(F_mpz_poly_t poly1, F_mpz_poly_t poly2)
   \brief  Efficiently swap poly1 and poly2
*/
void F_mpz_poly_swap(F_mpz_poly_t poly1, F_mpz_poly_t poly2);

/*===============================================================================

	Comparison

================================================================================*/

/** 
   \fn     int F_mpz_poly_equal(const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Returns 1 if poly1 and poly2 are equal (arithmetically), otherwise
	        returns 0.
*/
int F_mpz_poly_equal(const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/*===============================================================================

	Coefficient sizes

================================================================================*/

/** 
   \fn     long F_mpz_poly_max_bits(F_mpz_poly_t poly)
   \brief  Computes the largest number of bits n that any coefficient has and returns
	        -n if a negative coefficient exists in poly, else it returns n. Zero is 
			  returned for the zero polynomial.
*/
long F_mpz_poly_max_bits(const F_mpz_poly_t poly);

/** 
   \fn     ulong F_mpz_poly_max_limbs(F_mpz_poly_t poly)
   \brief  Returns the largest number of limbs required to store the absolute value
	        of coefficients of poly. Zero is returned for the zero polynomial.
*/
ulong F_mpz_poly_max_limbs(const F_mpz_poly_t poly);

/*===============================================================================

	Reverse

================================================================================*/

/** 
   \fn     void F_mpz_poly_reverse(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong length)
   \brief  Treats poly as though it were the given length (with leading zeroes if 
	        necessary) and sets res to the reverse polynomial.
*/
void F_mpz_poly_reverse(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong length);

/*===============================================================================

	Negation

================================================================================*/

/** 
   \fn     void F_mpz_poly_neg(F_mpz_poly_t res, const F_mpz_poly_t poly)
   \brief  Set res to the negative of poly.
*/
void F_mpz_poly_neg(F_mpz_poly_t res, const F_mpz_poly_t poly);

/*===============================================================================

	Addition/subtraction

================================================================================*/

/** 
   \fn     void F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Sets res to the sum of poly1 and poly2.
*/
void F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Sets res to the difference of poly1 and poly2.
*/
void F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/*===============================================================================

	Shifting

================================================================================*/

/** 
   \fn     void F_mpz_poly_left_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n)
   \brief  Multiplies poly by x^n and sets res to the result.
*/
void F_mpz_poly_left_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n);

/** 
   \fn     void F_mpz_poly_right_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n)
   \brief  Divides poly by x^n, discarding any remainder, and sets res to the result.
*/
void F_mpz_poly_right_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n);

/*===============================================================================

	Scalar multiplication

================================================================================*/

/** 
   \fn     void F_mpz_poly_scalar_mul_ui(F_mpz_poly_t poly1, F_mpz_poly_t poly2, ulong x)
   \brief  Multiply poly2 by the unsigned long x and set poly1 to the result.
*/
void F_mpz_poly_scalar_mul_ui(F_mpz_poly_t poly1, F_mpz_poly_t poly2, ulong x);

/** 
   \fn     void F_mpz_poly_scalar_mul_si(F_mpz_poly_t poly1, F_mpz_poly_t poly2, long x)
   \brief  Multiply poly2 by the signed long x and set poly1 to the result.
*/
void F_mpz_poly_scalar_mul_si(F_mpz_poly_t poly1, F_mpz_poly_t poly2, long x);

/** 
   \fn     void F_mpz_poly_scalar_mul(F_mpz_poly_t poly1, F_mpz_poly_t poly2, F_mpz_t x)
   \brief  Multiply poly2 by the F_mpz_t x and set poly1 to the result.
*/
void F_mpz_poly_scalar_mul(F_mpz_poly_t poly1, F_mpz_poly_t poly2, F_mpz_t x);

/*===============================================================================

	Multiplication

================================================================================*/

/** 
   \fn     void F_mpz_poly_mul_classical(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result, using the classical 
	        algorithm.
*/
void F_mpz_poly_mul_classical(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void _F_mpz_poly_mul_kara_recursive(F_mpz_poly_t out, ulong ostart, F_mpz_poly_t in1, ulong istart1, 
											              ulong len1, F_mpz_poly_t in2, ulong istart2, ulong len2, 
											              F_mpz_poly_t scratch, ulong sstart, ulong skip, ulong crossover)
   \brief  Recursive portion of karatsuba multiplication.

           Input polys are in1 and in2, each staggered by skip. We specify a starting point 
           in the coefficients for the output and inputs and scratch space and a length, 
           len1 and len2 for each of the intputs. Then out will be of length len1 + len2 - 1, 
           also staggered by skip.

           The scratch buffer should be length len1 + len2, also staggered by skip.

           All input/output/scratch polys should be initialised, and shouldn't overlap.

           Must have 1 <= len1 <= len2.

           If len1*len2 <= crossover, we use the classical multiplication algorithm. 
           The crossover parameter is passed down recursively to subproducts.
*/
void _F_mpz_poly_mul_kara_recursive(F_mpz_poly_t out, ulong ostart, F_mpz_poly_t in1, ulong istart1, 
											ulong len1, F_mpz_poly_t in2, ulong istart2, ulong len2, 
											F_mpz_poly_t scratch, ulong sstart, ulong skip, ulong crossover);

/** 
   \fn     void F_mpz_poly_mul_karatsuba(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result, using the karatsuba method.
*/
void F_mpz_poly_mul_karatsuba(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2);

/*===============================================================================

	Bit/byte/limb packing

================================================================================*/

/** 
   \fn     void F_mpz_poly_bit_pack(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                                           const long bitwidth, const ulong length, const long negate)

   \brief  Pack length coefficients into the given array of limbs, with each coefficient 
	        packed into a bitfield with the given bitwidth. If negate is -1L then each 
			  coefficient will be negated before being packed. The array is sign extended to 
			  the end of the array. If bitwidth is negative this means the coefficients 
			  are expected to be signed. Each coefficient is stored in bit packed format in 
			  twos complement format and if negative, 1 is borrowed from the next coefficient
			  before it is packed. The number of limbs in the array is n.
*/
void F_mpz_poly_bit_pack(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                         const long bitwidth, const ulong length, const long negate);


/** 
   \fn     void F_mpz_poly_bit_pack_unsigned(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                                                                    const ulong bits, const ulong length)

   \brief  Pack bundle coefficients at a time into the given array with each coefficient 
	        packed into a bitfield with the given bitwidth. The number of limbs in the 
			  array is n.
*/
void F_mpz_poly_bit_pack_unsigned(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                                                            const ulong bits, const ulong length);
/** 
   \fn     void F_mpz_poly_bit_pack2(mp_limb_t * array, mp_limb_t * array2, ulong n, const F_mpz_poly_t poly_F_mpz, 
                                            const long bitwidth, const ulong length, const long negate, long negate2)

   \brief  As for F_mpz_poly_bit_pack, except that two arrays are packed. The second
	        is packed from the same data, except the signs are alternated before packing, starting
			  with the sign given by negative2 for the constant coefficient.

*/
void F_mpz_poly_bit_pack2(mp_limb_t * array, mp_limb_t * array2, ulong n, const F_mpz_poly_t poly_F_mpz, 
                                 const long bitwidth, const ulong length, const long negate, long negate2);

/** 
   \fn     void F_mpz_poly_bit_unpack(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, ulong n,
                                                    const ulong bundle, const ulong bits);
   \brief  Coefficients are unpacked from the array from fields of the given number
	        of bits in width. The coefficients are assumed to be signed. If a negative 
			  coefficient is unpacked, the next coefficient gets 1 added to it for the borrow
			  that was effectively made by that previous coefficient. The output coefficients 
			  are written to the coefficients of the F_mpz_poly_t poly_F_mpz, and a total of
			  length coefficients will be written.
*/
void F_mpz_poly_bit_unpack(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, ulong n,
                                                    const ulong length, const ulong bits);

/** 
   \fn     void F_mpz_poly_bit_unpack_unsigned(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, ulong n, 
                                                                     const ulong length, const ulong bits)

   \brief  Coefficients are unpacked from the array from fields of the given number
	        of bits in width. The coefficients are assumed to be unsigned. The output coefficients 
			  are written to the coefficients of the F_mpz_poly_t poly_F_mpz. A total of length
			  coefficients will be written.
*/
void F_mpz_poly_bit_unpack_unsigned(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, ulong n, 
                                                             const ulong length, const ulong bits);

/** 
   \fn     F_mpz_poly_mul_KS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result, using the Kronecker method.
*/
void F_mpz_poly_mul_KS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

void F_mpz_poly_mul_KS2(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
