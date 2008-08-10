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
#include "zmod_poly.h"
#include "flint.h"

/*==============================================================================

   F_mpz_poly_t
   -----------

F_mpz_poly_t represents a dense polynomial in Z[x] 

"coeffs" is an array of limbs, one for each coefficient

There are two things each entry in this array can represent:

1) If the next to most significant bit is set, the entry represents
an index into the array "mpz_coeffs", and the mpz_t in that array
contains the coefficient.

2) If the next to msb is not set, the entry represents a signed coefficient
whose absolute value is no more than FLINT_BIT - 2 bits in length. The
coefficient is stored in sign/absolute value format with the sign being 
indicated by the most significant bit (1 for -ve) and the remaining bits 
corresponding to the absolute value.

----------------------------------------------
| sign | mpz | abs | abs | abs | ...... | abs |
----------------------------------------------

"mpz_coeffs" is an array of coefficients in mpz_t format (actually an array of
__mpz_struct's). Only coefficients which do not fit into FLINT_BITS - 2 bits are 
stored in an mpz_t.

"alloc" is the number of allocated coefficients. Obviously we always have
alloc >= length.

"length" is the length of the polynomial. If length == 0, this is the zero 
polynomial. All functions normalise so that the (length-1)-th coefficient 
is non-zero.

"mpz_alloc" is the number of allocated mpz_coefficients

"mpz_length" is the number of mpz_t's actually used

================================================================================*/
 
typedef struct
{
   mp_limb_t * coeffs;
   __mpz_struct * mpz_coeffs;
   ulong alloc;
   ulong length;
   ulong mpz_alloc;
   ulong mpz_length;
} F_mpz_poly_struct;

// fmpz_poly_t allows reference-like semantics for F_mpz_poly_struct
typedef F_mpz_poly_struct F_mpz_poly_t[1];

#define MPZ_BLOCK 16 // number of additional mpz_t's to initialise at a time

#define FMPZ_UVAL_MASK ((1L<<(FLINT_BITS-2))-1L) // mask with bottom FLINT_BITS - 2 bits set

#define FMPZ_MPZ_MASK (1L<<(FLINT_BITS - 2)) // mask with just the mpz bit set

#define FMPZ_SIGN_MASK (1L<<(FLINT_BITS - 1)) // mask with just the sign bit set

#define FMPZ_UVAL(xxx) ((xxx) & FMPZ_UVAL_MASK) // the bottom FLINT_BITS - 2 bits of xxx

#define FMPZ_MPZ(xxx) ((xxx) | FMPZ_MPZ_MASK) // xxx + the mpz bit set

#define FMPZ_IS_MPZ(xxx) ((xxx) & FMPZ_MPZ_MASK) // is the mpz bit set pf xxx set

/*===============================================================================

	Memory management

================================================================================*/

/** 
   \fn     void _F_mpz_poly_mpz_coeffs_new(F_mpz_poly_t poly)
   \brief  Add a new mpz_coeff. The mpz_t's are allocated and initialised
	        in blocks size MPZ_BLOCK. poly->mpz_length is incremented and 
			  poly->mpz_alloc is updated if new mpz_t coefficients were added.
*/
void _F_mpz_poly_mpz_coeffs_new(F_mpz_poly_t poly);

/** 
   \fn     void _F_mpz_poly_mpz_coeffs_clear(F_mpz_poly_t poly)
   \brief  Clear poly's array of mpz_t coefficients.
*/
void _F_mpz_poly_mpz_coeffs_clear(F_mpz_poly_t poly);

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
   \fn     void _F_mpz_poly_mpz_coeffs_clear(F_mpz_poly_t poly)
   \brief  Normalise poly so that the leading coefficient is nonzero or the 
	        polynomial has length zero
*/
void _F_mpz_poly_normalise(F_mpz_poly_t poly);

/*===============================================================================

	Coefficient operations

================================================================================*/

/** 
   \fn     void _F_mpz_promote(F_mpz_poly_t const poly, ulong coeff);
   \brief  Promote the given coefficient of poly to an mpz_t coefficient. The value
	        of the coefficient is not preserved.
*/
static inline
void _F_mpz_promote(F_mpz_poly_t poly, const ulong coeff)
{
   _F_mpz_poly_mpz_coeffs_new(poly);
	
	poly->coeffs[coeff] = FMPZ_MPZ(poly->mpz_length - 1);
}

/** 
   \fn     _F_mpz_demote_val(F_mpz_poly_t poly, const ulong coeff);
   \brief  If the given coefficient (which is assumed to be an mpz_t) will fit into
	        FLINT_BIT - 2 bits, it is demoted to a limb instead of an mpz_t, otherwise
			  nothing happens.
*/
void _F_mpz_demote_val(F_mpz_poly_t poly, const ulong coeff);

/** 
   \fn     _F_mpz_set_si(F_mpz_poly_t poly, ulong coeff, const long val)
   \brief  Set the given coefficient of poly to a signed value val
*/
void _F_mpz_set_si(F_mpz_poly_t poly, ulong coeff, const long val);

/** 
   \fn     _F_mpz_get_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong coeff)
   \brief  Returns the given coefficient of poly as an mpz_t
*/
void _F_mpz_get_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong coeff);

/** 
   \fn     _F_mpz_set_mpz(F_mpz_poly_t poly, ulong coeff, const mpz_t x)
   \brief  Sets the given coefficient to the given mpz_t
*/
void _F_mpz_set_mpz(F_mpz_poly_t poly, ulong coeff, const mpz_t x);

/** 
   \fn     _F_mpz_add(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
					                                       const F_mpz_poly_t poly2, const ulong coeff2)
   \brief  Add the given coefficients of poly1 and poly2 and set the given coefficient 
	        of res to the result. Assumes the coefficient of res is distinct from the 
			  other two coefficients.

*/
void _F_mpz_add(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
					                                 const F_mpz_poly_t poly2, const ulong coeff2);

/*===============================================================================

	Truncation

================================================================================*/

/** 
   \fn     void F_mpz_poly_truncate(F_mpz_poly_t poly, const ulong length)
   \brief  Truncate the polynomial to the given length. It is permissible for
	        length to be greater than the current length of the polynomial, in 
			  which case nothing will happen.
*/
static inline
void F_mpz_poly_truncate(F_mpz_poly_t poly, const ulong length)
{
	if (poly->length > length)
   {
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

#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
