/****************************************************************************

ZmodFpoly.h

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS_PER_LIMB. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_ZMODFPOLY_H
#define FLINT_ZMODFPOLY_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "memory-manager.h"
#include "mpn_extras.h"
#include "fmpz_poly.h"
#include "ZmodF.h"


/****************************************************************************

   ZmodFpoly_t
   -----------

ZmodFpoly_t represents a polynomial with coefficients in Z/pZ, where
p = B^n + 1, B = 2^FLINT_BITS_PER_LIMB. Coefficients are represented in the
format described in ZmodF.h.

Each polynomial has a fixed transform length 2^depth, specified at creation
time, where depth >= 0.

A polynomial may be in either "coefficient representation" (list of
coefficients of the polynomial), or "fourier representation" (list of
fourier coefficients). The polynomial does not keep track of which form it
is in, this is just a conceptual distinction.

x.length indicates how many coefficients contain meaningful data. If x is in
coefficient representation, the remaining coefficients are assumed to be
*zero*. If x is in fourier representation, the remaining coefficients are not
necessarily zero, they are simply *unknown*.

Always 0 <= length <= 2^depth.

Each polynomial carries a number of additional scratch buffers. The number of
scratch buffers is set at creation time. Various routines require a certain
number of scratch buffers to be present. The scratch buffers and coefficient
buffers are allocated as one large block, and routines may *permute* them,
so that outputs may well end up in what was originally a scratch buffer.

*/

typedef struct
{
   unsigned long depth;
   unsigned long n;
   unsigned long length;

   // Single chunk of memory where all coefficients live.
   mp_limb_t* storage;

   // Array of pointers to coefficients (length 2^depth).
   ZmodF_t* coeffs;

   // Array of pointers to scratch buffers (length scratch_count).
   unsigned long scratch_count;
   ZmodF_t* scratch;
   
} ZmodFpoly_struct;

// ZmodFpoly_t allows reference-like semantics for ZpolyFPoly_struct:
typedef ZmodFpoly_struct ZmodFpoly_t[1];

typedef ZmodFpoly_struct * ZmodFpoly_p;



/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

/*
   Initialises a ZmodFpoly_t with supplied parameters, and length = 0.
   Coefficients are not zeroed out.
*/
void ZmodFpoly_init(ZmodFpoly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count);
                    
void ZmodFpoly_stack_init(ZmodFpoly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count);


/*
   Frees resources for the given polynomial.
*/
void ZmodFpoly_clear(ZmodFpoly_t poly);

void ZmodFpoly_stack_clear(ZmodFpoly_t poly);


/* 
   Decrease the number of limbs n that are meaningful in a ZmodFpoly_t.
   The actual number of limbs allocated remains the same, only the field
   n is adjusted.
*/
static inline
void ZmodFpoly_decrease_n(ZmodFpoly_t poly, unsigned long n)
{
   FLINT_ASSERT(n <= poly->n);
   poly->n = n;
}


/****************************************************************************

   Conversion Routines
   
****************************************************************************/

/* 
   Converts fmpz_poly_t "poly_mpn" to a ZmodFpoly.
   
   Each coefficient of poly_mpn is assumed to fit into a coefficient 
   of poly_f.
   
   The maximum number of *bits* that any coefficient takes is returned
   and made negative if any of the coefficients was negative.
*/

long ZmodFpoly_convert_in_mpn(ZmodFpoly_t poly_f, fmpz_poly_t poly_mpn);


/* 
   Normalise and converts ZmodFpoly "poly_f" to a fmpz_poly_t. 
   
   Each coefficient of poly_f is assumed to fit into a coefficient 
   of poly_mpn. 
   
   The normalisation ensures that this function is the inverse of 
   ZmodFpoly_convert_in_mpn.
*/

void ZmodFpoly_convert_out_mpn(fmpz_poly_t poly_mpn, ZmodFpoly_t poly_f, long sign);

/* 
   Packs bundle coefficients, each padded out to the given number of limbs, into 
   the first coefficient of poly_f.
*/
void ZmodFpoly_limb_pack_mpn(ZmodFpoly_t poly_f, fmpz_poly_t poly_mpn,
                                           unsigned long bundle, long limbs);

/* 
   Unpacks bundle coefficients from the first coefficient of poly_f, each 
   assumed to be stored in a field of the given number of limbs.
*/
void ZmodFpoly_limb_unpack_mpn(fmpz_poly_t poly_mpn, ZmodFpoly_t poly_f, 
                                  unsigned long bundle, unsigned long limbs);

/* 
   Unpacks bundle coefficients from the first coefficient of poly_f, each 
   assumed to be stored in a field of the given number of limbs. Assumes the
   coefficients are unsigned.
*/
void ZmodFpoly_limb_unpack_unsigned_mpn(fmpz_poly_t poly_mpn, ZmodFpoly_t poly_f, 
                                  unsigned long bundle, unsigned long limbs);


/*
   Packs poly_mpn down to the bit into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it. Each of the original
   coefficients is packed into a bitfield "bits" bits wide including one 
   bit for a sign bit. 
   
   "bits" is assumed to be less than FLINT_BITS_PER_LIMB. If bits is 
   negative, the input poly is assumed to have signed coefficients.
*/ 
   
void ZmodFpoly_bit_pack_mpn(ZmodFpoly_t poly_f, fmpz_poly_t poly_mpn,
     unsigned long bundle, long bits);


/*
   Unpacks poly_f into poly_mpn. This is the inverse of ZmodFpoly_bitpack_mpn, so
   long as the final coefficient in the polynomial is positive.
   Each coeff of poly_f is assumed to contain "bundle" coefficients, each stored 
   in a bitfield "bits" bits wide with the most significant bit being reserved for
   the sign. 
   
   The total number of coefficients to be unpacked is given by the length of 
   poly_mpn. One must ensure each of the coefficients of poly_mpn are set to zero
   before calling this function for the first time since it adds to existing 
   coefficients of poly_mpn, rather than overwriting them.
   
   "bits" is assumed to be less than FLINT_BITS_PER_LIMB. 
*/ 
   
void ZmodFpoly_bit_unpack_mpn(fmpz_poly_t poly_mpn, ZmodFpoly_t poly_f, 
                              unsigned long bundle, unsigned long bits);
void ZmodFpoly_bit_unpack_unsigned_mpn(fmpz_poly_t poly_mpn, ZmodFpoly_t poly_f, 
                              unsigned long bundle, unsigned long bits);


     
/*
   Packs poly_mpn down to the byte into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it, each packed into a field
   "bytes" bytes wide.
   
   "coeff_bytes" is assumed to be at least FLINT_BITS_PER_LIMB/8, i.e. the
   coefficients are assumed to be at least a limb wide.
*/ 
   
void ZmodFpoly_byte_pack_mpn(ZmodFpoly_t poly_f, fmpz_poly_t poly_mpn,
                             unsigned long bundle, unsigned long coeff_bytes);

     
/*
   Unpacks poly_f into poly_mpn. Each coefficient of poly_f will have "bundle" 
   coefficients, each packed into a field "bytes" bytes wide.
   
   The total number of coefficients to be unpacked is given by the length of 
   poly_mpn.
   
   "coeff_bytes" is assumed to be at least FLINT_BITS_PER_LIMB/8, i.e. the
   coefficients are assumed to be at least a limb wide.
*/ 
   
void ZmodFpoly_byte_unpack_unsigned_mpn(fmpz_poly_t poly_m, mp_limb_t* array,
                               unsigned long bundle, unsigned long coeff_bytes);

void ZmodFpoly_byte_unpack_mpn(fmpz_poly_t poly_m, mp_limb_t* array,
                               unsigned long bundle, unsigned long coeff_bytes);


     
/*
   Splits each coefficient of poly_mpn into pieces "limbs" limbs long and 
   stores each piece into bundle coefficients of poly_f. 
*/
 
void ZmodFpoly_split_mpn(ZmodFpoly_t poly_f, fmpz_poly_t poly_mpn,
     unsigned long bundle, unsigned long limbs);

     
/*
   Combines each "bundle" coefficients of poly_f, each taken to be "limbs" 
   limbs long, into a coefficient of poly_mpn. 
   
   This function is used for testing purposed only, and is the exact inverse
   of ZmodFpoly_split_mpn.
   
   The number of coefficients extracted is given by the length of poly_mpn.
*/
 
void ZmodFpoly_unsplit_mpn(ZmodFpoly_t poly_f, fmpz_poly_t poly_mpn,
     unsigned long bundle, unsigned long limbs);



/****************************************************************************

   Basic Arithmetic Routines
   
****************************************************************************/

/*
   Sets x := y.

   Only y.length coefficients are copied.

   PRECONDITIONS:
      x and y must have compatible dimensions.
*/
void ZmodFpoly_set(ZmodFpoly_t x, ZmodFpoly_t y);


/*
   Sets res := pointwise product of x and y mod p.

   Only coefficients up to x.length are multiplied.

   PRECONDITIONS:
      Any combination of aliasing among res, x, y is allowed.
      x, y, res must have compatible dimensions.
      x and y must have the same length.

   NOTE:
      This function normalises the coefficients before multiplying.
*/
void ZmodFpoly_pointwise_mul(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y);


/*
   Sets res := x + y mod p.

   Only coefficients up to x.length are added.

   PRECONDITIONS:
      Any combination of aliasing among res, x, y is allowed.
      x, y, res must have compatible dimensions.
      x and y must have the same length.

   NOTE:
      This function does *not* normalise before subtracting. Be careful
      with the overflow limb.
*/
void ZmodFpoly_add(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y);


/*
   Sets res := x - y mod p.

   Only coefficients up to x.length are subtracted.

   PRECONDITIONS:
      Any combination of aliasing among res, x, y is allowed.
      x, y, res must have compatible dimensions.
      x and y must have the same length.
      
   NOTE:
      This function does *not* normalise before subtracting. Be careful
      with the overflow limb.
*/
void ZmodFpoly_sub(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y);


/*
   Normalises all coefficients (up to x.length) to be in the range [0, p).
*/
void ZmodFpoly_normalise(ZmodFpoly_t poly);


/*
   Divides all coefficients by 2^depth mod p. This should be used after
   running an inverse fourier transform.
*/
void ZmodFpoly_rescale(ZmodFpoly_t poly);


/****************************************************************************

   Fourier Transform Routines
   
For the following routines, 2^depth must divide 4*n*FLINT_BITS_PER_LIMB. This
ensures that Z/pZ has enough roots of unity.
   
****************************************************************************/


/*
This is the threshold for switching from a plain iterative FFT to an FFT
factoring algorithm. It should be set to about the number of limbs in L1 cache.
*/
//#define ZMODFPOLY_FFT_FACTOR_THRESHOLD 7500
#define ZMODFPOLY_FFT_FACTOR_THRESHOLD 7000


/*
   Converts from coefficient representation to fourier representation.

   "length" is the desired number of fourier coefficients; x.length is set
   to length when finished.

   Output is inplace. (Note that in general *all* 2^depth coefficients will
   get overwritten in intermediate steps.)

   PRECONDITIONS:
      0 <= length <= 2^poly.depth
      poly.scratch_count >= 1
*/
void ZmodFpoly_FFT(ZmodFpoly_t poly, unsigned long length);


/*
   Converts from fourier representation to coefficient representation.

   It *assumes* that the supplied fourier coefficients are actually the fourier
   transform of a polynomial whose coefficients beyond x.length are all zero.

   Result is inplace, x.length is not modified. (Note: after it's finished, the
   coefficients beyond x.length will contain garbage.)

   The output will be a factor of 2^depth too big. See ZmodFpoly_rescale().

   PRECONDITIONS:
      poly.scratch_count >= 1
*/
void ZmodFpoly_IFFT(ZmodFpoly_t poly);


/*
   Computes convolution of x and y, places result in res.

   The resulting length will be x.length + y.length - 1. If this is more
   than 2^depth, then the resulting length is 2^depth, and the convolution is
   actually cyclic of length 2^depth.

   PRECONDITIONS:
      Any combination of aliasing among res, x, y is allowed.
      x, y, res must have compatible dimensions.

   NOTE:
      x and y will both be converted to fourier representation.
      If you don't like it, make a copy first.

   PRECONDITIONS:
      x.scratch_count >= 1
      y.scratch_count >= 1
      res.scratch_count >= 1
*/
void ZmodFpoly_convolution(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y);


// internal functions

void _ZmodFpoly_FFT_iterative(
            ZmodF_t* x, unsigned long depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch);


void _ZmodFpoly_FFT_factor(
            ZmodF_t* x, unsigned long rows_depth, unsigned long cols_depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch);


void _ZmodFpoly_FFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                    unsigned long nonzero, unsigned long length,
                    unsigned long twist, unsigned long n,
                    ZmodF_t* scratch);


void _ZmodFpoly_IFFT_recursive(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long nonzero, unsigned long length, int extra,
               unsigned long twist, unsigned long n, ZmodF_t* scratch);


void _ZmodFpoly_IFFT_iterative(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long twist, unsigned long n, ZmodF_t* scratch);


void _ZmodFpoly_IFFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                     unsigned long nonzero, unsigned long length, int extra,
                     unsigned long twist, unsigned long n,
                     ZmodF_t* scratch);


/****************************************************************************

   Negacyclic Fourier Transform Routines
   
For the following routines, 2^(depth+1) must divide 4*n*FLINT_BITS_PER_LIMB.
This ensures that Z/pZ has enough roots of unity.

These routines are exactly the same as those listed in the previous section,
except that they evaluate at w^(2k+1), where w is a 2^(depth+1)-th root of
unity.
   
****************************************************************************/


void ZmodFpoly_negacyclic_FFT(ZmodFpoly_t poly);

void ZmodFpoly_negacyclic_IFFT(ZmodFpoly_t poly);

void ZmodFpoly_negacyclic_convolution(ZmodFpoly_t res,
                                      ZmodFpoly_t x, ZmodFpoly_t y);



#endif

// end of file ****************************************************************
