/****************************************************************************

ZmodFpoly.c

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS_PER_LIMB. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <sys/types.h>
#include "flint.h"
#include "flint-manager.h"
#include "ZmodFpoly.h"
#include "Zpoly_mpn.h"
#include "extras.h"


/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodFpoly_init(ZmodFpoly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->depth = depth;
   poly->scratch_count = scratch_count;
   poly->length = 0;
   
   unsigned long bufs = (1 << depth) + scratch_count;

   // todo: change these to use the FLINT heap allocator
   
   poly->storage = (mp_limb_t*) malloc(bufs * (n+1) * sizeof(mp_limb_t));

   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) malloc(bufs * sizeof(ZmodF_t*));
   poly->scratch = poly->coeffs + (1 << depth);
   
   poly->coeffs[0] = poly->storage;
   for (unsigned long i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}


void ZmodFpoly_clear(ZmodFpoly_t poly)
{
   free(poly->coeffs);
   free(poly->storage);
}



/****************************************************************************

   Conversion Routines
   
****************************************************************************/

void ZmodFpoly_convert_in_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   unsigned long size_f = poly_f->n + 1;
   unsigned long size_m = poly_mpn->limbs+1;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;
   
   long size_j;
   
   for (unsigned long i = 0, j = 0; i < poly_mpn->length; i++, j += size_m)
   {
      if ((size_j = coeffs_m[j]) < 0)
      {
         negate_limbs(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         set_limbs(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      } else
      {
         copy_limbs(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         clear_limbs(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      }
   }
   poly_f->length = poly_mpn->length;   
}

void ZmodFpoly_convert_out_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   unsigned long n = poly_f->n;
   unsigned long size_m = poly_mpn->limbs+1;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;

   for (unsigned long i = 0, j = 0; i < poly_f->length; i++, j += size_m)
   {
      ZmodF_normalise(coeffs_f[i], n);
      if (coeffs_f[i][n-1]>>(FLINT_BITS_PER_LIMB-1))
      {
         negate_limbs(coeffs_m + j + 1, coeffs_f[i], n);
         mpn_add_1(coeffs_m + j + 1, coeffs_m + j + 1, n, 1L);
         coeffs_m[j] = -n;
         NORM(coeffs_m + j);
      } else
      {
         copy_limbs(coeffs_m + j + 1, coeffs_f[i], n);
         coeffs_m[j] = n;
         NORM(coeffs_m + j);         
      }
   }
}

static inline long __get_next_coeff(mp_limb_t * coeff_m, long * borrow, long mask)
{ 
   long coeff;
   
   if ((long) coeff_m[0] >= 0) coeff = coeff_m[1] - *borrow;
   else coeff = (-coeff_m[1] - *borrow);
   *borrow = 0UL;
   if (coeff < 0) 
   {
      *borrow = 1UL;
      coeff&=mask;
   }
   
   return coeff;
}

void ZmodFpoly_bit_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                            unsigned long bundle, unsigned long bits)
{   
   unsigned long i, k, skip;
   unsigned long n = poly_f->n;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * array, * next_point;
    
   unsigned long temp;
   half_ulong lower;
   long coeff;
   long borrow;

   const unsigned long mask = (1UL<<bits)-1;
   
   poly_f->length = 0;
      
   while (coeff_m < poly_mpn->coeffs + 2*poly_mpn->length)
   {
      k=0; skip=0;
      coeff = 0UL; borrow = 0L; temp = 0;
      array = poly_f->coeffs[i];
      i++;
   
      next_point = coeff_m + 2*bundle;
      if (next_point > poly_mpn->coeffs + 2*poly_mpn->length) 
         next_point = poly_mpn->coeffs + 2*poly_mpn->length;
         
      while (coeff_m < next_point)
      {
         // k is guaranteed to be less than FLINT_BITS_PER_LIMB at this point
         while ((k<HALF_FLINT_BITS_PER_LIMB)&&(coeff_m < next_point))
         {
            temp+=(__get_next_coeff(coeff_m, &borrow, mask) << k);
            coeff_m+=2; k+=bits;
         }
         // k may exceed FLINT_BITS_PER_LIMB at this point but is less than 96

         if (k>FLINT_BITS_PER_LIMB)
         {
            // if k > FLINT_BITS_PER_LIMB write out a whole limb and read in remaining bits of coeff
            array[skip] = temp;
            skip++;
            temp=(coeff>>(bits+FLINT_BITS_PER_LIMB-k));
            k=(k-FLINT_BITS_PER_LIMB);
            // k < HALF_FLINT_BITS_PER_LIMB
         } else
         {
            // k <= FLINT_BITS_PER_LIMB
            if (k >= HALF_FLINT_BITS_PER_LIMB)
            {
               // if k >= HALF_FLINT_BITS_PER_LIMB store bottom HALF_FLINT_BITS_PER_LIMB bits
               lower = (half_ulong)temp;
               k-=HALF_FLINT_BITS_PER_LIMB;
               temp>>=HALF_FLINT_BITS_PER_LIMB;
               // k is now <= HALF_FLINT_BITS_PER_LIMB

               while ((k<HALF_FLINT_BITS_PER_LIMB)&&(coeff_m < next_point))
               {
                  temp+=(__get_next_coeff(coeff_m, &borrow, mask) << k);
                  coeff_m+=2; k+=bits;
               }
               // k may again exceed FLINT_BITS_PER_LIMB bits but is less than 96
               if (k>FLINT_BITS_PER_LIMB)
               {
                  // if k > FLINT_BITS_PER_LIMB, write out bottom HALF_FLINT_BITS_PER_LIMB bits (along with HALF_FLINT_BITS_PER_LIMB bits from lower)
                  // read remaining bits from coeff and reduce k by HALF_FLINT_BITS_PER_LIMB
                  array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+(unsigned long)lower;
                  skip++;
                  temp>>=HALF_FLINT_BITS_PER_LIMB;
                  temp+=((coeff>>(bits+FLINT_BITS_PER_LIMB-k))<<HALF_FLINT_BITS_PER_LIMB);
                  k=(k-HALF_FLINT_BITS_PER_LIMB);
                  // k < FLINT_BITS_PER_LIMB and we are ready to read next coefficient if there is one
               } else if (k >= HALF_FLINT_BITS_PER_LIMB) 
               {
                  // k <= FLINT_BITS_PER_LIMB
                  // if k >= HALF_FLINT_BITS_PER_LIMB write out bottom HALF_FLINT_BITS_PER_LIMB bits (along with lower)
                  // and reduce k by HALF_FLINT_BITS_PER_LIMB
                  k-=HALF_FLINT_BITS_PER_LIMB;
                  array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
                  temp>>=HALF_FLINT_BITS_PER_LIMB;
                  skip++;
                  // k is now less than or equal to HALF_FLINT_BITS_PER_LIMB and we are now ready to read 
                  // the next coefficient if there is one
               } else
               {
                  // k < HALF_FLINT_BITS_PER_LIMB
                  // there isn't enough to write out a whole FLINT_BITS_PER_LIMB bits, so put it all 
                  // together in temp
                  temp = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
                  k+=HALF_FLINT_BITS_PER_LIMB;
                  // k is now guaranteed to be less than FLINT_BITS_PER_LIMB and we are ready for the
                  // next coefficient if there is one
               }
            } // if
         } // else
         poly_f->length++;
      } // while

      // sign extend the last FLINT_BITS_PER_LIMB bits we write out
      if (skip < n)
      {
        if (borrow) temp+=(-1UL << k);
        array[skip] = temp;
        skip++;
      } 
      // sign extend the remainder of the array, reducing modulo p 
      if (borrow)
      {
         while (skip < n) 
         {
            array[skip] = -1UL;
            skip++;
         }
         mpn_add_1(array, array, n+1, 1UL);
      }
   } // for
}


void ZmodFpoly_bit_unpack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                              unsigned long bundle, unsigned long bits)
{
   unsigned long k, l, skip;

   unsigned long temp2;
   unsigned long temp;
   unsigned long full_limb;
   unsigned long carry;
    
   mp_limb_t* array;
    
   const unsigned long mask = (1UL<<bits)-1;
   const unsigned long sign_mask = (1UL<<(bits-1));

   unsigned long s;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * next_point;
   unsigned long size_m = poly_mpn->limbs+1;
    
   for (unsigned long i = 0; coeff_m < poly_mpn->coeffs + poly_mpn->length*size_m; i++)
   {
      array = poly_f->coeffs[i];
      
      k=0; skip=0; carry = 0UL; temp2 = 0;
      next_point = coeff_m + size_m*bundle;
      if (next_point > poly_mpn->coeffs + poly_mpn->length*size_m) next_point = poly_mpn->coeffs + poly_mpn->length*size_m;
      
      while (coeff_m < next_point)
      {
         // read in a full limb
         full_limb = array[skip];
         temp2 += l_shift(full_limb,k);
         s=FLINT_BITS_PER_LIMB-k;
         k+=s;
         while ((k >= bits)&&(coeff_m < next_point))
         {
            if (!(temp2&sign_mask)) 
            {
               __Zpoly_mpn_add_coeff_ui(coeff_m, (temp2&mask)+carry);
               carry = 0UL;
            }  
            else
            {
               temp = ((-temp2)&mask)-carry;
               __Zpoly_mpn_sub_coeff_ui(coeff_m, temp);
               carry = 1UL;
            }
            coeff_m += size_m;
            temp2>>=bits;
            k-=bits;
         }
         // k is now less than bits
         // read in remainder of full_limb
         temp2 += l_shift(r_shift(full_limb,s),k);
         k+=(FLINT_BITS_PER_LIMB-s);
       
         while ((k >= bits)&&(coeff_m < next_point))
         {
            if (!(temp2&sign_mask)) 
            {
               __Zpoly_mpn_add_coeff_ui(coeff_m, (temp2&mask)+carry);
               carry = 0UL;
            }
            else
            {
               temp = ((-temp2)&mask)-carry;
               __Zpoly_mpn_sub_coeff_ui(coeff_m, temp);
               carry = 1UL;
            }
            coeff_m += size_m;
            temp2>>=bits;
            l++;
            k-=bits;
         }
         // k is now less than bits
         skip++;
      }
   }
}
     
void ZmodFpoly_byte_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                             unsigned long bundle, unsigned long bytes)
{
   abort();
}
     
void ZmodFpoly_byte_unpack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                               unsigned long bundle, unsigned long bytes)
{
   abort();
}

     
void ZmodFpoly_split_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                         unsigned long bundle, unsigned long limbs)
{
   abort();
}
     
void ZmodFpoly_unsplit_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                           unsigned long bundle, unsigned long limbs)
{
   abort();
}
     

/****************************************************************************

   Basic Arithmetic Routines
   
****************************************************************************/


void ZmodFpoly_set(ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->n == y->n);

   for (unsigned long i = 0; i < y->length; i++)
      ZmodF_set(x->coeffs[i], y->coeffs[i], x->n);
   x->length = y->length;
}


void ZmodFpoly_mul(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   // todo: use FLINT memory allocation
   
   mp_limb_t* scratch = (mp_limb_t*) malloc((2 * x->n) * sizeof(mp_limb_t));

   if (x != y)
      for (unsigned long i = 0; i < x->length; i++)
         ZmodF_mul(res->coeffs[i], x->coeffs[i], y->coeffs[i], scratch, x->n);
   else
      for (unsigned long i = 0; i < x->length; i++)
         ZmodF_sqr(res->coeffs[i], x->coeffs[i], scratch, x->n);
   
   free(scratch);
}


void ZmodFpoly_add(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   for (unsigned long i = 0; i < x->length; i++)
      ZmodF_add(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodFpoly_sub(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   for (unsigned long i = 0; i < x->length; i++)
      ZmodF_sub(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodFpoly_normalise(ZmodFpoly_t poly)
{
   for (unsigned long i = 0; i < poly->length; i++)
      ZmodF_normalise(poly->coeffs[i], poly->n);
}


void ZmodFpoly_rescale(ZmodFpoly_t poly)
{
   if (poly->depth == 0)
      return;

   for (unsigned long i = 0; i < poly->length; i++)
      ZmodF_short_div_2exp(poly->coeffs[i], poly->coeffs[i],
                           poly->depth, poly->n);
}


/****************************************************************************

   Fourier Transform Routines

****************************************************************************/


/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of buffers assumed to be nonzero
length = number of fourier coefficients requested
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_FFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                    unsigned long nonzero, unsigned long length,
                    unsigned long twist, unsigned long n,
                    ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);
   
   // ========================
   // base cases
   
   if (depth == 0)
      return;
      
   if (depth == 1)
   {
      if (length == 1)
      {
         if (nonzero == 2)
            ZmodF_add(x[0], x[0], x[skip], n);
      }
      else   // length == 2
      {
         if (nonzero == 1)
            ZmodF_mul_sqrt2exp(x[skip], x[0], twist, n);
         else  // nonzero == 2
            ZmodF_forward_butterfly_sqrt2exp(x, x+skip, scratch, twist, n);
      }
      
      return;
   }
   
   // ========================
   // factoring case

   unsigned long rows_depth = depth >> 1;
   unsigned long cols_depth = depth - rows_depth;
   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long length_whole_rows = length_cols ?
                                     (length_rows + 1) : length_rows;
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // column transforms
   for (i = 0, y = x, j = twist; i < nonzero_cols; i++, y += skip, j += root)
      _ZmodFpoly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                     length_whole_rows, j, n, scratch);

   if (nonzero_rows)
   {
      for (; i < cols; i++, y += skip, j += root)
         _ZmodFpoly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                        length_whole_rows, j, n, scratch);
      nonzero_cols = cols;
   }
   
   // row transforms
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodFpoly_FFT(y, cols_depth, skip, nonzero_cols, cols,
                     twist << cols_depth, n, scratch);

   if (length_cols)
      // The relevant portion of the last row:
      _ZmodFpoly_FFT(y, cols_depth, skip, nonzero_cols, length_cols,
                     twist << cols_depth, n, scratch);
}



/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of *output* buffers assumed to be nonzero
length = number of untransformed coefficients requested
extra = indicates whether an extra *forward* coefficient should be computed
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_IFFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                     unsigned long nonzero, unsigned long length, int extra,
                     unsigned long twist, unsigned long n,
                     ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);
   
   // ========================
   // base cases
   
   if (depth == 0)
      return;
      
   if (depth == 1)
   {
      if (length == 0)
      {
         if (nonzero == 2)
            ZmodF_add(x[0], x[0], x[skip], n);
         ZmodF_short_div_2exp(x[0], x[0], 1, n);
      }
      else if (length == 1)
      {
         if (nonzero == 1)
         {
            if (extra)
               ZmodF_mul_sqrt2exp(x[skip], x[0], twist, n);
            ZmodF_add(x[0], x[0], x[0], n);
         }
         else  // nonzero == 2
         {
            if (extra)
            {
               ZmodF_sub(scratch[0], x[0], x[skip], n);
               ZmodF_add(x[0], x[0], scratch[0], n);
               ZmodF_mul_sqrt2exp(x[skip], scratch[0], twist, n);
            }
            else
            {
               ZmodF_add(x[0], x[0], x[0], n);
               ZmodF_sub(x[0], x[0], x[skip], n);
            }
         }
      }
      else  // length == 2
         ZmodF_inverse_butterfly_sqrt2exp(x, x + skip, scratch, twist, n);

      return;
   }
   
   // ========================
   // factoring case

   unsigned long rows_depth = depth >> 1;
   unsigned long cols_depth = depth - rows_depth;
   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // row transforms for the rows where we have all fourier coefficients
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodFpoly_IFFT(y, cols_depth, skip, cols, cols, 0,
                      twist << rows_depth, n, scratch);

   // column transforms where we have enough information
   for (i = length_cols, y = x + (skip * length_cols),
        j = twist + (root*length_cols);
        i < nonzero_cols; i++, y += skip, j += root)
   {
      _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                      length_rows, length_cols ? 1 : extra, j, n, scratch);
   }
   if (nonzero_rows)
      for (; i < cols; i++, y += skip, j += root)
         _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                         length_rows, length_cols ? 1 : extra, j, n, scratch);

   if (length_cols)
   {
      // a single switcheroo row transform
      _ZmodFpoly_IFFT(x + length_rows * (skip << cols_depth), cols_depth,
                      skip, (nonzero_rows ? cols : nonzero_cols),
                      length_cols, extra, twist << rows_depth, n, scratch);

      // remaining column transforms
      for (i = 0, y = x, j = twist; i < length_cols && i < nonzero_cols;
           i++, y += skip, j += root)
      {
         _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                         length_rows + 1, 0, j, n, scratch);
      }
      if (nonzero_rows)
      {
         for (; i < length_cols; i++, y += skip, j += root)
            _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                            length_rows + 1, 0, j, n, scratch);
      }
   }
   else if (extra)
   {
      // need one extra trivial fourier coefficient
      x += length_rows * (skip << cols_depth);
      for (i = 1, y = x + skip; i < (nonzero_rows ? cols : nonzero_cols);
           i++, y += skip)
      {
         ZmodF_add(x[0], x[0], y[0], n);
      }
      ZmodF_short_div_2exp(x[0], x[0], cols_depth, n);
   }
}


void ZmodFpoly_FFT(ZmodFpoly_t poly, unsigned long length)
{
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (length != 0)
   {
      if (poly->length == 0)
      {
         // input is zero, so output is zero too
         for (unsigned long i = 0; i < length; i++)
            ZmodF_zero(poly->coeffs[i], poly->n);
      }
      else
      {
         _ZmodFpoly_FFT(poly->coeffs, poly->depth, 1, poly->length,
                        length, 0, poly->n, poly->scratch);
      }
   }

   poly->length = length;
}


void ZmodFpoly_IFFT(ZmodFpoly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (poly->length != 0)
      _ZmodFpoly_IFFT(poly->coeffs, poly->depth, 1, poly->length,
                      poly->length, 0, 0, poly->n, poly->scratch);
}


void ZmodFpoly_convolution(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long length = x->length + y->length - 1;
   unsigned long size = 1UL << res->depth;
   if (length > size)
      length = size;
   
   ZmodFpoly_FFT(x, length);
   if (x != y)    // take care of aliasing
      ZmodFpoly_FFT(y, length);
      
   ZmodFpoly_mul(res, x, y);
   ZmodFpoly_IFFT(res);
}


/****************************************************************************

   Negacyclic Fourier Transform Routines
   
****************************************************************************/


void ZmodFpoly_negacyclic_FFT(ZmodFpoly_t poly, unsigned long length)
{
   // todo: I *think* this code is wrong
   abort();
      
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (length != 0)
   {
      if (poly->length == 0)
      {
         // input is zero, so output is zero too
         for (unsigned long i = 0; i < length; i++)
            ZmodF_zero(poly->coeffs[i], poly->n);
      }
      else
      {
         _ZmodFpoly_FFT(poly->coeffs, poly->depth, 1, poly->length, length,
                        (2 * poly->n * FLINT_BITS_PER_LIMB) >> poly->depth,
                        poly->n, poly->scratch);
      }
   }

   poly->length = length;
}


void ZmodFpoly_negacyclic_IFFT(ZmodFpoly_t poly)
{
   // todo: I *think* this code is wrong
   abort();
      
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (poly->length != 0)
   {
      _ZmodFpoly_IFFT(poly->coeffs, poly->depth, 1, poly->length, poly->length,
                      0, (2 * poly->n * FLINT_BITS_PER_LIMB) >> poly->depth,
                      poly->n, poly->scratch);
   }
}


void ZmodFpoly_negacyclic_convolution(ZmodFpoly_t res,
                                      ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long length = x->length + y->length - 1;
   unsigned long size = 1UL << res->depth;
   if (length > size)
      length = size;
   
   ZmodFpoly_negacyclic_FFT(x, length);
   if (x != y)    // take care of aliasing
      ZmodFpoly_negacyclic_FFT(y, length);
      
   ZmodFpoly_mul(res, x, y);
   ZmodFpoly_negacyclic_IFFT(res);
}


// end of file ****************************************************************
