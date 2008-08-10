/*============================================================================

    F_mpz_poly.c: Polynomials over Z (FLINT 2.0 polynomials)

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

===============================================================================*/

#include <stdint.h>
#include <string.h>
#include <math.h>

#include "mpz_poly.h"
#include "flint.h"
#include "F_mpz_poly.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "long_extras.h"
#include "zmod_poly.h"

/*===============================================================================

	Memory management

================================================================================*/

void _F_mpz_poly_mpz_coeffs_new(F_mpz_poly_t poly)
{
	if (poly->mpz_length == poly->mpz_alloc) // time to allocate some more mpz_t's
	{
	   if (poly->mpz_alloc) 
		   poly->mpz_coeffs = (__mpz_struct*) flint_heap_realloc_bytes(poly->mpz_coeffs, (poly->mpz_alloc + MPZ_BLOCK)*sizeof(__mpz_struct));
		else
		   poly->mpz_coeffs = (__mpz_struct*) flint_heap_alloc_bytes(MPZ_BLOCK*sizeof(__mpz_struct));
		
		// initialise the new mpz_t's
		for (ulong i = poly->mpz_alloc; i < poly->mpz_alloc + MPZ_BLOCK; i++)
			mpz_init(poly->mpz_coeffs + i);
		poly->mpz_alloc += MPZ_BLOCK;
	}

	poly->mpz_length++;
}

void _F_mpz_poly_mpz_coeffs_clear(F_mpz_poly_t poly)
{
   for (ulong i = 0; i < poly->mpz_alloc; i++)
	   mpz_clear(poly->mpz_coeffs + i);

	flint_heap_free(poly->mpz_coeffs);
	poly->mpz_coeffs = NULL;
	poly->mpz_alloc = 0;
	poly->mpz_length = 0;
}

void F_mpz_poly_init(F_mpz_poly_t poly)
{
   poly->coeffs = NULL;
   poly->mpz_coeffs = NULL;
   
   poly->alloc = 0;
   poly->length = 0;
   poly->mpz_alloc = 0;
   poly->mpz_length = 0;
}

void F_mpz_poly_init2(F_mpz_poly_t poly, const ulong alloc)
{
   if (alloc)
   {
      poly->coeffs = (mp_limb_t *) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
   }
   else poly->coeffs = NULL;
   
   poly->alloc = alloc;
   poly->length = 0;
   poly->mpz_alloc = 0;
   poly->mpz_length = 0;
}

void F_mpz_poly_realloc(F_mpz_poly_t poly, const ulong alloc)
{
   if (!alloc) // clear up
   {
         F_mpz_poly_clear(poly);
			return;
   }  
   
	if (poly->alloc) 
	{
		poly->coeffs = (mp_limb_t*) flint_heap_realloc(poly->coeffs, alloc);
		if (alloc > poly->alloc)
		   F_mpn_clear(poly->coeffs + poly->alloc, alloc - poly->alloc);
	}
   else 
	{
		poly->coeffs = (mp_limb_t*) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
	}
   
   poly->alloc = alloc;
   
   F_mpz_poly_truncate(poly, alloc); // truncate actual data if necessary   
}

void F_mpz_poly_fit_length(F_mpz_poly_t poly, const ulong length)
{
   ulong alloc = length;
   
	if (alloc <= poly->alloc) return;

   if (alloc < 2*poly->alloc) alloc = 2*poly->alloc;
   
   F_mpz_poly_realloc(poly, alloc);
}

void F_mpz_poly_clear(F_mpz_poly_t poly)
{
   if (poly->coeffs) flint_heap_free(poly->coeffs);
   if (poly->mpz_coeffs) _F_mpz_poly_mpz_coeffs_clear(poly);
   poly->coeffs = NULL;
	poly->alloc = 0;
	poly->length = 0;
}

/*===============================================================================

	Normalisation

================================================================================*/

void _F_mpz_poly_normalise(F_mpz_poly_t poly)
{
   ulong length = poly->length;
	
	while ((length) && (!poly->coeffs[length - 1])) length--;

	poly->length = length;
}

/*===============================================================================

	Coefficient operations

================================================================================*/

void _F_mpz_demote_val(F_mpz_poly_t poly, const ulong coeff)
{
   __mpz_struct * mpz_ptr = poly->mpz_coeffs + COEFF_TO_OFF(poly->coeffs[coeff]);

	long size = mpz_ptr->_mp_size;
	
	if (size == 0L)
	{
		poly->coeffs[coeff] = 0;
	} else if (size == 1L)
	{
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = uval;
	} else if (size == -1L)
   {
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = -uval;
	}
}

void _F_mpz_set_si(F_mpz_poly_t poly, ulong coeff, const long val)
{
   ulong uval = FLINT_ABS(val);
	
	if (uval > COEFF_MAX) // val is large
	{
		__mpz_struct * mpz_coeff = _F_mpz_promote(poly, coeff);
		mpz_set_si(mpz_coeff, val);
	} else poly->coeffs[coeff] = val; // val is small
}

void _F_mpz_get_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong coeff)
{
	mp_limb_t c = poly->coeffs[coeff];

	if (!COEFF_IS_MPZ(c)) mpz_set_si(x, c);
	else mpz_set(x, poly->mpz_coeffs + COEFF_TO_OFF(c));
}

void _F_mpz_set_mpz(F_mpz_poly_t poly, ulong coeff, const mpz_t x)
{
   long size = x->_mp_size;
	
	if (size == 0L)
	{
		poly->coeffs[coeff] = 0;
	} else if (size == 1L)
	{
	   ulong uval = mpz_get_ui(x);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = uval;
		else 
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(poly, coeff);
			mpz_set_ui(mpz_ptr, uval);
		}
	} else if (size == -1L)
   {
	   ulong uval = mpz_get_ui(x);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = -uval;
		else 
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(poly, coeff);
			mpz_set_ui(mpz_ptr, uval);
			mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // x is large
	{
		__mpz_struct * mpz_ptr = _F_mpz_promote(poly, coeff);
		mpz_set(mpz_ptr, x);
	}			
}

void _F_mpz_set(F_mpz_poly_t poly1, ulong coeff1, const F_mpz_poly_t poly2, const ulong coeff2)
{
   ulong c = poly2->coeffs[coeff2];

	if (!COEFF_IS_MPZ(c)) // coeff is small
	{
		poly1->coeffs[coeff1] = c;
	} else // coeff is large
	{
	   __mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
		mpz_set(mpz_ptr, poly2->mpz_coeffs + COEFF_TO_OFF(c));
	}
}

void _F_mpz_add(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
					                                 const F_mpz_poly_t poly2, const ulong coeff2) 
{
	mp_limb_t c1 = poly1->coeffs[coeff1];
	mp_limb_t c2 = poly2->coeffs[coeff2];
	
	if (!COEFF_IS_MPZ(c1))
	{
	   if (!COEFF_IS_MPZ(c2)) // both coefficients are small
		{
			_F_mpz_set_si(res, coeff3, c1 + c2);
		} else // c1 is small, c2 is large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + COEFF_TO_OFF(c2);
			if ((long) c1 < 0L) mpz_sub_ui(mpz3, mpz2, -c1);	
		   else mpz_add_ui(mpz3, mpz2, c1);
			_F_mpz_demote_val(res, coeff3);
		}
	} else
	{
		if (!COEFF_IS_MPZ(c2)) // c1 is large, c2 is small
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + COEFF_TO_OFF(c1);
			if ((long) c2 < 0L) mpz_sub_ui(mpz3, mpz1, -c2);	
			else mpz_add_ui(mpz3, mpz1, c2);
			_F_mpz_demote_val(res, coeff3);
		} else // c1 and c2 are large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + COEFF_TO_OFF(c1);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + COEFF_TO_OFF(c2);
			mpz_add(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(res, coeff3);
		}
	}
}

/*===============================================================================

	Conversions

================================================================================*/

void mpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const mpz_poly_t m_poly)
{
	F_mpz_poly_fit_length(F_poly, m_poly->length);

	F_poly->length = m_poly->length;
   
	for (ulong i = 0; i < m_poly->length; i++)
		_F_mpz_set_mpz(F_poly, i, m_poly->coeffs[i]);
}

void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly)
{
	mpz_poly_ensure_alloc(m_poly, F_poly->length);

   m_poly->length = F_poly->length;
   
   for (ulong i = 0; i < F_poly->length; i++)
	   _F_mpz_get_mpz(m_poly->coeffs[i], F_poly, i);
}

/*===============================================================================

	Addition/subtraction

================================================================================*/

void F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);
	ulong shorter = FLINT_MIN(poly1->length, poly2->length);

   F_mpz_poly_fit_length(res, longer);
	
   for (ulong i = 0; i < shorter; i++) // add up to the length of the shorter poly
      _F_mpz_add(res, i, poly1, i, poly2, i);   
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (ulong i = shorter; i < poly1->length; i++)
         _F_mpz_set(res, i, poly1, i);

   if (poly2 != res) // copy any remaining coefficients from poly2
      for (ulong i = shorter; i < poly2->length; i++)
         _F_mpz_set(res, i, poly2, i);
   
   if (poly1->length == poly2->length)
   {
      res->length = poly1->length;
      _F_mpz_poly_normalise(res); // there may have been cancellation
   } else
      res->length = longer;
}


