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

NTL-interface.h: Header file for NTL-interface.cpp

Copyright (C) 2007, William Hart

*****************************************************************************/

#ifndef FLINT_NTL_INT_H
#define FLINT_NTL_INT_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

NTL_CLIENT

/*
   Returns the number of limbs taken up by an NTL ZZ
*/

unsigned long ZZ_limbs(ZZ z);

/* 
   Convert an NTL ZZ to an fmpz_t
   Assumes the fmpz_t has already been allocated to have sufficient space
*/

void ZZ_to_fmpz(fmpz_t output, ZZ z);

/*
   Convert an fmpz_t to an NTL ZZ
*/

void fmpz_to_ZZ(ZZ& output, fmpz_t z);

#endif
