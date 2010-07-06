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
/******************************************************************************

mpq_mat.c: Matrices over Q, implemented as an array of mpq_t's
           Not intended to be efficient

Copyright (C) 2010 Andy Novocin, Max Flander

******************************************************************************/

#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "mpq_mat.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/


void mpq_mat_init(mpq_mat_t mat, ulong r, ulong c)
{
   if ((r) && (c)) mat->entries = (mpq_t*) flint_heap_alloc_bytes(r*c*sizeof(mpq_t));
	else mat->entries = NULL;

   long i;
   for (i = 0; i < r*c; i++)
	{
	   mpq_init(mat->entries[i]);
	}
	mat->r = r;
	mat->c = c;
}

void mpq_mat_clear(mpq_mat_t mat)
{
   long i;
   for (i = 0; i < mat->r*mat->c; i++)
      mpq_clear(mat->entries[i]);

   if (mat->entries) flint_heap_free(mat->entries);
	mat->entries = NULL;
	mat->r = 0;
	mat->c = 0;
}

// *************** end of file
